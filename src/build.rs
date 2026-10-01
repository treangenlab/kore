use crate::cli::*;
use crate::util::*;
use crate::consts::*;
use crate::lcb::*;
use crate::replicates::*;

use log::*;

use std::path::{Path, PathBuf};
use std::process::Command;
use std::fs;
use std::fs::File;
use std::io;
use std::thread;
use std::time::Duration;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use flate2::read::MultiGzDecoder;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;

fn check_args(args: &BuildArgs) -> Vec<String> {
    // Log file lives inside the output folder and mirrors the terminal.
    let log_path = Path::new(&args.output).join(LOG_FILE);
    init_logging(log_level(args.verbose, args.debug), &log_path);
    print_banner();

    let genomes = collect_genomes(&args.genomes, &args.file_input);
    if genomes.len() < 2 {
        error!("At least 2 genomes are required to build a core genome (got {})", genomes.len());
        std::process::exit(1);
    }

    if args.n_genomes < 2 {
        error!("--n-genomes must be at least 2");
        std::process::exit(1);
    }

    if args.replicates == 0 {
        error!("--replicates must be at least 1");
        std::process::exit(1);
    }

    if !(0.0..=1.0).contains(&args.max_diversity) {
        error!("--max-diversity must be between 0 and 1");
        std::process::exit(1);
    }

    for (name, value) in [("--replicate-identity", args.replicate_identity), ("--replicate-cov", args.replicate_cov)] {
        if !(0.0..=100.0).contains(&value) {
            error!("{} must be between 0 and 100", name);
            std::process::exit(1);
        }
    }

    if args.threads == 0 {
        error!("--threads must be at least 1");
        std::process::exit(1);
    }

    for bin in [PARSNP_BIN, BLASTN_BIN, MAKEBLASTDB_BIN] {
        if !in_path(bin) {
            error!("Unable to find {} on PATH, make sure it is installed (e.g. in the kore conda env)", bin);
            std::process::exit(1);
        }
    }

    genomes
}

fn shuffle_genomes(mut genomes: Vec<String>, seed: u64) -> Vec<String> {
    genomes.shuffle(&mut StdRng::seed_from_u64(seed));
    genomes
}

//split genomes into replicate sets of n, disjoint if possible, otherwise evenly overlapping
fn partition_genomes(genomes: &[String], n: usize, replicates: usize) -> Vec<Vec<String>> {
    let total = genomes.len();
    if total <= n {
        return vec![genomes.to_vec()];
    }

    let replicates = replicates.min(total - n + 1);
    let step = if replicates == 1 { 0 } else { n.min((total - n) / (replicates - 1)) };

    (0..replicates)
        .map(|i| genomes[i * step..i * step + n].to_vec())
        .collect()
}

//parsnp can't read gzipped fasta, so decompress any .gz genomes into genome_dir
fn decompress_genomes(sets: &[Vec<String>], genome_dir: &Path) -> Vec<Vec<String>> {
    fs::create_dir_all(genome_dir).unwrap_or_else(|e| {
        error!("Unable to create directory {}: {}", genome_dir.display(), e);
        std::process::exit(1);
    });

    sets.iter()
        .map(|set| set.iter().map(|g| decompress_genome(g, genome_dir)).collect())
        .collect()
}

fn decompress_genome(genome: &str, genome_dir: &Path) -> String {
    let Some(name) = Path::new(genome).file_name().unwrap().to_str().unwrap().strip_suffix(".gz") else {
        return genome.to_string();
    };

    let out_path = genome_dir.join(name);
    if !out_path.exists() {
        let result = File::open(genome)
            .and_then(|f| io::copy(&mut MultiGzDecoder::new(f), &mut File::create(&out_path)?));
        if let Err(e) = result {
            error!("Unable to decompress {}: {}", genome, e);
            std::process::exit(1);
        }
    }
    out_path.to_string_lossy().into_owned()
}

//runs parsnp, killing it early if another replicate has failed
fn absolute_path(path: &Path) -> PathBuf {
    std::path::absolute(path).unwrap_or_else(|e| {
        error!("Unable to resolve path {}: {}", path.display(), e);
        std::process::exit(1);
    })
}

fn run_parsnp(genomes: &[String], out_dir: &Path, threads: usize, failed: &AtomicBool) -> bool {
    debug!("Running parsnp on {} genomes -> {}", genomes.len(), out_dir.display());

    let log_path = out_dir.with_extension("log");
    let log = File::create(&log_path).unwrap_or_else(|e| {
        error!("Unable to create parsnp log {}: {}", log_path.display(), e);
        std::process::exit(1);
    });

    //parsnp writes allmums.out etc to its cwd, so each run gets its own
    let work_dir = out_dir.with_extension("work");
    fs::create_dir_all(&work_dir).unwrap_or_else(|e| {
        error!("Unable to create directory {}: {}", work_dir.display(), e);
        std::process::exit(1);
    });
    let genomes: Vec<PathBuf> = genomes.iter().map(|g| absolute_path(Path::new(g))).collect();

    let mut child = Command::new(PARSNP_BIN)
        .current_dir(&work_dir)
        .arg("-d").args(&genomes)
        .args(["-r", "!", "-c", "--skip-phylogeny", "--force-overwrite"])
        .arg("-p").arg(threads.to_string())
        .arg("-o").arg(absolute_path(out_dir))
        .stdout(log.try_clone().unwrap())
        .stderr(log)
        .spawn()
        .unwrap_or_else(|e| {
            error!("Unable to run {}: {}", PARSNP_BIN, e);
            std::process::exit(1);
        });

    let status = loop {
        if let Some(status) = child.try_wait().unwrap() {
            break status;
        }
        if failed.load(Ordering::SeqCst) {
            let _ = child.kill();
            let _ = child.wait();
            return false;
        }
        thread::sleep(Duration::from_millis(500));
    };

    if !status.success() {
        error!("parsnp failed for {} ({}), see {}", out_dir.display(), status, log_path.display());
        return false;
    }
    let xmfa = out_dir.join("parsnp.xmfa");
    if !xmfa.exists() {
        error!("parsnp produced no parsnp.xmfa for {}, see {}", out_dir.display(), log_path.display());
        return false;
    }

    //parsnp skips bad genomes but still exits 0
    let n_aligned = xmfa_genome_count(&xmfa);
    if n_aligned != genomes.len() {
        error!("parsnp aligned {} of {} genomes for {}, see {}", n_aligned, genomes.len(), out_dir.display(), log_path.display());
        return false;
    }
    true
}

fn xmfa_genome_count(xmfa: &Path) -> usize {
    xmfa_sequence_files(xmfa).iter().filter(|f| !f.ends_with(".ref")).count()
}

fn run_replicates(sets: &[Vec<String>], parsnp_dir: &Path, threads: usize) -> Vec<PathBuf> {
    fs::create_dir_all(parsnp_dir).unwrap_or_else(|e| {
        error!("Unable to create directory {}: {}", parsnp_dir.display(), e);
        std::process::exit(1);
    });

    let out_dirs: Vec<PathBuf> = (0..sets.len())
        .map(|i| parsnp_dir.join(format!("rep_{}", i)))
        .collect();
    let n_parallel = sets.len().min(MAX_PARALLEL_PARSNP);
    let threads_per = (threads / n_parallel).max(1);
    info!("Running {} parsnp job(s) at a time with {} thread(s) each", n_parallel, threads_per);

    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    thread::scope(|s| {
        for _ in 0..n_parallel {
            s.spawn(|| {
                loop {
                    let i = next.fetch_add(1, Ordering::SeqCst);
                    if i >= sets.len() || failed.load(Ordering::SeqCst) {
                        break;
                    }
                    if !run_parsnp(&sets[i], &out_dirs[i], threads_per, &failed) {
                        failed.store(true, Ordering::SeqCst);
                    }
                }
            });
        }
    });

    if failed.load(Ordering::SeqCst) {
        error!("A parsnp replicate failed, exiting");
        std::process::exit(1);
    }

    out_dirs
}

//majority vote consensus of each LCB, one fasta per replicate
fn build_consensus(replicate_dirs: &[PathBuf], consensus_dir: &Path, min_length: usize, max_diversity: f64) -> Vec<PathBuf> {
    fs::create_dir_all(consensus_dir).unwrap_or_else(|e| {
        error!("Unable to create directory {}: {}", consensus_dir.display(), e);
        std::process::exit(1);
    });

    let fastas: Vec<PathBuf> = (0..replicate_dirs.len())
        .map(|i| consensus_dir.join(format!("rep_{}.fasta", i)))
        .collect();

    thread::scope(|s| {
        for (i, (rep_dir, fasta)) in replicate_dirs.iter().zip(&fastas).enumerate() {
            s.spawn(move || {
                let stats = write_consensus_fasta(&rep_dir.join("parsnp.xmfa"), fasta, &format!("rep_{}", i), min_length, max_diversity);
                info!("rep_{}: kept {} consensus LCBs ({} bp), filtered {} too short and {} too diverse",
                    i, stats.written, stats.total_bp, stats.too_short, stats.too_diverse);
            });
        }
    });

    fastas
}

fn remove_intermediate(dir: &Path) {
    if let Err(e) = fs::remove_dir_all(dir) {
        warn!("Unable to remove intermediate files in {}: {}", dir.display(), e);
    }
}

pub fn build(args: BuildArgs) {
    let genomes = check_args(&args);
    info!("kore v{} building db from {} genomes", KORE_VERSION, genomes.len());

    let out_dir = PathBuf::from(&args.output);
    let tmp_dir = out_dir.join(INTERMEDIATE_DIR);
    let genomes = shuffle_genomes(genomes, args.seed);

    let sets = partition_genomes(&genomes, args.n_genomes, args.replicates);
    info!("Running {} parsnp replicate(s) of {} genomes", sets.len(), sets[0].len());

    let sets = decompress_genomes(&sets, &tmp_dir.join("genomes"));
    let replicate_dirs = run_replicates(&sets, &tmp_dir.join("parsnp"), args.threads);
    info!("Finished parsnp replicates");

    let consensus_fastas = build_consensus(&replicate_dirs, &tmp_dir.join("consensus"), args.min_length, args.max_diversity);

    let mut core = merge_replicates(
        &consensus_fastas,
        &tmp_dir.join("blast"),
        args.replicate_identity,
        args.replicate_cov,
        args.min_length,
        args.threads,
    );
    rename_core(&mut core);
    write_fasta(&core.iter().collect::<Vec<_>>(), &out_dir.join(CORE_FASTA));
    write_core_csv(&core, &out_dir.join(CORE_CSV));
    info!("Final core: {} LCBs ({} bp) written to {}", core.len(), core.iter().map(|r| r.seq.len()).sum::<usize>(), out_dir.join(CORE_FASTA).display());

    if !args.keep_intermediate {
        remove_intermediate(&tmp_dir);
    }
}
