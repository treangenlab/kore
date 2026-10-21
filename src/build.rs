use crate::cli::*;
use crate::util::*;
use crate::consts::*;
use crate::lcb::*;
use crate::replicates::*;
use crate::cluster::*;

use log::*;

use std::collections::{HashMap, HashSet};
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

    if !(0.0..=100.0).contains(&args.replicate_identity) {
        error!("--replicate-identity must be between 0 and 100");
        std::process::exit(1);
    }

    if args.threads == 0 {
        error!("--threads must be at least 1");
        std::process::exit(1);
    }

    if !(0.0 < args.cluster_ani && args.cluster_ani <= 100.0) {
        error!("--cluster-ani must be above 0 and at most 100");
        std::process::exit(1);
    }

    if let Some(p) = &args.skani_ani {
        if !Path::new(p).is_file() {
            error!("--skani-ani file {} does not exist", p);
            std::process::exit(1);
        }
    }

    let mut bins = vec![PARSNP_BIN, BLASTN_BIN, MAKEBLASTDB_BIN];
    if args.sampling == Sampling::Cluster && args.skani_ani.is_none() {
        bins.push(SKANI_BIN);
    }
    for bin in bins {
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

//genomes in sampling order with their cluster (None when not clustered)
fn order_genomes(genomes: Vec<String>, args: &BuildArgs, out_dir: &Path, skani_dir: &Path) -> Vec<(String, Option<usize>)> {
    match args.sampling {
        Sampling::Cluster => cluster_order(&genomes, args.cluster_ani, args.skani_ani.as_deref().map(Path::new), out_dir, skani_dir,
                                           args.threads, args.seed, args.replicate_grouping == Grouping::Clade)
            .into_iter().map(|(g, c)| (g, Some(c))).collect(),
        Sampling::Random => shuffle_genomes(genomes, args.seed).into_iter().map(|g| (g, None)).collect(),
    }
}

//replicate, genome and cluster of every genome used
fn write_sampled(sets: &[Vec<String>], cluster_of: &HashMap<String, Option<usize>>, path: &Path) {
    let mut out = String::from("replicate\tgenome\tcluster\n");
    for (i, set) in sets.iter().enumerate() {
        for g in set {
            let c = cluster_of.get(g).copied().flatten().map_or("NA".to_string(), |c| c.to_string());
            out.push_str(&format!("rep_{}\t{}\t{}\n", i, g, c));
        }
        let n_clusters = set.iter().filter_map(|g| cluster_of.get(g).copied().flatten()).collect::<HashSet<_>>().len();
        if n_clusters > 0 {
            info!("rep_{}: {} genomes from {} clusters", i, set.len(), n_clusters);
        }
    }
    fs::write(path, out).unwrap_or_else(|e| {
        error!("Unable to write {}: {}", path.display(), e);
        std::process::exit(1);
    });
}

//replicate sets of n from genomes in cluster order: genomes are dealt to the replicates in turn so
//each gets an equally diverse share; when there are too few for disjoint sets, each replicate is
//topped up with the earliest (most diverse) genomes it does not already have
fn deal_genomes(ordered: &[String], n: usize, replicates: usize) -> Vec<Vec<String>> {
    if ordered.len() <= n {
        return vec![ordered.to_vec()];
    }
    let mut sets: Vec<Vec<String>> = vec![Vec::with_capacity(n); replicates];
    for (k, g) in ordered.iter().enumerate() {
        let set = &mut sets[k % replicates];
        if set.len() < n {
            set.push(g.clone());
        }
    }
    for set in sets.iter_mut() {
        let have: HashSet<String> = set.iter().cloned().collect();
        let extra: Vec<String> = ordered.iter().filter(|g| !have.contains(*g)).take(n - set.len()).cloned().collect();
        set.extend(extra);
    }
    sets
}

//replicate sets of n from genomes in cluster order where cluster ids follow the ANI tree: clusters are split
//into consecutive blocks (one per replicate, about equal numbers of clusters), and each replicate takes one
//genome per cluster of its block in turn until it has n (or its block runs out)
fn clade_genomes(ordered: &[(String, usize)], n: usize, replicates: usize) -> Vec<Vec<String>> {
    if ordered.len() <= n {
        return vec![ordered.iter().map(|(g, _)| g.clone()).collect()];
    }
    let n_clusters = ordered.iter().map(|(_, c)| c + 1).max().unwrap_or(0);
    let block_of = |c: usize| c * replicates / n_clusters.max(1);
    let mut sets: Vec<Vec<String>> = vec![Vec::with_capacity(n); replicates];
    for (g, c) in ordered {  //ordered is round-robin over clusters, so each block is filled one genome per cluster per round
        let set = &mut sets[block_of(*c)];
        if set.len() < n {
            set.push(g.clone());
        }
    }
    sets
}

//replicates with cluster sampling: enough replicates of n that every cluster contributes a genome
//(clusters / n, rounded up), or more if requested
fn cluster_replicates(cluster_of: &HashMap<String, Option<usize>>, n: usize, requested: usize) -> usize {
    let clusters = cluster_of.values().flatten().collect::<HashSet<_>>().len();
    let needed = clusters.div_ceil(n.max(1));
    if needed > requested {
        info!("{} clusters for replicates of {} genomes: running {} replicates so every cluster is sampled", clusters, n, needed);
    }
    needed.max(requested)
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

//majority vote consensus of each LCB, one fasta per replicate. LCBs are dropped when too short or
//more diverse than max(max_diversity, DIVERSITY_PERCENTILE of the diversity of every LCB long enough)
fn build_consensus(replicate_dirs: &[PathBuf], consensus_dir: &Path, min_length: usize, max_diversity: f64) -> Vec<PathBuf> {
    fs::create_dir_all(consensus_dir).unwrap_or_else(|e| {
        error!("Unable to create directory {}: {}", consensus_dir.display(), e);
        std::process::exit(1);
    });

    let per_rep: Vec<(Vec<ConsensusLcb>, usize)> = thread::scope(|s| {
        let handles: Vec<_> = replicate_dirs.iter().enumerate()
            .map(|(i, rep_dir)| s.spawn(move || consensus_lcbs(&rep_dir.join("parsnp.xmfa"), &format!("rep_{}", i), min_length)))
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    let mut diversities: Vec<f64> = per_rep.iter().flat_map(|(lcbs, _)| lcbs.iter().map(|l| l.diversity)).collect();
    let percentile_div = percentile(&mut diversities, DIVERSITY_PERCENTILE);
    let cutoff = max_diversity.max(percentile_div);
    info!("LCB diversity cutoff: {:.4} (max of --max-diversity {} and percentile_{}={:.4} of {} LCBs >= {} bp)",
        cutoff, max_diversity, DIVERSITY_PERCENTILE, percentile_div, diversities.len(), min_length);

    per_rep.iter().enumerate().map(|(i, (lcbs, too_short))| {
        let kept: Vec<&FastaRecord> = lcbs.iter().filter(|l| l.diversity <= cutoff).map(|l| &l.record).collect();
        let fasta = consensus_dir.join(format!("rep_{}.fasta", i));
        write_fasta(&kept, &fasta);
        info!("rep_{}: kept {} consensus LCBs ({} bp), filtered {} too short and {} too diverse",
            i, kept.len(), kept.iter().map(|r| r.seq.len()).sum::<usize>(), too_short, lcbs.len() - kept.len());
        fasta
    }).collect()
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
    let ordered = order_genomes(genomes, &args, &out_dir, &tmp_dir.join("skani"));
    let ordered_clusters: Vec<(String, usize)> = ordered.iter().filter_map(|(g, c)| c.map(|c| (g.clone(), c))).collect();
    let cluster_of: HashMap<String, Option<usize>> = ordered.iter().cloned().collect();
    let genomes: Vec<String> = ordered.into_iter().map(|(g, _)| g).collect();

    let sets = match args.sampling {
        Sampling::Cluster => {
            let r = cluster_replicates(&cluster_of, args.n_genomes, args.replicates);
            match args.replicate_grouping {
                Grouping::Dealt => deal_genomes(&genomes, args.n_genomes, r),
                Grouping::Clade => clade_genomes(&ordered_clusters, args.n_genomes, r),
            }
        }
        Sampling::Random => partition_genomes(&genomes, args.n_genomes, args.replicates),
    };
    info!("Running {} parsnp replicate(s) of {} genomes", sets.len(), sets[0].len());
    write_sampled(&sets, &cluster_of, &out_dir.join(SAMPLED_GENOMES));

    let sets = decompress_genomes(&sets, &tmp_dir.join("genomes"));
    let replicate_dirs = run_replicates(&sets, &tmp_dir.join("parsnp"), args.threads);
    info!("Finished parsnp replicates");

    let consensus_fastas = build_consensus(&replicate_dirs, &tmp_dir.join("consensus"), args.min_length, args.max_diversity);

    let mut core = merge_replicates(
        &consensus_fastas,
        &tmp_dir.join("blast"),
        args.replicate_identity,
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

#[cfg(test)]
mod tests {
    use super::*;

    fn genomes(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("g{}", i)).collect()
    }

    #[test]
    fn deal_disjoint_and_balanced() {
        let sets = deal_genomes(&genomes(10), 3, 3);
        assert_eq!(sets, vec![vec!["g0", "g3", "g6"], vec!["g1", "g4", "g7"], vec!["g2", "g5", "g8"]]);
    }

    #[test]
    fn deal_tops_up_when_too_few() {
        let sets = deal_genomes(&genomes(5), 3, 2);
        assert_eq!(sets[0], vec!["g0", "g2", "g4"]);
        assert_eq!(sets[1], vec!["g1", "g3", "g0"]);
        assert!(sets.iter().all(|s| s.len() == 3 && s.iter().collect::<HashSet<_>>().len() == 3));
    }

    #[test]
    fn deal_single_set_when_all_fit() {
        assert_eq!(deal_genomes(&genomes(3), 5, 3), vec![genomes(3)]);
    }
}

#[cfg(test)]
mod cluster_replicate_tests {
    use super::*;

    //genomes in round-robin cluster order (one per cluster, then a second from each, ...) and their clusters
    fn ordered(sizes: &[usize]) -> (Vec<String>, HashMap<String, Option<usize>>) {
        let mut out = Vec::new();
        for round in 0..*sizes.iter().max().unwrap() {
            for (c, &n) in sizes.iter().enumerate() {
                if round < n {
                    out.push((format!("c{}_{}", c, round), c));
                }
            }
        }
        (out.iter().map(|(g, _)| g.clone()).collect(), out.into_iter().map(|(g, c)| (g, Some(c))).collect())
    }

    #[test]
    fn one_replicate_when_clusters_fit() {
        let (_, of) = ordered(&[3, 3, 3]);
        assert_eq!(cluster_replicates(&of, 4, 1), 1);
        assert_eq!(cluster_replicates(&of, 4, 3), 3);
    }

    #[test]
    fn extra_replicates_cover_every_cluster() {
        //101 clusters, 102 genomes, replicates of 100 (the example in the design)
        let mut sizes = vec![1; 101];
        sizes[0] = 2;
        let (genomes, of) = ordered(&sizes);
        let r = cluster_replicates(&of, 100, 1);
        assert_eq!(r, 2);
        let sets = deal_genomes(&genomes, 100, r);
        assert!(sets.iter().all(|s| s.len() == 100 && s.iter().collect::<HashSet<_>>().len() == 100));
        let covered: HashSet<_> = sets.iter().flatten().filter_map(|g| of[g]).collect();
        assert_eq!(covered.len(), 101);
    }

    #[test]
    fn clade_replicates_take_consecutive_clusters() {
        //4 clusters along the tree (ids 0..3) of 3 genomes each, round-robin order, 2 replicates of 4
        let mut ordered = Vec::new();
        for round in 0..3 {
            for c in 0..4 {
                ordered.push((format!("c{}_{}", c, round), c));
            }
        }
        let sets = clade_genomes(&ordered, 4, 2);
        assert!(sets[0].iter().all(|g| g.starts_with("c0") || g.starts_with("c1")));
        assert!(sets[1].iter().all(|g| g.starts_with("c2") || g.starts_with("c3")));
        assert!(sets.iter().all(|s| s.len() == 4));
    }
}
