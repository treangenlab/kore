use crate::cli::*;
use crate::util::*;
use crate::consts::*;
use crate::lcb::*;
use crate::replicates::*;
use crate::cluster::*;
use crate::taxonomy::{Lineage, Taxonomy, slug};
use crate::aggregate::{aggregate, checksum};

use log::*;

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::fs;
use std::fs::File;
use std::io;
use std::thread;
use std::time::Duration;
use std::io::{BufRead, BufReader};
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use flate2::read::MultiGzDecoder;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;

//checks of the build parameters (before any species is built)
fn check_args(args: &BuildArgs) {
    let fail = |m: &str| { error!("{}", m); std::process::exit(1); };
    if args.n_genomes < 2 { fail("--n-genomes must be at least 2"); }
    if args.replicates == 0 { fail("--replicates must be at least 1"); }
    if !(0.0..=1.0).contains(&args.max_diversity) { fail("--max-diversity must be between 0 and 1"); }
    if !(0.0..=100.0).contains(&args.replicate_identity) { fail("--replicate-identity must be between 0 and 100"); }
    if args.threads == 0 { fail("--threads must be at least 1"); }
    if args.cores_at_once == 0 { fail("--cores-at-once must be at least 1"); }
    if !(0.0 < args.cluster_ani && args.cluster_ani <= 100.0) { fail("--cluster-ani must be above 0 and at most 100"); }
    if let Some(p) = &args.skani_ani {
        if !Path::new(p).is_file() { fail(&format!("--skani-ani file {} does not exist", p)); }
    }
}

fn check_tools(args: &BuildArgs) {
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
}

fn shuffle_genomes(mut genomes: Vec<String>, seed: u64) -> Vec<String> {
    genomes.shuffle(&mut StdRng::seed_from_u64(seed));
    genomes
}

//genomes in sampling order with their cluster (None when not clustered), and the cluster counts at the
//REPORT_CLUSTER_ANI cuts (empty when not clustered)
fn order_genomes(genomes: Vec<String>, args: &BuildArgs, out_dir: &Path, skani_dir: &Path) -> (Vec<(String, Option<usize>)>, Vec<(f64, usize)>) {
    match args.sampling {
        Sampling::Cluster => {
            let (ordered, counts) = cluster_order(&genomes, args.cluster_ani, args.skani_ani.as_deref().map(Path::new), out_dir,
                                                  skani_dir, args.threads, args.seed, args.replicate_grouping == Grouping::Clade);
            (ordered.into_iter().map(|(g, c)| (g, Some(c))).collect(), counts)
        }
        Sampling::Random => (shuffle_genomes(genomes, args.seed).into_iter().map(|g| (g, None)).collect(), Vec::new()),
    }
}

//info.tsv key for the cluster count at an ANI cut, e.g. clusters_99.9 (also the species.tsv column)
pub fn clusters_key(ani: f64) -> String {
    format!("clusters_{}", ani)
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

//one species to build: its folder name, GTDB lineage (None: unlisted) and genomes
struct Job {
    name: String,
    lineage: Option<Lineage>,
    genomes: Vec<String>,
}

fn absolute(p: &str) -> String {
    std::path::absolute(p).map(|x| x.to_string_lossy().into_owned()).unwrap_or_else(|_| p.to_string())
}

//--species-table: genome path <TAB> species per line (a header line is skipped)
fn read_species_table(path: &str) -> Vec<(String, String)> {
    let file = File::open(path).unwrap_or_else(|e| {
        error!("Unable to open --species-table {}: {}", path, e);
        std::process::exit(1);
    });
    let base = Path::new(path).parent().unwrap_or(Path::new("."));
    let mut rows = Vec::new();
    let mut missing = Vec::new();
    for (i, line) in BufReader::new(file).lines().map_while(Result::ok).enumerate() {
        let f: Vec<&str> = line.split('\t').map(|x| x.trim()).collect();
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        if f.len() < 2 || f[1].is_empty() {
            error!("--species-table line {}: expected <genome path> <TAB> <species>: {}", i + 1, line);
            std::process::exit(1);
        }
        let g = if Path::new(f[0]).exists() { f[0].to_string() } else { base.join(f[0]).to_string_lossy().into_owned() };
        if !Path::new(&g).exists() {
            if i == 0 { continue; } //header
            missing.push(f[0].to_string());
            continue;
        }
        rows.push((absolute(&g), f[1].to_string()));
    }
    if !missing.is_empty() {
        error!("{} genome file(s) in --species-table do not exist, e.g. {}", missing.len(), missing[..missing.len().min(3)].join(", "));
        std::process::exit(1);
    }
    rows
}

//species to build, names checked against GTDB (all problems reported at once)
fn collect_jobs(args: &BuildArgs, tax: &Taxonomy) -> Vec<Job> {
    let pairs: Vec<(String, String)> = match &args.species_table {
        Some(t) => read_species_table(t),
        None => {
            let sp = args.species.clone().unwrap_or_else(|| {
                error!("Give the genomes' species with --species (or many species with --species-table)");
                std::process::exit(1);
            });
            collect_genomes(&args.genomes, &args.file_input).into_iter().map(|g| (absolute(&g), sp.clone())).collect()
        }
    };
    let mut names: Vec<&String> = pairs.iter().map(|(_, s)| s).collect::<HashSet<_>>().into_iter().collect();
    names.sort();
    let mut resolved: HashMap<String, (String, Option<Lineage>)> = HashMap::new();
    let mut problems = Vec::new();
    for n in names {
        match tax.resolve(n) {
            Ok(l) => { resolved.insert(n.clone(), (slug(&l.name), Some(l.clone()))); }
            Err(e) if args.allow_unlisted => {
                warn!("{}; building it as '{}' with no relatives (--allow-unlisted)", e, slug(n));
                resolved.insert(n.clone(), (slug(n), None));
            }
            Err(e) => problems.push(e),
        }
    }
    if !problems.is_empty() {
        for p in &problems {
            error!("{}", p);
        }
        error!("{} species name(s) not found in GTDB {}: fix them (or use --allow-unlisted for species not in GTDB)",
               problems.len(), GTDB_RELEASE);
        std::process::exit(1);
    }
    let mut jobs: HashMap<String, Job> = HashMap::new();
    for (g, n) in pairs {
        let (name, lineage) = resolved[&n].clone();
        jobs.entry(name.clone()).or_insert(Job { name, lineage, genomes: Vec::new() }).genomes.push(g);
    }
    let min = if args.species_table.is_some() { args.min_genomes.max(2) } else { 2 };
    let mut jobs: Vec<Job> = jobs.into_values().filter_map(|mut j| {
        j.genomes.sort();
        j.genomes.dedup();
        if j.genomes.len() < min {
            warn!("{}: {} genome(s), fewer than {}; skipped", j.name, j.genomes.len(), min);
            return None;
        }
        Some(j)
    }).collect();
    jobs.sort_by(|a, b| a.name.cmp(&b.name));
    if jobs.is_empty() {
        error!("No species with enough genomes to build");
        std::process::exit(1);
    }
    jobs
}

//writes <db>/cores/<sp>/genomes.txt and info.tsv; false if the core is already built from the same genomes
fn prepare(job: &Job, dir: &Path) -> bool {
    fs::create_dir_all(dir).unwrap();
    let list = job.genomes.join("\n") + "\n";
    let same = fs::read_to_string(dir.join(SPECIES_GENOMES)).is_ok_and(|old| old == list);
    if same && dir.join(CORE_FASTA).exists() {
        return false;
    }
    if dir.join(CORE_FASTA).exists() {
        warn!("{}: genome list changed, rebuilding its core", job.name);
        for f in [CORE_FASTA, CORE_CSV, CALIBRATION_LCBS, CALIBRATION_GENOMES] {
            let _ = fs::remove_file(dir.join(f));
        }
    }
    fs::write(dir.join(SPECIES_GENOMES), list).unwrap();
    let mut info = format!("species\t{}\nbuild_genomes\t{}\n", job.name, job.genomes.len());
    if let Some(l) = &job.lineage {
        for (k, v) in [("gtdb_species", &l.name), ("domain", &l.domain), ("phylum", &l.phylum), ("class", &l.class),
                       ("order", &l.order), ("family", &l.family), ("genus", &l.genus)] {
            info.push_str(&format!("{}\t{}\n", k, v));
        }
        info.push_str(&format!("gtdb_release\t{}\n", GTDB_RELEASE));
    }
    fs::write(dir.join(SPECIES_INFO), info).unwrap();
    true
}

//the build options passed on to each species' build (kore build --species-dir <dir>)
fn species_command(args: &BuildArgs, dir: &Path, threads: usize) -> Command {
    let mut c = Command::new(std::env::current_exe().unwrap());
    c.arg("build").arg("--species-dir").arg(dir).arg("--threads").arg(threads.to_string())
        .arg("--n-genomes").arg(args.n_genomes.to_string()).arg("--replicates").arg(args.replicates.to_string())
        .arg("--seed").arg(args.seed.to_string())
        .arg("--sampling").arg(match args.sampling { Sampling::Cluster => "cluster", Sampling::Random => "random" })
        .arg("--cluster-ani").arg(args.cluster_ani.to_string())
        .arg("--replicate-grouping").arg(match args.replicate_grouping { Grouping::Dealt => "dealt", Grouping::Clade => "clade" })
        .arg("--min-length").arg(args.min_length.to_string()).arg("--max-diversity").arg(args.max_diversity.to_string())
        .arg("--replicate-identity").arg(args.replicate_identity.to_string());
    if let Some(s) = &args.skani_ani { c.arg("--skani-ani").arg(s); }
    if args.keep_intermediate { c.arg("--keep-intermediate"); }
    if args.debug { c.arg("--debug"); }
    if args.verbose { c.arg("--verbose"); }
    //the child logs to <species>/kore.log; its terminal output also goes to a file, so a failure before logging
    //starts (or a crash) still leaves a trace
    match File::create(dir.join(BUILD_OUTPUT)) {
        Ok(f) => { c.stdout(f.try_clone().map_or(std::process::Stdio::null(), Into::into)).stderr(f); }
        Err(_) => { c.stdout(std::process::Stdio::null()).stderr(std::process::Stdio::null()); }
    }
    c
}

//last lines of a failed species build's terminal output, for the error message
fn tail(path: &Path, n: usize) -> String {
    let text = fs::read_to_string(path).unwrap_or_default();
    let lines: Vec<&str> = text.lines().filter(|l| !l.trim().is_empty()).collect();
    lines[lines.len().saturating_sub(n)..].iter().map(|l| format!("\n    {}", l)).collect()
}

pub fn build(args: BuildArgs) {
    if let Some(dir) = args.species_dir.clone() {
        build_species(&args, Path::new(&dir));
        return;
    }
    let db = PathBuf::from(&args.output);
    //never build over an existing db by accident (checked before the log file, which would be replaced)
    let exists = db.join(DB_SPECIES).exists()
        || fs::read_dir(db.join(DB_CORES_DIR)).is_ok_and(|mut r| r.next().is_some());
    if exists && !args.force && !args.aggregate_only {
        eprintln!("kore db already exists at {}, please use --force / -f (species already built from the same genomes are \
                   kept; new or changed species are built and the db is re-aggregated)", db.display());
        std::process::exit(1);
    }
    fs::create_dir_all(db.join(DB_CORES_DIR)).unwrap_or_else(|e| {
        error!("Unable to create {}: {}", db.display(), e);
        std::process::exit(1);
    });
    init_logging(log_level(args.verbose, args.debug), &db.join(LOG_FILE));
    print_banner();
    check_args(&args);
    if args.aggregate_only {
        aggregate(&db, args.threads, args.mask_genomes, args.mask_identity, args.seed);
        return;
    }
    check_tools(&args);
    let tax = Taxonomy::load();
    let jobs = collect_jobs(&args, &tax);
    let todo: Vec<&Job> = jobs.iter().filter(|j| prepare(j, &db.join(DB_CORES_DIR).join(&j.name))).collect();
    info!("kore v{}: {} species ({} genomes); {} to build, {} already built in {}", KORE_VERSION, jobs.len(),
          jobs.iter().map(|j| j.genomes.len()).sum::<usize>(), todo.len(), jobs.len() - todo.len(), db.display());

    let at_once = args.cores_at_once.min(todo.len()).max(1);
    let threads_per = (args.threads / at_once).max(1);
    if !todo.is_empty() {
        info!("Building {} core(s) at a time with {} thread(s) each; each species logs to cores/<species>/{}", at_once, threads_per, LOG_FILE);
    }
    let next = AtomicUsize::new(0);
    let failed: Mutex<Vec<String>> = Mutex::new(Vec::new());
    let done = AtomicUsize::new(0);
    thread::scope(|s| {
        for _ in 0..at_once {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= todo.len() {
                    break;
                }
                let dir = db.join(DB_CORES_DIR).join(&todo[i].name);
                let ok = species_command(&args, &dir, threads_per).status().is_ok_and(|st| st.success())
                    && dir.join(CORE_FASTA).exists();
                let d = done.fetch_add(1, Ordering::SeqCst) + 1;
                if ok {
                    let _ = fs::remove_file(dir.join(BUILD_OUTPUT));
                    info!("[{}/{}] {}: core built", d, todo.len(), todo[i].name);
                } else {
                    error!("[{}/{}] {}: build failed, see {} and {}; last output:{}", d, todo.len(), todo[i].name,
                           dir.join(LOG_FILE).display(), dir.join(BUILD_OUTPUT).display(), tail(&dir.join(BUILD_OUTPUT), 5));
                    failed.lock().unwrap().push(todo[i].name.clone());
                }
            });
        }
    });
    let failed = failed.into_inner().unwrap();
    if !failed.is_empty() {
        warn!("{} species failed to build and are left out: {}", failed.len(), failed.join(", "));
        if failed.len() == todo.len() && jobs.len() == todo.len() {
            error!("No core was built");
            std::process::exit(1);
        }
    }
    if args.skip_aggregate {
        info!("Cores built; aggregate them with: kore build --aggregate-only -o {}", db.display());
    } else {
        aggregate(&db, args.threads, args.mask_genomes, args.mask_identity, args.seed);
    }
}

//(shortest, longest) of all build genomes: the length range for report.csv kore_length_flag
pub fn genome_size_range(genomes: &[String]) -> (f64, f64) {
    let l: Vec<f64> = genomes.iter().map(|g| read_fasta(Path::new(g)).iter().map(|c| c.seq.len()).sum::<usize>() as f64).collect();
    (l.iter().cloned().fold(f64::INFINITY, f64::min), l.iter().cloned().fold(0.0, f64::max))
}

//median genome length of the genomes used for the core
pub fn genome_size(genomes: &[String]) -> f64 {
    let mut l: Vec<f64> = genomes.iter().map(|g| read_fasta(Path::new(g)).iter().map(|c| c.seq.len()).sum::<usize>() as f64).collect();
    percentile(&mut l, 50.0)
}

//calibration: up to CALIBRATION_N cluster-diverse build genomes mapped back onto the core; per LCB, how identical its
//best hit is in the species' own genomes (for a per-core species-match cutoff), and each genome's match score
pub fn calibrate(core: &[FastaRecord], genomes: &[String], dir: &Path, threads: usize) {
    let n = genomes.len().min(threads).max(1);
    let per = (threads / n).max(1);
    let next = AtomicUsize::new(0);
    let out: Mutex<Vec<(String, Vec<Option<f64>>)>> = Mutex::new(Vec::new());
    thread::scope(|s| {
        for _ in 0..n {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= genomes.len() {
                    break;
                }
                let hits = crate::estimate::map_genome(&genomes[i], core, per);
                out.lock().unwrap().push((genomes[i].clone(), crate::estimate::best_lcb_identities(core, &hits)));
            });
        }
    });
    let mut per_genome = out.into_inner().unwrap();
    per_genome.sort_by(|a, b| a.0.cmp(&b.0));
    let mut w = String::from("lcb\tgenomes\tfound\tidentity_min\tidentity_p5\tidentity_median\tidentity_mean\n");
    for (k, lcb) in core.iter().enumerate() {
        let mut ids: Vec<f64> = per_genome.iter().filter_map(|(_, v)| v[k]).collect();
        let empty = ids.is_empty();
        let mean = ids.iter().sum::<f64>() / ids.len().max(1) as f64;
        let min = ids.iter().copied().fold(f64::INFINITY, f64::min);
        let (p5, med) = if empty { (0.0, 0.0) } else { (percentile(&mut ids, 5.0), percentile(&mut ids, 50.0)) };
        let f = |x: f64| if empty { "NA".to_string() } else { format!("{:.3}", x) };
        w.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", lcb.id(), per_genome.len(), ids.len(), f(min), f(p5), f(med), f(mean)));
    }
    fs::write(dir.join(CALIBRATION_LCBS), w).unwrap();
    let mut g = String::from("genome\tfound_lcbs\tspecies_match_lcbs_percent\n");
    for (name, v) in &per_genome {
        let ids: Vec<f64> = v.iter().flatten().copied().collect();
        let pct = 100.0 * ids.iter().filter(|&&x| x >= SPECIES_MATCH_IDENTITY).count() as f64 / ids.len().max(1) as f64;
        g.push_str(&format!("{}\t{}\t{:.2}\n", name, ids.len(), pct));
    }
    fs::write(dir.join(CALIBRATION_GENOMES), g).unwrap();
    let mut pcts: Vec<f64> = per_genome.iter().map(|(_, v)| {
        let ids: Vec<f64> = v.iter().flatten().copied().collect();
        100.0 * ids.iter().filter(|&&x| x >= SPECIES_MATCH_IDENTITY).count() as f64 / ids.len().max(1) as f64
    }).collect();
    info!("Calibration: {} build genomes mapped back; species-match score min {:.1}%, median {:.1}%", per_genome.len(),
          pcts.iter().copied().fold(f64::INFINITY, f64::min), percentile(&mut pcts, 50.0));
}

//one species' core, built in its own folder of the db (run as its own process by build())
fn build_species(args: &BuildArgs, out_dir: &Path) {
    init_logging(log_level(args.verbose, args.debug), &out_dir.join(LOG_FILE));
    print_banner();
    check_args(args);
    check_tools(args);
    let genomes = read_file_list(out_dir.join(SPECIES_GENOMES).to_str().unwrap());
    if genomes.len() < 2 {
        error!("At least 2 genomes are required to build a core genome (got {})", genomes.len());
        std::process::exit(1);
    }
    info!("kore v{} building the core of {} from {} genomes", KORE_VERSION,
          out_dir.file_name().unwrap().to_string_lossy(), genomes.len());

    let tmp_dir = out_dir.join(INTERMEDIATE_DIR);
    let (ordered, cluster_counts) = order_genomes(genomes, args, out_dir, &tmp_dir.join("skani"));
    if !cluster_counts.is_empty() {
        info!("Clusters at {}", cluster_counts.iter().map(|(a, k)| format!("{}% ANI: {}", a, k)).collect::<Vec<_>>().join(", "));
    }
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
    //cluster-diverse genomes in sampling order, for the genome size and the calibration
    let mut sampled: Vec<String> = Vec::new();
    for k in 0..sets.iter().map(|s| s.len()).max().unwrap_or(0) {
        for set in &sets {
            if let Some(g) = set.get(k) && !sampled.contains(g) {
                sampled.push(g.clone());
            }
        }
    }

    let unzipped = decompress_genomes(&sets, &tmp_dir.join("genomes"));
    let replicate_dirs = run_replicates(&unzipped, &tmp_dir.join("parsnp"), args.threads);
    info!("Finished parsnp replicates");

    let consensus_fastas = build_consensus(&replicate_dirs, &tmp_dir.join("consensus"), args.min_length, args.max_diversity);
    let mut core = merge_replicates(&consensus_fastas, &tmp_dir.join("blast"), args.replicate_identity, args.min_length, args.threads);
    rename_core(&mut core);
    write_core_csv(&core, &out_dir.join(CORE_CSV));
    info!("Final core: {} LCBs ({} bp)", core.len(), core.iter().map(|r| r.seq.len()).sum::<usize>());

    let size = genome_size(&sampled);
    let (size_min, size_max) = genome_size_range(&genomes);
    let mut info = fs::read_to_string(out_dir.join(SPECIES_INFO)).unwrap_or_default();
    info.push_str(&format!("genome_size_min\t{:.0}\ngenome_size_max\t{:.0}\n", size_min, size_max));
    info.push_str(&format!("genome_size\t{:.0}\ncore_lcbs\t{}\ncore_bp\t{}\ngenomes_checksum\t{}\n", size, core.len(),
                           core.iter().map(|r| r.seq.len()).sum::<usize>(), checksum(&out_dir.join(SPECIES_GENOMES))));
    for (a, k) in &cluster_counts {
        info.push_str(&format!("{}\t{}\n", clusters_key(*a), k));
    }
    fs::write(out_dir.join(SPECIES_INFO), info).unwrap();
    let calib: Vec<String> = sampled.iter().take(CALIBRATION_N).cloned().collect();
    calibrate(&core, &calib, out_dir, args.threads);

    if !args.keep_intermediate {
        remove_intermediate(&tmp_dir);
    }
    //written last: its presence marks the species as built
    write_fasta(&core.iter().collect::<Vec<_>>(), &out_dir.join(CORE_FASTA));
    info!("Wrote {}", out_dir.join(CORE_FASTA).display());
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
