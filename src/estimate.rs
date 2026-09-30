use crate::cli::*;
use crate::util::*;
use crate::consts::*;

use log::*;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use minimap2::{Aligner, Built};

fn check_args(args: &EstimateArgs) -> Vec<String> {
    // Log file lives inside the output folder and mirrors the terminal.
    let log_path = Path::new(&args.output).join(LOG_FILE);
    init_logging(log_level(args.verbose, args.debug), &log_path);
    print_banner();

    let genomes = collect_genomes(&args.genomes, &args.file_input);

    if !Path::new(&args.core).exists() {
        error!("kore core does not exist: {}", args.core);
        std::process::exit(1);
    }

    if args.threads == 0 {
        error!("--threads must be at least 1");
        std::process::exit(1);
    }

    genomes
}

//"GCF_000005845.2_genomic.fna.gz" -> "GCF_000005845.2_genomic"
fn genome_name(path: &str) -> String {
    let name = Path::new(path).file_name().unwrap().to_string_lossy();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    Path::new(name).file_stem().unwrap().to_string_lossy().into_owned()
}

//minimap2 -c -x asm20 --secondary=yes -p 0.5 -N 20 -m 30 -s 40, indexed on the genome
fn build_aligner(genome: &str, threads: usize) -> Aligner<Built> {
    let mut aligner = Aligner::builder().asm20().with_cigar().with_index_threads(threads);
    aligner.mapopt.pri_ratio = MM_SECONDARY_RATIO;
    aligner.mapopt.best_n = MM_MAX_SECONDARY;
    aligner.mapopt.min_chain_score = MM_MIN_CHAIN_SCORE;
    aligner.mapopt.min_dp_max = MM_MIN_DP_SCORE;

    aligner.with_index(genome, None).unwrap_or_else(|e| {
        error!("Unable to build minimap2 index for {}: {}", genome, e);
        std::process::exit(1);
    })
}

pub struct Hit {
    pub contig: String,
    pub start: usize, //1-based, inclusive, on the genome
    pub end: usize,
    pub lcb_start: usize, //1-based, inclusive, on the LCB
    pub lcb_end: usize,
    pub identity: f64,
}

//alignments of an LCB to the genome passing the identity/length filters
fn map_lcb(aligner: &Aligner<Built>, lcb: &FastaRecord) -> Vec<Hit> {
    let mappings = aligner.map(&lcb.seq, false, false, None, None, None).unwrap_or_else(|e| {
        error!("minimap2 failed to map {}: {}", lcb.id(), e);
        std::process::exit(1);
    });

    mappings.iter()
        .filter(|m| m.block_len >= MIN_ALN_LEN && 100.0 * m.match_len as f64 / m.block_len as f64 >= MIN_ALN_IDENTITY)
        .map(|m| Hit {
            contig: m.target_name.as_ref().map_or(String::new(), |n| n.to_string()),
            start: m.target_start as usize + 1,
            end: m.target_end as usize,
            lcb_start: m.query_start as usize + 1,
            lcb_end: m.query_end as usize,
            identity: 100.0 * m.match_len as f64 / m.block_len as f64,
        })
        .collect()
}

//hits for each LCB in the core, in core order
fn map_genome(genome: &str, core: &[FastaRecord], threads: usize) -> Vec<Vec<Hit>> {
    let aligner = build_aligner(genome, threads);
    core.iter().map(|lcb| map_lcb(&aligner, lcb)).collect()
}

//percent of core bp covered by the hits
fn completeness(core: &[FastaRecord], hits: &[Vec<Hit>]) -> f64 {
    let core_bp: usize = core.iter().map(|lcb| lcb.seq.len()).sum();
    let covered: usize = hits.iter()
        .map(|h| covered_bp(&h.iter().map(|x| (x.lcb_start, x.lcb_end)).collect::<Vec<_>>()))
        .sum();
    100.0 * covered as f64 / core_bp as f64
}

fn create_csv(path: &Path) -> BufWriter<File> {
    let file = File::create(path).unwrap_or_else(|e| {
        error!("Unable to create {}: {}", path.display(), e);
        std::process::exit(1);
    });
    BufWriter::new(file)
}

fn write_completeness(genomes: &[String], core: &[FastaRecord], hits: &[Vec<Vec<Hit>>], path: &Path) {
    let mut writer = create_csv(path);
    writeln!(writer, "genome,completeness").unwrap();
    for (genome, h) in genomes.iter().zip(hits) {
        writeln!(writer, "{},{:.4}", genome_name(genome), completeness(core, h)).unwrap();
    }
    writer.flush().unwrap();
}

//one row per alignment, LCBs with no alignment get a row with 0 coverage
fn write_alignments(genomes: &[String], core: &[FastaRecord], hits: &[Vec<Vec<Hit>>], path: &Path) {
    let mut writer = create_csv(path);
    writeln!(writer, "genome,lcb,contig,start,end,coverage,identity").unwrap();
    for (genome, genome_hits) in genomes.iter().zip(hits) {
        let name = genome_name(genome);
        for (lcb, lcb_hits) in core.iter().zip(genome_hits) {
            if lcb_hits.is_empty() {
                writeln!(writer, "{},{},,,,0.0000,0.0000", name, lcb.id()).unwrap();
            }
            for hit in lcb_hits {
                let coverage = 100.0 * (hit.lcb_end - hit.lcb_start + 1) as f64 / lcb.seq.len() as f64;
                writeln!(writer, "{},{},{},{},{},{:.4},{:.4}", name, lcb.id(), hit.contig, hit.start, hit.end, coverage, hit.identity).unwrap();
            }
        }
    }
    writer.flush().unwrap();
}

pub fn estimate(args: EstimateArgs) {
    let genomes = check_args(&args);
    info!("kore v{} estimating completeness for {} genome(s)", KORE_VERSION, genomes.len());

    let core = read_fasta(Path::new(&args.core));
    info!("Loaded core with {} LCBs ({} bp)", core.len(), core.iter().map(|lcb| lcb.seq.len()).sum::<usize>());

    let n_parallel = genomes.len().min(args.threads);
    let threads_per = (args.threads / n_parallel).max(1);

    let next = AtomicUsize::new(0);
    let results: Mutex<Vec<Vec<Vec<Hit>>>> = Mutex::new((0..genomes.len()).map(|_| Vec::new()).collect());
    thread::scope(|s| {
        for _ in 0..n_parallel {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= genomes.len() {
                    break;
                }
                let hits = map_genome(&genomes[i], &core, threads_per);
                info!("{}: {:.2}% complete", genome_name(&genomes[i]), completeness(&core, &hits));
                results.lock().unwrap()[i] = hits;
            });
        }
    });

    let hits = results.into_inner().unwrap();
    let out_dir = PathBuf::from(&args.output);
    write_completeness(&genomes, &core, &hits, &out_dir.join(COMPLETENESS_CSV));
    write_alignments(&genomes, &core, &hits, &out_dir.join(ALIGNMENTS_CSV));
    info!("Wrote {} and {} to {}", COMPLETENESS_CSV, ALIGNMENTS_CSV, out_dir.display());
}
