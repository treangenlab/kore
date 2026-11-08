use crate::cli::*;
use crate::util::*;
use crate::consts::*;

use log::*;

use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use minimap2::{Aligner, Built};

//what estimate runs against
enum Target {
    Db(PathBuf),                       //whole kore db: species detection + contamination
    Species(PathBuf, String),          //one species' core from the db
    Core(PathBuf),                     //the user's own core (kore build)
}

fn check_args(args: &EstimateArgs) -> (Vec<String>, Target) {
    // Log file lives inside the output folder and mirrors the terminal.
    let log_path = Path::new(&args.output).join(LOG_FILE);
    init_logging(log_level(args.verbose, args.debug), &log_path);
    print_banner();

    let genomes = collect_genomes(&args.genomes, &args.file_input);
    if args.threads == 0 {
        error!("--threads must be at least 1");
        std::process::exit(1);
    }

    if let Some(core) = &args.core {
        if !Path::new(core).exists() {
            error!("Core does not exist: {}", core);
            std::process::exit(1);
        }
        return (genomes, Target::Core(PathBuf::from(core)));
    }
    let db = match args.db.clone().or_else(|| std::env::var(DB_ENV).ok()) {
        Some(d) => PathBuf::from(d),
        None => {
            error!("No kore db given: use --db <folder> or set ${} (or --core for your own core)", DB_ENV);
            std::process::exit(1);
        }
    };
    if !db.join(DB_SPECIES).exists() {
        error!("Not a kore db (no {} in {})", DB_SPECIES, db.display());
        std::process::exit(1);
    }
    if let Some(sp) = &args.species {
        let core = crate::aggregate::core_path(&db, &crate::taxonomy::slug(sp));
        if !core.exists() {
            error!("Species {} is not in the kore db {} (see {})", sp, db.display(), db.join(DB_SPECIES).display());
            std::process::exit(1);
        }
        return (genomes, Target::Species(core, crate::taxonomy::slug(sp)));
    }
    if !crate::aggregate::index_path(&db).exists() {
        error!("kore db index missing ({} in {})", AGGREGATE_INDEX, db.display());
        std::process::exit(1);
    }
    (genomes, Target::Db(db))
}

//"GCF_000005845.2_genomic.fna.gz" -> "GCF_000005845.2_genomic"
pub fn genome_name(path: &str) -> String {
    let name = Path::new(path).file_name().unwrap().to_string_lossy();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    Path::new(name).file_stem().unwrap().to_string_lossy().into_owned()
}

//minimap2 -c -x asm20 --secondary=yes -p 0.5 -N 20 -m 30 -s 40, indexed on the genome
pub fn build_aligner(genome: &str, threads: usize) -> Aligner<Built> {
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
pub fn map_lcb(aligner: &Aligner<Built>, lcb: &FastaRecord) -> Vec<Hit> {
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
pub fn map_genome(genome: &str, core: &[FastaRecord], threads: usize) -> Vec<Vec<Hit>> {
    let aligner = build_aligner(genome, threads);
    core.iter().map(|lcb| map_lcb(&aligner, lcb)).collect()
}

//percent of core bp covered by the hits
pub fn nucleotide_completeness(core: &[FastaRecord], hits: &[Vec<Hit>]) -> f64 {
    let core_bp: usize = core.iter().map(|lcb| lcb.seq.len()).sum();
    let covered: usize = hits.iter()
        .map(|h| covered_bp(&h.iter().map(|x| (x.lcb_start, x.lcb_end)).collect::<Vec<_>>()))
        .sum();
    100.0 * covered as f64 / core_bp as f64
}

//LCB is present if >= LCB_PRESENT_COV of it is covered by alignments >= LCB_PRESENT_IDENTITY
fn lcb_present(lcb: &FastaRecord, hits: &[Hit]) -> bool {
    let intervals: Vec<(usize, usize)> = hits.iter()
        .filter(|h| h.identity >= LCB_PRESENT_IDENTITY)
        .map(|h| (h.lcb_start, h.lcb_end))
        .collect();
    100.0 * covered_bp(&intervals) as f64 / lcb.seq.len() as f64 >= LCB_PRESENT_COV
}

//percent of core LCBs present
pub fn lcb_completeness(core: &[FastaRecord], hits: &[Vec<Hit>]) -> f64 {
    let present = core.iter().zip(hits).filter(|(lcb, h)| lcb_present(lcb, h)).count();
    100.0 * present as f64 / core.len() as f64
}

//Does the assembly belong to the core's species? Per LCB the best alignment (largest aligned bp x identity) is taken;
//the LCB is found if it covers >= SPECIES_MATCH_MIN_COVERAGE of the LCB. Matched when >= SPECIES_MATCH_FRACTION of
//found LCBs pass their bar (match_bars: calibrated per LCB, at most SPECIES_MATCH_IDENTITY; >= SPECIES_MATCH_MIN_LCBS found).
pub struct SpeciesMatch {
    pub found: usize,
    pub percent: f64, //% of found LCBs that pass their bar
}

impl SpeciesMatch {
    pub fn matched(&self) -> bool {
        self.found >= SPECIES_MATCH_MIN_LCBS && self.percent >= SPECIES_MATCH_FRACTION
    }
}

//per LCB: identity of its best alignment (largest aligned bp x identity) when that alignment covers
//>= SPECIES_MATCH_MIN_COVERAGE of the LCB (the LCB is found), else None
pub fn best_lcb_identities(core: &[FastaRecord], hits: &[Vec<Hit>]) -> Vec<Option<f64>> {
    core.iter().zip(hits).map(|(lcb, h)| {
        h.iter().max_by(|a, b| {
            let sa = (a.lcb_end - a.lcb_start + 1) as f64 * a.identity;
            let sb = (b.lcb_end - b.lcb_start + 1) as f64 * b.identity;
            sa.total_cmp(&sb)
        }).filter(|b| 100.0 * (b.lcb_end - b.lcb_start + 1) as f64 / lcb.seq.len() as f64 >= SPECIES_MATCH_MIN_COVERAGE)
          .map(|b| b.identity)
    }).collect()
}

//per-LCB species-match bars, in core order: min(SPECIES_MATCH_IDENTITY, own-genome minimum - SPECIES_MATCH_MARGIN) from
//the calibration.tsv next to the core; SPECIES_MATCH_IDENTITY for LCBs (or cores) without calibration
//calibration.tsv next to the core: lcb -> (identity_min, identity_median) in the species' own build genomes
fn read_calibration(core_path: &Path) -> std::collections::HashMap<String, (f64, f64)> {
    let cal = core_path.parent().map(|d| d.join(CALIBRATION_LCBS)).unwrap_or_default();
    let mut out = std::collections::HashMap::new();
    if let Ok(text) = std::fs::read_to_string(&cal) {
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        let col = |n: &str| header.iter().position(|h| *h == n);
        if let (Some(cl), Some(cm), Some(cd)) = (col("lcb"), col("identity_min"), col("identity_median")) {
            for l in lines {
                let f: Vec<&str> = l.split('\t').collect();
                if let (Some(lcb), Some(Ok(mn)), Some(Ok(md))) = (f.get(cl), f.get(cm).map(|x| x.parse::<f64>()),
                                                                 f.get(cd).map(|x| x.parse::<f64>())) {
                    out.insert(lcb.to_string(), (mn, md));
                }
            }
        }
    }
    out
}

//per LCB, in core order: species-match bar = min(SPECIES_MATCH_IDENTITY, own-genome minimum - SPECIES_MATCH_MARGIN)
//(SPECIES_MATCH_IDENTITY without calibration), and the LCB's median identity in the species' own genomes (None without)
pub fn calibration(core_path: &Path, core: &[FastaRecord]) -> (Vec<f64>, Vec<Option<f64>>) {
    let cal = read_calibration(core_path);
    if cal.is_empty() {
        debug!("No calibration for {}: species-match bar {}% for every LCB", core_path.display(), SPECIES_MATCH_IDENTITY);
    }
    core.iter().map(|l| match cal.get(l.id()) {
        Some((mn, md)) => (SPECIES_MATCH_IDENTITY.min(mn - SPECIES_MATCH_MARGIN), Some(*md)),
        None => (SPECIES_MATCH_IDENTITY, None),
    }).unzip()
}

pub fn match_bars(core_path: &Path, core: &[FastaRecord]) -> Vec<f64> {
    calibration(core_path, core).0
}

pub fn species_match(core: &[FastaRecord], hits: &[Vec<Hit>], bars: &[f64]) -> SpeciesMatch {
    let ids = best_lcb_identities(core, hits);
    let found: Vec<(f64, f64)> = ids.iter().zip(bars).filter_map(|(id, bar)| id.map(|x| (x, *bar))).collect();
    let pass = found.iter().filter(|(id, bar)| id >= bar).count();
    SpeciesMatch { found: found.len(), percent: if found.is_empty() { 0.0 } else { 100.0 * pass as f64 / found.len() as f64 } }
}

pub fn create_csv(path: &Path) -> BufWriter<File> {
    let file = File::create(path).unwrap_or_else(|e| {
        error!("Unable to create {}: {}", path.display(), e);
        std::process::exit(1);
    });
    BufWriter::new(file)
}

pub fn estimate(args: EstimateArgs) {
    let (genomes, target) = check_args(&args);
    let out = PathBuf::from(&args.output);
    let (core_path, name, source) = match target {
        Target::Db(db) => {
            crate::contam::estimate_db(&args, &genomes, &db);
            return;
        }
        Target::Species(path, sp) => (path, sp, "given"),
        Target::Core(path) => {
            //name the species after the core's build folder (core.fasta in <name>/), else the file name
            let p = Path::new(&path);
            let name = if p.file_stem().is_some_and(|s| s == "core") {
                p.parent().and_then(|d| d.file_name()).map(|s| s.to_string_lossy().into_owned())
            } else {
                p.file_stem().map(|s| s.to_string_lossy().into_owned())
            }.unwrap_or_else(|| "core".into());
            (path, name, "core")
        }
    };
    let core = std::sync::Arc::new(read_fasta(&core_path));
    let bars = match_bars(&core_path, &core);
    info!("kore v{} estimating completeness for {} genome(s) against the {} core ({} LCBs, {} bp)", KORE_VERSION, genomes.len(),
          name, core.len(), core.iter().map(|lcb| lcb.seq.len()).sum::<usize>());
    warn!("Only the {} core is used: contamination from other species is not assessed (contamination = same-species \
           redundancy only). Run without --species/--core to estimate contamination with the whole kore db", name);
    let results = crate::contam::run_all(&genomes, args.threads,
        |g, t| Some(crate::contam::analyse_single(&core, &bars, &name, source, g, t)));
    crate::contam::write_report(&results, &out);
}
