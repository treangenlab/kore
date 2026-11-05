//kore panel: combine species cores into one database for species detection and contamination estimates.
//
//Output folder:
//  panel.fasta   every core's LCBs, renamed <species>|<lcb>
//  panel.mmi     minimap2 (asm20) index of panel.fasta
//  cores/        <species>.fasta, each species' own core (completeness is estimated against it)
//  species.tsv   species, group, lcbs, core_bp, genome_size, mask_genomes, relatives, masked_bp
//  masks/        <species>.tsv (lcb, start, end): for species S, the parts of its relatives' cores that S's own
//                build genomes carry at >= --mask-identity. Relatives = other species of the same group.
//                Computed in kore's usual direction: each relative LCB is mapped onto each of S's build genomes
//                (secondaries kept), so relatives never compete for the same genome region and short LCBs are
//                not dropped next to a sister species' long ones.
//
//Why the mask: a relative's core can contain sequence most genomes of S carry outside S's own core (a diverse
//species such as E. coli has a small core; its GTDB sister species have large cores that include genes most
//E. coli carry). Without the mask a clean genome of S covers much of its relatives' cores and they look like
//contaminants. The mask is computed once, from the build genomes; only the intervals are stored.
use crate::cli::*;
use crate::consts::*;
use crate::util::*;

use log::*;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;

use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use minimap2::{Aligner, Built};

pub struct PanelSpecies {
    pub name: String,
    pub group: String,
    pub lcbs: usize,
    pub core_bp: usize,
    pub genome_size: f64, //median length of the genomes used for the mask
}

struct InputRow {
    species: String,
    core: String,
    genomes: Vec<String>,
    group: String,
}

fn read_input(path: &str) -> Vec<InputRow> {
    let file = File::open(path).unwrap_or_else(|e| {
        error!("Unable to open {}: {}", path, e);
        std::process::exit(1);
    });
    let mut lines = BufReader::new(file).lines().map_while(Result::ok);
    let header: Vec<String> = lines.next().unwrap_or_default().split('\t').map(|s| s.trim().to_string()).collect();
    let col = |name: &str| header.iter().position(|h| h == name).unwrap_or_else(|| {
        error!("{} has no '{}' column (needs species, core, genomes, group)", path, name);
        std::process::exit(1);
    });
    let (cs, cc, cg, cgr) = (col("species"), col("core"), col("genomes"), col("group"));

    let mut rows = Vec::new();
    let mut seen = HashSet::new();
    for line in lines.filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = line.split('\t').map(|s| s.trim()).collect();
        let species = f[cs].to_string();
        if species.contains(PANEL_SEP) || species.contains(char::is_whitespace) || !seen.insert(species.clone()) {
            error!("Species names must be unique and contain no whitespace or '{}': {}", PANEL_SEP, species);
            std::process::exit(1);
        }
        if !Path::new(f[cc]).exists() {
            error!("{}: core does not exist: {}", species, f[cc]);
            std::process::exit(1);
        }
        rows.push(InputRow { species, core: f[cc].to_string(), genomes: read_file_list(f[cg]), group: f[cgr].to_string() });
    }
    rows
}

//minimap2 -c -x asm20 --secondary=yes -p PANEL_SECONDARY_RATIO -N PANEL_MAX_SECONDARY -m 30 -s 40 against the panel
pub fn panel_aligner(index: &Path, threads: usize) -> Aligner<Built> {
    let mut aligner = Aligner::builder().asm20().with_cigar().with_index_threads(threads);
    aligner.mapopt.pri_ratio = PANEL_SECONDARY_RATIO;
    aligner.mapopt.best_n = PANEL_MAX_SECONDARY;
    aligner.mapopt.min_chain_score = MM_MIN_CHAIN_SCORE;
    aligner.mapopt.min_dp_max = MM_MIN_DP_SCORE;
    aligner.with_index(index, None).unwrap_or_else(|e| {
        error!("Unable to load minimap2 index {}: {}", index.display(), e);
        std::process::exit(1);
    })
}

pub struct PanelHit {
    pub contig: usize, //index of the contig in the assembly
    pub species: usize, //index in the panel species list
    pub lcb: String, //full panel name, <species>|<lcb>
    pub lcb_start: usize, //1-based, inclusive, on the LCB
    pub lcb_end: usize,
    pub block: usize,
    pub matches: usize,
}

impl PanelHit {
    pub fn identity(&self) -> f64 {
        100.0 * self.matches as f64 / self.block as f64
    }
}

//alignments of every contig to the panel (blocks >= MIN_ALN_LEN, identity >= MIN_ALN_IDENTITY)
pub fn map_to_panel(aligner: &Aligner<Built>, contigs: &[FastaRecord], species_index: &HashMap<String, usize>) -> Vec<PanelHit> {
    let mut hits = Vec::new();
    for (ci, contig) in contigs.iter().enumerate() {
        if contig.seq.is_empty() {
            continue;
        }
        let mappings = aligner.map(&contig.seq, false, false, None, None, None).unwrap_or_else(|e| {
            error!("minimap2 failed to map {}: {}", contig.id(), e);
            std::process::exit(1);
        });
        for m in mappings {
            if m.block_len < MIN_ALN_LEN || 100.0 * m.match_len as f64 / (m.block_len as f64) < MIN_ALN_IDENTITY {
                continue;
            }
            let lcb = m.target_name.as_ref().map_or(String::new(), |n| n.to_string());
            let Some(sp) = lcb.split(PANEL_SEP).next().and_then(|s| species_index.get(s)) else { continue };
            hits.push(PanelHit {
                contig: ci, species: *sp, lcb, lcb_start: m.target_start as usize + 1, lcb_end: m.target_end as usize,
                block: m.block_len as usize, matches: m.match_len as usize,
            });
        }
    }
    hits
}

pub fn read_species(dir: &Path) -> Vec<PanelSpecies> {
    let path = dir.join(PANEL_SPECIES);
    let file = File::open(&path).unwrap_or_else(|e| {
        error!("Unable to open {} (not a kore panel?): {}", path.display(), e);
        std::process::exit(1);
    });
    BufReader::new(file).lines().map_while(Result::ok).skip(1).filter(|l| !l.is_empty()).map(|l| {
        let f: Vec<&str> = l.split('\t').collect();
        PanelSpecies { name: f[0].to_string(), group: f[1].to_string(), lcbs: f[2].parse().unwrap(),
                       core_bp: f[3].parse().unwrap(), genome_size: f[4].parse().unwrap() }
    }).collect()
}

//masks/<species>.tsv -> merged intervals per relative LCB
pub fn read_mask(dir: &Path, species: &str) -> HashMap<String, Vec<(usize, usize)>> {
    let path = dir.join(PANEL_MASKS_DIR).join(format!("{}.tsv", species));
    let mut mask: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    if let Ok(file) = File::open(&path) {
        for l in BufReader::new(file).lines().map_while(Result::ok).skip(1) {
            let f: Vec<&str> = l.split('\t').collect();
            if f.len() == 3 {
                mask.entry(f[0].to_string()).or_default().push((f[1].parse().unwrap(), f[2].parse().unwrap()));
            }
        }
    }
    mask
}

fn mask_genomes(genomes: &[String], max: usize, seed: u64) -> Vec<String> {
    let mut g: Vec<String> = genomes.to_vec();
    g.sort();
    if g.len() > max {
        g.shuffle(&mut StdRng::seed_from_u64(seed));
        g.truncate(max);
        g.sort();
    }
    g
}

pub fn panel(args: PanelArgs) {
    let out = PathBuf::from(&args.output);
    init_logging(log_level(args.verbose, args.debug), &out.join(LOG_FILE));
    print_banner();
    let rows = read_input(&args.input);
    info!("kore v{} building a panel of {} species", KORE_VERSION, rows.len());
    fs::create_dir_all(out.join(PANEL_CORES_DIR)).unwrap();
    fs::create_dir_all(out.join(PANEL_MASKS_DIR)).unwrap();

    //panel.fasta and cores/
    let mut species: Vec<PanelSpecies> = Vec::new();
    let mut cores: Vec<Vec<FastaRecord>> = Vec::new();
    let mut panel_out = BufWriter::new(File::create(out.join(PANEL_FASTA)).unwrap());
    for r in &rows {
        let core = read_fasta(Path::new(&r.core));
        let refs: Vec<&FastaRecord> = core.iter().collect();
        write_fasta(&refs, &out.join(PANEL_CORES_DIR).join(format!("{}.fasta", r.species)));
        for lcb in &core {
            writeln!(panel_out, ">{}{}{}", r.species, PANEL_SEP, lcb.id()).unwrap();
            for chunk in lcb.seq.chunks(80) {
                panel_out.write_all(chunk).unwrap();
                panel_out.write_all(b"\n").unwrap();
            }
        }
        species.push(PanelSpecies { name: r.species.clone(), group: r.group.clone(), lcbs: core.len(),
                                    core_bp: core.iter().map(|l| l.seq.len()).sum(), genome_size: 0.0 });
        cores.push(core);
    }
    panel_out.flush().unwrap();
    drop(panel_out);
    info!("Wrote {} ({} LCBs, {} bp)", PANEL_FASTA, species.iter().map(|s| s.lcbs).sum::<usize>(),
          species.iter().map(|s| s.core_bp).sum::<usize>());

    let index = out.join(PANEL_INDEX);
    Aligner::builder().asm20().with_index_threads(args.threads)
        .with_index(out.join(PANEL_FASTA), Some(index.to_str().unwrap()))
        .unwrap_or_else(|e| { error!("Unable to build the panel index: {}", e); std::process::exit(1); });
    info!("Wrote {}", PANEL_INDEX);

    //masks: every relative LCB mapped onto each of the species' build genomes, kept at >= mask identity
    let mut jobs: Vec<(usize, String)> = Vec::new();
    let relatives: Vec<Vec<usize>> = rows.iter().map(|r| (0..rows.len())
        .filter(|&o| rows[o].group == r.group && rows[o].species != r.species).collect()).collect();
    let n_relatives: Vec<usize> = relatives.iter().map(|r| r.len()).collect();
    for (i, r) in rows.iter().enumerate() {
        for g in mask_genomes(&r.genomes, args.mask_genomes, args.seed) {
            jobs.push((i, g));
        }
    }
    info!("Masking: {} genomes from {} species ({} with relatives in the panel)", jobs.len(), rows.len(),
          n_relatives.iter().filter(|&&n| n > 0).count());

    let masks: Mutex<Vec<HashMap<String, Vec<(usize, usize)>>>> = Mutex::new((0..rows.len()).map(|_| HashMap::new()).collect());
    let lengths: Mutex<Vec<Vec<usize>>> = Mutex::new(vec![Vec::new(); rows.len()]);
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    thread::scope(|s| {
        for _ in 0..args.threads.max(1) {
            s.spawn(|| loop {
                let j = next.fetch_add(1, Ordering::SeqCst);
                if j >= jobs.len() {
                    break;
                }
                let (i, genome) = &jobs[j];
                let len: usize = read_fasta(Path::new(genome)).iter().map(|c| c.seq.len()).sum();
                lengths.lock().unwrap()[*i].push(len);
                if n_relatives[*i] > 0 {
                    let genome_aligner = crate::estimate::build_aligner(genome, 1);
                    let mut local: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
                    for &rel in &relatives[*i] {
                        for lcb in &cores[rel] {
                            for h in crate::estimate::map_lcb(&genome_aligner, lcb) {
                                if h.identity >= args.mask_identity {
                                    local.entry(format!("{}{}{}", species[rel].name, PANEL_SEP, lcb.id()))
                                        .or_default().push((h.lcb_start, h.lcb_end));
                                }
                            }
                        }
                    }
                    let mut m = masks.lock().unwrap();
                    for (lcb, iv) in local {
                        m[*i].entry(lcb).or_default().extend(iv);
                    }
                }
                let d = done.fetch_add(1, Ordering::SeqCst) + 1;
                if d % 200 == 0 || d == jobs.len() {
                    info!("Masking: {}/{} genomes", d, jobs.len());
                }
            });
        }
    });

    let masks = masks.into_inner().unwrap();
    let mut lengths = lengths.into_inner().unwrap();
    let mut table = BufWriter::new(File::create(out.join(PANEL_SPECIES)).unwrap());
    writeln!(table, "species\tgroup\tlcbs\tcore_bp\tgenome_size\tmask_genomes\trelatives\tmasked_bp").unwrap();
    for (i, sp) in species.iter_mut().enumerate() {
        let mut l: Vec<f64> = lengths[i].drain(..).map(|x| x as f64).collect();
        let n_genomes = l.len();
        sp.genome_size = percentile(&mut l, 50.0);
        let mut w = BufWriter::new(File::create(out.join(PANEL_MASKS_DIR).join(format!("{}.tsv", sp.name))).unwrap());
        writeln!(w, "lcb\tstart\tend").unwrap();
        let mut lcbs: Vec<&String> = masks[i].keys().collect();
        lcbs.sort();
        let mut masked = 0;
        for lcb in lcbs {
            for (a, b) in merge_intervals(&masks[i][lcb]) {
                writeln!(w, "{}\t{}\t{}", lcb, a, b).unwrap();
                masked += b - a + 1;
            }
        }
        w.flush().unwrap();
        writeln!(table, "{}\t{}\t{}\t{}\t{:.0}\t{}\t{}\t{}", sp.name, sp.group, sp.lcbs, sp.core_bp, sp.genome_size, n_genomes,
                 n_relatives[i], masked).unwrap();
    }
    table.flush().unwrap();
    info!("Wrote {}, {}/ and {}/ to {}", PANEL_SPECIES, PANEL_CORES_DIR, PANEL_MASKS_DIR, out.display());
}
