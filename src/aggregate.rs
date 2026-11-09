//Aggregate step of kore build: every species core in <db>/cores/ combined into the kore db used by kore estimate.
//
//  <db>/aggregate.fasta   every core's LCBs, renamed <species>|<lcb>
//  <db>/aggregate.mmi     minimap2 (asm20) index of aggregate.fasta (species detection, contamination)
//  <db>/species.tsv       one row per species: GTDB lineage, mask group (family), core size, genome size, masks
//  <db>/masks/<sp>.tsv    (lcb, start, end): the parts of its relatives' cores that species <sp>'s own build genomes
//                         carry at >= MASK_IDENTITY. Relatives = other species of the same GTDB family. Mapped in
//                         kore's usual direction (each relative LCB onto each build genome, secondaries kept), so
//                         relatives never compete for the same genome region and short LCBs are not dropped.
//  <db>/masks/<sp>.relatives.tsv  which relatives (and core checksums) the mask covers: rerunning the aggregate
//                         only maps relatives that are new or whose core changed, so adding a species is cheap.
//
//Why the mask: a relative's core can contain sequence most genomes of a species carry outside that species' own
//core (a conserved gene, or an accessory gene that is core in the relative). Without it a clean genome of the
//species covers part of its relatives' cores and they are wrongly called as contaminants.
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

pub struct DbSpecies {
    pub name: String,
    pub group: String,
    pub lcbs: usize,
    pub core_bp: usize,
    pub genome_size: f64, //median length of the species' build genomes
}

//core of a species in a db (older dbs: cores/<sp>.fasta)
pub fn core_path(db: &Path, species: &str) -> PathBuf {
    let p = db.join(DB_CORES_DIR).join(species).join(CORE_FASTA);
    if p.exists() { p } else { db.join(DB_CORES_DIR).join(format!("{}.fasta", species)) }
}

//aggregate index of a db (older dbs: panel.mmi)
pub fn index_path(db: &Path) -> PathBuf {
    let p = db.join(AGGREGATE_INDEX);
    if p.exists() { p } else { db.join("panel.mmi") }
}

//minimap2 -c -x asm20 --secondary=yes -p AGGREGATE_SECONDARY_RATIO -N AGGREGATE_MAX_SECONDARY -m 30 -s 40
pub fn aggregate_aligner(index: &Path, threads: usize) -> Aligner<Built> {
    let mut aligner = Aligner::builder().asm20().with_cigar().with_index_threads(threads);
    aligner.mapopt.pri_ratio = AGGREGATE_SECONDARY_RATIO;
    aligner.mapopt.best_n = AGGREGATE_MAX_SECONDARY;
    aligner.mapopt.min_chain_score = MM_MIN_CHAIN_SCORE;
    aligner.mapopt.min_dp_max = MM_MIN_DP_SCORE;
    aligner.with_index(index, None).unwrap_or_else(|e| {
        error!("Unable to load minimap2 index {}: {}", index.display(), e);
        std::process::exit(1);
    })
}

pub struct AggregateHit {
    pub contig: usize, //index of the contig in the assembly
    pub species: usize, //index in the db species list
    pub lcb: String, //full aggregate name, <species>|<lcb>
    pub lcb_start: usize, //1-based, inclusive, on the LCB
    pub lcb_end: usize,
    pub block: usize,
    pub matches: usize,
}

impl AggregateHit {
    pub fn identity(&self) -> f64 {
        100.0 * self.matches as f64 / self.block as f64
    }
}

//alignments of every contig to the aggregate (blocks >= MIN_ALN_LEN, identity >= MIN_ALN_IDENTITY)
pub fn map_to_aggregate(aligner: &Aligner<Built>, contigs: &[FastaRecord], species_index: &HashMap<String, usize>) -> Vec<AggregateHit> {
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
            let Some(sp) = lcb.split(DB_SEP).next().and_then(|s| species_index.get(s)) else { continue };
            hits.push(AggregateHit {
                contig: ci, species: *sp, lcb, lcb_start: m.target_start as usize + 1, lcb_end: m.target_end as usize,
                block: m.block_len as usize, matches: m.match_len as usize,
            });
        }
    }
    hits
}

//<db>/species.tsv, by column name
pub fn read_species(db: &Path) -> Vec<DbSpecies> {
    let path = db.join(DB_SPECIES);
    let file = File::open(&path).unwrap_or_else(|e| {
        error!("Unable to open {} (not a kore db?): {}", path.display(), e);
        std::process::exit(1);
    });
    let mut lines = BufReader::new(file).lines().map_while(Result::ok);
    let header: Vec<String> = lines.next().unwrap_or_default().split('\t').map(String::from).collect();
    let col = |n: &str| header.iter().position(|h| h == n).unwrap_or_else(|| {
        error!("{} has no '{}' column", path.display(), n);
        std::process::exit(1);
    });
    let (cs, cg, cl, cb, cz) = (col("species"), col("group"), col("lcbs"), col("core_bp"), col("genome_size"));
    lines.filter(|l| !l.is_empty()).map(|l| {
        let f: Vec<&str> = l.split('\t').collect();
        DbSpecies { name: f[cs].to_string(), group: f[cg].to_string(), lcbs: f[cl].parse().unwrap(),
                    core_bp: f[cb].parse().unwrap(), genome_size: f[cz].parse().unwrap() }
    }).collect()
}

//masks/<species>.tsv -> merged intervals per relative LCB
pub fn read_mask(db: &Path, species: &str) -> HashMap<String, Vec<(usize, usize)>> {
    let path = db.join(DB_MASKS_DIR).join(format!("{}.tsv", species));
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

//key -> value from <species>/info.tsv
pub fn read_info(dir: &Path) -> HashMap<String, String> {
    fs::read_to_string(dir.join(SPECIES_INFO)).unwrap_or_default().lines()
        .filter_map(|l| l.split_once('\t').map(|(k, v)| (k.to_string(), v.to_string()))).collect()
}

//FNV-1a of a file: detects a changed core or genome list between aggregate runs
pub fn checksum(path: &Path) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for b in fs::read(path).unwrap_or_default() {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{:016x}", h)
}

//genomes sampled into the core, cluster-diverse first (round-robin over the parsnp replicates of sampled_genomes.tsv);
//genomes.txt when there is no sampled_genomes.tsv
pub fn sampled_order(dir: &Path) -> Vec<String> {
    let mut reps: Vec<(String, Vec<String>)> = Vec::new();
    for l in fs::read_to_string(dir.join(SAMPLED_GENOMES)).unwrap_or_default().lines().skip(1) {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 2 { continue; }
        match reps.iter_mut().find(|(r, _)| r == f[0]) {
            Some((_, g)) => g.push(f[1].to_string()),
            None => reps.push((f[0].to_string(), vec![f[1].to_string()])),
        }
    }
    if reps.is_empty() {
        return read_file_list(dir.join(SPECIES_GENOMES).to_str().unwrap());
    }
    let mut out: Vec<String> = Vec::new();
    for k in 0..reps.iter().map(|(_, g)| g.len()).max().unwrap_or(0) {
        for (_, g) in &reps {
            if let Some(x) = g.get(k) && !out.contains(x) {
                out.push(x.clone());
            }
        }
    }
    out
}

//cores built by older kore or copied in: add what the aggregate and estimate need but is missing from <species>/
//(GTDB lineage in info.tsv, genome size, calibration), several species at a time
fn fill_missing(sps: &mut [Sp], threads: usize) {
    let need_lineage = sps.iter().any(|s| !s.info.contains_key("domain") && !s.info.contains_key("unlisted"));
    let tax = need_lineage.then(crate::taxonomy::Taxonomy::load);
    let need_calib: Vec<usize> = (0..sps.len()).filter(|&i| !sps[i].dir.join(CALIBRATION_LCBS).exists()).collect();
    let need_size: Vec<usize> = (0..sps.len()).filter(|&i| !sps[i].info.contains_key("genome_size")).collect();
    if tax.is_none() && need_calib.is_empty() && need_size.is_empty() {
        return;
    }
    info!("Completing species folders: {} without lineage, {} without genome size, {} without calibration",
          sps.iter().filter(|s| !s.info.contains_key("domain")).count(), need_size.len(), need_calib.len());
    for s in sps.iter_mut() {
        if let Some(t) = &tax && !s.info.contains_key("domain") {
            let mut add = String::new();
            match t.resolve(&s.name) {
                Ok(l) => for (k, v) in [("gtdb_species", &l.name), ("domain", &l.domain), ("phylum", &l.phylum), ("class", &l.class),
                                        ("order", &l.order), ("family", &l.family), ("genus", &l.genus)] {
                    add.push_str(&format!("{}\t{}\n", k, v));
                    s.info.insert(k.to_string(), v.clone());
                },
                Err(e) => {
                    warn!("{}: {}; no relatives (no masks)", s.name, e);
                    add.push_str("unlisted\ttrue\n");
                    s.info.insert("unlisted".into(), "true".into());
                }
            }
            let mut info = fs::read_to_string(s.dir.join(SPECIES_INFO)).unwrap_or_default();
            if !info.contains("species\t") { info.push_str(&format!("species\t{}\n", s.name)); }
            fs::write(s.dir.join(SPECIES_INFO), info + &add).unwrap();
        }
    }
    let work: Vec<usize> = (0..sps.len()).filter(|i| need_calib.contains(i) || need_size.contains(i)).collect();
    let per = CALIBRATION_N.min(threads).max(1);
    let at_once = (threads / per).max(1).min(work.len().max(1));
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    let sizes: Mutex<Vec<(usize, f64)>> = Mutex::new(Vec::new());
    let sps_ref: &[Sp] = sps;
    thread::scope(|sc| {
        for _ in 0..at_once {
            sc.spawn(|| loop {
                let k = next.fetch_add(1, Ordering::SeqCst);
                if k >= work.len() { break; }
                let i = work[k];
                let s = &sps_ref[i];
                let sampled = sampled_order(&s.dir);
                if need_size.contains(&i) {
                    let g: Vec<String> = sampled.iter().take(DEFAULT_MASK_GENOMES).cloned().collect();
                    sizes.lock().unwrap().push((i, crate::build::genome_size(&g)));
                }
                if need_calib.contains(&i) {
                    let g: Vec<String> = sampled.iter().take(CALIBRATION_N).cloned().collect();
                    crate::build::calibrate(&s.core, &g, &s.dir, per);
                }
                let d = done.fetch_add(1, Ordering::SeqCst) + 1;
                if d % 20 == 0 || d == work.len() {
                    info!("Completing species folders: {}/{}", d, work.len());
                }
            });
        }
    });
    for (i, size) in sizes.into_inner().unwrap() {
        let mut info = fs::read_to_string(sps[i].dir.join(SPECIES_INFO)).unwrap_or_default();
        info.push_str(&format!("genome_size\t{:.0}\n", size));
        fs::write(sps[i].dir.join(SPECIES_INFO), info).unwrap();
        sps[i].info.insert("genome_size".into(), format!("{:.0}", size));
    }
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

struct Sp {
    name: String,
    dir: PathBuf,
    info: HashMap<String, String>,
    genomes: Vec<String>,
    core: Vec<FastaRecord>,
    core_sum: String,
}

pub fn aggregate(db: &Path, threads: usize, n_mask_genomes: usize, mask_identity: f64, seed: u64) {
    let cores_dir = db.join(DB_CORES_DIR);
    let mut names: Vec<String> = fs::read_dir(&cores_dir).map(|rd| rd.filter_map(|e| e.ok())
        .filter(|e| e.path().join(CORE_FASTA).exists()).map(|e| e.file_name().to_string_lossy().into_owned()).collect())
        .unwrap_or_default();
    names.sort();
    if names.is_empty() {
        error!("No species cores in {} to aggregate", cores_dir.display());
        std::process::exit(1);
    }
    let mut sps: Vec<Sp> = names.iter().map(|n| {
        let dir = cores_dir.join(n);
        Sp { name: n.clone(), info: read_info(&dir), genomes: read_file_list(dir.join(SPECIES_GENOMES).to_str().unwrap()),
             core: read_fasta(&dir.join(CORE_FASTA)), core_sum: checksum(&dir.join(CORE_FASTA)), dir }
    }).collect();
    fill_missing(&mut sps, threads);
    let group = |s: &Sp| s.info.get("family").filter(|f| !f.is_empty()).cloned()
        .unwrap_or_else(|| format!("{}:{}", UNLISTED_GROUP, s.name));
    info!("Aggregating {} species cores into {}", sps.len(), db.display());
    fs::create_dir_all(db.join(DB_MASKS_DIR)).unwrap();

    //aggregate.fasta and its index
    let mut fasta = BufWriter::new(File::create(db.join(AGGREGATE_FASTA)).unwrap());
    for s in &sps {
        for lcb in &s.core {
            writeln!(fasta, ">{}{}{}", s.name, DB_SEP, lcb.id()).unwrap();
            for chunk in lcb.seq.chunks(80) {
                fasta.write_all(chunk).unwrap();
                fasta.write_all(b"\n").unwrap();
            }
        }
    }
    fasta.flush().unwrap();
    drop(fasta);
    Aligner::builder().asm20().with_index_threads(threads)
        .with_index(db.join(AGGREGATE_FASTA), Some(db.join(AGGREGATE_INDEX).to_str().unwrap()))
        .unwrap_or_else(|e| { error!("Unable to build the aggregate index: {}", e); std::process::exit(1); });
    info!("Wrote {} and {} ({} LCBs, {} bp)", AGGREGATE_FASTA, AGGREGATE_INDEX, sps.iter().map(|s| s.core.len()).sum::<usize>(),
          sps.iter().map(|s| s.core.iter().map(|l| l.seq.len()).sum::<usize>()).sum::<usize>());

    //masks: which relatives each species needs mapped (new, or core changed); the rest is kept from the last run
    let relatives: Vec<Vec<usize>> = (0..sps.len()).map(|i| (0..sps.len())
        .filter(|&o| o != i && group(&sps[o]) == group(&sps[i])).collect()).collect();
    let mut kept: Vec<HashMap<String, Vec<(usize, usize)>>> = Vec::new();
    let mut todo: Vec<Vec<usize>> = Vec::new();
    for (i, s) in sps.iter().enumerate() {
        let rel_file = db.join(DB_MASKS_DIR).join(format!("{}{}", s.name, MASK_RELATIVES_SUFFIX));
        let previous: HashMap<String, String> = fs::read_to_string(&rel_file).unwrap_or_default().lines()
            .filter_map(|l| l.split_once('\t').map(|(a, b)| (a.to_string(), b.to_string()))).collect();
        let own = format!("{}:{}:{}", checksum(&s.dir.join(SPECIES_GENOMES)), n_mask_genomes, mask_identity);
        let same_self = previous.get("self") == Some(&own);
        //relatives whose mask rows are still valid: same own genomes and settings, same relative core
        let valid: HashSet<&str> = relatives[i].iter()
            .filter(|&&r| same_self && previous.get(&sps[r].name) == Some(&sps[r].core_sum))
            .map(|&r| sps[r].name.as_str()).collect();
        let old = if same_self { read_mask(db, &s.name) } else { HashMap::new() };
        kept.push(old.into_iter().filter(|(lcb, _)| lcb.split(DB_SEP).next().is_some_and(|r| valid.contains(r))).collect());
        todo.push(relatives[i].iter().copied().filter(|&r| !valid.contains(sps[r].name.as_str())).collect());
    }
    let jobs: Vec<(usize, String)> = (0..sps.len()).filter(|&i| !todo[i].is_empty())
        .flat_map(|i| mask_genomes(&sps[i].genomes, n_mask_genomes, seed).into_iter().map(move |g| (i, g))).collect();
    info!("Masks: {} species with relatives, {} need (re)mapping: {} genome x relative-set jobs",
          relatives.iter().filter(|r| !r.is_empty()).count(), todo.iter().filter(|t| !t.is_empty()).count(), jobs.len());

    let masks: Mutex<Vec<HashMap<String, Vec<(usize, usize)>>>> = Mutex::new(kept);
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    thread::scope(|s| {
        for _ in 0..threads.max(1).min(jobs.len().max(1)) {
            s.spawn(|| loop {
                let j = next.fetch_add(1, Ordering::SeqCst);
                if j >= jobs.len() {
                    break;
                }
                let (i, genome) = &jobs[j];
                let genome_aligner = crate::estimate::build_aligner(genome, 1);
                let mut local: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
                for &rel in &todo[*i] {
                    for lcb in &sps[rel].core {
                        for h in crate::estimate::map_lcb(&genome_aligner, lcb) {
                            if h.identity >= mask_identity {
                                local.entry(format!("{}{}{}", sps[rel].name, DB_SEP, lcb.id())).or_default().push((h.lcb_start, h.lcb_end));
                            }
                        }
                    }
                }
                let mut m = masks.lock().unwrap();
                for (lcb, iv) in local {
                    m[*i].entry(lcb).or_default().extend(iv);
                }
                let d = done.fetch_add(1, Ordering::SeqCst) + 1;
                if d % 200 == 0 || d == jobs.len() {
                    info!("Masks: {}/{} genomes", d, jobs.len());
                }
            });
        }
    });

    //masks/, species.tsv
    let masks = masks.into_inner().unwrap();
    let mut table = BufWriter::new(File::create(db.join(DB_SPECIES)).unwrap());
    let lineage = ["gtdb_species", "domain", "phylum", "class", "order", "family", "genus"];
    let ckeys: Vec<String> = REPORT_CLUSTER_ANI.iter().map(|&a| crate::build::clusters_key(a)).collect();
    writeln!(table, "species\t{}\tgroup\tlcbs\tcore_bp\tgenome_size\tbuild_genomes\t{}\tmask_genomes\trelatives\tmasked_bp", lineage.join("\t"),
             ckeys.join("\t")).unwrap();
    for (i, s) in sps.iter().enumerate() {
        let mut w = BufWriter::new(File::create(db.join(DB_MASKS_DIR).join(format!("{}.tsv", s.name))).unwrap());
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
        let mut rel = format!("self\t{}:{}:{}\n", checksum(&s.dir.join(SPECIES_GENOMES)), n_mask_genomes, mask_identity);
        for &r in &relatives[i] {
            rel.push_str(&format!("{}\t{}\n", sps[r].name, sps[r].core_sum));
        }
        fs::write(db.join(DB_MASKS_DIR).join(format!("{}{}", s.name, MASK_RELATIVES_SUFFIX)), rel).unwrap();
        let lin: Vec<String> = lineage.iter().map(|k| s.info.get(*k).cloned().unwrap_or_default()).collect();
        //cluster counts: NA for cores built without clustering or by older kore
        let clusters: Vec<&str> = ckeys.iter().map(|k| s.info.get(k).map_or("NA", |v| v.as_str())).collect();
        writeln!(table, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", s.name, lin.join("\t"), group(s), s.core.len(),
                 s.core.iter().map(|l| l.seq.len()).sum::<usize>(), s.info.get("genome_size").map_or("0", |v| v),
                 s.genomes.len(), clusters.join("\t"), s.genomes.len().min(n_mask_genomes), relatives[i].len(), masked).unwrap();
    }
    table.flush().unwrap();
    //masks of species no longer in the db
    let present: HashSet<&str> = sps.iter().map(|s| s.name.as_str()).collect();
    for e in fs::read_dir(db.join(DB_MASKS_DIR)).unwrap().filter_map(|e| e.ok()) {
        let f = e.file_name().to_string_lossy().into_owned();
        let sp = f.trim_end_matches(MASK_RELATIVES_SUFFIX).trim_end_matches(".tsv");
        if !present.contains(sp) {
            let _ = fs::remove_file(e.path());
        }
    }
    info!("Wrote {}, {}/ to {}", DB_SPECIES, DB_MASKS_DIR, db.display());
}
