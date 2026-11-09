//kore estimate: species detection, completeness and contamination against the kore db (all species cores), or
//completeness and same-species redundancy against one core (--species / --core: analyse_single).
//
//Per assembly, whole db:
// 1. Every contig is mapped to all cores in the kore db (the aggregate) (secondaries kept, so every species' match is seen).
// 2. Species: candidates are the species whose cores the aggregate scan covers most (fraction of core bp
//    at >= SPECIES_IDENTITY). Each candidate's core is then mapped onto the assembly in kore's usual direction; the
//    species is the candidate passing the species check whose found LCBs are most identical (length-weighted identity
//    of each LCB's best hit), independent of how much of the core is present. (The aggregate scan alone
//    is biased: minimap2 drops secondary hits scoring < 20 % of the best, so a species with short LCBs, such as
//    E. coli, loses to a sister species with long LCBs over the same sequence.)
// 3. Contig assignment: each contig goes to the species whose core it matches at the highest bp-weighted
//    identity over >= CONTIG_MIN_ALIGNED aligned bp (the target species competes too). A contig therefore counts
//    for one species only, and a contaminant is not split across its sister species.
// 4. Foreign species S (S != target): breadth = % of S's diagnostic core (S's core minus the target's mask:
//    the parts the target's own genomes carry) covered at >= FOREIGN_IDENTITY by contigs assigned to S.
//    S is called present when breadth >= --call-breadth and >= --call-lcbs LCBs have >= CALL_LCB_BP covered.
//    Its contamination = breadth x genome size(S) / genome size(target).
// 5. Completeness and redundancy against the target core, in kore's usual direction (core -> assembly, second
//    copies kept), ignoring contigs assigned to a called foreign species:
//      completeness = % of core bp covered;  redundancy = % of core bp covered twice or more
//    Redundancy estimates same-species (strain) contamination; the identity of the second copy tells a
//    strain mix (e.g. 98.5 %) from a duplication (~100 %).
// 6. contamination = redundancy + sum of called foreign species.
// 7. Divergent sequence (reported separately, not in the total): contigs >= CONTIG_MIN_ALIGNED with no
//    target-core alignment whose best core match is < DIVERGENT_IDENTITY: most likely a relative with no core
//    in the kore db.
// 8. Species match (estimate::species_match): if the assembly is not the target species (too few of the target's
//    found LCBs pass their calibrated identity bar), every estimate is reported as NA with a warning.
use crate::cli::*;
use crate::consts::*;
use crate::estimate::{Hit, SpeciesMatch, best_lcb_identities, build_aligner, create_csv, genome_name, lcb_completeness, map_genome, map_lcb,
                      nucleotide_completeness, species_match, calibration};
use crate::aggregate::*;
use crate::util::*;

use log::*;

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

pub struct Foreign {
    pub species: String,
    pub breadth: f64,
    pub lcbs: usize,
    pub estimate: f64,
    pub contigs: usize,
    pub contig_bp: usize,
    pub called: bool,
}

pub struct Result {
    pub genome: String,
    pub species: String,
    pub species_source: &'static str,
    pub species_coverage: f64,
    pub species_match: SpeciesMatch,
    pub foreign_assessed: bool, //false for --species / --core: other species are not looked for
    pub completeness: f64,
    pub lcb_completeness: f64,
    pub redundancy: f64,
    pub second_copy_identity: Option<f64>,
    pub foreign: Vec<Foreign>,
    pub divergent_bp: usize,
    pub divergent: f64,
    pub assembly_bp: usize,
    pub core: Arc<Vec<FastaRecord>>,
    pub hits: Vec<Vec<Hit>>,
}

impl Result {
    pub fn foreign_contamination(&self) -> f64 {
        self.foreign.iter().filter(|f| f.called).map(|f| f.estimate).sum::<f64>() + 0.0 //+0.0: no "-0" when empty
    }
    pub fn contamination(&self) -> f64 {
        self.redundancy + self.foreign_contamination()
    }
}

//core bp at depth >= 1 and >= 2 (1-based inclusive hits on LCB coordinates), and the median (bp-weighted) identity
//of the lower copy where depth >= 2
pub fn depth_stats(hits: &[Vec<Hit>]) -> (usize, usize, Option<f64>) {
    let (mut cov1, mut cov2) = (0usize, 0usize);
    let mut second: Vec<(usize, f64)> = Vec::new();
    for lcb_hits in hits {
        let mut ev: Vec<(usize, i32, f64)> = Vec::new(); //(position, +1 start / -1 end, identity), half-open
        for h in lcb_hits {
            ev.push((h.lcb_start - 1, 1, h.identity));
            ev.push((h.lcb_end, -1, h.identity));
        }
        ev.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut active: Vec<f64> = Vec::new();
        let mut last = 0usize;
        for (pos, kind, id) in ev {
            if pos > last && !active.is_empty() {
                cov1 += pos - last;
                if active.len() >= 2 {
                    cov2 += pos - last;
                    let mut a = active.clone();
                    a.sort_by(|x, y| y.total_cmp(x));
                    second.push((pos - last, a[1]));
                }
            }
            if kind == 1 {
                active.push(id);
            } else if let Some(i) = active.iter().position(|&x| x == id) {
                active.remove(i);
            }
            last = pos;
        }
    }
    let median = if second.is_empty() {
        None
    } else {
        second.sort_by(|a, b| a.1.total_cmp(&b.1));
        let half = second.iter().map(|s| s.0).sum::<usize>() as f64 / 2.0;
        let mut acc = 0.0;
        second.iter().find(|s| { acc += s.0 as f64; acc >= half }).map(|s| s.1)
    };
    (cov1, cov2, median)
}

fn subtract(iv: &[(usize, usize)], mask: &[(usize, usize)]) -> usize {
    //bp of merged 1-based inclusive intervals iv not covered by merged mask
    let mut total = 0;
    for &(a, b) in iv {
        let mut covered = 0;
        for &(ma, mb) in mask {
            let (s, e) = (a.max(ma), b.min(mb));
            if s <= e {
                covered += e - s + 1;
            }
        }
        total += (b - a + 1).saturating_sub(covered);
    }
    total
}

struct Db {
    dir: PathBuf,
    species: Vec<DbSpecies>,
    index: HashMap<String, usize>,
    aligner: minimap2::Aligner<minimap2::Built>,
    cores: Mutex<HashMap<usize, Arc<Vec<FastaRecord>>>>,
    cal: Mutex<HashMap<usize, Arc<(Vec<f64>, Vec<Option<f64>>)>>>,
    masks: Mutex<HashMap<usize, Arc<HashMap<String, Vec<(usize, usize)>>>>>,
}

impl Db {
    fn core(&self, sp: usize) -> Arc<Vec<FastaRecord>> {
        let mut c = self.cores.lock().unwrap();
        c.entry(sp).or_insert_with(|| Arc::new(read_fasta(
            &core_path(&self.dir, &self.species[sp].name)))).clone()
    }
    //(species-match bars, own-genome median identity) per LCB
    fn cal(&self, sp: usize) -> Arc<(Vec<f64>, Vec<Option<f64>>)> {
        let core = self.core(sp);
        let mut c = self.cal.lock().unwrap();
        c.entry(sp).or_insert_with(|| Arc::new(calibration(&core_path(&self.dir, &self.species[sp].name), &core))).clone()
    }
    fn mask(&self, sp: usize) -> Arc<HashMap<String, Vec<(usize, usize)>>> {
        let mut m = self.masks.lock().unwrap();
        m.entry(sp).or_insert_with(|| {
            let raw = read_mask(&self.dir, &self.species[sp].name);
            Arc::new(raw.into_iter().map(|(k, v)| (k, merge_intervals(&v))).collect())
        }).clone()
    }
}

//given: the species (--species), the only candidate; otherwise the species is detected
fn analyse(db: &Db, genome: &str, threads: usize, args: &EstimateArgs, given: Option<usize>) -> Option<Result> {
    let contigs = read_fasta(Path::new(genome));
    let assembly_bp: usize = contigs.iter().map(|c| c.seq.len()).sum();
    let hits = map_to_aggregate(&db.aligner, &contigs, &db.index);

    //species coverage at >= SPECIES_IDENTITY
    let mut by_lcb: HashMap<(usize, &str), Vec<(usize, usize)>> = HashMap::new();
    for h in hits.iter().filter(|h| h.identity() >= SPECIES_IDENTITY) {
        by_lcb.entry((h.species, h.lcb.as_str())).or_default().push((h.lcb_start, h.lcb_end));
    }
    let mut covered = vec![0usize; db.species.len()];
    for ((sp, _), iv) in &by_lcb {
        covered[*sp] += covered_bp(iv);
    }
    let coverage = |sp: usize| 100.0 * covered[sp] as f64 / db.species[sp].core_bp.max(1) as f64;
    let mut order: Vec<usize> = (0..db.species.len()).filter(|&s| covered[s] > 0).collect();
    order.sort_by(|&a, &b| coverage(b).total_cmp(&coverage(a)));
    let source = if given.is_some() { "given" } else { "detected" };
    if order.is_empty() && given.is_none() {
        warn!("{}: does not match any species in the kore db (no core aligned at >= {}% identity); estimates reported as {}",
              genome_name(genome), SPECIES_IDENTITY, NO_MATCH);
        return Some(Result {
            genome: genome_name(genome), species: "none".into(), species_source: "detected", species_coverage: 0.0,
            species_match: SpeciesMatch { found: 0, percent: 0.0 }, foreign_assessed: true, completeness: 0.0,
            lcb_completeness: 0.0, redundancy: 0.0, second_copy_identity: None, foreign: Vec::new(), divergent_bp: 0,
            divergent: 0.0, assembly_bp, core: Arc::new(Vec::new()), hits: Vec::new(),
        });
    }
    let candidates: Vec<usize> = match given {
        Some(sp) => vec![sp],
        None => {
            let top = coverage(order[0]);
            order.into_iter().enumerate()
                .take_while(|(i, s)| *i < SPECIES_CANDIDATES || (*i < 2 * SPECIES_CANDIDATES && coverage(*s) >= SPECIES_CANDIDATE_FRACTION * top))
                .map(|(_, s)| s).collect()
        }
    };
    //verify candidates core -> assembly. The candidates competing for the species are those whose alignments fall on the
    //same assembly bases as the best-covered candidate's (SPECIES_SHARED_SEQUENCE: sister species do, a contaminant does
    //not). Among them, the species is the one passing the species check whose found LCBs are most identical: identity of
    //each LCB's best hit, weighted by LCB length (an ANI over the core). Measured on found LCBs only, so incompleteness
    //does not matter. (Identity relative to each species' own median identity was tested and favoured diverse species.)
    let genome_aligner = build_aligner(genome, threads);
    struct Cand { sp: usize, identity: f64, covered: f64, m: SpeciesMatch, hits: Vec<Vec<Hit>>, on: HashMap<String, Vec<(usize, usize)>> }
    let mut cands: Vec<Cand> = Vec::new();
    for &c in &candidates {
        let core_c = db.core(c);
        let cal = db.cal(c);
        let hits_c: Vec<Vec<Hit>> = core_c.iter().map(|l| map_lcb(&genome_aligner, l)).collect();
        let (mut num, mut den) = (0.0, 0.0);
        for (lcb, id) in core_c.iter().zip(best_lcb_identities(&core_c, &hits_c)) {
            if let Some(id) = id {
                num += lcb.seq.len() as f64 * id;
                den += lcb.seq.len() as f64;
            }
        }
        let identity = if den > 0.0 { num / den } else { f64::NEG_INFINITY };
        let mut on: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
        for h in hits_c.iter().flatten() {
            on.entry(h.contig.clone()).or_default().push((h.start, h.end));
        }
        for iv in on.values_mut() {
            *iv = merge_intervals(iv);
        }
        let covered = nucleotide_completeness(&core_c, &hits_c);
        let m = species_match(&core_c, &hits_c, &cal.0);
        debug!("{}: candidate {} scan {:.2}% covered {:.2}% identityity {:.3} over {} LCBs, match {:.1}% ({})", genome_name(genome),
               db.species[c].name, coverage(c), covered, identity, m.found, m.percent, m.matched());
        cands.push(Cand { sp: c, identity, covered, m, hits: hits_c, on });
    }
    let top = cands.iter().enumerate().max_by(|a, b| a.1.covered.total_cmp(&b.1.covered)).map(|(i, _)| i)?;
    let shared = |c: &Cand| -> f64 {
        let own: usize = c.on.values().map(|iv| covered_bp(iv)).sum();
        let both: usize = c.on.iter().map(|(ctg, iv)| cands[top].on.get(ctg)
            .map_or(0, |t| covered_bp(&intersect_intervals(iv, t)))).sum();
        let top_own: usize = cands[top].on.values().map(|iv| covered_bp(iv)).sum();
        both as f64 / own.min(top_own).max(1) as f64 //relative to the smaller footprint: a small core is a subset of a large one
    };
    let keep: Vec<bool> = cands.iter().enumerate().map(|(i, c)| i == top || shared(c) >= SPECIES_SHARED_SEQUENCE).collect();
    let mut cands: Vec<Cand> = cands.into_iter().zip(keep).filter(|(_, k)| *k).map(|(c, _)| c).collect();
    cands.sort_by(|a, b| b.m.matched().cmp(&a.m.matched()).then(b.identity.total_cmp(&a.identity)));
    let best = cands.into_iter().next()?;
    let (target, target_match, all_hits) = (best.sp, best.m, best.hits);
    if !target_match.matched() {
        let what = if given.is_some() { format!("the species {}", db.species[target].name) }
                   else { format!("any species in the kore db (closest: {})", db.species[target].name) };
        warn_no_match(genome, &what, &target_match);
    }

    //contig assignment: highest bp-weighted identity over >= CONTIG_MIN_ALIGNED aligned bp
    let mut agg: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
    for h in &hits {
        let e = agg.entry((h.contig, h.species)).or_insert((0, 0));
        e.0 += h.block;
        e.1 += h.matches;
    }
    let mut winner: HashMap<usize, (usize, f64, usize)> = HashMap::new(); //contig -> (species, identity, bp)
    for (&(c, sp), &(bp, m)) in &agg {
        if bp < CONTIG_MIN_ALIGNED {
            continue;
        }
        let id = 100.0 * m as f64 / bp as f64;
        let better = match winner.get(&c) {
            None => true,
            Some(&(_, wid, wbp)) => id > wid || (id == wid && bp > wbp),
        };
        if better {
            winner.insert(c, (sp, id, bp));
        }
    }

    //foreign species on their diagnostic core
    let mask = db.mask(target);
    let mut masked_bp = vec![0usize; db.species.len()];
    for (lcb, iv) in mask.iter() {
        if let Some(sp) = lcb.split(DB_SEP).next().and_then(|s| db.index.get(s)) {
            masked_bp[*sp] += iv.iter().map(|(a, b)| b - a + 1).sum::<usize>();
        }
    }
    let mut foreign_iv: HashMap<usize, HashMap<&str, Vec<(usize, usize)>>> = HashMap::new();
    for h in hits.iter().filter(|h| h.species != target && h.identity() >= FOREIGN_IDENTITY) {
        if winner.get(&h.contig).is_some_and(|w| w.0 == h.species) {
            foreign_iv.entry(h.species).or_default().entry(h.lcb.as_str()).or_default().push((h.lcb_start, h.lcb_end));
        }
    }
    let g_target = db.species[target].genome_size.max(1.0);
    let mut foreign: Vec<Foreign> = Vec::new();
    for (sp, lcbs) in &foreign_iv {
        let diagnostic = db.species[*sp].core_bp.saturating_sub(masked_bp[*sp]);
        if diagnostic < MIN_DIAGNOSTIC_BP {
            continue; //(almost) all of this core is carried by the target's genomes: cannot be told apart
        }
        let empty = Vec::new();
        let mut cov = 0;
        let mut n_lcbs = 0;
        for (lcb, iv) in lcbs {
            let d = subtract(&merge_intervals(iv), mask.get(*lcb).unwrap_or(&empty));
            cov += d;
            n_lcbs += (d >= CALL_LCB_BP) as usize;
        }
        if cov == 0 {
            continue;
        }
        let breadth = 100.0 * cov as f64 / diagnostic as f64;
        let won: Vec<usize> = winner.iter().filter(|(_, w)| w.0 == *sp).map(|(c, _)| *c).collect();
        foreign.push(Foreign {
            species: db.species[*sp].name.clone(), breadth, lcbs: n_lcbs,
            estimate: breadth * db.species[*sp].genome_size / g_target,
            contigs: won.len(), contig_bp: won.iter().map(|&c| contigs[c].seq.len()).sum(),
            called: breadth >= args.call_breadth && n_lcbs >= args.call_lcbs,
        });
    }
    foreign.sort_by(|a, b| b.breadth.total_cmp(&a.breadth));
    let called: HashSet<&str> = foreign.iter().filter(|f| f.called).map(|f| f.species.as_str()).collect();
    let claimed: HashSet<&str> = winner.iter()
        .filter(|(_, w)| called.contains(db.species[w.0].name.as_str()))
        .map(|(c, _)| contigs[*c].id())
        .collect();

    //completeness and redundancy against the target core, without contigs of called foreign species
    let core = db.core(target);
    let kept: Vec<Vec<Hit>> = all_hits.into_iter()
        .map(|lh| lh.into_iter().filter(|h| !claimed.contains(h.contig.as_str())).collect())
        .collect();
    let core_bp: usize = core.iter().map(|l| l.seq.len()).sum();
    let (_, cov2, second) = depth_stats(&kept);

    //divergent contigs
    let on_target: HashSet<&str> = kept.iter().flatten().map(|h| h.contig.as_str()).collect();
    let divergent_bp: usize = winner.iter()
        .filter(|(c, w)| w.1 < DIVERGENT_IDENTITY && !on_target.contains(contigs[**c].id())
                && !claimed.contains(contigs[**c].id()) && contigs[**c].seq.len() >= CONTIG_MIN_ALIGNED)
        .map(|(c, _)| contigs[*c].seq.len())
        .sum();

    Some(Result {
        genome: genome_name(genome), species: db.species[target].name.clone(), species_source: source,
        species_coverage: coverage(target), species_match: target_match, foreign_assessed: true, completeness: nucleotide_completeness(&core, &kept),
        lcb_completeness: lcb_completeness(&core, &kept), redundancy: 100.0 * cov2 as f64 / core_bp.max(1) as f64,
        second_copy_identity: second, foreign, divergent_bp, divergent: 100.0 * divergent_bp as f64 / g_target,
        assembly_bp, core, hits: kept,
    })
}

fn warn_no_match(genome: &str, what: &str, m: &SpeciesMatch) {
    warn!("{}: does not match {} ({:.1}% of {} found LCBs as identical as in the species' own genomes; {}% needed); \
           estimates reported as {}", genome_name(genome), what, m.percent, m.found, SPECIES_MATCH_FRACTION, NO_MATCH);
}

//--species / --core: one core only. Completeness and same-species redundancy; no species detection and no search for
//other species (foreign contamination and divergent sequence are not assessed).
pub fn analyse_single(core: &Arc<Vec<FastaRecord>>, bars: &[f64], species: &str, source: &'static str, genome: &str, threads: usize) -> Result {
    let assembly_bp: usize = read_fasta(Path::new(genome)).iter().map(|c| c.seq.len()).sum();
    let hits = map_genome(genome, core, threads);
    let m = species_match(core, &hits, bars);
    if !m.matched() {
        warn_no_match(genome, &format!("the species {}", species), &m);
    }
    let core_bp: usize = core.iter().map(|l| l.seq.len()).sum();
    let (cov1, cov2, second) = depth_stats(&hits);
    Result {
        genome: genome_name(genome), species: species.to_string(), species_source: source,
        species_coverage: 100.0 * cov1 as f64 / core_bp.max(1) as f64, species_match: m, foreign_assessed: false,
        completeness: nucleotide_completeness(core, &hits), lcb_completeness: lcb_completeness(core, &hits),
        redundancy: 100.0 * cov2 as f64 / core_bp.max(1) as f64, second_copy_identity: second, foreign: Vec::new(),
        divergent_bp: 0, divergent: 0.0, assembly_bp, core: core.clone(), hits,
    }
}

//runs f on every genome, n at a time with the threads split between them, keeping input order
pub fn run_all<F>(genomes: &[String], threads: usize, f: F) -> Vec<Result>
where F: Fn(&str, usize) -> Option<Result> + Sync {
    let n_parallel = genomes.len().min(threads).max(1);
    let threads_per = (threads / n_parallel).max(1);
    let next = AtomicUsize::new(0);
    let results: Mutex<Vec<Option<Result>>> = Mutex::new((0..genomes.len()).map(|_| None).collect());
    thread::scope(|s| {
        for _ in 0..n_parallel {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= genomes.len() {
                    break;
                }
                if let Some(r) = f(&genomes[i], threads_per) {
                    if r.species_match.matched() {
                        let other = if r.foreign_assessed {
                            format!(", foreign {:.2}%{}, divergent {:.2}%", r.foreign_contamination(),
                                    r.foreign.iter().filter(|x| x.called).map(|x| format!(" {}", x.species)).collect::<String>(), r.divergent)
                        } else { String::new() };
                        info!("{}: {} ({}), {:.2}% complete, {:.2}% contamination (redundancy {:.2}%{})",
                              r.genome, r.species, r.species_source, r.completeness, r.contamination(), r.redundancy, other);
                    }
                    results.lock().unwrap()[i] = Some(r);
                }
            });
        }
    });
    results.into_inner().unwrap().into_iter().flatten().collect()
}

//report.csv (one row per assembly, every estimate), contaminants.csv (whole db only) and alignments.csv
pub fn write_report(results: &[Result], out: &Path) {
    let na = NO_MATCH;
    let mut r_out = create_csv(&out.join(REPORT_CSV));
    writeln!(r_out, "genome,species,species_source,species_match,species_match_lcbs_percent,closest_species,completeness,\
                     lcb_completeness,contamination,redundancy,second_copy_identity,foreign_contamination,foreign_species,\
                     divergent_bp,divergent_percent,species_core_coverage,assembly_bp").unwrap();
    let whole_db = results.iter().any(|r| r.foreign_assessed);
    let mut f = whole_db.then(|| {
        let mut f = create_csv(&out.join(CONTAMINANTS_CSV));
        writeln!(f, "genome,target_species,species,diagnostic_breadth,lcbs,estimate,contigs,contig_bp,called").unwrap();
        f
    });
    let mut a = create_csv(&out.join(ALIGNMENTS_CSV));
    writeln!(a, "genome,lcb,contig,start,end,coverage,identity").unwrap();
    for r in results {
        let pct = r.species_match.percent;
        if !r.species_match.matched() {
            //not the given species / no species in the db: no estimates
            writeln!(r_out, "{},{na},{},false,{:.2},{},{na},{na},{na},{na},,{na},,{na},{na},{:.4},{}", r.genome, r.species_source,
                     pct, r.species, r.species_coverage, r.assembly_bp).unwrap();
            continue;
        }
        let second = r.second_copy_identity.map_or(String::new(), |x| format!("{:.4}", x));
        let (foreign, named, div_bp, div) = if r.foreign_assessed {
            let named: Vec<String> = r.foreign.iter().filter(|x| x.called).map(|x| format!("{}:{:.2}", x.species, x.estimate)).collect();
            (format!("{:.4}", r.foreign_contamination()), named.join(";"), r.divergent_bp.to_string(), format!("{:.4}", r.divergent))
        } else {
            (na.into(), na.into(), na.into(), na.into())
        };
        writeln!(r_out, "{},{},{},true,{:.2},{},{:.4},{:.4},{:.4},{:.4},{},{},{},{},{},{:.4},{}", r.genome, r.species, r.species_source,
                 pct, r.species, r.completeness, r.lcb_completeness, r.contamination(), r.redundancy, second, foreign, named,
                 div_bp, div, r.species_coverage, r.assembly_bp).unwrap();
        if let Some(f) = f.as_mut() {
            for x in &r.foreign {
                writeln!(f, "{},{},{},{:.4},{},{:.4},{},{},{}", r.genome, r.species, x.species, x.breadth, x.lcbs, x.estimate,
                         x.contigs, x.contig_bp, x.called).unwrap();
            }
        }
        for (lcb, lh) in r.core.iter().zip(&r.hits) {
            if lh.is_empty() {
                writeln!(a, "{},{},,,,0.0000,0.0000", r.genome, lcb.id()).unwrap();
            }
            for h in lh {
                let cov = 100.0 * (h.lcb_end - h.lcb_start + 1) as f64 / lcb.seq.len() as f64;
                writeln!(a, "{},{},{},{},{},{:.4},{:.4}", r.genome, lcb.id(), h.contig, h.start, h.end, cov, h.identity).unwrap();
            }
        }
    }
    r_out.flush().unwrap();
    a.flush().unwrap();
    if let Some(mut f) = f {
        f.flush().unwrap();
    }
    let files = if whole_db { format!("{}, {} and {}", REPORT_CSV, CONTAMINANTS_CSV, ALIGNMENTS_CSV) }
                else { format!("{} and {}", REPORT_CSV, ALIGNMENTS_CSV) };
    info!("Wrote {} to {}", files, out.display());
}

//whole db: species detection (unless given) and contamination against every core
pub fn estimate_db(args: &EstimateArgs, genomes: &[String], db: &Path, given: Option<&str>) {
    let species = read_species(db);
    let index: HashMap<String, usize> = species.iter().enumerate().map(|(i, s)| (s.name.clone(), i)).collect();
    info!("kore v{} estimating completeness and contamination for {} genome(s) with the kore db ({} species)",
          KORE_VERSION, genomes.len(), species.len());
    let aligner = aggregate_aligner(&index_path(db), args.threads);
    info!("Loaded {}", index_path(db).display());
    let given = given.map(|sp| *index.get(sp).unwrap_or_else(|| {
        error!("Species {} is not in {}", sp, db.join(DB_SPECIES).display());
        std::process::exit(1);
    }));
    if let Some(sp) = given {
        info!("Species given: {} (no species detection)", species[sp].name);
    }
    let kdb = Db { dir: db.to_path_buf(), species, index, aligner, cores: Mutex::new(HashMap::new()), cal: Mutex::new(HashMap::new()), masks: Mutex::new(HashMap::new()) };
    let results = run_all(genomes, args.threads, |g, t| analyse(&kdb, g, t, args, given));
    write_report(&results, Path::new(&args.output));
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(s: usize, e: usize, id: f64) -> Hit {
        Hit { contig: String::new(), start: 0, end: 0, lcb_start: s, lcb_end: e, identity: id }
    }

    #[test]
    fn depth_counts_second_copies() {
        //LCB covered 1-100 once, 51-100 twice (second copy at 98.5 %)
        let (c1, c2, sec) = depth_stats(&[vec![hit(1, 100, 100.0), hit(51, 100, 98.5)]]);
        assert_eq!((c1, c2), (100, 50));
        assert_eq!(sec, Some(98.5));
    }

    fn lcb(len: usize) -> FastaRecord {
        FastaRecord { header: "x".into(), seq: vec![b'A'; len] }
    }

    #[test]
    fn species_match_rule() {
        //10 LCBs of 1 kb: 6 at 99 % and 4 at 90 % -> 60 % same-species -> match; 5/10 -> no match
        let core: Vec<FastaRecord> = (0..10).map(|_| lcb(1000)).collect();
        let at = |n_same: usize| -> Vec<Vec<Hit>> {
            (0..10).map(|i| vec![hit(1, 1000, if i < n_same { 99.0 } else { 90.0 })]).collect()
        };
        let fixed = vec![SPECIES_MATCH_IDENTITY; 10];
        assert!(species_match(&core, &at(6), &fixed).matched());
        assert!(!species_match(&core, &at(5), &fixed).matched());
        //calibrated: LCBs whose own genomes reach only 90 % get a lower bar, so 90 % hits pass there
        let mut cal = fixed.clone();
        for b in cal.iter_mut().skip(5) { *b = 90.0 - SPECIES_MATCH_MARGIN; }
        assert!(species_match(&core, &at(5), &cal).matched());
        //a short 100 % hit on a conserved stretch (< 50 % of the LCB) does not make an LCB found
        let short: Vec<Vec<Hit>> = (0..10).map(|_| vec![hit(1, 300, 100.0)]).collect();
        assert_eq!(species_match(&core, &short, &fixed).found, 0);
        //fewer than SPECIES_MATCH_MIN_LCBS found LCBs: no match even if all are identical
        let few: Vec<Vec<Hit>> = (0..10).map(|i| if i < 5 { vec![hit(1, 1000, 100.0)] } else { vec![] }).collect();
        assert!(!species_match(&core, &few, &fixed).matched());
    }

    #[test]
    fn subtract_mask() {
        assert_eq!(subtract(&[(1, 100)], &[(11, 20), (91, 200)]), 80);
        assert_eq!(subtract(&[(1, 100)], &[]), 100);
    }
}
