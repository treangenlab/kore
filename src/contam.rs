//kore estimate: species detection, completeness and contamination against the kore db (all species cores), or
//completeness and same-species redundancy against one core (--species / --core: analyse_single).
//
//Per assembly, whole db:
// 1. Every contig is mapped to all cores in the kore db (the aggregate) (secondaries kept, so every species' match is seen).
// 2. Species: candidates are the species whose cores the aggregate scan covers most (fraction of core bp at >=
//    SPECIES_IDENTITY). Each candidate's core is then mapped onto the assembly in kore's usual direction; among candidates
//    passing the species check, the species is the one with the highest match percent (share of found LCBs passing their
//    calibrated identity bars), then the most identical found LCBs. Without a confident match (>= SPECIES_CONFIDENT_MATCH),
//    up to SPECIES_GENUS_CANDIDATES more species of the candidates' genera are verified too: the aggregate scan alone is
//    biased (minimap2 drops secondary hits scoring < 20 % of the best, so a species with short LCBs, such as E. coli, can
//    lose to a sister species with long LCBs over the same sequence).
// 3. Contig assignment, target first: a contig whose alignments to the target core mostly meet the target's calibrated
//    per-LCB bars belongs to the target; any other contig goes to the species whose core it matches at the highest
//    bp-weighted identity over >= CONTIG_MIN_ALIGNED aligned bp (the target species competes too). A contig therefore counts
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
// 7. LCBs found (kore_lcbs_found = found/core LCBs: best alignment, >= 80 % identity, covers >= 50 % of the LCB).
//    Core divergence: per found target LCB, 100 - identity of its best alignment; median and sd over found LCBs (how far
//    the assembly is from the species' core). Length: assembly length relative to the species' build genomes (kore_length_
//    ratio = assembly / median; kore_length_flag 'long' above the longest build genome + LENGTH_FLAG_SLACK x median): sequence
//    the contamination estimate cannot see (species outside the db, non-bacterial DNA) shows up here. Core fraction: core bp
//    found / assembly bp, against the species' core bp / median genome size, with the
//    range over its build genomes: 'median [core/longest-core/shortest]' (kore_core_fraction_expected); kore_core_
//    fraction_flag 'low' below core / (longest + LENGTH_FLAG_SLACK x median), NA below
//    CORE_FRACTION_FLAG_MIN_COMPLETENESS. kore_divergence_flag 'high' when the median divergence exceeds the species'
//    own build genomes' (estimate::divergence_reference + DIVERGENCE_FLAG_MARGIN; DIVERGENCE_FLAG_FALLBACK without
//    calibration): a clade the build genomes did not sample, completeness may read low. Percentages to 2 decimals, fractions to 4; NA for every
//    empty or not-applicable value.
// 8. Species match (estimate::species_match): if the assembly is not the target species (too few of the target's
//    found LCBs pass their calibrated identity bar), every estimate is reported as NA with a warning.
use crate::cli::*;
use crate::consts::*;
use crate::estimate::{Hit, SpeciesMatch, best_lcb_identities, build_aligner, create_csv, genome_name, map_genome, map_lcb,
                      nucleotide_completeness, species_match, calibration, divergence_reference};
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
    pub redundancy: f64,
    pub second_copy_identity: Option<f64>,
    pub foreign: Vec<Foreign>,
    pub lcbs_found: usize, //core LCBs whose best alignment (>= MIN_ALN_IDENTITY) covers >= SPECIES_MATCH_MIN_COVERAGE of the LCB
    pub divergence: Option<(f64, f64)>, //(median, sd) over found core LCBs of 100 - identity of the LCB's best alignment
    pub assembly_bp: usize,
    pub genome_size: f64, //median genome size of the species' build genomes (NaN when unknown)
    pub genome_size_min: f64, //shortest build genome of the species (NaN when unknown)
    pub divergence_ref: Option<f64>, //species divergence reference (estimate::divergence_reference), None without calibration
    pub genome_size_max: f64, //longest build genome of the species (NaN when unknown)
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
//how far the assembly is from the core: per found LCB (best alignment covering >= SPECIES_MATCH_MIN_COVERAGE), 100 - identity of
//its best alignment; (median, population sd) over found LCBs. Small for a typical member of the species; a higher median for a
//divergent lineage, a higher sd for a mosaic (recombinant or mixed) assembly.
fn found_count(core: &[FastaRecord], hits: &[Vec<Hit>]) -> usize {
    best_lcb_identities(core, hits).iter().filter(|x| x.is_some()).count()
}

fn core_divergence(core: &[FastaRecord], hits: &[Vec<Hit>]) -> Option<(f64, f64)> {
    let mut d: Vec<f64> = best_lcb_identities(core, hits).into_iter().flatten().map(|id| 100.0 - id).collect();
    if d.is_empty() {
        return None;
    }
    let mean = d.iter().sum::<f64>() / d.len() as f64;
    let sd = (d.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / d.len() as f64).sqrt();
    Some((percentile(&mut d, 50.0), sd))
}

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
    divref: Mutex<HashMap<usize, Option<f64>>>,
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
    fn div_ref(&self, sp: usize) -> Option<f64> {
        let mut c = self.divref.lock().unwrap();
        *c.entry(sp).or_insert_with(|| divergence_reference(&core_path(&self.dir, &self.species[sp].name)))
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
            redundancy: 0.0, second_copy_identity: None, foreign: Vec::new(), lcbs_found: 0, divergence: None,
            assembly_bp, genome_size: f64::NAN, genome_size_min: f64::NAN, divergence_ref: None, genome_size_max: f64::NAN, core: Arc::new(Vec::new()), hits: Vec::new(),
        });
    }
    //(candidates, sister species of the best-covered one verified only without a confident match among the candidates)
    let (candidates, extra): (Vec<usize>, Vec<usize>) = match given {
        Some(sp) => (vec![sp], Vec::new()),
        None => {
            let top = coverage(order[0]);
            let c: Vec<usize> = order.iter().copied().enumerate()
                .take_while(|(i, s)| *i < SPECIES_CANDIDATES || (*i < 2 * SPECIES_CANDIDATES && coverage(*s) >= SPECIES_CANDIDATE_FRACTION * top))
                .map(|(_, s)| s).collect();
            //same genus as any candidate, most scan coverage first, then those the scan missed (a GTDB split can put the
            //best-covered sister in another genus, e.g. ECMA0423 sp047199055 for E. coli)
            let genera: HashSet<&str> = c.iter().map(|&s| db.species[s].genus.as_str()).filter(|g| !g.is_empty()).collect();
            let in_scan: HashSet<usize> = order.iter().copied().collect();
            let extra: Vec<usize> = order.iter().copied().chain((0..db.species.len()).filter(|s| !in_scan.contains(s)))
                .filter(|&s| genera.contains(db.species[s].genus.as_str()) && !c.contains(&s))
                .take(SPECIES_GENUS_CANDIDATES).collect();
            (c, extra)
        }
    };
    //verify candidates core -> assembly. The candidates competing for the species are those whose alignments fall on the
    //same assembly bases as the best-covered candidate's (SPECIES_SHARED_SEQUENCE: sister species do, a contaminant does
    //not). Among them, the species is the one passing the species check with the highest match percent (share of found
    //LCBs passing their calibrated bars), then the most identical found LCBs: identity of each LCB's best hit, weighted by
    //LCB length (an ANI over the core). Measured on found LCBs only, so incompleteness does not matter. Identity alone
    //favoured sisters whose core is the conserved backbone of the genus (e.g. S. enterica over S. diarizonae); identity
    //relative to each species' own median identity favoured diverse species.
    let genome_aligner = build_aligner(genome, threads);
    struct Cand { sp: usize, identity: f64, covered: f64, m: SpeciesMatch, hits: Vec<Vec<Hit>>, on: HashMap<String, Vec<(usize, usize)>> }
    let verify = |c: usize| -> Cand {
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
        Cand { sp: c, identity, covered, m, hits: hits_c, on }
    };
    let mut cands: Vec<Cand> = candidates.iter().map(|&c| verify(c)).collect();
    let confident = cands.iter().any(|c| c.m.matched() && c.m.percent >= SPECIES_CONFIDENT_MATCH);
    if !confident && !extra.is_empty() {
        debug!("{}: no candidate matched at >= {}%, verifying {} more of the genus", genome_name(genome), SPECIES_CONFIDENT_MATCH, extra.len());
        cands.extend(extra.iter().map(|&c| verify(c)));
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
    //species passing the check: the one whose LCBs pass their calibrated bars most often, then the most identical. Identity
    //alone (over each one's own found LCBs) favours a sister whose core is only the conserved backbone of the genus.
    cands.sort_by(|a, b| b.m.matched().cmp(&a.m.matched()).then(b.m.percent.total_cmp(&a.m.percent))
        .then(b.identity.total_cmp(&a.identity)));
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

    //target first: a contig stays with the target if >= half of its target-aligned bp is in alignments meeting the
    //target's calibrated per-LCB species-match bars, i.e. it lies within the target species' own variation; only
    //contigs the target cannot explain are left to the foreign species they match best
    {
        let core_t = db.core(target);
        let cal_t = db.cal(target);
        let idx: HashMap<&str, usize> = core_t.iter().enumerate().map(|(i, r)| (r.id(), i)).collect();
        let mut pass: HashMap<usize, (usize, usize)> = HashMap::new(); //contig -> (bp meeting bars, bp aligned to target)
        for h in hits.iter().filter(|h| h.species == target) {
            let name = h.lcb.split(DB_SEP).nth(1).unwrap_or("");
            if let Some(&i) = idx.get(name) {
                let e = pass.entry(h.contig).or_insert((0, 0));
                e.1 += h.block;
                if h.identity() >= cal_t.0[i] {
                    e.0 += h.block;
                }
            }
        }
        for (&c, w) in winner.iter_mut() {
            if w.0 != target {
                if let Some(&(p, t)) = pass.get(&c) {
                    if t >= CONTIG_MIN_ALIGNED && 2 * p >= t {
                        if let Some(&(bp, m)) = agg.get(&(c, target)) {
                            *w = (target, 100.0 * m as f64 / bp as f64, bp);
                        }
                    }
                }
            }
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

    Some(Result {
        genome: genome_name(genome), species: db.species[target].name.clone(), species_source: source,
        species_coverage: coverage(target), species_match: target_match, foreign_assessed: true, completeness: nucleotide_completeness(&core, &kept),
        redundancy: 100.0 * cov2 as f64 / core_bp.max(1) as f64,
        second_copy_identity: second, foreign, lcbs_found: found_count(&core, &kept), divergence: core_divergence(&core, &kept), assembly_bp,
        genome_size: db.species[target].genome_size, divergence_ref: db.div_ref(target), genome_size_min: db.species[target].genome_size_min, genome_size_max: db.species[target].genome_size_max, core, hits: kept,
    })
}

fn warn_no_match(genome: &str, what: &str, m: &SpeciesMatch) {
    warn!("{}: does not match {} ({:.1}% of {} found LCBs as identical as in the species' own genomes; {}% needed); \
           estimates reported as {}", genome_name(genome), what, m.percent, m.found,
           if m.found < SPECIES_MATCH_FEW_LCBS { SPECIES_MATCH_FEW_FRACTION } else { SPECIES_MATCH_FRACTION }, NO_MATCH);
}

//--species / --core: one core only. Completeness and same-species redundancy; no species detection and no search for
//other species (foreign contamination is not assessed). sizes: the species' (median, shortest, longest) genome length, if known.
pub fn analyse_single(core: &Arc<Vec<FastaRecord>>, bars: &[f64], species: &str, source: &'static str, sizes: (f64, f64, f64), div_ref: Option<f64>, genome: &str,
                      threads: usize) -> Result {
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
        completeness: nucleotide_completeness(core, &hits),
        redundancy: 100.0 * cov2 as f64 / core_bp.max(1) as f64, second_copy_identity: second, foreign: Vec::new(),
        lcbs_found: found_count(core, &hits), divergence: core_divergence(core, &hits), assembly_bp, genome_size: sizes.0, divergence_ref: div_ref, genome_size_min: sizes.1, genome_size_max: sizes.2, core: core.clone(), hits,
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
                            format!(", foreign {:.2}%{}", r.foreign_contamination(),
                                    r.foreign.iter().filter(|x| x.called).map(|x| format!(" {}", x.species)).collect::<String>())
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

//report.csv (one row per assembly), contaminants.csv (whole db only) and, with alignments, alignments.csv
pub fn write_report(results: &[Result], out: &Path, alignments: bool) {
    let na = NO_MATCH;
    let f4 = |x: f64| format!("{:.4}", x); //fractions and ratios
    let f2 = |x: f64| format!("{:.2}", x); //percentages
    let mut r_out = create_csv(&out.join(REPORT_CSV));
    writeln!(r_out, "genome,species,species_source,kore_species_match_percent,closest_species,kore_completeness,kore_lcbs_found,\
                     kore_contamination,kore_contamination_redundant,kore_contamination_redundant_identity,kore_contamination_foreign,\
                     kore_contamination_foreign_species,kore_divergence_median,kore_divergence_sd,kore_divergence_flag,assembly_bp,kore_length_ratio,\
                     kore_length_flag,core_bp,kore_core_fraction,kore_core_fraction_expected,kore_core_fraction_flag").unwrap();
    let whole_db = results.iter().any(|r| r.foreign_assessed);
    let mut f = whole_db.then(|| {
        let mut f = create_csv(&out.join(CONTAMINANTS_CSV));
        writeln!(f, "genome,target_species,species,diagnostic_breadth,lcbs,estimate,contigs,contig_bp,called").unwrap();
        f
    });
    let mut a = alignments.then(|| {
        let mut a = create_csv(&out.join(ALIGNMENTS_CSV));
        writeln!(a, "genome,lcb,contig,start,end,coverage,identity").unwrap();
        a
    });
    for r in results {
        let pct = r.species_match.percent;
        if !r.species_match.matched() {
            //not the given species / no species in the db: no estimates; closest_species = the best candidate
            writeln!(r_out, "{},{na},{},{:.2},{},{na},{na},{na},{na},{na},{na},{na},{na},{na},{na},{},{na},{na},{na},{na},{na},{na}", r.genome,
                     r.species_source, pct, r.species, r.assembly_bp).unwrap();
            continue;
        }
        let red_id = r.second_copy_identity.filter(|_| r.redundancy >= REDUNDANT_IDENTITY_MIN).map_or(na.to_string(), f2);
        let (foreign, named) = if r.foreign_assessed {
            let named: Vec<String> = r.foreign.iter().filter(|x| x.called).map(|x| format!("{}:{:.2}", x.species, x.estimate)).collect();
            (f2(r.foreign_contamination()), if named.is_empty() { na.to_string() } else { named.join(";") })
        } else {
            (na.into(), na.into())
        };
        let (div_med, div_sd) = r.divergence.map_or((na.to_string(), na.to_string()), |(m, s)| (f2(m), f2(s)));
        //'high': far from the species core (a clade the build genomes did not sample), completeness may read low
        let div_bar = r.divergence_ref.map_or(DIVERGENCE_FLAG_FALLBACK, |d| d + DIVERGENCE_FLAG_MARGIN);
        let div_flag = r.divergence.map_or(na, |(m, _)| if m > div_bar { "high" } else { "ok" });
        let known = |x: f64| x.is_finite() && x > 0.0;
        let core_bp: usize = r.core.iter().map(|l| l.seq.len()).sum();
        //share of the assembly that is the species' core (core bp found / assembly bp), and the same for the species'
        //typical genome (core bp / median genome size): extra sequence (contamination the estimate cannot see, plasmids)
        //lowers the first below the second; lost non-core sequence raises it
        let core_found_bp = r.completeness / 100.0 * core_bp as f64;
        let core_fraction = if r.assembly_bp > 0 { f4(core_found_bp / r.assembly_bp as f64) } else { na.to_string() };
        //expected value 'median [min-max]': core bp over the median, longest and shortest build genome
        let core_fraction_expected = if !known(r.genome_size) { na.to_string() }
            else if known(r.genome_size_min) && known(r.genome_size_max) {
                format!("{} [{}-{}]", f4(core_bp as f64 / r.genome_size), f4(core_bp as f64 / r.genome_size_max), f4(core_bp as f64 / r.genome_size_min))
            } else { f4(core_bp as f64 / r.genome_size) };
        //'low' below the shortest expected fraction, with the length flag's slack: core bp / (longest build genome +
        //LENGTH_FLAG_SLACK x median); a complete assembly is 'low' about when it is 'long', but extra sequence still shows
        //in an incomplete one
        let core_fraction_flag = if known(r.genome_size) && known(r.genome_size_max) && r.assembly_bp > 0
                                    && r.completeness >= CORE_FRACTION_FLAG_MIN_COMPLETENESS {
            let lowest = core_bp as f64 / (r.genome_size_max + LENGTH_FLAG_SLACK * r.genome_size);
            if core_found_bp / (r.assembly_bp as f64) < lowest { "low" } else { "ok" }
        } else { na };
        let length_ratio = if known(r.genome_size) { f4(r.assembly_bp as f64 / r.genome_size) } else { na.to_string() };
        let length_flag = if known(r.genome_size) && known(r.genome_size_max) {
            if r.assembly_bp as f64 > r.genome_size_max + LENGTH_FLAG_SLACK * r.genome_size { "long" } else { "ok" }
        } else { na };
        writeln!(r_out, "{},{},{},{:.2},{na},{:.2},{}/{},{:.2},{:.2},{},{},{},{},{},{},{},{},{},{},{},{},{}", r.genome, r.species,
                 r.species_source, pct, r.completeness, r.lcbs_found, r.core.len(), r.contamination(), r.redundancy, red_id, foreign, named, div_med, div_sd, div_flag,
                 r.assembly_bp, length_ratio, length_flag, core_bp, core_fraction, core_fraction_expected, core_fraction_flag).unwrap();
        if let Some(f) = f.as_mut() {
            for x in &r.foreign {
                writeln!(f, "{},{},{},{:.4},{},{:.4},{},{},{}", r.genome, r.species, x.species, x.breadth, x.lcbs, x.estimate,
                         x.contigs, x.contig_bp, x.called).unwrap();
            }
        }
        if let Some(a) = a.as_mut() {
            for (lcb, lh) in r.core.iter().zip(&r.hits) {
                for h in lh {
                    let cov = 100.0 * (h.lcb_end - h.lcb_start + 1) as f64 / lcb.seq.len() as f64;
                    writeln!(a, "{},{},{},{},{},{:.4},{:.4}", r.genome, lcb.id(), h.contig, h.start, h.end, cov, h.identity).unwrap();
                }
            }
        }
    }
    r_out.flush().unwrap();
    if let Some(mut a) = a {
        a.flush().unwrap();
    }
    if let Some(mut f) = f {
        f.flush().unwrap();
    }
    let mut files = vec![REPORT_CSV];
    if whole_db { files.push(CONTAMINANTS_CSV); }
    if alignments { files.push(ALIGNMENTS_CSV); }
    info!("Wrote {} to {}", files.join(", "), out.display());
}

//whole db: species detection (unless given) and contamination against every core.
//given: one species per genome, in genome order (--species or --species-table)
pub fn estimate_db(args: &EstimateArgs, genomes: &[String], db: &Path, given: Option<&[String]>) {
    let species = read_species(db);
    let index: HashMap<String, usize> = species.iter().enumerate().map(|(i, s)| (s.name.clone(), i)).collect();
    info!("kore v{} estimating completeness and contamination for {} genome(s) with the kore db ({} species)",
          KORE_VERSION, genomes.len(), species.len());
    let aligner = aggregate_aligner(&index_path(db), args.threads);
    info!("Loaded {}", index_path(db).display());
    //genome path -> given species index (a path listed twice must name the same species)
    let mut given_of: HashMap<&str, usize> = HashMap::new();
    if let Some(given) = given {
        for (g, sp) in genomes.iter().zip(given) {
            let i = *index.get(sp).unwrap_or_else(|| {
                error!("Species {} is not in {}", sp, db.join(DB_SPECIES).display());
                std::process::exit(1);
            });
            if given_of.insert(g.as_str(), i).is_some_and(|prev| prev != i) {
                error!("{} is listed with more than one species", g);
                std::process::exit(1);
            }
        }
        let mut names: Vec<&str> = given_of.values().map(|&i| species[i].name.as_str()).collect();
        names.sort();
        names.dedup();
        if names.len() == 1 {
            info!("Species given: {} (no species detection)", names[0]);
        } else {
            info!("Species given for every genome: {} species (no species detection)", names.len());
        }
    }
    let kdb = Db { dir: db.to_path_buf(), species, index, aligner, cores: Mutex::new(HashMap::new()), cal: Mutex::new(HashMap::new()), divref: Mutex::new(HashMap::new()), masks: Mutex::new(HashMap::new()) };
    let results = run_all(genomes, args.threads, |g, t| analyse(&kdb, g, t, args, given_of.get(g).copied()));
    write_report(&results, Path::new(&args.output), args.alignments);
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
        //60 LCBs of 1 kb (>= SPECIES_MATCH_FEW_LCBS found): 60 % at 99 %, the rest at 90 % -> match; 50 % -> no match
        let n = 60;
        let core: Vec<FastaRecord> = (0..n).map(|_| lcb(1000)).collect();
        let at = |n_same: usize| -> Vec<Vec<Hit>> {
            (0..n).map(|i| vec![hit(1, 1000, if i < n_same { 99.0 } else { 90.0 })]).collect()
        };
        let fixed = vec![SPECIES_MATCH_IDENTITY; n];
        assert!(species_match(&core, &at(36), &fixed).matched());
        assert!(!species_match(&core, &at(30), &fixed).matched());
        //calibrated: LCBs whose own genomes reach only 90 % get a lower bar, so 90 % hits pass there
        let mut cal = fixed.clone();
        for b in cal.iter_mut().skip(30) { *b = 90.0 - SPECIES_MATCH_MARGIN; }
        assert!(species_match(&core, &at(30), &cal).matched());
        //few found LCBs (< SPECIES_MATCH_FEW_LCBS): SPECIES_MATCH_FEW_FRACTION of them must pass (another species aligns
        //only its conserved LCBs): 10 found, 8 pass -> no match; 9 pass -> match
        let few_at = |n_found: usize, n_same: usize| -> Vec<Vec<Hit>> {
            (0..n).map(|i| if i < n_found { vec![hit(1, 1000, if i < n_same { 99.0 } else { 90.0 })] } else { vec![] }).collect()
        };
        assert!(!species_match(&core, &few_at(10, 8), &fixed).matched());
        assert!(species_match(&core, &few_at(10, 9), &fixed).matched());
        //a short 100 % hit on a conserved stretch (< 50 % of the LCB) does not make an LCB found
        let short: Vec<Vec<Hit>> = (0..n).map(|_| vec![hit(1, 300, 100.0)]).collect();
        assert_eq!(species_match(&core, &short, &fixed).found, 0);
        //fewer than SPECIES_MATCH_MIN_LCBS found LCBs: no match even if all are identical
        assert!(!species_match(&core, &few_at(5, 5), &fixed).matched());
    }

    #[test]
    fn subtract_mask() {
        assert_eq!(subtract(&[(1, 100)], &[(11, 20), (91, 200)]), 80);
        assert_eq!(subtract(&[(1, 100)], &[]), 100);
    }
}
