//kore estimate --panel: species detection, completeness and contamination against a kore panel.
//
//Per assembly:
// 1. Every contig is mapped to all cores in the panel (secondaries kept, so every species' match is seen).
// 2. Species (or --species): candidates are the species whose cores the panel scan covers most (fraction of core bp
//    at >= SPECIES_IDENTITY). Each candidate's core is then mapped onto the assembly in kore's usual direction and
//    scored as sum over LCBs of covered bp x identity of the LCB's best hit, over core bp; the best is the species. (The panel scan alone
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
//    in the panel.
use crate::cli::*;
use crate::consts::*;
use crate::estimate::{Hit, build_aligner, create_csv, genome_name, lcb_completeness, map_lcb, nucleotide_completeness};
use crate::panel::*;
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

struct Panel {
    dir: PathBuf,
    species: Vec<PanelSpecies>,
    index: HashMap<String, usize>,
    aligner: minimap2::Aligner<minimap2::Built>,
    cores: Mutex<HashMap<usize, Arc<Vec<FastaRecord>>>>,
    masks: Mutex<HashMap<usize, Arc<HashMap<String, Vec<(usize, usize)>>>>>,
}

impl Panel {
    fn core(&self, sp: usize) -> Arc<Vec<FastaRecord>> {
        let mut c = self.cores.lock().unwrap();
        c.entry(sp).or_insert_with(|| Arc::new(read_fasta(
            &self.dir.join(PANEL_CORES_DIR).join(format!("{}.fasta", self.species[sp].name))))).clone()
    }
    fn mask(&self, sp: usize) -> Arc<HashMap<String, Vec<(usize, usize)>>> {
        let mut m = self.masks.lock().unwrap();
        m.entry(sp).or_insert_with(|| {
            let raw = read_mask(&self.dir, &self.species[sp].name);
            Arc::new(raw.into_iter().map(|(k, v)| (k, merge_intervals(&v))).collect())
        }).clone()
    }
}

fn analyse(panel: &Panel, genome: &str, forced: Option<usize>, threads: usize, args: &EstimateArgs) -> Option<Result> {
    let contigs = read_fasta(Path::new(genome));
    let assembly_bp: usize = contigs.iter().map(|c| c.seq.len()).sum();
    let hits = map_to_panel(&panel.aligner, &contigs, &panel.index);

    //species coverage at >= SPECIES_IDENTITY
    let mut by_lcb: HashMap<(usize, &str), Vec<(usize, usize)>> = HashMap::new();
    for h in hits.iter().filter(|h| h.identity() >= SPECIES_IDENTITY) {
        by_lcb.entry((h.species, h.lcb.as_str())).or_default().push((h.lcb_start, h.lcb_end));
    }
    let mut covered = vec![0usize; panel.species.len()];
    for ((sp, _), iv) in &by_lcb {
        covered[*sp] += covered_bp(iv);
    }
    let coverage = |sp: usize| 100.0 * covered[sp] as f64 / panel.species[sp].core_bp.max(1) as f64;
    let candidates: Vec<usize> = match forced {
        Some(t) => vec![t],
        None => {
            let mut order: Vec<usize> = (0..panel.species.len()).filter(|&s| covered[s] > 0).collect();
            order.sort_by(|&a, &b| coverage(b).total_cmp(&coverage(a)));
            if order.is_empty() {
                warn!("{}: no species in the panel matched; skipped", genome_name(genome));
                return None;
            }
            let top = coverage(order[0]);
            order.into_iter().enumerate()
                .take_while(|(i, s)| *i < SPECIES_CANDIDATES || (*i < 2 * SPECIES_CANDIDATES && coverage(*s) >= SPECIES_CANDIDATE_FRACTION * top))
                .map(|(_, s)| s).collect()
        }
    };
    //verify candidates core -> assembly; score = covered core fraction x bp-weighted identity
    let genome_aligner = build_aligner(genome, threads);
    let mut best: Option<(usize, f64, Vec<Vec<Hit>>)> = None;
    for &c in &candidates {
        let core_c = panel.core(c);
        let hits_c: Vec<Vec<Hit>> = core_c.iter().map(|l| map_lcb(&genome_aligner, l)).collect();
        let cbp: usize = core_c.iter().map(|l| l.seq.len()).sum();
        //per LCB: covered bp x identity of its best hit (second copies from a contaminant must not lower the score)
        let mut matched = 0.0;
        for h in &hits_c {
            if let Some(best_id) = h.iter().map(|x| x.identity).max_by(|a, b| a.total_cmp(b)) {
                matched += covered_bp(&h.iter().map(|x| (x.lcb_start, x.lcb_end)).collect::<Vec<_>>()) as f64 * best_id / 100.0;
            }
        }
        let score = matched / cbp.max(1) as f64;
        debug!("{}: candidate {} scan {:.2}% verified score {:.3}", genome_name(genome), panel.species[c].name, coverage(c), score);
        if best.as_ref().is_none_or(|b| score > b.1) {
            best = Some((c, score, hits_c));
        }
    }
    let (target, _, all_hits) = best?;
    let source = if forced.is_some() { "given" } else { "detected" };

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
    let mask = panel.mask(target);
    let mut masked_bp = vec![0usize; panel.species.len()];
    for (lcb, iv) in mask.iter() {
        if let Some(sp) = lcb.split(PANEL_SEP).next().and_then(|s| panel.index.get(s)) {
            masked_bp[*sp] += iv.iter().map(|(a, b)| b - a + 1).sum::<usize>();
        }
    }
    let mut foreign_iv: HashMap<usize, HashMap<&str, Vec<(usize, usize)>>> = HashMap::new();
    for h in hits.iter().filter(|h| h.species != target && h.identity() >= FOREIGN_IDENTITY) {
        if winner.get(&h.contig).is_some_and(|w| w.0 == h.species) {
            foreign_iv.entry(h.species).or_default().entry(h.lcb.as_str()).or_default().push((h.lcb_start, h.lcb_end));
        }
    }
    let g_target = panel.species[target].genome_size.max(1.0);
    let mut foreign: Vec<Foreign> = Vec::new();
    for (sp, lcbs) in &foreign_iv {
        let diagnostic = panel.species[*sp].core_bp.saturating_sub(masked_bp[*sp]);
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
            species: panel.species[*sp].name.clone(), breadth, lcbs: n_lcbs,
            estimate: breadth * panel.species[*sp].genome_size / g_target,
            contigs: won.len(), contig_bp: won.iter().map(|&c| contigs[c].seq.len()).sum(),
            called: breadth >= args.call_breadth && n_lcbs >= args.call_lcbs,
        });
    }
    foreign.sort_by(|a, b| b.breadth.total_cmp(&a.breadth));
    let called: HashSet<&str> = foreign.iter().filter(|f| f.called).map(|f| f.species.as_str()).collect();
    let claimed: HashSet<&str> = winner.iter()
        .filter(|(_, w)| called.contains(panel.species[w.0].name.as_str()))
        .map(|(c, _)| contigs[*c].id())
        .collect();

    //completeness and redundancy against the target core, without contigs of called foreign species
    let core = panel.core(target);
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
        genome: genome_name(genome), species: panel.species[target].name.clone(), species_source: source,
        species_coverage: coverage(target), completeness: nucleotide_completeness(&core, &kept),
        lcb_completeness: lcb_completeness(&core, &kept), redundancy: 100.0 * cov2 as f64 / core_bp.max(1) as f64,
        second_copy_identity: second, foreign, divergent_bp, divergent: 100.0 * divergent_bp as f64 / g_target,
        assembly_bp, core, hits: kept,
    })
}

fn write_outputs(results: &[Result], out: &Path) {
    let mut c = create_csv(&out.join(COMPLETENESS_CSV));
    writeln!(c, "genome,species,nucleotide_completeness,lcb_completeness").unwrap();
    let mut t = create_csv(&out.join(CONTAMINATION_CSV));
    writeln!(t, "genome,species,species_source,species_core_coverage,contamination,redundancy,second_copy_identity,\
                 foreign_contamination,foreign_species,divergent_bp,divergent_percent,assembly_bp").unwrap();
    let mut f = create_csv(&out.join(CONTAMINANTS_CSV));
    writeln!(f, "genome,target_species,species,diagnostic_breadth,lcbs,estimate,contigs,contig_bp,called").unwrap();
    let mut a = create_csv(&out.join(ALIGNMENTS_CSV));
    writeln!(a, "genome,lcb,contig,start,end,coverage,identity").unwrap();
    for r in results {
        writeln!(c, "{},{},{:.4},{:.4}", r.genome, r.species, r.completeness, r.lcb_completeness).unwrap();
        let named: Vec<String> = r.foreign.iter().filter(|x| x.called).map(|x| format!("{}:{:.2}", x.species, x.estimate)).collect();
        writeln!(t, "{},{},{},{:.4},{:.4},{:.4},{},{:.4},{},{},{:.4},{}", r.genome, r.species, r.species_source, r.species_coverage,
                 r.contamination(), r.redundancy, r.second_copy_identity.map_or(String::new(), |x| format!("{:.4}", x)),
                 r.foreign_contamination(), named.join(";"), r.divergent_bp, r.divergent, r.assembly_bp).unwrap();
        for x in &r.foreign {
            writeln!(f, "{},{},{},{:.4},{},{:.4},{},{},{}", r.genome, r.species, x.species, x.breadth, x.lcbs, x.estimate,
                     x.contigs, x.contig_bp, x.called).unwrap();
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
    for w in [&mut c, &mut t, &mut f, &mut a] {
        w.flush().unwrap();
    }
}

pub fn estimate_panel(args: &EstimateArgs, genomes: &[String]) {
    let dir = PathBuf::from(args.panel.as_ref().unwrap());
    let species = read_species(&dir);
    let index: HashMap<String, usize> = species.iter().enumerate().map(|(i, s)| (s.name.clone(), i)).collect();
    let forced = args.species.as_ref().map(|s| *index.get(s).unwrap_or_else(|| {
        error!("Species {} is not in the panel", s);
        std::process::exit(1);
    }));
    info!("kore v{} estimating completeness and contamination for {} genome(s) against a panel of {} species",
          KORE_VERSION, genomes.len(), species.len());
    let aligner = panel_aligner(&dir.join(PANEL_INDEX), args.threads);
    info!("Loaded {}", dir.join(PANEL_INDEX).display());
    let panel = Panel { dir, species, index, aligner, cores: Mutex::new(HashMap::new()), masks: Mutex::new(HashMap::new()) };

    let n_parallel = genomes.len().min(args.threads).max(1);
    let threads_per = (args.threads / n_parallel).max(1);
    let next = AtomicUsize::new(0);
    let results: Mutex<Vec<Option<Result>>> = Mutex::new((0..genomes.len()).map(|_| None).collect());
    thread::scope(|s| {
        for _ in 0..n_parallel {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= genomes.len() {
                    break;
                }
                if let Some(r) = analyse(&panel, &genomes[i], forced, threads_per, args) {
                    info!("{}: {} ({}), {:.2}% complete, {:.2}% contamination (redundancy {:.2}%, foreign {:.2}%{}), divergent {:.2}%",
                          r.genome, r.species, r.species_source, r.completeness, r.contamination(), r.redundancy,
                          r.foreign_contamination(),
                          r.foreign.iter().filter(|x| x.called).map(|x| format!(" {}", x.species)).collect::<String>(),
                          r.divergent);
                    results.lock().unwrap()[i] = Some(r);
                }
            });
        }
    });
    let results: Vec<Result> = results.into_inner().unwrap().into_iter().flatten().collect();
    let out = PathBuf::from(&args.output);
    write_outputs(&results, &out);
    info!("Wrote {}, {}, {} and {} to {}", COMPLETENESS_CSV, CONTAMINATION_CSV, CONTAMINANTS_CSV, ALIGNMENTS_CSV, out.display());
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

    #[test]
    fn subtract_mask() {
        assert_eq!(subtract(&[(1, 100)], &[(11, 20), (91, 200)]), 80);
        assert_eq!(subtract(&[(1, 100)], &[]), 100);
    }
}
