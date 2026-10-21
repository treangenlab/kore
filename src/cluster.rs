use crate::consts::*;

use log::*;

use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

use std::collections::{HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

//all-vs-all skani ANI of the genomes (sparse output: one row per reported pair)
fn run_skani(genomes: &[String], skani_dir: &Path, threads: usize) -> PathBuf {
    fs::create_dir_all(skani_dir).unwrap_or_else(|e| {
        error!("Unable to create directory {}: {}", skani_dir.display(), e);
        std::process::exit(1);
    });
    let list = skani_dir.join("genomes.txt");
    fs::write(&list, genomes.iter().map(|g| format!("{}\n", g)).collect::<String>()).unwrap_or_else(|e| {
        error!("Unable to write {}: {}", list.display(), e);
        std::process::exit(1);
    });

    let out = skani_dir.join("ani.tsv");
    let output = Command::new(SKANI_BIN)
        .arg("triangle")
        .arg("-l").arg(&list)
        .arg("--sparse")
        .arg("-t").arg(threads.to_string())
        .arg("-o").arg(&out)
        .output()
        .unwrap_or_else(|e| {
            error!("Unable to run {}: {}", SKANI_BIN, e);
            std::process::exit(1);
        });
    if !output.status.success() {
        error!("{} failed:\n{}", SKANI_BIN, String::from_utf8_lossy(&output.stderr));
        std::process::exit(1);
    }
    out
}

fn open_text(path: &Path) -> Box<dyn BufRead> {
    let file = File::open(path).unwrap_or_else(|e| {
        error!("Unable to open {}: {}", path.display(), e);
        std::process::exit(1);
    });
    let reader: Box<dyn Read> = if path.extension().is_some_and(|e| e == "gz") { Box::new(MultiGzDecoder::new(file)) } else { Box::new(file) };
    Box::new(BufReader::new(reader))
}

fn file_name(path: &str) -> &str {
    Path::new(path).file_name().and_then(|f| f.to_str()).unwrap_or(path)
}

//keeps the first SKANI_KEEP_COLUMNS columns of skani's sparse output (paths, ANI, aligned fractions), gzipped
fn save_ani(raw: &Path, saved: &Path, list: &Path, genomes: &[String]) {
    let write = || -> std::io::Result<()> {
        let mut out = GzEncoder::new(BufWriter::new(File::create(saved)?), Compression::default());
        for line in open_text(raw).lines() {
            let line = line?;
            let cols: Vec<&str> = line.split('\t').take(SKANI_KEEP_COLUMNS).collect();
            writeln!(out, "{}", cols.join("\t"))?;
        }
        out.finish()?.flush()?;
        let mut sorted: Vec<&String> = genomes.iter().collect();
        sorted.sort();
        fs::write(list, sorted.iter().map(|g| format!("{}\n", g)).collect::<String>())
    };
    write().unwrap_or_else(|e| {
        error!("Unable to save skani ANI to {}: {}", saved.display(), e);
        std::process::exit(1);
    });
}

//saved ANI covers exactly these genomes
fn same_genomes(list: &Path, genomes: &[String]) -> bool {
    match fs::read_to_string(list) {
        Ok(text) => text.lines().filter(|l| !l.is_empty()).map(String::from).collect::<HashSet<_>>() == genomes.iter().cloned().collect(),
        Err(_) => false,
    }
}

//skani ANI for the genomes: --skani-ani if given, else the ANI saved in the output folder by an
//earlier build of the same genomes, else run skani and save its output there for next time
fn ani_source(genomes: &[String], skani_ani: Option<&Path>, out_dir: &Path, skani_dir: &Path, threads: usize) -> PathBuf {
    if let Some(p) = skani_ani {
        info!("Using skani ANI from {}", p.display());
        return p.to_path_buf();
    }
    let saved = out_dir.join(SKANI_ANI);
    let list = out_dir.join(SKANI_GENOMES);
    if saved.exists() && same_genomes(&list, genomes) {
        info!("Reusing skani ANI from {} (same {} genomes)", saved.display(), genomes.len());
        return saved;
    }
    if saved.exists() {
        info!("Saved skani ANI in {} is for a different genome set, rerunning skani", out_dir.display());
    }
    let raw = run_skani(genomes, skani_dir, threads);
    save_ani(&raw, &saved, &list, genomes);
    info!("Saved skani ANI to {}", saved.display());
    saved
}

//n x n ANI distance (100 - ANI); pairs skani did not report get its reporting floor. Genomes are
//matched to the ANI table by path, or by file name when the paths differ (e.g. a moved genome folder)
fn distance_matrix(ani: &Path, genomes: &[String]) -> Vec<f32> {
    let n = genomes.len();
    let index: HashMap<&str, usize> = genomes.iter().enumerate().map(|(i, g)| (g.as_str(), i)).collect();
    let by_name: HashMap<&str, usize> = genomes.iter().enumerate().map(|(i, g)| (file_name(g), i)).collect();
    let lookup = |p: &str| index.get(p).or_else(|| by_name.get(file_name(p))).copied();
    let mut dist = vec![(100.0 - SKANI_MISSING_ANI) as f32; n * n];
    for i in 0..n {
        dist[i * n + i] = 0.0;
    }

    let mut seen = vec![false; n];
    for line in open_text(ani).lines().map_while(Result::ok).skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        if let (Some(i), Some(j), Ok(a)) = (lookup(f[0]), lookup(f[1]), f[2].parse::<f32>()) {
            dist[i * n + j] = 100.0 - a;
            dist[j * n + i] = 100.0 - a;
            seen[i] = true;
            seen[j] = true;
        }
    }
    let unseen = seen.iter().filter(|&&s| !s).count();
    if n > 1 && unseen > 0 {
        warn!("{} of {} genomes have no pairs in {}; they are treated as {}% ANI to everything", unseen, n, ani.display(), SKANI_MISSING_ANI);
    }
    dist
}

fn find(parent: &mut [usize], x: usize) -> usize {
    let mut root = x;
    while parent[root] != root {
        root = parent[root];
    }
    let mut x = x;
    while parent[x] != root {
        let next = parent[x];
        parent[x] = root;
        x = next;
    }
    root
}

//average linkage (UPGMA) by the nearest-neighbour chain algorithm, O(n^2), then cut so genomes
//joined at <= max_dist share a cluster (same flat clusters as scipy fcluster criterion="distance").
//Also returns the dendrogram's leaf order, so neighbouring leaves are close relatives.
fn upgma_clusters(mut dist: Vec<f32>, n: usize, max_dist: f32) -> (Vec<usize>, Vec<usize>) {
    let mut size = vec![1usize; n];
    let mut active = vec![true; n];
    let mut merges: Vec<(usize, usize, f32)> = Vec::with_capacity(n.saturating_sub(1));
    let mut chain: Vec<usize> = Vec::new();
    let mut remaining = n;

    while remaining > 1 {
        if chain.is_empty() {
            chain.push((0..n).find(|&i| active[i]).unwrap());
        }
        loop {
            let a = *chain.last().unwrap();
            let prev = if chain.len() >= 2 { Some(chain[chain.len() - 2]) } else { None };
            //nearest active neighbour of a; prefer the previous chain element on ties so the chain terminates
            let mut best = usize::MAX;
            let mut best_d = f32::INFINITY;
            for k in 0..n {
                if k == a || !active[k] {
                    continue;
                }
                let d = dist[a * n + k];
                if d < best_d || (d == best_d && Some(k) == prev) {
                    best = k;
                    best_d = d;
                }
            }
            if Some(best) == prev {
                chain.pop();
                chain.pop();
                //merge best into a: Lance-Williams update for average linkage
                let (sa, sb) = (size[a] as f32, size[best] as f32);
                for k in 0..n {
                    if active[k] && k != a && k != best {
                        let d = (sa * dist[a * n + k] + sb * dist[best * n + k]) / (sa + sb);
                        dist[a * n + k] = d;
                        dist[k * n + a] = d;
                    }
                }
                size[a] += size[best];
                active[best] = false;
                merges.push((a, best, best_d));
                remaining -= 1;
                break;
            }
            chain.push(best);
        }
    }

    //leaf order: rebuild the tree from the merges (node n + k is the k-th merge) and walk it depth first
    let mut node_of: Vec<usize> = (0..n).collect();
    let mut children: Vec<(usize, usize)> = Vec::with_capacity(merges.len());
    for (k, &(a, b, _)) in merges.iter().enumerate() {
        children.push((node_of[a], node_of[b]));
        node_of[a] = n + k;
    }
    let mut leaf_order = Vec::with_capacity(n);
    let mut stack = vec![if merges.is_empty() { 0 } else { n + merges.len() - 1 }];
    while let Some(node) = stack.pop() {
        if node < n {
            leaf_order.push(node);
        } else {
            let (l, r) = children[node - n];
            stack.push(r);
            stack.push(l);
        }
    }

    let mut parent: Vec<usize> = (0..n).collect();
    for &(a, b, d) in &merges {
        if d <= max_dist {
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            parent[rb] = ra;
        }
    }
    let mut label: HashMap<usize, usize> = HashMap::new();
    let labels = (0..n).map(|i| {
        let r = find(&mut parent, i);
        let next = label.len();
        *label.entry(r).or_insert(next)
    }).collect();
    (labels, leaf_order)
}

//genomes ordered for sampling: clusters at cluster_ani, then one genome from each cluster per round
//(seeded order of genomes within a cluster) until every cluster is exhausted. Clusters are in seeded
//random order, or with by_tree in dendrogram order (cluster ids then follow the tree, so consecutive ids
//are close relatives). Returns the ordered genomes and each one's cluster.
pub fn cluster_order(genomes: &[String], cluster_ani: f64, skani_ani: Option<&Path>, out_dir: &Path, skani_dir: &Path,
                     threads: usize, seed: u64, by_tree: bool) -> Vec<(String, usize)> {
    let ani = ani_source(genomes, skani_ani, out_dir, skani_dir, threads);
    let n = genomes.len();
    let (labels, leaf_order) = upgma_clusters(distance_matrix(&ani, genomes), n, (100.0 - cluster_ani) as f32 + 1e-6);

    let mut clusters: Vec<Vec<usize>> = vec![Vec::new(); labels.iter().max().map_or(0, |m| m + 1)];
    for (i, &c) in labels.iter().enumerate() {
        clusters[c].push(i);
    }
    let mut rng = StdRng::seed_from_u64(seed);
    for c in clusters.iter_mut() {
        c.shuffle(&mut rng);
    }
    if by_tree {
        //each cluster's position = its first leaf in the dendrogram
        let mut first = vec![usize::MAX; clusters.len()];
        for (pos, &leaf) in leaf_order.iter().enumerate() {
            first[labels[leaf]] = first[labels[leaf]].min(pos);
        }
        let mut idx: Vec<usize> = (0..clusters.len()).collect();
        idx.sort_by_key(|&c| first[c]);
        clusters = idx.into_iter().map(|c| std::mem::take(&mut clusters[c])).collect();
    } else {
        clusters.shuffle(&mut rng);
    }

    let mut sizes: Vec<usize> = clusters.iter().map(|c| c.len()).collect();
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    info!("Clustered {} genomes into {} clusters at {}% ANI (largest {}, {} singletons)",
        n, clusters.len(), cluster_ani, sizes[0], sizes.iter().filter(|&&s| s == 1).count());

    let rounds = sizes[0];
    let mut ordered = Vec::with_capacity(n);
    for r in 0..rounds {
        for (c, members) in clusters.iter().enumerate() {
            if let Some(&i) = members.get(r) {
                ordered.push((genomes[i].clone(), c));
            }
        }
    }
    ordered
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn upgma_cuts_two_groups() {
        //two tight pairs far apart: {0,1} and {2,3}
        let d = vec![
            0.0, 0.1, 2.0, 2.0,
            0.1, 0.0, 2.0, 2.0,
            2.0, 2.0, 0.0, 0.2,
            2.0, 2.0, 0.2, 0.0,
        ];
        let (l, _) = upgma_clusters(d.clone(), 4, 0.5);
        assert_eq!(l[0], l[1]);
        assert_eq!(l[2], l[3]);
        assert_ne!(l[0], l[2]);
        assert!(upgma_clusters(d.clone(), 4, 0.05).0.iter().collect::<std::collections::HashSet<_>>().len() == 4);
        assert!(upgma_clusters(d, 4, 5.0).0.iter().all(|&x| x == 0));
    }

    #[test]
    fn upgma_uses_average_not_single_linkage() {
        //0-1 close, 2 close to 1 but far from 0: average distance of 2 to {0,1} is 1.05 > 0.6 cut
        let d = vec![
            0.0, 0.5, 1.6,
            0.5, 0.0, 0.5,
            1.6, 0.5, 0.0,
        ];
        let (l, _) = upgma_clusters(d, 3, 0.6);
        assert_eq!(l[0], l[1]);
        assert_ne!(l[0], l[2]);
    }

    #[test]
    fn leaf_order_keeps_clades_together() {
        //{0,2} close, {1,3} close, far apart: the two clades must be contiguous in the leaf order
        let d = vec![
            0.0, 2.0, 0.1, 2.0,
            2.0, 0.0, 2.0, 0.2,
            0.1, 2.0, 0.0, 2.0,
            2.0, 0.2, 2.0, 0.0,
        ];
        let (_, order) = upgma_clusters(d, 4, 0.5);
        let pos = |x: usize| order.iter().position(|&o| o == x).unwrap();
        assert_eq!(order.len(), 4);
        assert_eq!(pos(0).abs_diff(pos(2)), 1);
        assert_eq!(pos(1).abs_diff(pos(3)), 1);
    }
}
