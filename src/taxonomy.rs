//GTDB taxonomy for kore build: downloaded once (bac120 + ar53 taxonomy of GTDB_RELEASE, ~10 MB) into $KORE_TAXONOMY
//or ~/.kore/gtdb_<release>, then used to check species names and to give each species its lineage. The family is
//the mask group: masks are computed between species of the same family.
//
//Species names are accepted as the GTDB name ("Escherichia coli", "Bacillus_A anthracis"), with or without "s__",
//or as kore's folder name (lowercase, spaces -> "_": "escherichia_coli", "bacillus_a_anthracis").
use crate::consts::*;
use crate::util::in_path;

use log::*;

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

use flate2::read::MultiGzDecoder;

#[derive(Clone)]
pub struct Lineage {
    pub name: String, //GTDB species name, without s__
    pub domain: String,
    pub phylum: String,
    pub class: String,
    pub order: String,
    pub family: String,
    pub genus: String,
}

pub struct Taxonomy {
    by_name: HashMap<String, Lineage>,     //GTDB species name -> lineage
    by_slug: HashMap<String, Vec<String>>, //kore folder name -> GTDB names (more than one: ambiguous)
}

//kore's folder name for a species: lowercase, no "s__", spaces -> "_"
pub fn slug(name: &str) -> String {
    name.trim().trim_start_matches("s__").to_lowercase().replace(' ', "_")
}

pub fn taxonomy_dir() -> PathBuf {
    if let Ok(d) = std::env::var(TAXONOMY_ENV) {
        return PathBuf::from(d);
    }
    let home = std::env::var("HOME").unwrap_or_else(|_| ".".into());
    Path::new(&home).join(".kore").join(format!("gtdb_{}", GTDB_RELEASE))
}

//downloads the GTDB taxonomy files that are missing from dir (curl, else wget)
fn ensure_downloaded(dir: &Path) {
    let missing: Vec<&str> = GTDB_FILES.iter().copied().filter(|f| !dir.join(f).exists()).collect();
    if missing.is_empty() {
        return;
    }
    fs::create_dir_all(dir).unwrap_or_else(|e| {
        error!("Unable to create the taxonomy folder {}: {} (set ${} to another folder)", dir.display(), e, TAXONOMY_ENV);
        std::process::exit(1);
    });
    let tool = if in_path("curl") { "curl" } else if in_path("wget") { "wget" } else {
        error!("GTDB taxonomy {} not found in {} and neither curl nor wget is available to download it", GTDB_RELEASE, dir.display());
        std::process::exit(1);
    };
    for f in missing {
        let url = format!("{}/{}", GTDB_URL, f);
        let part = dir.join(format!("{}.part", f));
        info!("Downloading GTDB {} taxonomy: {}", GTDB_RELEASE, url);
        let ok = if tool == "curl" {
            Command::new("curl").args(["-fsSL", "--retry", "3", "-o"]).arg(&part).arg(&url).status()
        } else {
            Command::new("wget").args(["-q", "-O"]).arg(&part).arg(&url).status()
        }.is_ok_and(|s| s.success());
        if !ok || fs::rename(&part, dir.join(f)).is_err() {
            let _ = fs::remove_file(&part);
            error!("Unable to download {} (no internet?). Download it into {} by hand or set ${}", url, dir.display(), TAXONOMY_ENV);
            std::process::exit(1);
        }
    }
}

impl Taxonomy {
    pub fn load() -> Taxonomy {
        let dir = taxonomy_dir();
        ensure_downloaded(&dir);
        let mut by_name: HashMap<String, Lineage> = HashMap::new();
        for f in GTDB_FILES {
            let file = File::open(dir.join(f)).unwrap_or_else(|e| {
                error!("Unable to open {}: {}", dir.join(f).display(), e);
                std::process::exit(1);
            });
            for line in BufReader::new(MultiGzDecoder::new(file)).lines().map_while(Result::ok) {
                let Some((_, tax)) = line.split_once('\t') else { continue };
                let r: Vec<&str> = tax.split(';').map(|x| x.get(3..).unwrap_or("")).collect();
                if r.len() != 7 || r[6].is_empty() || by_name.contains_key(r[6]) {
                    continue;
                }
                by_name.insert(r[6].to_string(), Lineage {
                    name: r[6].into(), domain: r[0].into(), phylum: r[1].into(), class: r[2].into(), order: r[3].into(),
                    family: r[4].into(), genus: r[5].into(),
                });
            }
        }
        let mut by_slug: HashMap<String, Vec<String>> = HashMap::new();
        for n in by_name.keys() {
            by_slug.entry(slug(n)).or_default().push(n.clone());
        }
        info!("GTDB {} taxonomy: {} species ({})", GTDB_RELEASE, by_name.len(), dir.display());
        Taxonomy { by_name, by_slug }
    }

    //the GTDB species a user's name refers to, or an error with suggestions
    pub fn resolve(&self, name: &str) -> Result<&Lineage, String> {
        let exact = name.trim().trim_start_matches("s__").replace('_', " ");
        if let Some(l) = self.by_name.get(name.trim().trim_start_matches("s__")) {
            return Ok(l);
        }
        match self.by_slug.get(&slug(name)).map(|v| v.as_slice()) {
            Some([one]) => return Ok(&self.by_name[one]),
            Some(many) if many.len() > 1 => return Err(format!("'{}' is ambiguous in GTDB {}: give the exact name ({})",
                                                                name, GTDB_RELEASE, many.join(" / "))),
            _ => {}
        }
        Err(format!("'{}' is not a GTDB {} species{}", name, GTDB_RELEASE, self.suggest(&exact)))
    }

    //up to 3 GTDB names closest to an unknown one (same genus first)
    fn suggest(&self, name: &str) -> String {
        let q = slug(name);
        let genus = q.split('_').next().unwrap_or("");
        let mut scored: Vec<(usize, &String)> = self.by_slug.keys()
            .filter(|s| genus.len() >= 3 && s.starts_with(genus))
            .map(|s| (edit_distance(&q, s), s)).collect();
        if scored.is_empty() {
            scored = self.by_slug.keys().map(|s| (edit_distance(&q, s), s)).collect();
        }
        scored.sort();
        let best: Vec<String> = scored.iter().take(3).filter(|(d, _)| *d <= q.len().max(4) / 2)
            .flat_map(|(_, s)| self.by_slug[*s].iter().cloned()).collect();
        if best.is_empty() { String::new() } else { format!(" (did you mean: {}?)", best.join(", ")) }
    }
}

fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b): (Vec<char>, Vec<char>) = (a.chars().collect(), b.chars().collect());
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for i in 1..=a.len() {
        let mut cur = vec![i; b.len() + 1];
        for j in 1..=b.len() {
            cur[j] = (prev[j] + 1).min(cur[j - 1] + 1).min(prev[j - 1] + (a[i - 1] != b[j - 1]) as usize);
        }
        prev = cur;
    }
    prev[b.len()]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slugs() {
        assert_eq!(slug("s__Bacillus_A anthracis"), "bacillus_a_anthracis");
        assert_eq!(slug("Escherichia coli"), "escherichia_coli");
    }

    #[test]
    fn distance() {
        assert_eq!(edit_distance("escherichia_coli", "escherichia_colli"), 1);
        assert_eq!(edit_distance("abc", "abc"), 0);
    }
}
