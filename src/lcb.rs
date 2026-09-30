use log::*;

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

//IUPAC code indexed by a bitmask of the tied bases: A=1, C=2, G=4, T=8
const IUPAC: &[u8; 16] = b"NACMGRSVTWYHKDBN";
const LINE_WRAP: usize = 80;

pub struct ConsensusStats {
    pub written: usize,
    pub total_bp: usize,
    pub too_short: usize,
    pub too_diverse: usize,
}

pub struct Lcb {
    pub start: usize, //coordinates on the reference
    pub end: usize,
    pub rows: Vec<Vec<u8>>, //aligned sequence for each genome
}

fn open_xmfa(xmfa: &Path) -> BufReader<File> {
    let file = File::open(xmfa).unwrap_or_else(|e| {
        error!("Unable to open {}: {}", xmfa.display(), e);
        std::process::exit(1);
    });
    BufReader::new(file)
}

//##SequenceFile entries from the xmfa header, in ##SequenceIndex order
pub fn xmfa_sequence_files(xmfa: &Path) -> Vec<String> {
    open_xmfa(xmfa).lines()
        .map_while(Result::ok)
        .take_while(|l| l.starts_with('#'))
        .filter_map(|l| l.strip_prefix("##SequenceFile ").map(String::from))
        .collect()
}

//parsnp adds a duplicate '<reference>.ref' entry
fn ref_duplicate_index(xmfa: &Path) -> Option<usize> {
    xmfa_sequence_files(xmfa).iter().position(|f| f.ends_with(".ref")).map(|i| i + 1)
}

//"> 1:22681-22997 + cluster27 s1:p22680" -> (1, 22681, 22997)
fn parse_block_header(line: &str) -> Option<(usize, usize, usize)> {
    let (idx, coords) = line[1..].split_whitespace().next()?.split_once(':')?;
    let (start, end) = coords.split_once('-')?;
    Some((idx.parse().ok()?, start.parse().ok()?, end.parse().ok()?))
}

//streams each LCB in the xmfa to f, skipping the parsnp .ref duplicate so the reference only gets one vote
pub fn for_each_lcb(xmfa: &Path, mut f: impl FnMut(Lcb)) {
    let skip_idx = ref_duplicate_index(xmfa);
    let mut rows: Vec<Vec<u8>> = Vec::new();
    let mut coords: Option<(usize, usize)> = None;
    let mut keep_row = false;

    for line in open_xmfa(xmfa).lines().map_while(Result::ok) {
        if line.starts_with('#') {
            continue;
        } else if line.starts_with('>') {
            let Some((idx, start, end)) = parse_block_header(&line) else {
                error!("Unable to parse xmfa header in {}: {}", xmfa.display(), line);
                std::process::exit(1);
            };
            coords.get_or_insert((start, end));
            keep_row = Some(idx) != skip_idx;
            if keep_row {
                rows.push(Vec::new());
            }
        } else if line.starts_with('=') {
            if let Some((start, end)) = coords.take() {
                f(Lcb { start, end, rows: std::mem::take(&mut rows) });
            }
        } else if keep_row {
            rows.last_mut().unwrap().extend(line.trim_end().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
}

//majority vote for one alignment column, None if gaps are the strict majority
fn column_consensus(rows: &[Vec<u8>], col: usize) -> Option<u8> {
    let mut counts = [0usize; 4];
    let mut gaps = 0;
    for row in rows {
        match row[col] {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            b'-' => gaps += 1,
            _ => {} //N/ambiguous: base present but unknown
        }
    }

    if gaps * 2 > rows.len() {
        return None;
    }

    let max = *counts.iter().max().unwrap();
    if max == 0 {
        return Some(b'N');
    }
    let mask = counts.iter().enumerate()
        .filter(|(_, c)| **c == max)
        .fold(0, |m, (i, _)| m | (1 << i));
    Some(IUPAC[mask])
}

pub fn lcb_consensus(lcb: &Lcb) -> Vec<u8> {
    let len = lcb.rows[0].len();
    if lcb.rows.iter().any(|r| r.len() != len) {
        error!("Aligned rows differ in length in LCB at {}-{}", lcb.start, lcb.end);
        std::process::exit(1);
    }
    (0..len).filter_map(|col| column_consensus(&lcb.rows, col)).collect()
}

//fraction of alignment columns where any genome differs (gaps and Ns included)
pub fn lcb_diversity(lcb: &Lcb) -> f64 {
    let len = lcb.rows[0].len();
    let variable = (0..len)
        .filter(|&col| lcb.rows.iter().any(|r| r[col] != lcb.rows[0][col]))
        .count();
    variable as f64 / len as f64
}

//one consensus record per LCB that passes the length and diversity filters
pub fn write_consensus_fasta(xmfa: &Path, out_fasta: &Path, prefix: &str, min_length: usize, max_diversity: f64) -> ConsensusStats {
    let file = File::create(out_fasta).unwrap_or_else(|e| {
        error!("Unable to create {}: {}", out_fasta.display(), e);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(file);

    let mut stats = ConsensusStats { written: 0, total_bp: 0, too_short: 0, too_diverse: 0 };
    let mut n_lcbs = 0;
    for_each_lcb(xmfa, |lcb| {
        n_lcbs += 1;
        let seq = lcb_consensus(&lcb);
        if seq.len() < min_length {
            stats.too_short += 1;
            return;
        }
        let diversity = lcb_diversity(&lcb);
        if diversity > max_diversity {
            stats.too_diverse += 1;
            return;
        }

        writeln!(writer, ">{}_lcb{} ref={}-{} n_genomes={} diversity={:.4}", prefix, n_lcbs, lcb.start, lcb.end, lcb.rows.len(), diversity).unwrap();
        for chunk in seq.chunks(LINE_WRAP) {
            writer.write_all(chunk).unwrap();
            writer.write_all(b"\n").unwrap();
        }
        stats.written += 1;
        stats.total_bp += seq.len();
    });
    writer.flush().unwrap();

    stats
}
