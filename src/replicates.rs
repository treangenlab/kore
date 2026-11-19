use crate::consts::*;
use crate::util::*;

use log::*;

use std::collections::{HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

//"rep_0_lcb12" -> "rep_0"
fn replicate_of(id: &str) -> &str {
    id.rsplit_once("_lcb").map(|(rep, _)| rep).unwrap_or(id)
}

fn run_command(cmd: &mut Command, name: &str) {
    let output = cmd.output().unwrap_or_else(|e| {
        error!("Unable to run {}: {}", name, e);
        std::process::exit(1);
    });
    if !output.status.success() {
        error!("{} failed:\n{}", name, String::from_utf8_lossy(&output.stderr));
        std::process::exit(1);
    }
}

//all-vs-all megablast of every replicate's consensus LCBs
fn run_blast(fastas: &[PathBuf], blast_dir: &Path, threads: usize) -> PathBuf {
    fs::create_dir_all(blast_dir).unwrap_or_else(|e| {
        error!("Unable to create directory {}: {}", blast_dir.display(), e);
        std::process::exit(1);
    });

    let all_fasta = blast_dir.join("all_lcbs.fasta");
    let records: Vec<FastaRecord> = fastas.iter().flat_map(|f| read_fasta(f)).collect();
    write_fasta(&records.iter().collect::<Vec<_>>(), &all_fasta);

    let db = blast_dir.join("all_lcbs");
    run_command(
        Command::new(MAKEBLASTDB_BIN)
            .arg("-in").arg(&all_fasta)
            .args(["-dbtype", "nucl"])
            .arg("-out").arg(&db),
        MAKEBLASTDB_BIN,
    );

    let hits = blast_dir.join("all_vs_all.tsv");
    debug!("Running blastn all-vs-all on {} LCBs -> {}", records.len(), hits.display());
    run_command(
        Command::new(BLASTN_BIN)
            .args(["-task", "megablast"])
            .arg("-query").arg(&all_fasta)
            .arg("-db").arg(&db)
            .args(["-evalue", BLAST_EVALUE, "-max_target_seqs", BLAST_MAX_TARGETS])
            .args(["-outfmt", "6 qseqid sseqid pident qstart qend"])
            .arg("-num_threads").arg(threads.to_string())
            .arg("-out").arg(&hits),
        BLASTN_BIN,
    );

    hits
}

//for each (LCB, other replicate), the regions of the LCB hit at >= min_identity
fn load_hit_intervals(hits: &Path, min_identity: f64) -> HashMap<(String, String), Vec<(usize, usize)>> {
    let file = File::open(hits).unwrap_or_else(|e| {
        error!("Unable to open {}: {}", hits.display(), e);
        std::process::exit(1);
    });

    let mut intervals: HashMap<(String, String), Vec<(usize, usize)>> = HashMap::new();
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        let f: Vec<&str> = line.split('\t').collect();
        let (query, subject) = (f[0], f[1]);
        let (q_rep, s_rep) = (replicate_of(query), replicate_of(subject));
        if q_rep == s_rep || f[2].parse::<f64>().unwrap() < min_identity {
            continue;
        }

        let (start, end): (usize, usize) = (f[3].parse().unwrap(), f[4].parse().unwrap());
        intervals.entry((query.to_string(), s_rep.to_string()))
            .or_default()
            .push((start.min(end), start.max(end)));
    }
    intervals
}

//stretches of an LCB that are covered in every other replicate
fn shared_segments(
    record: &FastaRecord,
    replicates: &[String],
    intervals: &HashMap<(String, String), Vec<(usize, usize)>>,
) -> Vec<(usize, usize)> {
    let own = replicate_of(record.id());
    replicates.iter().filter(|r| r.as_str() != own).fold(vec![(1, record.seq.len())], |shared, other| {
        let hits = intervals.get(&(record.id().to_string(), other.clone())).map_or(vec![], |h| merge_intervals(h));
        intersect_intervals(&shared, &hits)
    })
}

fn segment_record(record: &FastaRecord, (start, end): (usize, usize)) -> FastaRecord {
    let rest = record.header.split_once(' ').map_or("", |(_, r)| r);
    FastaRecord {
        header: format!("{}:{}-{} lcb_len={} {}", record.id(), start, end, record.seq.len(), rest),
        seq: record.seq[start - 1..end].to_vec(),
    }
}

//trims every LCB of the first replicate to the stretches present in all replicates
pub fn merge_replicates(
    fastas: &[PathBuf],
    blast_dir: &Path,
    min_identity: f64,
    min_length: usize,
    threads: usize,
) -> Vec<FastaRecord> {
    let per_rep: Vec<Vec<FastaRecord>> = fastas.iter().map(|f| read_fasta(f)).collect();

    if fastas.len() == 1 {
        per_rep[0].iter().map(|r| segment_record(r, (1, r.seq.len()))).collect()
    } else {
        let hits = run_blast(fastas, blast_dir, threads);
        let intervals = load_hit_intervals(&hits, min_identity);
        let replicates: Vec<String> = (0..fastas.len()).map(|i| format!("rep_{}", i)).collect();

        let segments: Vec<FastaRecord> = per_rep[0].iter()
            .flat_map(|r| shared_segments(r, &replicates, &intervals).into_iter().map(move |seg| segment_record(r, seg)))
            .filter(|r| r.seq.len() >= min_length)
            .collect();

        let total_bp = |records: &[FastaRecord]| records.iter().map(|r| r.seq.len()).sum::<usize>();
        let n_lcbs = segments.iter().map(|r| r.id().split_once(':').map_or(r.id(), |(lcb, _)| lcb)).collect::<HashSet<_>>().len();
        info!("Trimmed rep_0 ({} LCBs, {} bp) to {} segments from {} LCBs ({} bp) present in all replicates",
            per_rep[0].len(), total_bp(&per_rep[0]), segments.len(), n_lcbs, total_bp(&segments));
        segments
    }
}

//renames the final core records to lcb1, lcb2, ...
pub fn rename_core(core: &mut [FastaRecord]) {
    for (i, record) in core.iter_mut().enumerate() {
        let rest = record.header.split_once(' ').map_or("", |(_, r)| r).to_string();
        record.header = format!("lcb{} {}", i + 1, rest);
    }
}

//"lcb1 ... n_genomes=25 diversity=0.0123" -> one tsv row
fn tsv_row(record: &FastaRecord) -> String {
    let field = |key: &str| record.header.split_whitespace()
        .find_map(|f| f.strip_prefix(key))
        .unwrap_or("");
    [record.id(), &record.seq.len().to_string(), field("diversity="), field("n_genomes=")].join("\t")
}

pub fn write_core_tsv(core: &[FastaRecord], path: &Path) {
    let file = File::create(path).unwrap_or_else(|e| {
        error!("Unable to create {}: {}", path.display(), e);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(file);
    writeln!(writer, "id\tlength\tdiversity\tn_genomes").unwrap();
    for record in core {
        writeln!(writer, "{}", tsv_row(record)).unwrap();
    }
    writer.flush().unwrap();
}
