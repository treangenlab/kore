use std::path::Path;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use log::*;
use simplelog::{CombinedLogger, TermLogger, WriteLogger, Config, TerminalMode, ColorChoice};

use crate::consts::KORE_VERSION;

/// Initialize logging so that every log is written both to the terminal and the log file
pub fn init_logging(level: log::LevelFilter, log_path: &Path) {
    // Make sure the directory that will hold the log file exists.
    if let Some(parent) = log_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).unwrap_or_else(|e| {
                eprintln!("Unable to create directory for log file {}: {}", parent.display(), e);
                std::process::exit(1);
            });
        }
    }

    let mut file = File::create(log_path).unwrap_or_else(|e| {
        eprintln!("Unable to create log file {}: {}", log_path.display(), e);
        std::process::exit(1);
    });

    // Record the exact command at the top of the log file only (not the terminal) for reproducibility.
    let command = std::env::args()
        .map(|a| {
            if a.is_empty() || a.chars().any(|c| c.is_whitespace() || c == '\'' || c == '"') {
                format!("'{}'", a.replace('\'', "'\\''"))
            } else {
                a
            }
        })
        .collect::<Vec<_>>()
        .join(" ");
    writeln!(file, "Command: {}\n", command).unwrap_or_else(|e| {
        eprintln!("Unable to write to log file {}: {}", log_path.display(), e);
        std::process::exit(1);
    });

    CombinedLogger::init(vec![
        TermLogger::new(level, Config::default(), TerminalMode::Stdout, ColorChoice::Auto),
        WriteLogger::new(level, Config::default(), file),
    ])
    .unwrap();
}

/// Pick the log level from the --verbose/--debug flags
pub fn log_level(verbose: bool, debug: bool) -> log::LevelFilter {
    if verbose {
        log::LevelFilter::Trace
    } else if debug {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    }
}

/// Emit the kore startup banner through the logger so it lands in the .log file
/// as well as the terminal. Must be called after `init_logging`.
pub fn print_banner() {
    info!(
        "kore v{}\nDeveloped by Ryan Doughty (Rice University)\nCorrespondence: rdd4@rice.edu\n",
        KORE_VERSION
    );
}

/// Check that a file has a fasta extension (optionally gzipped)
pub fn check_fasta(file: &str) -> bool {
    let file = file.strip_suffix(".gz").unwrap_or(file);
    [".fasta", ".fa", ".fna", ".fas"].iter().any(|ext| file.ends_with(ext))
}

/// Read a .txt file containing one path per line (blank lines skipped)
pub fn read_file_list(path: &str) -> Vec<String> {
    let contents = fs::read_to_string(path).unwrap_or_else(|e| {
        error!("Unable to read file list {}: {}", path, e);
        std::process::exit(1);
    });
    contents
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty())
        .map(String::from)
        .collect()
}

/// Gather genome paths from the positional/flag list and an optional --file-input list,
/// exiting if none are provided or any are missing/not fasta.
pub fn collect_genomes(genomes: &[String], file_input: &Option<String>) -> Vec<String> {
    let mut all: Vec<String> = genomes.to_vec();
    if let Some(list) = file_input {
        all.extend(read_file_list(list));
    }

    if all.is_empty() {
        error!("No genomes provided (use -g/--genomes or --file-input)");
        std::process::exit(1);
    }

    for g in &all {
        if !Path::new(g).exists() {
            error!("Genome file does not exist: {}", g);
            std::process::exit(1);
        }
        if !check_fasta(g) {
            error!("Genome file is not a fasta(.gz) file: {}", g);
            std::process::exit(1);
        }
    }
    all
}

/// Check whether an executable is available on PATH (without running it)
pub fn in_path(bin: &str) -> bool {
    std::env::var_os("PATH")
        .map(|paths| std::env::split_paths(&paths).any(|dir| dir.join(bin).is_file()))
        .unwrap_or(false)
}

const LINE_WRAP: usize = 80;

pub struct FastaRecord {
    pub header: String, //full header line without '>'
    pub seq: Vec<u8>,
}

impl FastaRecord {
    pub fn id(&self) -> &str {
        self.header.split_whitespace().next().unwrap()
    }
}

pub fn read_fasta(path: &Path) -> Vec<FastaRecord> {
    let file = File::open(path).unwrap_or_else(|e| {
        error!("Unable to open {}: {}", path.display(), e);
        std::process::exit(1);
    });

    let mut records: Vec<FastaRecord> = Vec::new();
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        if let Some(header) = line.strip_prefix('>') {
            records.push(FastaRecord { header: header.to_string(), seq: Vec::new() });
        } else if let Some(record) = records.last_mut() {
            record.seq.extend(line.trim_end().bytes());
        }
    }
    records
}

pub fn write_fasta(records: &[&FastaRecord], path: &Path) {
    let file = File::create(path).unwrap_or_else(|e| {
        error!("Unable to create {}: {}", path.display(), e);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(file);
    for record in records {
        writeln!(writer, ">{}", record.header).unwrap();
        for chunk in record.seq.chunks(LINE_WRAP) {
            writer.write_all(chunk).unwrap();
            writer.write_all(b"\n").unwrap();
        }
    }
    writer.flush().unwrap();
}

//sorted union of 1-based inclusive intervals
pub fn merge_intervals(intervals: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut sorted = intervals.to_vec();
    sorted.sort();

    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (start, end) in sorted {
        match merged.last_mut() {
            Some((_, e)) if start <= *e + 1 => *e = (*e).max(end),
            _ => merged.push((start, end)),
        }
    }
    merged
}

pub fn covered_bp(intervals: &[(usize, usize)]) -> usize {
    merge_intervals(intervals).iter().map(|(s, e)| e - s + 1).sum()
}

//overlap of two sorted, merged interval lists
pub fn intersect_intervals(a: &[(usize, usize)], b: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        let start = a[i].0.max(b[j].0);
        let end = a[i].1.min(b[j].1);
        if start <= end {
            out.push((start, end));
        }
        if a[i].1 < b[j].1 { i += 1 } else { j += 1 }
    }
    out
}

//p-th percentile (0-100) with linear interpolation between ranks, 0 if empty
pub fn percentile(values: &mut [f64], p: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let rank = p / 100.0 * (values.len() - 1) as f64;
    let (lo, hi) = (rank.floor() as usize, rank.ceil() as usize);
    values[lo] + (values[hi] - values[lo]) * (rank - lo as f64)
}
