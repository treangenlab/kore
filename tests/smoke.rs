//Smoke tests: kore runs to completion on the small real-genome data in tests/data (GTDB r232 genomes from the paper's
//benchmark, see tests/data/make_test_data.py and accessions.tsv) and writes its output files.
use std::fs;
use std::io::Read;
use std::path::PathBuf;
use std::process::Command;

const SPECIES: [&str; 3] = ["mycoplasma_capricolum", "mycoplasma_feriruminatoris", "chlamydia_trachomatis"];

fn data() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

//fresh output folder under target/
fn tmp(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

//kore with the fake GTDB taxonomy in tests/data (no download)
fn kore(args: &[&str]) {
    let out = Command::new(env!("CARGO_BIN_EXE_kore")).args(args).env("KORE_TAXONOMY", data().join("gtdb"))
        .env_remove("KORE_DB").output().unwrap();
    assert!(out.status.success(), "kore {} failed:\n{}{}", args.join(" "),
            String::from_utf8_lossy(&out.stdout), String::from_utf8_lossy(&out.stderr));
}

fn has_rows(p: PathBuf) {
    let n = fs::read_to_string(&p).unwrap_or_else(|_| panic!("missing {}", p.display())).lines().count();
    assert!(n > 1, "{} is empty", p.display());
}

//estimate every test assembly against db, check the outputs
fn estimate(db: &PathBuf, out: &PathBuf) {
    let mut args = vec!["estimate".to_string(), "-d".into(), db.display().to_string(), "-o".into(), out.display().to_string(),
                        "-t".into(), "4".into(), "--alignments".into(), "-g".into()];
    for e in fs::read_dir(data().join("assemblies")).unwrap() {
        let p = e.unwrap().path();
        if p.to_string_lossy().ends_with(".fna.gz") {
            args.push(p.display().to_string());
        }
    }
    kore(&args.iter().map(String::as_str).collect::<Vec<_>>());
    for f in ["report.csv", "contaminants.csv", "alignments.csv", "kore.log"] {
        has_rows(out.join(f));
    }
}

//kore estimate on a db made from the ready-made cores in tests/data/cores (no external tools needed)
#[test]
fn estimate_runs() {
    let db = tmp("db_from_cores");
    for sp in SPECIES {
        let (from, dir) = (data().join("cores").join(sp), db.join("cores").join(sp));
        fs::create_dir_all(&dir).unwrap();
        fs::copy(from.join("info.tsv"), dir.join("info.tsv")).unwrap();
        let mut core = Vec::new();
        flate2::read::GzDecoder::new(fs::File::open(from.join("core.fasta.gz")).unwrap()).read_to_end(&mut core).unwrap();
        fs::write(dir.join("core.fasta"), core).unwrap();
        let genomes: Vec<String> = fs::read_dir(data().join("genomes").join(sp)).unwrap()
            .map(|e| e.unwrap().path().display().to_string()).collect();
        fs::write(dir.join("genomes.txt"), genomes.join("\n")).unwrap();
    }
    kore(&["build", "--aggregate-only", "-t", "4", "-o", db.to_str().unwrap()]);
    for f in ["species.tsv", "aggregate.mmi", "masks/mycoplasma_capricolum.tsv", "cores/mycoplasma_capricolum/calibration.tsv"] {
        assert!(db.join(f).exists(), "missing {}", f);
    }
    estimate(&db, &tmp("estimate_out"));
}

//full kore build (parsnp2, skani, blast), then estimate. Skipped when the tools are not on PATH, unless
//$KORE_TEST_REQUIRE_TOOLS is set (so CI with the conda env cannot skip it silently)
#[test]
fn build_runs() {
    let have = ["parsnp", "skani", "blastn", "makeblastdb"].iter()
        .all(|b| std::env::var_os("PATH").is_some_and(|p| std::env::split_paths(&p).any(|d| d.join(b).is_file())));
    if !have {
        assert!(std::env::var_os("KORE_TEST_REQUIRE_TOOLS").is_none(), "parsnp/skani/blastn/makeblastdb not on PATH");
        eprintln!("build_runs skipped: parsnp, skani and blast are not on PATH");
        return;
    }
    let db = tmp("db_built");
    let table = data().join("species_table.tsv");
    kore(&["build", "--species-table", table.to_str().unwrap(), "-o", db.to_str().unwrap(), "-t", "4", "--cores-at-once", "3"]);
    for sp in SPECIES {
        has_rows(db.join("cores").join(sp).join("core.fasta"));
    }
    has_rows(db.join("species.tsv"));
    estimate(&db, &tmp("build_estimate_out"));
}
