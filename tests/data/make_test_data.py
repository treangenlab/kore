#!/usr/bin/env python3
"""Write the test data in tests/data from the kore paper's 380-species benchmark (real GTDB r232 genomes).

    python3 tests/data/make_test_data.py

Only runs where the benchmark lives (/dodo/rdd4/kore_paper_experiments); the files it writes are committed, and the
accessions are listed in accessions.tsv. Species (all ~1 Mb genomes, so the data stays small and tests fast):
  Mycoplasma capricolum        target species
  Mycoplasma feriruminatoris   its sister species (same genus): the contaminant
  Chlamydia trachomatis        unrelated species (other phylum)
  Mycoplasmoides pneumoniae    not in the test db: an assembly kore cannot assign (NA)

Written:
  gtdb/                      GTDB r232 taxonomy of the three db species only (tests set $KORE_TAXONOMY: no download)
  genomes/<species>/*.fna.gz 5 build genomes per species (the first 5 of the benchmark core's cluster-diverse sample):
                             input for kore build (species_table.tsv)
  cores/<species>/           the benchmark db's core (core.fasta.gz) and lineage (info.tsv): a ready-made core, so the
                             estimate test needs no parsnp (the test aggregates them with kore build --aggregate-only)
  assemblies/*.fna.gz        held-out test genomes (released after the benchmark's 2024-01-01 cutoff, not build genomes)
                             cut into 5-50 kb contigs, and mixtures of them; true values in assemblies/truth.tsv
"""
import csv
import gzip
import random
import shutil
from pathlib import Path

SEED = 7
HERE = Path(__file__).resolve().parent
BENCH = Path("/dodo/rdd4/kore_paper_experiments")
DB = BENCH / "kore_db_380_rel"
TEST_GENOMES = BENCH / "fig3" / "manifests" / "test_genomes.tsv"
SPECIES = ["mycoplasma_capricolum", "mycoplasma_feriruminatoris", "chlamydia_trachomatis"]
BUILD_GENOMES = 5
#held-out test genomes (fig3/manifests/test_genomes.tsv)
TARGET = "GCA_051053445.1"       # Mycoplasma capricolum
STRAIN = "GCA_051052485.1"       # another M. capricolum: same-species contamination
SISTER = "GCA_047715995.1"       # Mycoplasma feriruminatoris: foreign contamination
OTHER = "GCA_039778415.1"        # Chlamydia trachomatis
UNKNOWN = "GCA_041931765.1"      # Mycoplasmoides pneumoniae (not in the test db)

rng = random.Random(SEED)


def read_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    seqs, name = {}, None
    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(f">{n}\n" + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n" for n, s in records)
    with gzip.GzipFile(path, "wb", mtime=0) as f:  # mtime 0: the same bytes on every run
        f.write(text.encode())


def contigs(genome, lo=5_000, hi=50_000):
    """every replicon cut into random 5-50 kb contigs"""
    out = []
    for seq in genome.values():
        i = 0
        while i < len(seq):
            out.append(seq[i:i + rng.randint(lo, hi)])
            i += len(out[-1])
        if len(out) > 1 and len(out[-1]) < 1_000:
            out[-2] += out.pop()
    return out


def take(cs, frac):
    """random contigs adding up to at most frac of the total length (as close to it as they fit)"""
    total, out = sum(map(len, cs)), []
    for c in rng.sample(cs, len(cs)):
        if sum(map(len, out)) + len(c) <= frac * total:
            out.append(c)
    return out


def main():
    for d in ["gtdb", "genomes", "cores", "assemblies"]:
        shutil.rmtree(HERE / d, ignore_errors=True)
    lineage = {r["species"]: r for r in csv.DictReader(open(DB / "species.tsv"), delimiter="\t")}
    tests = {r["accession"]: r for r in csv.DictReader(open(TEST_GENOMES), delimiter="\t")}
    accessions = [("accession", "species", "use")]

    # GTDB taxonomy: the three db species (and one archaeon so ar53 is not empty)
    (HERE / "gtdb").mkdir(parents=True)
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "gtdb_species"]
    tax = "".join(f"RS_GCF_00000000{i}.1\t" + ";".join(f"{k[0] if k != 'gtdb_species' else 's'}__{lineage[s][k]}" for k in ranks) + "\n"
                  for i, s in enumerate(SPECIES))
    for name, text in [("bac120_taxonomy_r232.tsv.gz", tax),
                       ("ar53_taxonomy_r232.tsv.gz", "RS_GCF_900000001.1\td__Archaea;p__A;c__A;o__A;f__A;g__A;s__A one\n")]:
        with gzip.GzipFile(HERE / "gtdb" / name, "wb", mtime=0) as f:
            f.write(text.encode())

    table = ["genome\tspecies"]
    for sp in SPECIES:
        # build genomes: the core's cluster-diverse sample order
        sampled = [l.split("\t")[1] for l in (DB / "cores" / sp / "sampled_genomes.tsv").read_text().splitlines()[1:]]
        for path in list(dict.fromkeys(sampled))[:BUILD_GENOMES]:
            acc = Path(path).parent.name
            out = HERE / "genomes" / sp / f"{acc}.fna.gz"
            write_fasta(out, read_fasta(path).items())
            table.append(f"genomes/{sp}/{out.name}\t{lineage[sp]['gtdb_species']}")
            accessions.append((acc, lineage[sp]["gtdb_species"], "build genome"))
        # ready-made core and its lineage
        (HERE / "cores" / sp).mkdir(parents=True)
        with open(DB / "cores" / sp / "core.fasta", "rb") as f, gzip.GzipFile(HERE / "cores" / sp / "core.fasta.gz", "wb", mtime=0) as g:
            g.write(f.read())
        info = {"species": sp, **{k: lineage[sp][k] for k in ranks}}
        (HERE / "cores" / sp / "info.tsv").write_text("".join(f"{k}\t{v}\n" for k, v in info.items()))
    (HERE / "species_table.tsv").write_text("\n".join(table) + "\n")

    # assemblies from held-out test genomes
    g = {a: contigs(read_fasta(tests[a]["fasta"])) for a in [TARGET, STRAIN, SISTER, OTHER, UNKNOWN]}
    for a, use in [(TARGET, "assemblies: target"), (STRAIN, "target_strain30: same-species contaminant"),
                   (SISTER, "target_sister30: foreign contaminant"), (OTHER, "other_complete"), (UNKNOWN, "unknown: not in the db")]:
        accessions.append((a, tests[a]["species"], use))
    bp = lambda cs: sum(map(len, cs))  # noqa: E731
    named = lambda p, cs: [(f"{p}_{i + 1}", c) for i, c in enumerate(cs)]  # noqa: E731
    half, sister, strain = take(g[TARGET], 0.5), take(g[SISTER], 0.3), take(g[STRAIN], 0.3)
    A = HERE / "assemblies"
    write_fasta(A / "target_complete.fna.gz", named("target", g[TARGET]))
    write_fasta(A / "target_half.fna.gz", named("target", half))
    write_fasta(A / "target_sister30.fna.gz", named("target", g[TARGET]) + named("sister", sister))
    write_fasta(A / "target_strain30.fna.gz", named("target", g[TARGET]) + named("strain", strain))
    write_fasta(A / "other_complete.fna.gz", named("other", g[OTHER]))
    write_fasta(A / "unknown.fna.gz", named("unknown", g[UNKNOWN]))
    # truth: completeness = share of the target genome kept, contamination = added bp / target genome bp
    t = bp(g[TARGET])
    truth = [("target_complete", "mycoplasma_capricolum", 100.0, 0.0),
             ("target_half", "mycoplasma_capricolum", 100.0 * bp(half) / t, 0.0),
             ("target_sister30", "mycoplasma_capricolum", 100.0, 100.0 * bp(sister) / t),
             ("target_strain30", "mycoplasma_capricolum", 100.0, 100.0 * bp(strain) / t),
             ("other_complete", "chlamydia_trachomatis", 100.0, 0.0), ("unknown", "NA", float("nan"), float("nan"))]
    (A / "truth.tsv").write_text("assembly\tspecies\tcompleteness\tcontamination\n"
                                 + "".join(f"{a}\t{s}\t{c:.2f}\t{x:.2f}\n" for a, s, c, x in truth))
    (HERE / "accessions.tsv").write_text("".join("\t".join(r) + "\n" for r in accessions))
    total = sum(f.stat().st_size for f in HERE.rglob("*") if f.is_file())
    print(f"wrote {HERE} ({total / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
