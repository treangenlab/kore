# *kore* - accurate and robust bacterial genome completeness and contamination estimation

## Introduction
*kore* is a bacterial genome completeness and contamination estimator that utilizes species specific core genomes for it's estimations. It is pretty much a stand-in replacement for BUSCO/CheckM/CheckM2.

*kore* compares each assembly directly against the core genome of its species, and names the species of any contamination it finds. It is intended for species with core genomes in our database, not for distant relatives or novel species. In that case, we still recommend using CheckM2.

### Why use *kore*
1. **Accuracy** -- *kore* is nearly 1-2 points better on MAE than existing tools, and had lower per-species bias for most species.
2. **Speed** -- *kore* is nearly 2-3x faster than existing approaches, and more if you run large batches of genomes altogether.
3. **Annotation free** -- *kore* does not require any annotations, it simply matches DNA to DNA.
4. **Robust to errors and assembly N50** -- *kore* is much less affected by common assembly quality issues (fragmentation, low N50) than existing tools.
5. **Names the contaminant** -- contamination from another species in the database is reported with the species it comes from.

### How does *kore* work
For each species, *kore* builds a core genome: the sequence (locally collinear blocks, or LCBs) shared by the species' high quality genomes, found with a Parsnp2 core-genome alignment. All species' cores are combined into a single database.

To estimate, *kore* first detects the species of each assembly from the database, then aligns that species' core to the assembly. Completeness is the percent of the core found in the assembly. Contamination is the sum of two parts:
- **same species** -- the percent of the core found more than once (e.g. a second strain)
- **other species** -- the core of each other species in the database found in the assembly, leaving out the sequence that genomes of the assembly's own species also carry (masks)

## Some notes before running
If you are working with rare genomes that have less sampled diversity (think biodiversity studies), then *kore* may not support the genomes you are working with. In this case, we suggest using other tools (although note that these estimates will suffer from the same biases as shown in our paper).

If an assembly does not match its species' core (too few of the core's LCBs are found, or too few of them are as identical as in the species' own genomes), *kore* does not guess: its estimates are reported as `NA` with a warning. This is usually a species that is not in the database, or a very divergent strain.

## Requirements
`kore estimate` has no requirements beyond *kore* itself.

`kore build` also needs [Parsnp2](https://github.com/marbl/parsnp) (v2.1 or newer), [skani](https://github.com/bluenote-1577/skani) and [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) on your `PATH`. The easiest way to get these is through conda, but these will be automatically instaled if you install *kore* directly. 

## Installation Instructions

#### Conda Installation (recommended)
*kore* is easily installable via bioconda using conda/mamba (miniforge available [here](https://github.com/conda-forge/miniforge))
```
mamba create -n kore kore
mamba activate kore
```

#### From source
*kore* is written in Rust. You will need Rust 1.88 or newer (https://rustup.rs) and a C compiler.
```
git clone https://github.com/treangenlab/kore.git
cd kore
cargo build --release
```

The binary is then at `target/release/kore`. You should be able to test if installed successfully by running
```
kore --help
```

To run the tests (the build test needs the tools above, otherwise it is skipped):
```
cargo test --release
```

## Basic Usage
For full details, please see the wiki.

### Getting the *kore* database
We provide a prebuilt database of core genomes for 1,023 species (GTDB r232): **[link coming soon]**. It is ~17 GB.

Once downloaded, you can point *kore* to it with `-d/--db`, or set it once so you can leave out `-d`:
```
export KORE_DB=/path/to/kore_db
```

### Estimating completeness and contamination
To estimate completeness and contamination, use the `kore estimate` command:
```
kore estimate -d /path/to/kore_db -g /path/to/assemblies/*.fna -o /path/to/output
```

This detects the species of each assembly, estimates its completeness against that species' core, and its contamination against every core in the database. Assemblies can be from any number of different species.

If you already know the species, you can skip species detection with `-s/--species` (contamination from other species is still estimated):
```
kore estimate -d /path/to/kore_db -s "Escherichia coli" -g /path/to/assemblies/*.fna -o /path/to/output
```

or for assemblies of different species, with a tab-separated file of assembly path and species (header optional):
```
kore estimate -d /path/to/kore_db --species-table assemblies.tsv -o /path/to/output
```

Species names can be written as in GTDB (`Escherichia coli`, `s__Bacillus_A anthracis`) or as *kore* folder names (`escherichia_coli`). A list of assemblies can also be given with `--file-input` (one path per line).

### Building a *kore* database
If your species is not in our database, you can build its core with `kore build`. You only give the genomes and their species; species names are checked against the GTDB taxonomy (downloaded automatically on first use into `~/.kore/gtdb_r232`, or `$KORE_TAXONOMY`), and each species' GTDB family decides which relatives it is masked against.

For a single species:
```
kore build -g /path/to/genomes/*.fna -s "Escherichia coli" -o kore_db
```

or for many species at once, with a tab-separated file of genome path and species (header optional):
```
kore build --species-table genomes.tsv -o kore_db --threads 64 --cores-at-once 8
```

Species with fewer than 5 genomes in a table are skipped (`--min-genomes`), and misspelled species names stop the build with suggestions. `--allow-unlisted` builds species missing from GTDB (without masks).

To add species to an existing database (including ours), use `-f/--force`. Species already built from the same genomes are kept, and masks are only recomputed for the new species' family. To build species as separate cluster jobs, use `-f --skip-aggregate` for each job and then `kore build --aggregate-only -o kore_db` at the end.

#### Other Parameters
Run `kore build --help` or `kore estimate --help` for every parameter. The defaults are what we used in the paper.

#### Description of Output Files
`kore estimate` writes the following tab-separated files to the output folder (`-o`):
- `report.tsv` -- one row per assembly (columns below)
- `contaminants.tsv` -- every other species considered for each assembly, and whether it was called
- `alignments.tsv` -- every core alignment in every assembly (only with `--alignments`)
- `kore.log` -- the log, including the command that was run


`kore build` writes a database folder:
```
kore_db/
  species.tsv          one row per species: GTDB lineage, core size, genome size, masks
  aggregate.fasta/.mmi every core together (species detection and contamination)
  masks/               per species, the parts of its relatives' cores its own genomes carry
  cores/<species>/     core.fasta, core_lcbs.tsv, genomes.txt, info.tsv, calibration.tsv, kore.log
```

## Issues
*kore* is still in early versions and we welcome any feedback, suggestions, or feature requests. If you run into any issues, please email Ryan Doughty at rdd4@rice.edu or raise an issue on github

## Citation
> [!NOTE]
> **Citation:** If you use *kore* in your research, please cite:
>
> **[preprint coming soon]**

## AI Usage Declaration
Portions of the code base were generated or edited using Claude. All outputs were reviewed and edited as appropriate by the authors. The authors take full responsibility for the accuracy, integrity, and originality of the code, analyses, and content presented in this repository.
