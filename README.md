# kore
*kore* is a bacterial genome completeness and contamination estimator using core-genomes. It is pretty much a stand-in replacement for BUSCO/CheckM/CheckM2. 

kore is intended for use for species with existing built core genomes in our database, not for distant relatives or novel species. In that case, we still recommend using CheckM2.  
## Usage

kore has three commands:

```
kore build     # build a core for one species from its genomes (species not in the kore db)
kore download  # download species cores or the whole kore db (not available yet)
kore estimate  # completeness and contamination of assemblies
```

### Estimate

```
# default: whole kore db - each assembly's species is detected, completeness is estimated against its core and
# contamination against every core (other species are named)
kore estimate --db kore_db -g assembly1.fna assembly2.fna.gz -o results

# one species from the db: completeness and same-species redundancy only (other species are not looked for)
kore estimate --db kore_db --species escherichia_coli -g assembly.fna -o results

# your own core from kore build: as --species
kore estimate --core my_species/core.fasta -g assembly.fna -o results
```

`--db` can be left out when `$KORE_DB` points to the db. With `--species` or `--core`, kore warns that contamination is
only same-species redundancy: a close relative's genome can still raise it, but is not named.

If an assembly is not the species (fewer than 60 % of its found core LCBs as identical as in the species' own genomes:
each LCB's bar is its lowest identity among the build genomes minus 0.5, at most 97 %; with the whole db: not
any species in it), its estimates are reported as `NA` with a warning.

Output (`-o`): `report.csv` (one row per assembly: species, completeness, contamination, redundancy, foreign species,
divergent sequence), `contaminants.csv` (whole db only: every foreign species considered), `alignments.csv` and
`kore.log`.

### Build

`kore build` makes a kore db: a core for each species (parsnp2 core-genome alignment of its genomes), aggregated so
`kore estimate` can detect species and contamination. You only give genomes and their species; species names are
checked against the GTDB taxonomy (downloaded automatically on first use into `~/.kore/gtdb_r232`, or `$KORE_TAXONOMY`),
and each species' GTDB family decides which relatives it is masked against.

```
# one species
kore build -g genomes/*.fna --species "Escherichia coli" -o kore_db

# many species: a tab-separated table of genome path and species (header optional)
kore build --species-table genomes.tsv -o kore_db --threads 64 --cores-at-once 8
```

Species names can be written as in GTDB (`Escherichia coli`, `s__Bacillus_A anthracis`) or as kore folder names
(`escherichia_coli`); misspelled names stop the build with suggestions. `--min-genomes` (default 5) skips species with
too few genomes in a table; `--allow-unlisted` builds species missing from GTDB (with no relatives).

Building into an existing db needs `--force` / `-f` (kore stops otherwise, so a db is never overwritten by accident); it
then extends the db: species already built from the same genomes are skipped, and masks are only
recomputed for the new species' family. To build species as separate cluster jobs, use `-f --skip-aggregate` per job and
then `kore build --aggregate-only -o kore_db`.

The db folder:

```
kore_db/
  species.tsv          one row per species: GTDB lineage, core size, genome size, masks
  aggregate.fasta/.mmi every core together (species detection and contamination)
  masks/               per species, the parts of its relatives' cores its own genomes carry
  cores/<species>/     core.fasta, core_lcbs.csv, genomes.txt, info.tsv, calibration.tsv, kore.log
```
