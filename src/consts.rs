pub const KORE_VERSION: &str = env!("CARGO_PKG_VERSION");
pub const DEFAULT_THREADS: usize = 30; //number of threads
pub const PARSNP_BIN: &str = "parsnp"; //parsnp2 executable (expected on PATH via the conda env)
pub const SKANI_BIN: &str = "skani"; //skani executable for ANI clustering (expected on PATH via the conda env)
pub const DEFAULT_CLUSTER_ANI: f64 = 99.5; //ANI (%) genomes are clustered at before sampling
pub const REPORT_CLUSTER_ANI: [f64; 4] = [99.9, 99.5, 99.0, 98.0]; //ANI (%) cuts whose cluster counts are recorded (info.tsv, species.tsv; information only)
pub const SKANI_MISSING_ANI: f64 = 80.0; //ANI given to pairs skani does not report (its reporting floor)
pub const SAMPLED_GENOMES: &str = "sampled_genomes.tsv"; //genomes used in each replicate, with their cluster
pub const SKANI_ANI: &str = "skani_ani.tsv.gz"; //all-vs-all skani ANI kept in the output folder and reused by later builds
pub const SKANI_GENOMES: &str = "skani_genomes.txt"; //genome set SKANI_ANI covers
pub const SKANI_KEEP_COLUMNS: usize = 5; //Ref_file, Query_file, ANI, Align_fraction_ref, Align_fraction_query
pub const BLASTN_BIN: &str = "blastn"; //blastn executable (expected on PATH via the conda env)
pub const MAKEBLASTDB_BIN: &str = "makeblastdb"; //makeblastdb executable (expected on PATH via the conda env)
pub const BLAST_EVALUE: &str = "1e-10"; //evalue for all-vs-all LCB blast
pub const BLAST_MAX_TARGETS: &str = "50"; //an LCB can hit several shorter LCBs per replicate
pub const DEFAULT_SEED: u64 = 42; //seed for shuffling genomes (within clusters when clustering)
pub const DEFAULT_N_GENOMES: usize = 100; //number of genomes per parsnp replicate
pub const DEFAULT_REPLICATES: usize = 1; //number of parsnp replicates
pub const MAX_PARALLEL_PARSNP: usize = 4; //max parsnp jobs running at once
pub const DEFAULT_MIN_LCB_LENGTH: usize = 150; //minimum consensus length of an LCB
pub const DEFAULT_MAX_LCB_DIVERSITY: f64 = 0.15; //max fraction of variable columns in an LCB (raised to DIVERSITY_PERCENTILE if higher)
pub const DIVERSITY_PERCENTILE: f64 = 95.0; //percentile of LCB diversity (LCBs >= min length, all replicates) used as the cutoff when above max diversity
pub const DEFAULT_REPLICATE_IDENTITY: f64 = 95.0; //min blast identity (%) for an LCB hit in another replicate
pub const DEFAULT_BUILD_OUTPUT: &str = "kore_build"; //folder for build outputs
pub const CORE_FASTA: &str = "core.fasta"; //final core LCBs (the kore db)
pub const CORE_CSV: &str = "core_lcbs.csv"; //info on each LCB in the core fasta
pub const LOG_FILE: &str = "kore.log"; //log file written in the output folder
pub const INTERMEDIATE_DIR: &str = "intermediate"; //parsnp, consensus and blast files (deleted unless --keep-intermediate)
pub const DEFAULT_OUT_FOLDER: &str = "kore_output"; //folder for estimate outputs
pub const REPORT_CSV: &str = "report.csv"; //per genome report: species, completeness, contamination
pub const DB_ENV: &str = "KORE_DB"; //environment variable holding the default kore db folder
pub const ALIGNMENTS_CSV: &str = "alignments.csv"; //every LCB alignment in every genome
pub const MM_SECONDARY_RATIO: f32 = 0.5; //minimap2 -p, looser than default so second copies are still reported
pub const MM_MAX_SECONDARY: i32 = 20; //minimap2 -N
pub const MM_MIN_CHAIN_SCORE: i32 = 30; //minimap2 -m, asm20 default drops short LCBs
pub const MM_MIN_DP_SCORE: i32 = 40; //minimap2 -s, asm20 default drops short LCBs
pub const MIN_ALN_IDENTITY: f64 = 80.0; //min identity (%) of an alignment to count towards coverage
pub const MIN_ALN_LEN: i32 = 100; //min alignment block length to count towards coverage
pub const LCB_PRESENT_COV: f64 = 80.0; //min coverage (%) of an LCB for it to count as present
pub const LCB_PRESENT_IDENTITY: f64 = 90.0; //min identity (%) of alignments counted for LCB presence

//SPECIES MATCH (is the assembly the core's species at all? otherwise completeness is NA)
//Per LCB, its best alignment (largest aligned bp x identity); the LCB is found if that alignment covers
//>= SPECIES_MATCH_MIN_COVERAGE of it. Each LCB has its own bar, calibrated on the species' own build genomes
//(<species>/calibration.tsv, identity_min): bar = min(SPECIES_MATCH_IDENTITY, identity_min - SPECIES_MATCH_MARGIN), so an
//LCB that varies within the species gets a lower bar while conserved LCBs keep 97 %; LCBs (or cores) without
//calibration use SPECIES_MATCH_IDENTITY. The assembly matches the species when >= SPECIES_MATCH_FRACTION of found LCBs
//pass their bar and >= SPECIES_MATCH_MIN_LCBS LCBs are found. Within a species nearly every LCB is at 98-100 %; in a
//sister species (~88-95 % ANI) most are below 97 % (only conserved ones, ~5-20 %, reach it).
//Tested on 380 species (lcb_divergence/gate_full/calib_eval.tsv): accepts 99.94 % of own complete genomes and 99.19 % of
//own 10 % subsets, 2.35 % of same-genus other species (mostly GTDB boundary splits; median completeness reported 98.1 %).
pub const SPECIES_MATCH_IDENTITY: f64 = 97.0; //highest bar (%) an LCB's best alignment needs to count as same-species
pub const SPECIES_MATCH_MARGIN: f64 = 0.5; //an LCB's bar sits this far below its lowest identity in the species' own genomes
pub const SPECIES_MATCH_FRACTION: f64 = 60.0; //min % of found LCBs that pass their bar
pub const SPECIES_MATCH_MIN_COVERAGE: f64 = 50.0; //min % of an LCB its best alignment covers for the LCB to be found
pub const SPECIES_MATCH_MIN_LCBS: usize = 10; //fewer found LCBs than this: no evidence for the species
pub const NO_MATCH: &str = "NA"; //reported in place of estimates when the assembly does not match the species

//KORE DB (kore build -o <db>): one folder per species core plus the aggregate of all of them
pub const DEFAULT_DB: &str = "kore_db"; //default db folder for kore build
pub const DB_CORES_DIR: &str = "cores"; //<db>/cores/<species>/: core.fasta, core_lcbs.csv, genomes.txt, info.tsv, calibration*.tsv, kore.log
pub const DB_SPECIES: &str = "species.tsv"; //<db>/species.tsv: one row per species (taxonomy, core size, genome size, masks)
pub const AGGREGATE_FASTA: &str = "aggregate.fasta"; //every core's LCBs, named <species>|<lcb>
pub const AGGREGATE_INDEX: &str = "aggregate.mmi"; //minimap2 index of AGGREGATE_FASTA (species detection, contamination)
pub const DB_MASKS_DIR: &str = "masks"; //<db>/masks/<species>.tsv: relatives' core intervals the species' own genomes carry
pub const MASK_RELATIVES_SUFFIX: &str = ".relatives.tsv"; //relatives (and their core checksums) a mask was computed for
pub const SPECIES_GENOMES: &str = "genomes.txt"; //<species>/genomes.txt: every build genome of the species
pub const BUILD_OUTPUT: &str = "build_output.log"; //<species>/: terminal output of its build process (kept only when it fails)
pub const SPECIES_INFO: &str = "info.tsv"; //<species>/info.tsv: GTDB taxonomy, genome size, genome counts
pub const CALIBRATION_LCBS: &str = "calibration.tsv"; //<species>/: per LCB, identity of its best hit in sampled build genomes
pub const CALIBRATION_GENOMES: &str = "calibration_genomes.tsv"; //<species>/: per sampled build genome, species-match score
pub const CALIBRATION_N: usize = 20; //build genomes (cluster-diverse, from the parsnp samples) mapped back for calibration
pub const DB_SEP: char = '|'; //separates species and LCB in aggregate sequence names
pub const AGGREGATE_SECONDARY_RATIO: f32 = 0.2; //minimap2 -p against the aggregate: every relative's core must still be reported
pub const AGGREGATE_MAX_SECONDARY: i32 = 100; //minimap2 -N against the aggregate
pub const DEFAULT_MASK_GENOMES: usize = 100; //max build genomes per species used for its mask (seeded sample)
pub const MASK_IDENTITY: f64 = 97.0; //relative-core sequence found at >= this identity in a species' own genomes is masked for it
pub const MIN_DIAGNOSTIC_BP: usize = 10_000; //a relative with less diagnostic core than this cannot be told apart from the target
pub const DEFAULT_MIN_GENOMES: usize = 5; //species with fewer genomes in a --species-table are skipped
pub const DEFAULT_CORES_AT_ONCE: usize = 1; //species cores built at the same time (threads are split between them)

//GTDB TAXONOMY (downloaded on first use; species names are checked against it and families become mask groups)
pub const GTDB_RELEASE: &str = "r232";
pub const GTDB_URL: &str = "https://data.gtdb.ecogenomic.org/releases/release232/232.0";
pub const GTDB_FILES: [&str; 2] = ["bac120_taxonomy_r232.tsv.gz", "ar53_taxonomy_r232.tsv.gz"];
pub const TAXONOMY_ENV: &str = "KORE_TAXONOMY"; //folder holding (or to download) the GTDB taxonomy; default ~/.kore/gtdb_<release>
pub const UNLISTED_GROUP: &str = "unlisted"; //group of species not in GTDB (built with --allow-unlisted): no relatives

//CONTAMINATION (kore estimate, whole db)
pub const SPECIES_IDENTITY: f64 = 95.0; //alignments >= this identity count towards a species' core coverage (aggregate scan)
pub const SPECIES_CANDIDATES: usize = 5; //species with the highest aggregate-scan coverage verified core -> assembly
pub const SPECIES_SHARED_SEQUENCE: f64 = 0.5; //candidates whose core alignments share >= this fraction (of the smaller
//footprint) of assembly bases with the best-covered candidate's compete for the species (sister species align to the same contigs; a contaminant
//aligns to other contigs, so it never competes); among them the most identical one passing the species check is the species
pub const SPECIES_CANDIDATE_FRACTION: f64 = 0.5; //...plus any with >= this fraction of the best scan coverage (max 2x)
pub const CONTIG_MIN_ALIGNED: usize = 1_000; //a contig is assigned to a species only on >= this many aligned bp
pub const FOREIGN_IDENTITY: f64 = 98.0; //alignments >= this identity count towards a foreign species' diagnostic core
pub const DEFAULT_CALL_BREADTH: f64 = 1.0; //min % of a foreign species' diagnostic core covered to call it present
pub const DEFAULT_CALL_LCBS: usize = 10; //min LCBs with >= CALL_LCB_BP diagnostic bp covered to call a foreign species
pub const CALL_LCB_BP: usize = 500; //diagnostic bp an LCB needs covered to count towards DEFAULT_CALL_LCBS
pub const DIVERGENT_IDENTITY: f64 = 95.0; //contigs whose best core match is below this are from a relative with no core
pub const CONTAMINANTS_CSV: &str = "contaminants.csv"; //every foreign species considered, per genome
