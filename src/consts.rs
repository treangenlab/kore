pub const KORE_VERSION: &str = env!("CARGO_PKG_VERSION");
pub const DEFAULT_THREADS: usize = 30; //number of threads
pub const PARSNP_BIN: &str = "parsnp"; //parsnp2 executable (expected on PATH via the conda env)
pub const SKANI_BIN: &str = "skani"; //skani executable for ANI clustering (expected on PATH via the conda env)
pub const DEFAULT_CLUSTER_ANI: f64 = 99.5; //ANI (%) genomes are clustered at before sampling
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
pub const COMPLETENESS_CSV: &str = "completeness.csv"; //per genome completeness estimates
pub const ALIGNMENTS_CSV: &str = "alignments.csv"; //every LCB alignment in every genome
pub const MM_SECONDARY_RATIO: f32 = 0.5; //minimap2 -p, looser than default so second copies are still reported
pub const MM_MAX_SECONDARY: i32 = 20; //minimap2 -N
pub const MM_MIN_CHAIN_SCORE: i32 = 30; //minimap2 -m, asm20 default drops short LCBs
pub const MM_MIN_DP_SCORE: i32 = 40; //minimap2 -s, asm20 default drops short LCBs
pub const MIN_ALN_IDENTITY: f64 = 80.0; //min identity (%) of an alignment to count towards coverage
pub const MIN_ALN_LEN: i32 = 100; //min alignment block length to count towards coverage
pub const LCB_PRESENT_COV: f64 = 80.0; //min coverage (%) of an LCB for it to count as present
pub const LCB_PRESENT_IDENTITY: f64 = 90.0; //min identity (%) of alignments counted for LCB presence

//PANEL (kore panel: every species' core in one database, for species detection and contamination)
pub const DEFAULT_PANEL_OUTPUT: &str = "kore_panel"; //folder for panel outputs
pub const PANEL_FASTA: &str = "panel.fasta"; //every core's LCBs, named <species>|<lcb>
pub const PANEL_INDEX: &str = "panel.mmi"; //minimap2 index of PANEL_FASTA
pub const PANEL_SPECIES: &str = "species.tsv"; //species, group, LCBs, core bp, genome size
pub const PANEL_CORES_DIR: &str = "cores"; //<species>.fasta: each species' core
pub const PANEL_MASKS_DIR: &str = "masks"; //<species>.tsv: relatives' core intervals the species' own genomes carry
pub const PANEL_SEP: char = '|'; //separates species and LCB in panel sequence names
pub const PANEL_SECONDARY_RATIO: f32 = 0.2; //minimap2 -p against the panel: every relative's core must still be reported
pub const PANEL_MAX_SECONDARY: i32 = 100; //minimap2 -N against the panel
pub const DEFAULT_MASK_GENOMES: usize = 100; //max build genomes per species used for its mask (seeded sample)
pub const MASK_IDENTITY: f64 = 97.0; //relative-core sequence found at >= this identity in a species' own genomes is masked for it
pub const MIN_DIAGNOSTIC_BP: usize = 10_000; //a relative with less diagnostic core than this cannot be told apart from the target

//CONTAMINATION (kore estimate --panel)
pub const SPECIES_IDENTITY: f64 = 95.0; //alignments >= this identity count towards a species' core coverage (panel scan)
pub const SPECIES_CANDIDATES: usize = 5; //species with the highest panel-scan coverage verified core -> assembly
pub const SPECIES_CANDIDATE_FRACTION: f64 = 0.5; //...plus any with >= this fraction of the best scan coverage (max 2x)
pub const CONTIG_MIN_ALIGNED: usize = 1_000; //a contig is assigned to a species only on >= this many aligned bp
pub const FOREIGN_IDENTITY: f64 = 98.0; //alignments >= this identity count towards a foreign species' diagnostic core
pub const DEFAULT_CALL_BREADTH: f64 = 1.0; //min % of a foreign species' diagnostic core covered to call it present
pub const DEFAULT_CALL_LCBS: usize = 10; //min LCBs with >= CALL_LCB_BP diagnostic bp covered to call a foreign species
pub const CALL_LCB_BP: usize = 500; //diagnostic bp an LCB needs covered to count towards DEFAULT_CALL_LCBS
pub const DIVERGENT_IDENTITY: f64 = 95.0; //contigs whose best core match is below this are from a relative with no core
pub const CONTAMINATION_CSV: &str = "contamination.csv"; //per genome contamination estimates
pub const CONTAMINANTS_CSV: &str = "contaminants.csv"; //every foreign species considered, per genome
