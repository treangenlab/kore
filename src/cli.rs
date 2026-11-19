use crate::consts::*;
use clap::{Args, Parser, Subcommand, ValueEnum};
use clap::builder::styling::{AnsiColor, Effects, Styles};


fn custom_styles() -> Styles {
    Styles::styled()
        .header(AnsiColor::Blue.on_default() | Effects::BOLD)
        .usage(AnsiColor::White.on_default() | Effects::BOLD)
        .literal(AnsiColor::White.on_default() | Effects::BOLD)
        .placeholder(AnsiColor::White.on_default())
}

#[derive(Parser)]
#[command(author, version, about, long_about = None, styles=custom_styles())]
#[command(propagate_version = true)]
pub struct Cli {
    #[clap(subcommand)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    Build(BuildArgs), //Build --> species cores (parsnp2 core genome alignment) aggregated into a kore db
    Download(DownloadArgs), //Download --> fetch species cores / the whole kore db (not implemented yet)
    Estimate(EstimateArgs), //Estimate --> completeness and contamination of assemblies against the kore db (or one core)
}

#[derive(Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Grouping {
    Dealt, //lineages dealt across replicates in turn: equally diverse replicates
    Clade, //consecutive lineages along the ANI tree: each replicate is a clade of close relatives
}

#[derive(Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Sampling {
    Cluster, //skani ANI clusters, one genome from each cluster in turn
    Random, //random shuffle of all genomes
}

#[derive(Args)]
#[clap(about="Build a kore db by generating species cores from high quality genomes", arg_required_else_help = true)]
pub struct BuildArgs {

    //INPUT
    #[clap(num_args=1.., short='g', long = "genomes", requires = "species", help_heading = "INPUT", help = "Genomes of one species (fasta/gzip); with --species")]
    pub genomes: Vec<String>,

    #[clap(long="file-input", requires = "species", help_heading="INPUT", help = "Path to .txt file with the paths of one species' genomes, one per line; with --species")]
    pub file_input: Option<String>,

    #[clap(short='s', long="species", help_heading="INPUT", help = "GTDB species of the genomes given with -g/--file-input (e.g. 'Escherichia coli' or escherichia_coli)")]
    pub species: Option<String>,

    #[clap(long="species-table", conflicts_with_all = ["species", "genomes", "file_input"], help_heading="INPUT", help = "Many species at once: tab-separated genome path and GTDB species per line (header optional). Every species with >= --min-genomes genomes gets a core")]
    pub species_table: Option<String>,

    #[clap(long="min-genomes", default_value_t = DEFAULT_MIN_GENOMES, help_heading="INPUT", help = "Species in --species-table with fewer genomes are skipped")]
    pub min_genomes: usize,

    #[clap(long="allow-unlisted", help_heading="INPUT", help = "Build species that are not in the GTDB taxonomy (otherwise an error): they get no relatives, so no masks")]
    pub allow_unlisted: bool,

    //DB
    #[clap(short='o', long="output", default_value = DEFAULT_DB, help_heading="DB", help = "kore db folder. Species already built in it are kept (an existing db can be extended); every core in it is aggregated at the end")]
    pub output: String,

    #[clap(short='f', long="force", help_heading="DB", help = "Build into an existing kore db: species already built from the same genomes are kept, the others are (re)built, and the db is re-aggregated")]
    pub force: bool,

    #[clap(long="cores-at-once", default_value_t = DEFAULT_CORES_AT_ONCE, help_heading="DB", help = "Species cores built at the same time; --threads is split between them")]
    pub cores_at_once: usize,

    #[clap(long="aggregate-only", conflicts_with_all = ["species", "species_table", "genomes", "file_input"], help_heading="DB", help = "Only (re)aggregate the cores already in the db (e.g. after building species in separate jobs with --skip-aggregate)")]
    pub aggregate_only: bool,

    #[clap(long="skip-aggregate", help_heading="DB", help = "Build the cores only (e.g. one job per species on a cluster); aggregate later with --aggregate-only")]
    pub skip_aggregate: bool,

    #[clap(long="mask-genomes", default_value_t = DEFAULT_MASK_GENOMES, help_heading="DB", help="Max build genomes per species used to find the parts of its relatives' cores it carries (seeded sample)")]
    pub mask_genomes: usize,

    #[clap(long="mask-identity", default_value_t = MASK_IDENTITY, help_heading="DB", help="Relative-core sequence found in a species' build genomes at >= this identity (%) is masked for that species")]
    pub mask_identity: f64,

    #[clap(long="species-dir", hide = true)]
    pub species_dir: Option<String>, //internal: build one species core in this folder (genomes.txt), no aggregate

    //PARSNP
    #[clap(short='n', long="n-genomes", default_value_t = DEFAULT_N_GENOMES, help_heading="PARSNP", help="Number of genomes used to build the core genome in each replicate")]
    pub n_genomes: usize,

    #[clap(short='r', long="replicates", default_value_t = DEFAULT_REPLICATES, help_heading="PARSNP", help="Number of parsnp2 replicates (disjoint genome sets when possible). With cluster sampling, raised to clusters / n-genomes (rounded up) so every cluster is sampled; replicates get equally diverse shares, topped up to n-genomes with reused genomes if needed")]
    pub replicates: usize,

    #[clap(long="seed", default_value_t = DEFAULT_SEED, help_heading="PARSNP", help="Random seed used to order genomes (clusters and genomes within them when clustering) before splitting into replicates")]
    pub seed: u64,

    //CLUSTERING
    #[clap(long="sampling", value_enum, default_value_t = Sampling::Cluster, help_heading="CLUSTERING", help="How genomes are chosen for the replicates: cluster = skani ANI clusters, one genome from each cluster in turn until they run out; random = random shuffle")]
    pub sampling: Sampling,

    #[clap(long="cluster-ani", default_value_t = DEFAULT_CLUSTER_ANI, help_heading="CLUSTERING", help="ANI (%) at which genomes are clustered (average linkage of skani ANI) before sampling")]
    pub cluster_ani: f64,

    #[clap(long="skani-ani", conflicts_with = "species_table", help_heading="CLUSTERING", help="Existing skani triangle --sparse output (plain or .gz) to cluster from instead of running skani (one species only). Without it, ANI saved in the output folder by an earlier build of the same genomes is reused")]
    pub skani_ani: Option<String>,

    #[clap(long="replicate-grouping", value_enum, default_value_t = Grouping::Dealt, help_heading="CLUSTERING", help="How lineages are split across replicates when there are more than --n-genomes: dealt = in turn, so replicates are equally diverse; clade = consecutive lineages along the ANI tree, so each replicate aligns close relatives")]
    pub replicate_grouping: Grouping,

    //LCB FILTERS
    #[clap(long="min-length", default_value_t = DEFAULT_MIN_LCB_LENGTH, help_heading="LCB FILTERS", help="Minimum consensus length of an LCB")]
    pub min_length: usize,

    #[clap(long="max-diversity", default_value_t = DEFAULT_MAX_LCB_DIVERSITY, help_heading="LCB FILTERS", help="Maximum fraction of alignment columns in an LCB where any genome differs (gaps included); raised to the 95th percentile of LCB diversity when that is higher")]
    pub max_diversity: f64,

    #[clap(long="replicate-identity", default_value_t = DEFAULT_REPLICATE_IDENTITY, help_heading="LCB FILTERS", help="Minimum blast identity (%) for an LCB to match LCBs in another replicate")]
    pub replicate_identity: f64,

    #[clap(long="keep-intermediate", help_heading="DB", help="Keep intermediate files (genomes, parsnp2, consensus and blast output) in each species' intermediate/ folder")]
    pub keep_intermediate: bool,

    // OTHER PARAMETERS
    //Number of threads
    #[clap(short, long="threads", default_value_t=DEFAULT_THREADS, help="Number of threads")]
    pub threads: usize,

    //Debug mode
    #[clap(long = "debug", help = "Debug output")]
    pub debug: bool,

    //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,
}

#[derive(Args)]
#[clap(about="Download species cores or the whole kore db (not available yet)", arg_required_else_help = true)]
pub struct DownloadArgs {
    #[clap(num_args=1.., short='s', long="species", required_unless_present = "all", conflicts_with = "all", help = "Species to download (GTDB names, e.g. escherichia_coli)")]
    pub species: Vec<String>,

    #[clap(long="all", help = "Download the whole kore db (every species core, the aggregate index and masks)")]
    pub all: bool,

    #[clap(short='d', long="db", help = "Folder to download into (default: $KORE_DB)")]
    pub db: Option<String>,
}

#[derive(Args)]
#[clap(about="Estimate completeness and contamination of assemblies with a kore db", arg_required_else_help = true)]
pub struct EstimateArgs {

    //ASSEMBLY INPUT
    #[clap(num_args=1.., short='g', long = "genomes", help_heading = "ASSEMBLY INPUT", help = "Assemblies to assess (fasta/gzip)")]
    pub genomes: Vec<String>,

    #[clap(long="file-input", help_heading="ASSEMBLY INPUT", help = "Path to .txt file containing paths to each assembly, one line per file")]
    pub file_input: Option<String>,

    //DATABASE
    #[clap(short='d', long="db", help_heading="DATABASE", help = "kore db folder (default: $KORE_DB). By default every species core is used: each assembly's species is detected and contamination is estimated against all cores")]
    pub db: Option<String>,

    #[clap(short='s', long="species", conflicts_with = "core", help_heading="DATABASE", help = "Use only this species' core from the db (no species detection; contamination from other species is not assessed, only same-species redundancy)")]
    pub species: Option<String>,

    #[clap(short='c', long="core", conflicts_with = "db", help_heading="DATABASE", help = "Use your own species core (core.fasta from kore build) instead of the db; as --species")]
    pub core: Option<String>,

    //CONTAMINATION
    #[clap(long="call-breadth", default_value_t = DEFAULT_CALL_BREADTH, help_heading="CONTAMINATION", help="Min % of a foreign species' diagnostic core covered to call it present (whole db only)")]
    pub call_breadth: f64,

    #[clap(long="call-lcbs", default_value_t = DEFAULT_CALL_LCBS, help_heading="CONTAMINATION", help="Min LCBs of a foreign species with diagnostic sequence covered to call it present (whole db only)")]
    pub call_lcbs: usize,

    //OUTPUT
    #[clap(short='o', long="output", default_value = DEFAULT_OUT_FOLDER, help_heading="OUTPUT", help="Output folder (report.csv, contaminants.csv, kore.log; alignments.csv with --alignments)")]
    pub output: String,
    #[clap(long="alignments", help_heading="OUTPUT", help="Also write alignments.csv: every core LCB alignment in every assembly (large: ~100 KB per assembly)")]
    pub alignments: bool,

    // OTHER PARAMETERS
    //Number of threads
    #[clap(short, long="threads", default_value_t=DEFAULT_THREADS, help="Number of threads")]
    pub threads: usize,

    //Debug mode
    #[clap(long = "debug", help = "Debug output")]
    pub debug: bool,

    //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,
}


pub fn parse_args() -> Cli {
    Cli::parse()
}
