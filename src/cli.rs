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
    #[command(hide = true)]
    Download(DownloadArgs), //Download --> fetch species cores / the whole kore db (not implemented yet, hidden)
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
#[clap(about="Build a kore db of species core genomes from high-quality genomes", arg_required_else_help = true)]
pub struct BuildArgs {

    //REFERENCE INPUT
    #[clap(num_args=1.., short='g', long = "genomes", requires = "species", help_heading = "REFERENCE INPUT", help = "Genomes of a single species (fasta/gzip), with --species")]
    pub genomes: Vec<String>,

    #[clap(long="file-input", requires = "species", help_heading="REFERENCE INPUT", help = "Path to .txt file containing paths to each genome, one line per file, with --species")]
    pub file_input: Option<String>,

    #[clap(short='s', long="species", help_heading="REFERENCE INPUT", help = "GTDB species of the input genomes (e.g. 'Escherichia coli')")]
    pub species: Option<String>,

    #[clap(long="species-table", conflicts_with_all = ["species", "genomes", "file_input"], help_heading="REFERENCE INPUT", help = "Tab-separated file of genome path and GTDB species, to build many species at once (alternative to -g and -s together)")]
    pub species_table: Option<String>,

    #[clap(long="min-genomes", default_value_t = DEFAULT_MIN_GENOMES, help_heading="REFERENCE INPUT", help = "Minimum genomes for a species in --species-table to be built")]
    pub min_genomes: usize,

    #[clap(long="allow-unlisted", help_heading="REFERENCE INPUT", help = "Allow species not in GTDB (built without masking to other species)")]
    pub allow_unlisted: bool,

    //DATABASE
    #[clap(long="cores-at-once", default_value_t = DEFAULT_CORES_AT_ONCE, help_heading="DATABASE", help = "Number of species built in parallel (threads are split between them)")]
    pub cores_at_once: usize,

    #[clap(long="skip-aggregate", help_heading="DATABASE", help = "Only build the species cores (aggregate later with --aggregate-only)")]
    pub skip_aggregate: bool,

    #[clap(long="aggregate-only", conflicts_with_all = ["species", "species_table", "genomes", "file_input"], help_heading="DATABASE", help = "Only aggregate the species cores already in the db")]
    pub aggregate_only: bool,

    //hidden: mask tuning (the paper db uses the defaults)
    #[clap(long="mask-genomes", default_value_t = DEFAULT_MASK_GENOMES, hide = true, help="Max genomes per species used for its masks")]
    pub mask_genomes: usize,

    #[clap(long="mask-identity", default_value_t = MASK_IDENTITY, hide = true, help="Min identity (%) for relative-core sequence to be masked")]
    pub mask_identity: f64,

    #[clap(long="species-dir", hide = true)]
    pub species_dir: Option<String>, //internal: build one species core in this folder (genomes.txt), no aggregate

    //CORE GENOME
    #[clap(short='n', long="n-genomes", default_value_t = DEFAULT_N_GENOMES, help_heading="CORE GENOME", help="Number of genomes in each parsnp2 replicate")]
    pub n_genomes: usize,

    #[clap(short='r', long="replicates", default_value_t = DEFAULT_REPLICATES, help_heading="CORE GENOME", help="Number of parsnp2 replicates (raised if needed so every cluster is sampled)")]
    pub replicates: usize,

    #[clap(long="cluster-ani", default_value_t = DEFAULT_CLUSTER_ANI, help_heading="CORE GENOME", help="ANI (%) used to cluster genomes before sampling")]
    pub cluster_ani: f64,

    #[clap(long="seed", default_value_t = DEFAULT_SEED, help_heading="CORE GENOME", help="Random seed for genome sampling")]
    pub seed: u64,

    //hidden: alternative sampling, kept for testing
    #[clap(long="sampling", value_enum, default_value_t = Sampling::Cluster, hide = true, help="Genome sampling: cluster (one genome per ANI cluster in turn) or random")]
    pub sampling: Sampling,

    #[clap(long="replicate-grouping", value_enum, default_value_t = Grouping::Dealt, hide = true, help="Lineages split across replicates: dealt (equally diverse) or clade (close relatives together)")]
    pub replicate_grouping: Grouping,

    #[clap(long="skani-ani", conflicts_with = "species_table", hide = true, help="Existing skani triangle --sparse output to cluster from (one species only)")]
    pub skani_ani: Option<String>,

    //LCB FILTERS
    #[clap(long="min-length", default_value_t = DEFAULT_MIN_LCB_LENGTH, help_heading="LCB FILTERS", help="Minimum LCB length (bp)")]
    pub min_length: usize,

    #[clap(long="max-diversity", default_value_t = DEFAULT_MAX_LCB_DIVERSITY, help_heading="LCB FILTERS", help="Maximum fraction of variable columns in an LCB")]
    pub max_diversity: f64,

    #[clap(long="replicate-identity", default_value_t = DEFAULT_REPLICATE_IDENTITY, help_heading="LCB FILTERS", help="Minimum identity (%) to match an LCB across replicates")]
    pub replicate_identity: f64,

    //OUTPUT
    #[clap(short='o', long="output", default_value = DEFAULT_DB, help_heading="OUTPUT", help = "kore db folder")]
    pub output: String,

    #[clap(short='f', long="force", help_heading="OUTPUT", help = "Add to an existing kore db (species already built are kept)")]
    pub force: bool,
    
    #[clap(long="keep-intermediate", help_heading="OUTPUT", help="Keep intermediate files (parsnp2, consensus and blast output)")]
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
    #[clap(num_args=1.., short='s', long="species", required_unless_present = "all", conflicts_with = "all", help = "Species to download")]
    pub species: Vec<String>,

    #[clap(long="all", help = "Download the whole kore db")]
    pub all: bool,

    #[clap(short='d', long="db", help = "Folder to download into (default: $KORE_DB)")]
    pub db: Option<String>,
}

#[derive(Args)]
#[clap(about="Estimate completeness and contamination of genome assemblies with a kore db", arg_required_else_help = true)]
pub struct EstimateArgs {

    //ASSEMBLY INPUT
    #[clap(num_args=1.., short='g', long = "genomes", help_heading = "ASSEMBLY INPUT", help = "Assemblies to assess (fasta/gzip)")]
    pub genomes: Vec<String>,

    #[clap(long="file-input", help_heading="ASSEMBLY INPUT", help = "Path to .txt file containing paths to each assembly, one line per file")]
    pub file_input: Option<String>,

    #[clap(long="species-table", conflicts_with_all = ["genomes", "file_input", "species", "core"], help_heading="ASSEMBLY INPUT", help = "Tab-separated file of assembly path and species (skips species detection)")]
    pub species_table: Option<String>,

    //DATABASE
    #[clap(short='d', long="db", help_heading="DATABASE", help = "kore db folder (default: $KORE_DB)")]
    pub db: Option<String>,

    #[clap(short='s', long="species", conflicts_with = "core", help_heading="DATABASE", help = "Species of every assembly (skips species detection)")]
    pub species: Option<String>,

    //hidden: estimate against a single core without a db (completeness and redundancy only)
    #[clap(short='c', long="core", conflicts_with = "db", hide = true, help = "Use a single core from kore build instead of a db (no foreign contamination)")]
    pub core: Option<String>,

    //CONTAMINATION
    #[clap(long="call-breadth", default_value_t = DEFAULT_CALL_BREADTH, help_heading="CONTAMINATION", help="Minimum % of a foreign species' core found to call it")]
    pub call_breadth: f64,

    #[clap(long="call-lcbs", default_value_t = DEFAULT_CALL_LCBS, help_heading="CONTAMINATION", help="Minimum LCBs of a foreign species found to call it")]
    pub call_lcbs: usize,

    //OUTPUT
    #[clap(short='o', long="output", default_value = DEFAULT_OUT_FOLDER, help_heading="OUTPUT", help="Folder to output all resulting files")]
    pub output: String,

    #[clap(long="alignments", help_heading="OUTPUT", help="Also output alignments.tsv of every core alignment (~100 KB per assembly)")]
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
