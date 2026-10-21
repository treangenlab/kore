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
    Build(BuildArgs), //Build --> core genome alignment (parsnp2) of a species' genomes + downstream analyses to build a kore db
    Estimate(EstimateArgs), //Estimate --> estimates completeness of assemblies against a kore db
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
#[clap(about="Build a kore database from a set of species reference genomes", arg_required_else_help = true)]
pub struct BuildArgs {

    //REFERENCE INPUT
    #[clap(num_args=1.., short='g', long = "genomes", help_heading = "REFERENCE INPUT", help = "Genome files to build the core genome from (fasta/gzip)")]
    pub genomes: Vec<String>,

    #[clap(long="file-input", help_heading="REFERENCE INPUT", help = "Path to .txt file containing paths to each file, one line per file, each containing a single genome")]
    pub file_input: Option<String>,

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

    #[clap(long="skani-ani", help_heading="CLUSTERING", help="Existing skani triangle --sparse output (plain or .gz) to cluster from instead of running skani. Without it, ANI saved in the output folder by an earlier build of the same genomes is reused")]
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

    //OUTPUT
    #[clap(short='o', long="output", default_value = DEFAULT_BUILD_OUTPUT, help_heading="OUTPUT", help="Output folder (core.fasta, core_lcbs.csv and kore.log)")]
    pub output: String,

    #[clap(long="keep-intermediate", help_heading="OUTPUT", help="Keep intermediate files (genomes, parsnp2, consensus and blast output) in <output>/intermediate")]
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
#[clap(about="Estimate completeness of genome assemblies using a kore database", arg_required_else_help = true)]
pub struct EstimateArgs {

    //ASSEMBLY INPUT
    #[clap(num_args=1.., short='g', long = "genomes", help_heading = "ASSEMBLY INPUT", help = "Assemblies to estimate completeness for (fasta/gzip)")]
    pub genomes: Vec<String>,

    #[clap(long="file-input", help_heading="ASSEMBLY INPUT", help = "Path to .txt file containing paths to each assembly, one line per file")]
    pub file_input: Option<String>,

    //DATABASE INPUT
    #[clap(short='c', long="core", required = true, help_heading="DATABASE INPUT", help = "kore db (core.fasta from kore build)")]
    pub core: String,

    //OUTPUT
    #[clap(short='o', long="output", default_value = DEFAULT_OUT_FOLDER, help_heading="OUTPUT", help="Folder to output all resulting files")]
    pub output: String,

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
