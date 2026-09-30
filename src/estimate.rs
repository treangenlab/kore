use crate::cli::*;
use crate::util::*;
use crate::consts::*;

use log::*;

use std::path::Path;

fn check_args(args: &EstimateArgs) -> Vec<String> {
    // Log file sits inside the output folder and mirrors the terminal.
    let log_path = Path::new(&args.output).join("kore.log");
    init_logging(log_level(args.verbose, args.debug), &log_path);
    print_banner();

    let genomes = collect_genomes(&args.genomes, &args.file_input);

    if !Path::new(&args.db).exists() {
        error!("kore db does not exist: {}", args.db);
        std::process::exit(1);
    }

    if args.threads == 0 {
        error!("--threads must be at least 1");
        std::process::exit(1);
    }

    genomes
}

pub fn estimate(args: EstimateArgs) {
    let genomes = check_args(&args);
    info!("kore v{} estimating completeness for {} assemblies", KORE_VERSION, genomes.len());

    // TODO: load db
    // TODO: estimate completeness for each assembly
    // TODO: write results to output folder
}
