use std::time::Instant;

pub mod cli;
use cli::*;

use crate::consts::KORE_VERSION;

pub mod consts;
pub mod build;
pub mod cluster;
pub mod estimate;
pub mod lcb;
pub mod replicates;
pub mod util;

fn main() {
    // banner is emitted in both build and estimate
    let start = Instant::now();

    let args = cli::parse_args();
    match args.mode {
        Mode::Build(build_args) => build::build(build_args),
        Mode::Estimate(estimate_args) => estimate::estimate(estimate_args),
    }

    let end = Instant::now();
    log::info!("kore v{} finished in {}s", KORE_VERSION, end.duration_since(start).as_secs_f32());
}
