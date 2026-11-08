//kore download: fetch species cores or the whole kore db. Not implemented yet: the hosting location of the db is
//still to be decided. The command and its arguments exist so the interface is fixed.
use crate::cli::DownloadArgs;
use crate::consts::*;

pub fn download(args: DownloadArgs) {
    let target = args.db.or_else(|| std::env::var(DB_ENV).ok()).unwrap_or_else(|| "(no --db and no $KORE_DB)".into());
    let what = if args.all { "the whole kore db".to_string() } else { args.species.join(", ") };
    eprintln!("kore download is not available yet (requested: {} into {}).", what, target);
    std::process::exit(2);
}
