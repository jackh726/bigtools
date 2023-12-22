use std::error::Error;

use bigtools::utils::cli::bigwigmerge::{bigwigmerge, BigWigMergeArgs};
use bigtools::utils::cli::compat_args;
use clap::Parser;

#[derive(Parser)]
#[command(about = "Merges multiple bigwigs.", long_about = None)]
struct Cli {
    #[command(flatten)]
    args: BigWigMergeArgs,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Cli::parse_from(compat_args(std::env::args_os()));

    bigwigmerge(matches.args)
}

#[test]
fn verify_cli_bigwigmerge() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
