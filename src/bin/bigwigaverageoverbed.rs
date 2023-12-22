use std::error::Error;

use bigtools::utils::cli::bigwigaverageoverbed::{bigwigaverageoverbed, BigWigAverageOverBedArgs};
use bigtools::utils::cli::compat_args;
use clap::Parser;

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    args: BigWigAverageOverBedArgs,
}

fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
    let matches = Cli::parse_from(compat_args(std::env::args_os()));

    bigwigaverageoverbed(matches.args)
}

#[test]
fn verify_cli_bigwigaverageoverbed() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
