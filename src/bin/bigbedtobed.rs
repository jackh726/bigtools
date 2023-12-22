use std::error::Error;

use bigtools::utils::cli::bigbedtobed::{bigbedtobed, BigBedToBedArgs};
use bigtools::utils::cli::compat_args;
use clap::Parser;

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    args: BigBedToBedArgs,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Cli::parse_from(compat_args(std::env::args_os()));
    bigbedtobed(matches.args)
}

#[test]
fn verify_cli_bigbedtobed() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
