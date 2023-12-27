use bigtools::utils::cli::bigbedinfo::{bigbedinfo, BigBedInfoArgs};
use bigtools::utils::cli::compat_args;
use clap::Parser;

use std::error::Error;

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    args: BigBedInfoArgs,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Cli::parse_from(compat_args(std::env::args_os()));

    bigbedinfo(matches.args)
}
