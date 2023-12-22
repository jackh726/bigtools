use bigtools::utils::cli::bigwiginfo::{bigwiginfo, BigWigInfoArgs};
use bigtools::utils::cli::compat_args;
use clap::Parser;

use std::error::Error;

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    args: BigWigInfoArgs,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Cli::parse_from(compat_args(std::env::args_os()));

    bigwiginfo(matches.args)
}
