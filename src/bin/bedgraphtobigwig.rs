use bigtools::utils::cli::{
    bedgraphtobigwig::{bedgraphtobigwig, BedGraphToBigWigArgs},
    compat_args,
};
use clap::Parser;
use std::error::Error;

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    args: BedGraphToBigWigArgs,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Cli::parse_from(compat_args(std::env::args_os()));

    bedgraphtobigwig(matches.args)
}
