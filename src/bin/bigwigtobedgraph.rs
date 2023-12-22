use std::error::Error;

use bigtools::utils::cli::bigwigtobedgraph::{bigwigtobedgraph, BigWigToBedGraphArgs};
use bigtools::utils::cli::compat_args;
use clap::Parser;

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    args: BigWigToBedGraphArgs,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Cli::parse_from(compat_args(std::env::args_os()));

    bigwigtobedgraph(matches.args)
}

#[test]
fn verify_cli_bigwigtobedgraph() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
