use std::error::Error;

use bigtools::utils::cli::bedtobigbed::{bedtobigbed, BedToBigBedArgs};
use bigtools::utils::cli::compat_args;
use clap::Parser;

#[derive(Parser)]
#[command(about = "Converts a bed to a bigBed.", long_about = None)]
struct Cli {
    #[command(flatten)]
    args: BedToBigBedArgs,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Cli::parse_from(compat_args(std::env::args_os()));

    bedtobigbed(matches.args)
}

#[test]
fn verify_cli_bedtobigbed() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
