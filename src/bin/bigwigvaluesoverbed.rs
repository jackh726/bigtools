use std::error::Error;

use bigtools::utils::cli::bigwigvaluesoverbed::{bigwigvaluesoverbed, BigWigValuesOverBedArgs};
use clap::Parser;

#[derive(Parser)]
#[command(about = "Gets statistics of a bigWig over a bed.", long_about = None)]
struct Cli {
    #[command(flatten)]
    args: BigWigValuesOverBedArgs,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Cli::parse();

    bigwigvaluesoverbed(matches.args)
}
