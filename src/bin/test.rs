use std::error::Error;

use clap::{App, Arg};

use bigtools::bbi::BigBedRead;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("Testing")
        .arg(
            Arg::new("bigbed")
                .help("the bigbed to get info for")
                .index(1)
                .required(true),
        )
        .get_matches();

    let bigwigpath = matches.value_of("bigbed").unwrap().to_owned();

    let mut bigwig = BigBedRead::from_file_and_attach(bigwigpath)?;
    println!("info: {:?}", bigwig.info);
    println!("Header: {:?}", bigwig.info.header);

    let intervals: Vec<_> = bigwig
        .get_interval("chr1", 12_244_400, 12_258_000)?
        .collect::<Result<_, _>>()?;
    println!("Intervals {:?}", intervals.len());

    Ok(())
}
