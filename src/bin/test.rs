use clap::{App, Arg};

use bigtools::bigwig::{BigBedRead, BigBedReadAttachError};
use bigtools::remote_file::*;

fn main() -> Result<(), BigBedReadAttachError> {
    let matches = App::new("Testing")
        .arg(
            Arg::with_name("bigbed")
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

    let f =
        RemoteFile::new("http://users.wenglab.org/hueyj/temp/ENCSR711VWL/hg38/GRCh38-EDGEs.bigBed");
    let mut remote = BigBedRead::from(f)?;
    println!("remote info: {:?}", remote.info);

    let remote_intervals: Vec<_> = remote
        .get_interval("chr1", 12_244_400, 12_248_000)?
        .collect::<Result<_, _>>()?;
    println!("Remote intervals {:?}", remote_intervals.len());
    println!("{:?}", remote_intervals);

    Ok(())
}
