use clap::{Arg, Command};

use bigtools::{BigWigRead, BigWigReadOpenError};

fn main() -> Result<(), BigWigReadOpenError> {
    let matches = Command::new("Verify")
        .about("Verifies different parts of a bigwig or bigbed file. By default, it only verifies the header is correct and that file offsets are valid.")
        .arg(Arg::new("input")
            .help("the bigwig or bigbed to verify")
            .index(1)
            .required(true)
        )
        .arg(Arg::new("index")
            .short('i')
            .help("If set, the entire index will be verified. The file offsets pointed by leaf nodes will not be checked.")
        )
        .arg(Arg::new("data")
            .short('d')
            .help("If set, every data block will be uncompressed to verify their valid. The contents will not be checked. This implies -i.")
        )
        .get_matches();

    let bigwigpath = matches.get_one::<String>("input").unwrap();
    let _index = matches.get_count("index") > 0;
    let _data = matches.get_count("data") > 0;

    let _bigwig = BigWigRead::open_file(bigwigpath)?;

    unimplemented!();
    //Ok(())
}
