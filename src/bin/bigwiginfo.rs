use clap::{App, Arg};

use bigtools::bigwig::{BigWigRead, BigWigReadAttachError};

fn main() -> Result<(), BigWigReadAttachError> {
    let matches = App::new("BigWigInfo")
        .arg(
            Arg::new("bigwig")
                .about("the bigwig to get info for")
                .index(1)
                .required(true),
        )
        .get_matches();

    let bigwigpath = matches.value_of("bigwig").unwrap();

    #[cfg(feature = "remote")]
    {
        if bigwigpath.starts_with("http") {
            use bigtools::utils::remote_file::RemoteFile;
            let f = RemoteFile::new(bigwigpath);
            let mut bigwig = BigWigRead::from(f)?;
            println!("Header: {:?}", bigwig.info.header);
            println!("Summary: {:?}", bigwig.get_summary());
        } else {
            let mut bigwig = BigWigRead::from_file_and_attach(bigwigpath)?;
            println!("Header: {:?}", bigwig.info.header);
            println!("Summary: {:?}", bigwig.get_summary());
        }
    }
    #[cfg(not(feature = "remote"))]
    {
        let mut bigwig = BigWigRead::from_file_and_attach(bigwigpath)?;
        println!("Header: {:?}", bigwig.info.header);
        println!("Summary: {:?}", bigwig.get_summary());
    }

    Ok(())
}
