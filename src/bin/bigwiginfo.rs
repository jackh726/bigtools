use clap::{Arg, Command};

use bigtools::{BBIRead, BigWigRead, BigWigReadAttachError};

fn main() -> Result<(), BigWigReadAttachError> {
    let matches = Command::new("BigWigInfo")
        .arg(
            Arg::new("bigwig")
                .help("the bigwig to get info for")
                .index(1)
                .required(true),
        )
        .get_matches();

    let bigwigpath = matches.get_one::<String>("bigwig").unwrap();

    #[cfg(feature = "remote")]
    {
        if bigwigpath.starts_with("http") {
            use bigtools::utils::remote_file::RemoteFile;
            let f = RemoteFile::new(bigwigpath);
            let mut bigwig = BigWigRead::open(f)?;
            println!("Header: {:#?}", bigwig.get_info().header);
            println!("Summary: {:#?}", bigwig.get_summary());
            println!("Zooms: {:#?}", bigwig.get_info().zoom_headers);
        } else {
            let mut bigwig = BigWigRead::open_file(bigwigpath)?;
            println!("Header: {:#?}", bigwig.get_info().header);
            println!("Summary: {:#?}", bigwig.get_summary());
            println!("Zooms: {:#?}", bigwig.get_info().zoom_headers);
        }
    }
    #[cfg(not(feature = "remote"))]
    {
        let mut bigwig = BigWigRead::open_file(bigwigpath)?;
        println!("Header: {:#?}", bigwig.get_info().header);
        println!("Summary: {:#?}", bigwig.get_summary());
        println!("Zooms: {:#?}", bigwig.get_info().zoom_headers);
    }

    Ok(())
}
