use clap::{Arg, Command};

use bigtools::{BBIFileRead, BigWigRead, BigWigReadOpenError};

fn main() -> Result<(), BigWigReadOpenError> {
    let matches = Command::new("BigWigInfo")
        .arg(
            Arg::new("bigwig")
                .help("the bigwig to get info for")
                .index(1)
                .required(true),
        )
        .get_matches();

    let bigwigpath = matches.get_one::<String>("bigwig").unwrap();

    fn print_info<R: BBIFileRead>(mut bigwig: BigWigRead<R>) {
        println!("Header: {:#?}", bigwig.info().header);
        println!("Summary: {:#?}", bigwig.get_summary());
        println!("Zooms: {:#?}", bigwig.info().zoom_headers);
    }

    #[cfg(feature = "remote")]
    {
        if bigwigpath.starts_with("http") {
            use bigtools::utils::remote_file::RemoteFile;
            let f = RemoteFile::new(bigwigpath);
            let bigwig = BigWigRead::open(f)?;
            print_info(bigwig);
            return Ok(());
        }
    }

    let bigwig = BigWigRead::open_file(bigwigpath)?;
    print_info(bigwig);

    Ok(())
}
