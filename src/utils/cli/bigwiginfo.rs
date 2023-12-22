use std::error::Error;

use clap::Parser;

use crate::BigWigRead;

#[derive(Debug, Parser)]
#[command(
    name = "bigwiginfo",
    about = "Converts a bed to a bigBed.",
    long_about = None,
)]
pub struct BigWigInfoArgs {
    /// The bigwig to get info for.
    bigwig: String,
}

pub fn bigwiginfo(args: BigWigInfoArgs) -> Result<(), Box<dyn Error>> {
    let bigwigpath = args.bigwig;

    #[cfg(feature = "remote")]
    {
        if bigwigpath.starts_with("http") {
            use crate::utils::remote_file::RemoteFile;
            let f = RemoteFile::new(&bigwigpath);
            let mut bigwig = BigWigRead::open(f)?;
            println!("Header: {:#?}", bigwig.info().header);
            println!("Summary: {:#?}", bigwig.get_summary());
            println!("Zooms: {:#?}", bigwig.info().zoom_headers);
        } else {
            let mut bigwig = BigWigRead::open_file(&bigwigpath)?;
            println!("Header: {:#?}", bigwig.info().header);
            println!("Summary: {:#?}", bigwig.get_summary());
            println!("Zooms: {:#?}", bigwig.info().zoom_headers);
        }
    }
    #[cfg(not(feature = "remote"))]
    {
        let mut bigwig = BigWigRead::open_file(bigwigpath)?;
        println!("Header: {:#?}", bigwig.info().header);
        println!("Summary: {:#?}", bigwig.get_summary());
        println!("Zooms: {:#?}", bigwig.info().zoom_headers);
    }

    Ok(())
}
