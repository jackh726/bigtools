use std::error::Error;

use clap::Parser;

use crate::{BBIFileRead, BigWigRead};

#[derive(Clone, Debug, PartialEq, Parser)]
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

    fn print_info<R: BBIFileRead>(mut bigwig: BigWigRead<R>) {
        println!("Header: {:#?}", bigwig.info().header);
        println!("Summary: {:#?}", bigwig.get_summary());
        println!("Zooms: {:#?}", bigwig.info().zoom_headers);
    }

    #[cfg(feature = "remote")]
    {
        if bigwigpath.starts_with("http") {
            use crate::utils::remote_file::RemoteFile;
            let f = RemoteFile::new(&bigwigpath);
            let bigwig = BigWigRead::open(f)?;
            print_info(bigwig);
            return Ok(());
        }
    }

    let bigwig = BigWigRead::open_file(&bigwigpath)?;
    print_info(bigwig);

    Ok(())
}
