use std::io;

use bigwig2::bigwig::BigWigRead;

fn main() -> io::Result<()> {
    let mut args = std::env::args();
    args.next();
    let bigwigpath = args.next().expect("Must pass a bigwig.");
    println!("Args: {:}", bigwigpath);

    let bigwig = BigWigRead::from_file_and_attach(bigwigpath)?;
    println!("Header: {:?}", bigwig.info.header);

    Ok(())
}
