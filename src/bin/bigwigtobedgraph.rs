#![feature(async_await)]

use std::fs::File;
use std::io::{self, Write};

use futures::future::FutureExt;

use bigwig2::bigwig::BigWigRead;
use bigwig2::tempfilebuffer::TempFileBuffer;

pub fn write_bg(bigwig: BigWigRead, mut out_file: File) -> std::io::Result<()> {
    let chrom_files: Vec<io::Result<(_, TempFileBuffer)>> = bigwig.get_chroms().into_iter().map(|chrom| {
        let bigwig = bigwig.clone();
        let (buf, file) = TempFileBuffer::new()?;
        let mut writer = std::io::BufWriter::new(file);
        let file_future = async move || -> io::Result<()> {
            for raw_val in bigwig.get_interval(&chrom.name, 1, chrom.length)? {
                let val = raw_val?;
                writer.write_fmt(format_args!("{}\t{}\t{}\t{}\n", chrom.name, val.start, val.end, val.value))?;
            }
            Ok(())
        };
        let (remote, handle) = file_future().remote_handle();
        std::thread::spawn(move || {
            futures::executor::block_on(remote);
        });
        Ok((handle,buf))
    }).collect::<Vec<_>>();

    for res in chrom_files {
        let (f, mut buf) = res.unwrap();
        buf.switch(out_file).unwrap();
        futures::executor::block_on(f).unwrap();
        out_file = buf.await_file();
    }

    Ok(())
}

fn main() -> io::Result<()> {
    let mut args = std::env::args();
    args.next();
    let bigwigpath = args.next().expect("Must pass a bigwig.");
    let bedgraphpath = args.next().expect("Must pass a bedgraph.");
    println!("Args: {:} {:}", bigwigpath, bedgraphpath);

    let bigwig = BigWigRead::from_file_and_attach(bigwigpath)?;
    let bedgraph = File::create(bedgraphpath)?;

    write_bg(bigwig, bedgraph)?;

    Ok(())
}
