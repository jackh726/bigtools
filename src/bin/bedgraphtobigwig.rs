#![feature(async_await)]

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

use futures::future::FutureExt;

use bigwig2::bigwig::{BigWigRead, BigWigWrite};
use bigwig2::bedgraphparser::{self, BedGraphParser};
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
    let bedgraphpath = args.next().expect("Must pass a bedgraph.");
    let chrom_map = args.next().expect("Must pass a chromosome sizes file.");
    let bigwigpath = args.next().expect("Must pass a bigwig.");
    println!("Args: {} {} {}", bedgraphpath, chrom_map, bigwigpath);


    let outb = BigWigWrite::create_file(bigwigpath)?;
    let chrom_map = BufReader::new(File::open(chrom_map)?)
        .lines()
        .filter(|l| match l { Ok(s) => !s.is_empty(), _ => true })
        .map(|l| {
            let words = l.expect("Split error");
            let mut split = words.split_whitespace();
            (split.next().expect("Missing chrom").to_owned(), split.next().expect("Missing size").parse::<u32>().unwrap())
        })
        .collect();

    let infile = File::open(bedgraphpath)?;
    let vals_iter = BedGraphParser::<BufReader<File>>::new(infile);
    let chsi = bedgraphparser::get_chromgroupstreamingiterator(vals_iter, outb.options.clone());
    outb.write_groups(chrom_map, chsi)?;

    Ok(())
}
