use std::fs::File;
use std::io::{self, Write};

use clap::{App, Arg};

use futures::task::SpawnExt;

use bigtools::bigwig::{BBIRead, BigWigRead, BigWigReadAttachError, ChromAndSize};
use bigtools::seekableread::{Reopen, SeekableRead};
use bigtools::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};

pub fn write_bg<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
    bigwig: BigWigRead<R, S>,
    mut out_file: File,
    nthreads: usize,
) -> std::io::Result<()> {
    /*
    // This is the simple single-threaded approach
    let mut chroms: Vec<ChromAndSize> = bigwig.get_chroms();
    chroms.sort_by(|a, b| a.name.cmp(&b.name));
    let mut writer = io::BufWriter::new(out_file);
    for chrom in chroms {
        let mut values = bigwig.get_interval(&chrom.name, 0, chrom.length)?;
        while let Some(raw_val) = values.next() {
            let val = raw_val?;
            writer.write_fmt(format_args!("{}\t{}\t{}\t{}\n", chrom.name, val.start, val.end, val.value))?;
        }
    }
    */

    let pool = futures::executor::ThreadPoolBuilder::new()
        .pool_size(nthreads)
        .create()
        .expect("Unable to create thread pool.");

    let chrom_files: Vec<io::Result<(_, TempFileBuffer<File>)>> = bigwig
        .get_chroms()
        .into_iter()
        .map(|chrom| {
            let bigwig = bigwig.clone();
            let (buf, file): (TempFileBuffer<File>, TempFileBufferWriter<File>) =
                TempFileBuffer::new()?;
            let writer = io::BufWriter::new(file);
            async fn file_future<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
                mut bigwig: BigWigRead<R, S>,
                chrom: ChromAndSize,
                mut writer: io::BufWriter<TempFileBufferWriter<File>>,
            ) -> io::Result<()> {
                for raw_val in bigwig.get_interval(&chrom.name, 0, chrom.length)? {
                    let val = raw_val?;
                    writer.write_fmt(format_args!(
                        "{}\t{}\t{}\t{}\n",
                        chrom.name, val.start, val.end, val.value
                    ))?;
                }
                Ok(())
            };
            let handle = pool
                .spawn_with_handle(file_future(bigwig, chrom, writer))
                .expect("Couldn't spawn.");
            Ok((handle, buf))
        })
        .collect::<Vec<_>>();

    for res in chrom_files {
        let (f, mut buf) = res.unwrap();
        buf.switch(out_file);
        futures::executor::block_on(f).unwrap();
        out_file = buf.await_real_file();
    }

    Ok(())
}

fn main() -> Result<(), BigWigReadAttachError> {
    let matches = App::new("BigWigToBedGraph")
        .about("Converts an input bigWig to a bedGraph. Can be multi-threaded for substantial speedups. Note for roughly each core, one temporary file will be opened.")
        .arg(Arg::with_name("bigwig")
            .help("the bigwig to get convert to bedgraph")
            .index(1)
            .required(true)
        )
        .arg(Arg::with_name("bedgraph")
            .help("the path of the bedgraph to output to")
            .index(2)
            .required(true)
        )
        .arg(Arg::with_name("nthreads")
            .short("t")
            .help("Set the number of threads to use. This tool will nearly always benefit from more cores (<= # chroms). Note: for parts of the runtime, the actual usage may be nthreads+1")
            .takes_value(true)
            .default_value("6"))
        .get_matches();

    let bigwigpath = matches.value_of("bigwig").unwrap();
    let bedgraphpath = matches.value_of("bedgraph").unwrap();

    let nthreads = {
        let nthreads = matches.value_of("nthreads").unwrap();
        let parsed = nthreads.parse();
        if parsed.is_err() {
            eprintln!("Invalid argument for `nthreads`: must be a positive number");
            return Ok(());
        }
        parsed.unwrap()
    };

    let bigwig = BigWigRead::from_file_and_attach(bigwigpath)?;
    let bedgraph = File::create(bedgraphpath)?;

    write_bg(bigwig, bedgraph, nthreads)?;

    Ok(())
}
