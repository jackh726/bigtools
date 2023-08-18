use std::error::Error;
use std::fs::File;
use std::io::{self, Write};

use clap::{Arg, Command};

use futures::task::SpawnExt;

use bigtools::bbi::{BBIRead, BigWigRead, ChromAndSize};
use bigtools::bbiread::BBIReadError;
use bigtools::utils::reopen::{Reopen, SeekableRead};
use bigtools::utils::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};
use ufmt::uwrite;

pub fn write_bg<R: Reopen + SeekableRead + Send + 'static>(
    bigwig: BigWigRead<R>,
    mut out_file: File,
    nthreads: usize,
) -> Result<(), BBIReadError> {
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
            let bigwig = bigwig.reopen()?;
            let (buf, file): (TempFileBuffer<File>, TempFileBufferWriter<File>) =
                TempFileBuffer::new()?;
            let writer = io::BufWriter::new(file);
            async fn file_future<R: Reopen + SeekableRead + 'static>(
                mut bigwig: BigWigRead<R>,
                chrom: ChromAndSize,
                mut writer: io::BufWriter<TempFileBufferWriter<File>>,
            ) -> Result<(), BBIReadError> {
                for raw_val in bigwig.get_interval(&chrom.name, 0, chrom.length)? {
                    let val = raw_val?;
                    let mut buf = String::with_capacity(50); // Estimate

                    // Using ryu for f32 to string conversion has a ~15% speedup
                    uwrite!(
                        &mut buf,
                        "{}\t{}\t{}\t{}\n",
                        chrom.name.as_str(),
                        val.start,
                        val.end,
                        ryu::Buffer::new().format(val.value)
                    )
                    .unwrap();
                    writer.write(buf.as_bytes())?;
                }
                Ok(())
            }
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

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("BigWigToBedGraph")
        .about("Converts an input bigWig to a bedGraph. Can be multi-threaded for substantial speedups. Note for roughly each core, one temporary file will be opened.")
        .arg(Arg::new("bigwig")
            .help("the bigwig to get convert to bedgraph")
            .index(1)
            .required(true)
        )
        .arg(Arg::new("bedgraph")
            .help("the path of the bedgraph to output to")
            .index(2)
            .required(true)
        )
        .arg(Arg::new("nthreads")
            .short('t')
            .help("Set the number of threads to use. This tool will nearly always benefit from more cores (<= # chroms). Note: for parts of the runtime, the actual usage may be nthreads+1")
            .num_args(1)
            .default_value("6")
            .value_parser(clap::value_parser!(usize)))
        .get_matches();

    let bigwigpath = matches.get_one::<String>("bigwig").unwrap();
    let bedgraphpath = matches.get_one::<String>("bedgraph").unwrap();

    let nthreads = *matches.get_one::<usize>("nthreads").unwrap();

    let bigwig = BigWigRead::open_file(bigwigpath)?;
    let bedgraph = File::create(bedgraphpath)?;

    write_bg(bigwig, bedgraph, nthreads)?;

    Ok(())
}
