use std::error::Error;
use std::fs::File;
use std::io::{self, Write};

use clap::{Arg, Command};

use futures::task::SpawnExt;

use bigtools::bbi::{BBIRead, BigBedRead, ChromAndSize};
use bigtools::bbiread::BBIReadError;
use bigtools::utils::reopen::{Reopen, SeekableRead};
use bigtools::utils::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};

pub fn write_bed<R: Reopen + SeekableRead + Send + 'static>(
    bigbed: BigBedRead<R>,
    mut out_file: File,
    nthreads: usize,
) -> Result<(), BBIReadError> {
    let pool = futures::executor::ThreadPoolBuilder::new()
        .pool_size(nthreads)
        .create()
        .expect("Unable to create thread pool.");

    let chrom_files: Vec<io::Result<(_, TempFileBuffer<File>)>> = bigbed
        .get_chroms()
        .into_iter()
        .map(|chrom| {
            let bigbed = bigbed.reopen()?;
            let (buf, file): (TempFileBuffer<File>, TempFileBufferWriter<File>) =
                TempFileBuffer::new()?;
            let writer = io::BufWriter::new(file);
            async fn file_future<R: SeekableRead + 'static>(
                mut bigbed: BigBedRead<R>,
                chrom: ChromAndSize,
                mut writer: io::BufWriter<TempFileBufferWriter<File>>,
            ) -> Result<(), BBIReadError> {
                for raw_val in bigbed.get_interval(&chrom.name, 0, chrom.length)? {
                    let val = raw_val?;
                    let end = if !val.rest.is_empty() {
                        format!("\t{}\n", val.rest)
                    } else {
                        "\n".to_string()
                    };
                    writer.write_fmt(format_args!(
                        "{}\t{}\t{}{}",
                        chrom.name, val.start, val.end, end
                    ))?;
                }
                Ok(())
            }
            let handle = pool
                .spawn_with_handle(file_future(bigbed, chrom, writer))
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
    let matches = Command::new("BigBedToBedGraph")
        .about("Converts an input bigBed to a bed. Can be multi-threaded for substantial speedups. Note for roughly each core, one temporary file will be opened.")
        .arg(Arg::new("bigbed")
            .help("the bigbed to get convert to bed")
            .index(1)
            .required(true)
        )
        .arg(Arg::new("bed")
            .help("the path of the bed to output to")
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

    let bigbedpath = matches.get_one::<String>("bigbed").unwrap().to_owned();
    let bedpath = matches.get_one::<String>("bed").unwrap().to_owned();

    let nthreads = *matches.get_one::<usize>("nthreads").unwrap();

    let bigbed = BigBedRead::open_file(&bigbedpath)?;
    let bed = File::create(bedpath)?;

    write_bed(bigbed, bed, nthreads)?;

    Ok(())
}
