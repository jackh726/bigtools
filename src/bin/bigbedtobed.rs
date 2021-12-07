use std::fs::File;
use std::io::{self, Write};

use clap::{App, Arg};

use futures::task::SpawnExt;

use bigtools::bigwig::{BBIRead, BigBedRead, BigBedReadAttachError, ChromAndSize};
use bigtools::seekableread::{Reopen, SeekableRead};
use bigtools::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};

pub fn write_bed<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
    bigbed: BigBedRead<R, S>,
    mut out_file: File,
    nthreads: usize,
) -> std::io::Result<()> {
    let pool = futures::executor::ThreadPoolBuilder::new()
        .pool_size(nthreads)
        .create()
        .expect("Unable to create thread pool.");

    let chrom_files: Vec<io::Result<(_, TempFileBuffer<File>)>> = bigbed
        .get_chroms()
        .into_iter()
        .map(|chrom| {
            let bigbed = bigbed.clone();
            let (buf, file): (TempFileBuffer<File>, TempFileBufferWriter<File>) =
                TempFileBuffer::new()?;
            let writer = io::BufWriter::new(file);
            async fn file_future<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
                mut bigbed: BigBedRead<R, S>,
                chrom: ChromAndSize,
                mut writer: io::BufWriter<TempFileBufferWriter<File>>,
            ) -> io::Result<()> {
                for raw_val in bigbed.get_interval(&chrom.name, 0, chrom.length)? {
                    let val = raw_val?;
                    let end = if val.rest.len() > 0 {
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

fn main() -> Result<(), BigBedReadAttachError> {
    let matches = App::new("BigBedToBedGraph")
        .about("Converts an input bigBed to a bed. Can be multi-threaded for substantial speedups. Note for roughly each core, one temporary file will be opened.")
        .arg(Arg::new("bigbed")
            .about("the bigbed to get convert to bed")
            .index(1)
            .required(true)
        )
        .arg(Arg::new("bed")
            .about("the path of the bed to output to")
            .index(2)
            .required(true)
        )
        .arg(Arg::new("nthreads")
            .short('t')
            .about("Set the number of threads to use. This tool will nearly always benefit from more cores (<= # chroms). Note: for parts of the runtime, the actual usage may be nthreads+1")
            .takes_value(true)
            .default_value("6"))
        .get_matches();

    let bigbedpath = matches.value_of("bigbed").unwrap().to_owned();
    let bedpath = matches.value_of("bed").unwrap().to_owned();

    let nthreads = {
        let nthreads = matches.value_of("nthreads").unwrap();
        let parsed = nthreads.parse();
        if parsed.is_err() {
            eprintln!("Invalid argument for `nthreads`: must be a positive number");
            return Ok(());
        }
        parsed.unwrap()
    };

    let bigbed = BigBedRead::from_file_and_attach(bigbedpath)?;
    let bed = File::create(bedpath)?;

    write_bed(bigbed, bed, nthreads)?;

    Ok(())
}
