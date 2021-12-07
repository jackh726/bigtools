use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use clap::{App, Arg};

use bigtools::bigwig::{BigWigRead, BigWigReadAttachError};
use bigtools::seekableread::{Reopen, SeekableRead};
use bigtools::streaming_linereader::StreamingLineReader;

struct Options {
    withnames: bool,
    delimiter: String,
}

fn write<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
    bedinpath: &Path,
    mut bigwigin: BigWigRead<R, S>,
    out: File,
    options: Options,
) -> io::Result<()> {
    let uniquenames = {
        if !options.withnames {
            true
        } else {
            let reader = BufReader::new(File::open(bedinpath)?);
            let mut lines = reader
                .lines()
                .take(10)
                .map(|line| -> io::Result<Option<String>> {
                    let l = line?;
                    let mut split = l.splitn(5, '\t');
                    let _chrom = split
                        .next()
                        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData));
                    let _start = split
                        .next()
                        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData));
                    let _end = split
                        .next()
                        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData));
                    let name = split.next();
                    Ok(name.map(|name| name.to_string()))
                })
                .collect::<io::Result<Vec<_>>>()?;
            lines.sort();
            lines.dedup();
            lines.len() == 10
        }
    };

    let bedin = File::open(bedinpath).unwrap();
    let mut bedstream = StreamingLineReader::new(BufReader::new(bedin));
    let mut outwriter = BufWriter::new(out);

    while let Some(line) = bedstream.read() {
        let line = line?;
        let mut split = line.trim().splitn(5, '\t');
        let chrom = split.next().expect("Missing chrom");
        let start = split.next().expect("Missing start").parse::<u32>().unwrap();
        let end = split.next().expect("Missing end").parse::<u32>().unwrap();
        let name = split.next();
        let interval = bigwigin
            .get_interval(chrom, start, end)?
            .collect::<Result<Vec<_>, _>>()?;
        let size = end - start;
        let mut vals: Vec<f32> = vec![0f32; size as usize];
        for val in interval {
            for i in val.start..val.end {
                vals[(i - start) as usize] = val.value;
            }
        }
        let vals_strings: Vec<String> = vals.into_iter().map(|v| v.to_string()).collect();
        let vals_string = &vals_strings[..].join(&options.delimiter);
        if options.withnames {
            let uniquename = if uniquenames {
                name.expect("Bad bed format (no name).").to_owned()
            } else {
                format!("{}:{}-{}", chrom, start, end)
            };
            outwriter.write_fmt(format_args!(
                "{}{}{}\n",
                uniquename, &options.delimiter, vals_string
            ))?;
        } else {
            outwriter.write_fmt(format_args!("{}\n", vals_string))?;
        }
    }
    Ok(())
}

fn main() -> Result<(), BigWigReadAttachError> {
    let matches = App::new("BigWigInfo")
        .arg(Arg::new("bigwig")
                .about("The input bigwig file")
                .index(1)
                .required(true)
            )
        .arg(Arg::new("bedin")
                .about("The input bed file")
                .index(2)
                .required(true)
            )
        .arg(Arg::new("output")
                .about("The output file")
                .index(3)
                .required(true)
            )
        .arg(Arg::new("names")
                .short('n')
                .about("If set, the output file will print the name of each bed entry (or `chrom:start-end` if names are not unique) in the first column of each output line.")
            )
        .arg(Arg::new("delimiter")
                .short('d')
                .takes_value(true)
                .about("Sets the delimiter to use for the output file. (Defaults to tab).")
            )
        .get_matches();

    let bigwigpath = matches.value_of("bigwig").unwrap();
    let bedinpath = matches.value_of("bedin").unwrap();
    let outputpath = matches.value_of("output").unwrap();

    let withnames = matches.is_present("names");
    let mut delimiter = matches.value_of("delimiter").unwrap_or("\t").to_owned();
    if delimiter == "\\t" {
        delimiter = String::from("\t");
    }

    let bedin = Path::new(&bedinpath);
    if !bedin.exists() {
        eprintln!("File does not exist: {}", bedin.display());
        return Ok(());
    }

    let out = File::create(outputpath)?;
    let options = Options {
        withnames,
        delimiter,
    };

    #[cfg(feature = "remote")]
    {
        if bigwigpath.starts_with("http") {
            use bigtools::remote_file::RemoteFile;
            let f = RemoteFile::new(bigwigpath);
            let inbigwig = BigWigRead::from(f)?;
            write(bedin, inbigwig, out, options)?;
        } else {
            let inbigwig = BigWigRead::from_file_and_attach(bigwigpath)?;
            write(bedin, inbigwig, out, options)?;
        }
    }
    #[cfg(not(feature = "remote"))]
    {
        let inbigwig = BigWigRead::from_file_and_attach(bigwigpath)?;
        write(&bedin, inbigwig, out, options)?;
    }

    Ok(())
}
