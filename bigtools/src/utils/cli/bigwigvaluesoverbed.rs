use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use clap::Parser;

use crate::utils::streaming_linereader::StreamingLineReader;
use crate::BigWigRead;
use crate::{BBIFileRead, BBIReadError};

#[derive(Clone, Debug, PartialEq, Parser)]
#[command(
    name = "bigwigvaluesoverbed",
    about = "For each bed interval, gets the base-resolution values from a bigWig.",
    long_about = None,
)]
pub struct BigWigValuesOverBedArgs {
    /// The input bigwig file
    pub bigwig: String,

    /// The input bed file
    pub bedin: String,

    /// The output bed file
    pub output: String,

    /// If set, the output file will print the name of each bed entry (or `chrom:start-end` if names are not unique) in the first column of each output line.
    #[arg(short = 'n', long)]
    pub names: bool,

    /// Sets the delimiter to use for the output file. (Defaults to tab).
    #[arg(short = 'd', long)]
    #[arg(default_value = "\t")]
    pub delimiter: String,
}

pub fn bigwigvaluesoverbed(args: BigWigValuesOverBedArgs) -> Result<(), Box<dyn Error>> {
    let bigwigpath = args.bigwig;
    let bedinpath = args.bedin;
    let outputpath = args.output;

    let withnames = args.names;
    let mut delimiter = args.delimiter;
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
            use crate::utils::remote_file::RemoteFile;
            let f = RemoteFile::new(&bigwigpath);
            let inbigwig = BigWigRead::open(f)?;
            write(bedin, inbigwig, out, options)?;
        } else {
            let inbigwig = BigWigRead::open_file(&bigwigpath)?.cached();
            write(bedin, inbigwig, out, options)?;
        }
    }
    #[cfg(not(feature = "remote"))]
    {
        let inbigwig = BigWigRead::open_file(&bigwigpath)?.cached();
        write(&bedin, inbigwig, out, options)?;
    }

    Ok(())
}

struct Options {
    withnames: bool,
    delimiter: String,
}

fn write<R: BBIFileRead>(
    bedinpath: &Path,
    mut bigwigin: BigWigRead<R>,
    out: File,
    options: Options,
) -> Result<(), BBIReadError> {
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
