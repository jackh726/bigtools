use std::error::Error;

use byteordered::Endianness;
use clap::Parser;

use crate::{BBIFileRead, BigBedRead};

#[derive(Clone, Debug, Parser, PartialEq)]
#[command(
    name = "bigbedinfo",
    about = "Gets information about a bigBed.",
    long_about = None,
)]
pub struct BigBedInfoArgs {
    /// The bigbed to get info for.
    pub bigbed: String,

    /// If set, will print out the list of chromosomes in the bigBed and their sizes.
    #[arg(long)]
    #[arg(default_value_t = false)]
    pub chroms: bool,

    /// If set, will print out the list of all zoom levels.
    #[arg(long)]
    #[arg(default_value_t = false)]
    pub zooms: bool,

    /// If set, will print out the autosql spec.
    #[arg(long)]
    #[arg(default_value_t = false)]
    pub autosql: bool,
}

pub fn bigbedinfo(args: BigBedInfoArgs) -> Result<(), Box<dyn Error>> {
    let bigbedpath = &args.bigbed;

    fn print_info<R: BBIFileRead>(
        mut bigbed: BigBedRead<R>,
        args: &BigBedInfoArgs,
    ) -> Result<(), Box<dyn Error>> {
        let header = bigbed.info().header;
        println!("version: {}", header.version);
        println!("fieldCount: {}", header.field_count);
        println!(
            "isCompressed: {}",
            (header.uncompress_buf_size > 0)
                .then(|| "yes")
                .unwrap_or("no")
        );
        println!(
            "isSwapped: {}",
            (matches!(header.endianness, Endianness::Big))
                .then(|| "1")
                .unwrap_or("0")
        );
        println!("itemCount: {}", bigbed.item_count()?);
        println!(
            "primaryDataSize: {}",
            num_with_commas(header.full_index_offset - header.full_data_offset)
        );
        let first_zoom_start = bigbed.info().zoom_headers.first().map(|z| z.data_offset);
        if let Some(first_zoom_start) = first_zoom_start {
            println!(
                "primaryIndexSize: {}",
                num_with_commas(first_zoom_start - header.full_index_offset)
            );
        }
        println!("zoomLevels: {}", bigbed.info().zoom_headers.len());
        if args.zooms {
            for zoom in bigbed.info().zoom_headers.iter() {
                println!(
                    "\t{}\t{}",
                    zoom.reduction_level,
                    zoom.index_offset - zoom.data_offset
                );
            }
        }
        println!("chromCount: {}", bigbed.info().chrom_info.len());
        if args.chroms {
            for chrom in bigbed.info().chrom_info.iter() {
                println!("\t{} {} {}", chrom.name, chrom.id, chrom.length);
            }
        }
        if args.autosql {
            let autosql = bigbed.autosql()?;
            if autosql.len() == 0 {
                println!("as:  n/a");
            } else {
                println!("as:");
                print!("{}", autosql);
            }
        }
        let summary = bigbed.get_summary()?;
        println!("basesCovered: {}", num_with_commas(summary.bases_covered));
        println!(
            "meanDepth: {:.6}",
            summary.sum / summary.bases_covered as f64
        );
        println!("minDepth: {:.6}", summary.min_val);
        println!("maxDepth: {:.6}", summary.max_val);
        let var = (summary.sum_squares
            - (summary.sum * summary.sum) / summary.bases_covered as f64)
            / (summary.bases_covered as f64 - 1.0);
        let std = var.sqrt();
        println!("std of depth: {:.6}", std);

        Ok(())
    }

    #[cfg(feature = "remote")]
    {
        if bigbedpath.starts_with("http") {
            use crate::utils::remote_file::RemoteFile;
            let f = RemoteFile::new(bigbedpath);
            let bigbed = BigBedRead::open(f)?;
            print_info(bigbed, &args)?;
            return Ok(());
        }
    }

    let bigbed = BigBedRead::open_file(bigbedpath)?;
    print_info(bigbed, &args)?;

    Ok(())
}

fn num_with_commas(mut num: u64) -> String {
    if num == 0 {
        return format!("0");
    }
    let mut formatted = String::new();
    loop {
        let remainder = num % 1000;

        if remainder != 0 {
            if num > 1000 && remainder < 10 {
                formatted = format!("00{remainder}") + &formatted;
            } else if num > 1000 && remainder < 100 {
                formatted = format!("0{remainder}") + &formatted;
            } else {
                formatted = format!("{remainder}") + &formatted;
            }
            num -= remainder;
        } else {
            formatted = format!("000") + &formatted;
        }

        num = num / 1000;

        if num > 0 {
            formatted = ",".to_string() + &formatted;
        }

        if num == 0 {
            break formatted;
        }
    }
}

#[test]
fn test_num_with_commas() {
    assert_eq!("0", num_with_commas(0));
    assert_eq!("987", num_with_commas(987));
    assert_eq!("1,000", num_with_commas(1000));
    assert_eq!("1,987", num_with_commas(1987));
    assert_eq!("12,000", num_with_commas(12000));
    assert_eq!("12,987", num_with_commas(12987));
    assert_eq!("123,987", num_with_commas(123987));
    assert_eq!("4,023,987", num_with_commas(4023987));
    assert_eq!("4,123,987", num_with_commas(4123987));
    assert_eq!("45,123,987", num_with_commas(45123987));
    assert_eq!("456,123,987", num_with_commas(456123987));
    assert_eq!("9,000,123,987", num_with_commas(9000123987));
    assert_eq!("9,456,000,987", num_with_commas(9456000987));
    assert_eq!("9,456,123,000", num_with_commas(9456123000));
    assert_eq!("9,456,123,987", num_with_commas(9456123987));
}
