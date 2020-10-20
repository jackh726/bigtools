use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use clap::{App, Arg};

use bigtools::bbiwrite::InputSortType;
use bigtools::bedparser::{self, BedParser};
use bigtools::bigwig::{BigBedWrite, WriteGroupsError};

fn main() -> Result<(), WriteGroupsError> {
    let matches = App::new("BedToBigBed")
        .arg(Arg::new("bed")
                .about("the n to convert to a bigbed")
                .index(1)
                .required(true)
            )
        .arg(Arg::new("chromsizes")
                .about("A chromosome sizes file. Each line should be have a chromosome and its size in bases, separated by whitespace.")
                .index(2)
                .required(true)
            )
        .arg(Arg::new("output")
                .about("The output bigbed path")
                .index(3)
                .required(true)
            )
        .arg(Arg::new("nthreads")
                .short('t')
                .about("Set the number of threads to use")
                .takes_value(true)
                .default_value("6"))
        .arg(Arg::new("nzooms")
                .short('z')
                .about("Set the maximum of zooms to create.")
                .takes_value(true)
                .default_value("10"))
        .arg(Arg::new("uncompressed")
                .short('u')
                .about("Don't use compression."))
        .arg(Arg::new("sorted")
                .short('s')
                .about("Sets whether the input is sorted. Can take `all`, `start`, or `none`. `all` means that the input bedGraph is sorted by chroms and start (`sort -k1,1 -k2,2n`). `start` means that the the chroms are out of order but the starts within a chrom is sorted. `none` means that the file is not sorted at all. `all` is default. `none` currently errors but may be supported in the future. Note that using a value other than `all` will not guarantee (though likely) support for third-party tools.")
                .takes_value(true)
                .default_value("all"))
        .arg(Arg::new("autosql")
                .short('a')
                .about("The path to an .as file containing the autosql that defines the fields in this bigBed")
                .takes_value(true))
        .get_matches();

    let bedpath = matches.value_of("bed").unwrap().to_owned();
    let chrom_map = matches.value_of("chromsizes").unwrap().to_owned();
    let bigwigpath = matches.value_of("output").unwrap().to_owned();
    let nthreads = {
        let nthreads = matches.value_of("nthreads").unwrap();
        let parsed = nthreads.parse();
        if parsed.is_err() {
            eprintln!("Invalid argument for `nthreads`: must be a positive number");
            return Ok(());
        }
        parsed.unwrap()
    };
    let nzooms = {
        let nzooms = matches.value_of("nzooms").unwrap();
        let parsed = nzooms.parse();
        if parsed.is_err() {
            eprintln!("Invalid argument for `nzooms`: must be a positive number");
            return Ok(());
        }
        parsed.unwrap()
    };
    let uncompressed = { matches.is_present("uncompressed") };
    let input_sort_type = match matches.value_of("sorted") {
        None => InputSortType::ALL,
        Some("all") => InputSortType::ALL,
        Some("start") => InputSortType::START,
        Some("none") => {
            eprintln!("Using completely unsorted input is not implemented yet.");
            return Ok(());
        }
        Some(sorted) => {
            eprintln!(
                "Invalid option for `sorted`: `{}`. Options are `all`, `start`, or `none`.",
                sorted
            );
            return Ok(());
        }
    };

    let mut outb = BigBedWrite::create_file(bigwigpath);
    outb.options.max_zooms = nzooms;
    outb.options.compress = !uncompressed;
    outb.options.input_sort_type = input_sort_type;
    let chrom_map: HashMap<String, u32> = BufReader::new(File::open(chrom_map)?)
        .lines()
        .filter(|l| match l {
            Ok(s) => !s.is_empty(),
            _ => true,
        })
        .map(|l| {
            let words = l.expect("Split error");
            let mut split = words.split_whitespace();
            (
                split.next().expect("Missing chrom").to_owned(),
                split.next().expect("Missing size").parse::<u32>().unwrap(),
            )
        })
        .collect();

    let pool = futures::executor::ThreadPoolBuilder::new()
        .pool_size(nthreads)
        .create()
        .expect("Unable to create thread pool.");

    let infile = File::open(bedpath)?;
    let mut vals_iter = BedParser::from_bed_file(infile);
    let options = outb.options.clone();

    let autosql = match matches.value_of("autosql") {
        None => {
            use bigtools::chromvalues::ChromGroups;
            use bigtools::chromvalues::ChromValues;
            let (_, mut group) = vals_iter.peek()?.unwrap();
            let first = group.peek().unwrap();
            bigtools::autosql::bed_autosql(&first.rest)
        }
        Some(file) => {
            std::fs::read_to_string(file)?
        }
    };
    outb.autosql = Some(autosql);

    let parse_fn = move |chrom, chrom_id, chrom_length, group| {
        BigBedWrite::begin_processing_chrom(
            chrom,
            chrom_id,
            chrom_length,
            group,
            pool.clone(),
            options.clone(),
        )
    };
    let allow_out_of_order_chroms = match &outb.options.input_sort_type {
        InputSortType::ALL => false,
        _ => true,
    };
    let chsi = bedparser::BedParserChromGroupStreamingIterator::new(
        vals_iter,
        chrom_map.clone(),
        Box::new(parse_fn),
        allow_out_of_order_chroms,
    );
    outb.write_groups(chrom_map, chsi)?;

    Ok(())
}
