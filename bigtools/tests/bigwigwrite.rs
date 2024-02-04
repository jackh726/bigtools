use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self};
use std::path::PathBuf;

use tempfile;

use bigtools::bed::bedparser::BedParser;
use bigtools::bedchromdata::BedParserStreamingIterator;
use bigtools::utils::chromvalues::ChromValues;
use bigtools::{BigWigRead, BigWigWrite, Value};
use tokio::runtime;

#[test]
fn test() -> Result<(), Box<dyn Error>> {
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    let mut single_chrom_bedgraph = dir.clone();
    single_chrom_bedgraph.push("single_chrom.bedGraph");

    let first = {
        let infile = File::open(single_chrom_bedgraph.clone())?;
        let mut vals_iter = BedParser::from_bedgraph_file(infile);
        let (_, mut group) = vals_iter.next_chrom().unwrap().unwrap();
        group.next().unwrap().unwrap()
    };

    let runtime = runtime::Builder::new_multi_thread()
        .worker_threads(6)
        .build()
        .expect("Unable to create runtime.");

    let infile = File::open(single_chrom_bedgraph)?;
    let tempfile = tempfile::NamedTempFile::new()?;
    let vals_iter = BedParser::from_bedgraph_file(infile);
    let outb = BigWigWrite::create_file(tempfile.path().to_string_lossy().to_string());

    let mut chrom_map = HashMap::new();
    chrom_map.insert("chr17".to_string(), 83257441);

    let chsi = BedParserStreamingIterator::new(vals_iter, false);
    outb.write(chrom_map, chsi, runtime).unwrap();

    let mut bwread = BigWigRead::open_file(&tempfile.path().to_string_lossy()).unwrap();

    let chroms = bwread.chroms();
    assert_eq!(chroms.len(), 1);
    assert_eq!(chroms[0].name, "chr17");
    assert_eq!(chroms[0].length, 83257441);

    let mut intervals = bwread.get_interval("chr17", 0, 83257441)?;
    let first_interval = intervals.next().unwrap().unwrap();
    assert_eq!(first.start, first_interval.start);
    assert_eq!(first.end, first_interval.end);
    assert_eq!(first.value, first_interval.value);

    Ok(())
}

#[test]
fn test_multi_pass() -> Result<(), Box<dyn Error>> {
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    let mut single_chrom_bedgraph = dir.clone();
    single_chrom_bedgraph.push("single_chrom.bedGraph");

    let first = {
        let infile = File::open(single_chrom_bedgraph.clone())?;
        let mut vals_iter = BedParser::from_bedgraph_file(infile);
        let (_, mut group) = vals_iter.next_chrom().unwrap().unwrap();
        group.next().unwrap().unwrap()
    };

    let runtime = runtime::Builder::new_multi_thread()
        .worker_threads(6)
        .build()
        .expect("Unable to create runtime.");

    let tempfile = tempfile::NamedTempFile::new()?;

    let outb = BigWigWrite::create_file(tempfile.path().to_string_lossy().to_string());

    let mut chrom_map = HashMap::new();
    chrom_map.insert("chr17".to_string(), 83257441);

    outb.write_multipass(
        || {
            let infile = File::open(single_chrom_bedgraph.clone())?;
            let vals_iter = BedParser::from_bedgraph_file(infile);
            let chsi = BedParserStreamingIterator::new(vals_iter, false);
            Ok(chsi)
        },
        chrom_map,
        runtime,
    )
    .unwrap();

    let mut bwread = BigWigRead::open_file(&tempfile.path().to_string_lossy()).unwrap();

    let chroms = bwread.chroms();
    assert_eq!(chroms.len(), 1);
    assert_eq!(chroms[0].name, "chr17");
    assert_eq!(chroms[0].length, 83257441);

    let mut intervals = bwread.get_interval("chr17", 0, 83257441)?;
    let first_interval = intervals.next().unwrap().unwrap();
    assert_eq!(first.start, first_interval.start);
    assert_eq!(first.end, first_interval.end);
    assert_eq!(first.value, first_interval.value);

    Ok(())
}

#[test]
fn test_multi_chrom() -> io::Result<()> {
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    let mut multi_chrom_bedgraph = dir.clone();
    multi_chrom_bedgraph.push("multi_chrom.bedGraph");

    let runtime = runtime::Builder::new_multi_thread()
        .worker_threads(6)
        .build()
        .expect("Unable to create runtime.");

    let infile = File::open(multi_chrom_bedgraph)?;
    let tempfile = tempfile::NamedTempFile::new()?;
    let vals_iter = BedParser::from_bedgraph_file(infile);
    let outb = BigWigWrite::create_file(tempfile.path().to_string_lossy().to_string());

    let mut chrom_map = HashMap::new();
    chrom_map.insert("chr1".to_string(), 248956422);
    chrom_map.insert("chr2".to_string(), 242193529);
    chrom_map.insert("chr3".to_string(), 198295559);
    chrom_map.insert("chr4".to_string(), 190214555);
    chrom_map.insert("chr5".to_string(), 181538259);
    chrom_map.insert("chr6".to_string(), 170805979);

    let chsi = BedParserStreamingIterator::new(vals_iter, false);
    outb.write(chrom_map, chsi, runtime).unwrap();

    let mut bwread = BigWigRead::open_file(&tempfile.path().to_string_lossy()).unwrap();

    let chroms = bwread.chroms();
    assert_eq!(chroms.len(), 6);

    assert_eq!(
        bwread.get_interval("chr1", 0, 248956422).unwrap().count(),
        200
    );
    assert_eq!(
        bwread.get_interval("chr6", 0, 170805979).unwrap().count(),
        2000
    );

    Ok(())
}

#[test]
fn test_iter() {
    let chroms = vec![
        "chrY", "chrY", "chrY", "chrY", "chrY", "chrY", "chrY", "chrY", "chrY", "chrY",
    ];
    let starts: Vec<u32> = vec![530, 538, 584, 713, 751, 860, 865, 873, 879, 902];
    let ends: Vec<u32> = vec![538, 584, 713, 751, 860, 865, 873, 879, 902, 923];

    let iter = chroms
        .into_iter()
        .zip(starts.into_iter().zip(ends.into_iter()));
    let iter = iter.map(|(chrom, (start, end))| {
        (
            chrom,
            Value {
                start,
                end,
                value: 0.0,
            },
        )
    });

    let vals_iter = BedParserStreamingIterator::new(BedParser::wrap_infallible_iter(iter), true);

    let chrom_map = HashMap::from([("chrY".to_string(), 57_227_415)]);

    let runtime = tokio::runtime::Builder::new_current_thread()
        .build()
        .expect("Unable to create runtime.");

    let tempfile = tempfile::NamedTempFile::new().unwrap();
    let outb = BigWigWrite::create_file(tempfile.path().to_string_lossy().to_string());
    outb.write_singlethreaded(chrom_map, vals_iter, runtime)
        .unwrap();
}
