use std::error::Error;

use bigtools::bedchromdata::BedParserStreamingIterator;
use tokio::runtime;

#[test]
fn bigbedwrite_test() -> Result<(), Box<dyn Error>> {
    use std::collections::HashMap;
    use std::fs::File;
    use std::path::PathBuf;

    use tempfile;

    use bigtools::bed::bedparser::BedParser;
    use bigtools::utils::chromvalues::ChromValues;
    use bigtools::{BigBedRead, BigBedWrite};

    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    let mut bed = dir.clone();
    bed.push("small.bed");

    let first = {
        let infile = File::open(bed.clone())?;
        let mut vals_iter = BedParser::from_bed_file(infile);
        let (_, mut group) = vals_iter.next_chrom().unwrap().unwrap();
        group.next().unwrap().unwrap()
    };

    let runtime = runtime::Builder::new_multi_thread()
        .worker_threads(6)
        .build()
        .expect("Unable to create runtime.");

    let infile = File::open(bed)?;
    let tempfile = tempfile::NamedTempFile::new()?;
    let mut vals_iter = BedParser::from_bed_file(infile);
    let mut outb = BigBedWrite::create_file(tempfile.path().to_string_lossy().to_string());
    outb.autosql = {
        let (_, mut group) = vals_iter.next_chrom().unwrap().unwrap();
        let first = group.peek().unwrap().unwrap();
        Some(bigtools::bed::autosql::bed_autosql(&first.rest))
    };
    outb.options.compress = false;

    let mut chrom_map = HashMap::new();
    chrom_map.insert("chr17".to_string(), 83257441);
    chrom_map.insert("chr18".to_string(), 80373285);
    chrom_map.insert("chr19".to_string(), 58617616);

    let chsi = BedParserStreamingIterator::new(vals_iter, false);
    outb.write(chrom_map, chsi, runtime).unwrap();

    let mut bwread = BigBedRead::open_file(&tempfile.path().to_string_lossy()).unwrap();

    let chroms = bwread.chroms();
    assert_eq!(chroms.len(), 3);
    assert_eq!(chroms[0].name, "chr17");
    assert_eq!(chroms[0].length, 83257441);

    assert_eq!(
        &bwread.autosql().unwrap(),
        "table bed\n\"Browser Extensible Data\"\n(\n    string chrom;       \"Reference sequence chromosome or scaffold\"\n    uint   chromStart;  \"Start position in chromosome\"\n    uint   chromEnd;    \"End position in chromosome\"\n   string name;        \"Name of item.\"\n   uint score;          \"Score (0-1000)\"\n)",
    );

    let mut intervals = bwread.get_interval("chr17", 0, 83257441)?;
    let first_interval = intervals.next().unwrap().unwrap();
    assert_eq!(first.start, first_interval.start);
    assert_eq!(first.end, first_interval.end);
    assert_eq!(first.rest, first_interval.rest);

    Ok(())
}
