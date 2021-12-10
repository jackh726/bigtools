use std::io::{self};

#[test]
fn bigbedwrite_test() -> io::Result<()> {
    use std::collections::HashMap;
    use std::fs::File;
    use std::path::PathBuf;

    use tempfile;

    use bigtools::bed::bedparser::{self, BedParser};
    use bigtools::bigwig::{BBIRead, BigBedRead, BigBedWrite};
    use bigtools::utils::chromvalues::ChromValues;

    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    let mut bed = dir.clone();
    bed.push("small.bed");

    let first = {
        let infile = File::open(bed.clone())?;
        let mut vals_iter = BedParser::from_bed_file(infile);
        let (_, mut group) = vals_iter.next().unwrap()?;
        group.next().unwrap().unwrap()
    };

    let pool = futures::executor::ThreadPoolBuilder::new()
        .pool_size(6)
        .create()
        .expect("Unable to create thread pool.");

    let infile = File::open(bed)?;
    let tempfile = tempfile::NamedTempFile::new()?;
    let mut vals_iter = BedParser::from_bed_file(infile);
    let mut outb = BigBedWrite::create_file(tempfile.path().to_string_lossy().to_string());
    outb.autosql = {
        let (_, mut group) = vals_iter.next().unwrap()?;
        let first = group.peek().unwrap();
        Some(bigtools::bed::autosql::bed_autosql(&first.rest))
    };
    outb.options.compress = false;

    let mut chrom_map = HashMap::new();
    chrom_map.insert("chr17".to_string(), 83257441);
    chrom_map.insert("chr18".to_string(), 80373285);
    chrom_map.insert("chr19".to_string(), 58617616);

    let chsi = bedparser::BedParserStreamingIterator::new(
        vals_iter,
        chrom_map.clone(),
        pool.clone(),
        false,
    );
    outb.write(chrom_map, chsi).unwrap();

    let mut bwread =
        BigBedRead::from_file_and_attach(tempfile.path().to_string_lossy().to_string()).unwrap();

    let chroms = bwread.get_chroms();
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
