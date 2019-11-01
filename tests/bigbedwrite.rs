use std::io::{self};

#[test]
fn test() -> io::Result<()> {
    use std::collections::HashMap;
    use std::fs::File;
    use std::path::PathBuf;

    use tempfile;

    use bigtools::bedparser::{self, BedParser};
    use bigtools::chromvalues::{ChromGroups, ChromValues};
    use bigtools::bigwig::{BBIRead, BigBedRead, BigBedWrite};

    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    let mut bed = dir.clone();
    bed.push("small.bed");
    
    let first = {
        let infile = File::open(bed.clone())?;
        let mut vals_iter = BedParser::from_bed_file(infile);
        let (_, mut group) = vals_iter.next()?.unwrap();
        group.next()?.unwrap()
    };

    let pool = futures::executor::ThreadPoolBuilder::new().pool_size(6).create().expect("Unable to create thread pool.");

    let infile = File::open(bed)?;
    let tempfile = tempfile::NamedTempFile::new()?;
    let vals_iter = BedParser::from_bed_file(infile);
    let outb = BigBedWrite::create_file(tempfile.path().to_string_lossy().to_string());

    let options = outb.options.clone();

    let mut chrom_map = HashMap::new();
    chrom_map.insert("chr17".to_string(), 83257441);
    chrom_map.insert("chr18".to_string(), 80373285);
    chrom_map.insert("chr19".to_string(), 58617616);

    let parse_fn = move |chrom, chrom_id, chrom_length, group| {
        BigBedWrite::begin_processing_chrom(chrom, chrom_id, chrom_length, group, pool.clone(), options.clone())
    };
    let chsi = bedparser::BedParserChromGroupStreamingIterator::new(vals_iter, chrom_map.clone(), Box::new(parse_fn));
    outb.write_groups(chrom_map, chsi).unwrap();

    let mut bwread = BigBedRead::from_file_and_attach(tempfile.path().to_string_lossy().to_string()).unwrap(); 

    let chroms = bwread.get_chroms();
    assert_eq!(chroms.len(), 3);
    assert_eq!(chroms[0].name, "chr17");
    assert_eq!(chroms[0].length, 83257441);

    let mut intervals = bwread.get_interval("chr17", 0, 83257441)?;
    let first_interval = intervals.next().unwrap().unwrap();
    assert_eq!(first.start, first_interval.start);
    assert_eq!(first.end, first_interval.end);
    assert_eq!(first.rest, first_interval.rest);

    Ok(())
}
