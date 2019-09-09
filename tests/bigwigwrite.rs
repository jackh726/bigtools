use std::io::{self};

#[test]
fn test() -> io::Result<()> {
    use std::collections::HashMap;
    use std::fs::File;
    use std::path::PathBuf;

    use tempfile;

    use bigwig2::bedgraphparser::{self, BedGraphParser};
    use bigwig2::chromvalues::{ChromGroups, ChromValues};
    use bigwig2::bigwig::{BBIRead, BigWigRead, BigWigWrite};

    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    let mut single_chrom_bedgraph = dir.clone();
    single_chrom_bedgraph.push("single_chrom.bedGraph");
    
    let first = 'f: loop {
        let infile = File::open(single_chrom_bedgraph.clone())?;
        let mut vals_iter = BedGraphParser::from_file(infile);
        let (_, mut group) = vals_iter.next()?.unwrap();
        break group.next()?.unwrap();
    };

    let infile = File::open(single_chrom_bedgraph)?;
    let tempfile = tempfile::NamedTempFile::new()?;
    let vals_iter = BedGraphParser::from_file(infile);
    let outb = BigWigWrite::create_file(tempfile.path().to_string_lossy().to_string());

    let mut chrom_map = HashMap::new();
    chrom_map.insert("chr17".to_string(), 83257441);
    let chsi = bedgraphparser::get_chromgroupstreamingiterator(vals_iter, outb.options.clone(), chrom_map.clone());
    outb.write_groups(chrom_map, chsi).unwrap();

    let mut bwread = BigWigRead::from_file_and_attach(tempfile.path().to_string_lossy().to_string()).unwrap(); 

    let chroms = bwread.get_chroms();
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