use std::error::Error;

#[test]
fn test_valid_read() -> Result<(), Box<dyn Error>> {
    use std::path::PathBuf;

    use bigtools::bigwig::{BBIRead, BigWigRead};

    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    let mut valid_bigwig = dir.clone();
    valid_bigwig.push("valid.bigWig");

    let mut bwread = BigWigRead::from_file_and_attach(&valid_bigwig.to_string_lossy()).unwrap();

    // Test that chrom tree parsing works
    let chroms = bwread.get_chroms();
    assert_eq!(chroms.len(), 1);
    // chr17
    assert_eq!(chroms[0].length, 83257441);

    // Test that header/summary is correct
    let summary = bwread.get_summary()?;
    assert_eq!(summary.bases_covered, 137894);
    assert_eq!(summary.max_val, 14254.0);

    // This tests simply reading an interval. Importantly, the provided interval actually splits an interval in two, so this also tests correct splitting on read.
    let first_interval = bwread
        .get_interval("chr17", 0, 59899)?
        .next()
        .unwrap()
        .unwrap();
    assert_eq!(first_interval.start, 59898);
    assert_eq!(first_interval.end, 59899);
    assert_eq!(first_interval.value, 0.06792);

    Ok(())
}

#[test]
fn test_values() -> Result<(), Box<dyn Error>> {
    use std::path::PathBuf;

    use bigtools::bigwig::BigWigRead;

    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    let mut valid_bigwig = dir.clone();
    valid_bigwig.push("valid.bigWig");

    let mut bwread = BigWigRead::from_file_and_attach(&valid_bigwig.to_string_lossy()).unwrap();

    let vals = bwread.values("chr17", 0, 59899)?;
    assert_eq!(vals.len(), 59899);
    assert_eq!(vals[59898], 0.06792);
    Ok(())
}
