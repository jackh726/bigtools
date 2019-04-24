use std::fs::File;
use std::io::Write;

use crate::bigwig::BigWigRead;

pub fn write_bg(bigwig: BigWigRead, out_file: File) -> std::io::Result<()> {
    let mut writer = std::io::BufWriter::new(out_file);
    for chrom in bigwig.get_chroms() {
        for val in bigwig.get_interval(&chrom.name, 1, chrom.length)? {
            writer.write_fmt(format_args!("{}\t{}\t{}\t{}\n", chrom.name, val.start, val.end, val.value))?;
        }
    }

    Ok(())
}