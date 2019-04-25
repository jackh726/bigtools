use std::fs::File;
use std::io::{self, Seek, Write};

use futures::future::FutureExt;

use crate::bigwig::BigWigRead;
use crate::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};

pub fn write_bg(bigwig: BigWigRead, mut out_file: File) -> std::io::Result<()> {
    let chrom_files: Vec<io::Result<(_, TempFileBuffer)>> = bigwig.get_chroms().into_iter().map(|chrom| {
        let bigwig = bigwig.clone();
        let (mut buf, mut file) = TempFileBuffer::new()?;
        let mut writer = std::io::BufWriter::new(file);
        let file_future = async move || -> io::Result<()> {
            for val in bigwig.get_interval(&chrom.name, 1, chrom.length)? {
                writer.write_fmt(format_args!("{}\t{}\t{}\t{}\n", chrom.name, val.start, val.end, val.value))?;
            }
            Ok(())
        };
        let (remote, handle) = file_future().remote_handle();
        std::thread::spawn(move || {
            futures::executor::block_on(remote);
        });
        Ok((handle,buf))
    }).collect::<Vec<_>>();

    for res in chrom_files {
        let (f, mut buf) = res.unwrap();
        buf.switch(out_file).unwrap();
        futures::executor::block_on(f).unwrap();
        out_file = buf.await();
    }

    Ok(())
}