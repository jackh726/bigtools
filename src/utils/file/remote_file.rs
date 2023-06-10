use std::fs::File;
use std::io::{self, Cursor, Read, Seek, SeekFrom, Write};

use tempfile;

use crate::utils::file::reopen::Reopen;

const READ_SIZE: usize = 10 * 1024; // 10 KB chunks

// Remote file reads are cached to a temporary file. The size of each block
// (with the exception of the last block) is `READ_SIZE + 1` bytes. The first
// byte of a block is `0` if the data hasn't been written yet, or `1` if it
// has (and has enough data to fill the block). A value of `2` signifies that
// there wasn't enough data to fill the block, which only should happen for the
// last block.

pub struct RemoteFile {
    url: String,
    current_position: u64,
    current: Option<(u64, Cursor<Vec<u8>>)>,
    cache: Option<File>,
}

impl RemoteFile {
    pub fn new(url: &str) -> RemoteFile {
        RemoteFile {
            url: url.to_string(),
            current_position: 0,
            current: None,
            cache: None,
        }
    }
}

impl RemoteFile {
    fn read_current_block(&mut self, read_size: u64) -> io::Result<u64> {
        let block = self.current_position / READ_SIZE as u64;
        let block_start = block * READ_SIZE as u64;
        let cache_block_start = block * (READ_SIZE as u64 + 1);
        let cache = match self.cache.as_mut() {
            None => {
                self.cache = Some(tempfile::tempfile()?);
                self.cache.as_mut().unwrap()
            }
            Some(cache) => cache,
        };
        use byteorder::ReadBytesExt;
        use byteorder::WriteBytesExt;
        cache.seek(SeekFrom::Start(cache_block_start))?;
        let status = cache.read_u8().unwrap_or(0);
        if status == 1 {
            let mut bytes = vec![0u8; READ_SIZE];
            cache.read_exact(&mut bytes)?;
            self.current = Some((block_start, Cursor::new(bytes.to_vec())));
            return Ok(READ_SIZE as u64);
        } else if status == 2 {
            let bytes_available = cache.read_u64::<byteorder::BigEndian>()?;
            let mut bytes = vec![0u8; bytes_available as usize];
            cache.read_exact(&mut bytes)?;
            self.current = Some((block_start, Cursor::new(bytes.to_vec())));
            return Ok(bytes_available);
        }

        let read_len = {
            let cur_pos = self.current_position;
            let block = cur_pos / (READ_SIZE as u64);
            let block_start = block * (READ_SIZE as u64);
            let blocks_to_read = (cur_pos - block_start + read_size - 1) / (READ_SIZE as u64) + 1;
            blocks_to_read * (READ_SIZE as u64)
        };

        let resp = attohttpc::get(&self.url)
            .header(
                "range",
                format!(
                    "bytes={}-{}",
                    block_start,
                    block_start + read_len as u64 - 1
                ),
            )
            .send()?;
        let bytes = if resp.is_success() {
            resp.bytes()?
        } else {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Unable to connect to server to receive file.".to_string(),
            ));
        };
        cache.seek(SeekFrom::Start(cache_block_start))?;
        let blocks_to_write = if bytes.len() == read_len as usize {
            bytes.len() / READ_SIZE
        } else {
            (bytes.len() + READ_SIZE - 1) / READ_SIZE
        };
        for start in 0..blocks_to_write {
            let begin = start * READ_SIZE;
            let end = ((start + 1) * READ_SIZE).min(bytes.len());
            let block_data = &bytes[begin..end];
            if block_data.len() == READ_SIZE {
                cache.write_u8(1)?;
            } else {
                cache.write_u8(2)?;
                cache.write_u64::<byteorder::BigEndian>(block_data.len() as u64)?;
            }
            cache.write_all(block_data)?;
        }
        let len = bytes.len() as u64;
        self.current = Some((block_start, Cursor::new(bytes)));
        Ok(len)
    }
}

impl Read for RemoteFile {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut remaining_buf = buf;
        let mut total_read = 0;
        loop {
            // At this point there's a few cases to consider:
            // 1) We have not read the current block + maybe extra into memory
            // 2) We have read the complete current block (+ maybe extra) into
            //    memory.
            // 3) Whatever is left in the current memory is leftover from a
            //    a previous read, but it's enough.
            // 4) Whatever is left in the current memory is leftover from a
            //    a previous read, and it's not enough.
            let reset_cursor = |this: &mut Self| -> io::Result<u64> {
                let cursor_start = (this.current_position / READ_SIZE as u64) * READ_SIZE as u64;
                let in_block = this.current_position - cursor_start;
                // If we not at the start of the block, then the length that we need
                // is longer than the length of the buf itself, since we have to
                // acount for how far into the block we are.
                let bytes_available = this.read_current_block(remaining_buf.len() as u64)?;
                // If we are not at the beginning of a block, then skip to where
                // we need to be.
                if in_block > 0 {
                    this.current
                        .as_mut()
                        .unwrap()
                        .1
                        .seek(SeekFrom::Start(in_block))?;
                }
                Ok(bytes_available - in_block.min(bytes_available))
            };
            let bytes_available = match self.current.as_ref() {
                None => reset_cursor(self)?,
                Some((_, cursor)) => {
                    let bytes_in_cursor = cursor.get_ref().len() as u64;
                    let cursor_position = cursor.position();
                    let bytes_available = bytes_in_cursor - cursor_position;
                    // We don't have enough bytes in the cursor. Let's reload this
                    // block just to ensure that we have the data loaded.
                    if bytes_available < remaining_buf.len() as u64 {
                        reset_cursor(self)?
                    } else {
                        bytes_available
                    }
                }
            };
            let read = self.current.as_mut().unwrap().1.read(remaining_buf)?;
            self.current_position += read as u64;
            total_read += read;
            if read == 0 || read == remaining_buf.len() || read == bytes_available as usize {
                break;
            }
            let cursor_start = (self.current_position / READ_SIZE as u64) * READ_SIZE as u64;
            let in_block = self.current_position - cursor_start;
            let remaining_in_block = READ_SIZE - in_block as usize;
            // If we didn't read everything, we *must* have at least read until
            // the end of the block
            assert!(read >= remaining_in_block);
            remaining_buf = &mut remaining_buf[read..];
        }
        Ok(total_read)
    }
}

impl Seek for RemoteFile {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.current_position = match pos {
            SeekFrom::Start(s) => s,
            SeekFrom::End(_) => unimplemented!(),
            SeekFrom::Current(s) => {
                if s >= 0 {
                    self.current_position + (s as u64)
                } else {
                    if self.current_position < s.checked_neg().unwrap() as u64 {
                        panic!("Seeked to <0");
                    }
                    self.current_position - s.checked_neg().unwrap() as u64
                }
            }
        };
        if let Some((cursor_start, cursor)) = self.current.as_mut() {
            let cursor_end = *cursor_start + READ_SIZE as u64;
            if *cursor_start <= self.current_position && self.current_position < cursor_end {
                let new_position = self.current_position - *cursor_start;
                cursor.set_position(new_position);
                return Ok(self.current_position);
            }
            self.current = None;
        }
        Ok(self.current_position)
    }
}

impl Clone for RemoteFile {
    fn clone(&self) -> Self {
        RemoteFile {
            url: self.url.clone(),
            current_position: 0,
            current: None,
            cache: None,
        }
    }
}

impl Reopen for RemoteFile {
    fn reopen(&self) -> io::Result<RemoteFile> {
        Ok(RemoteFile {
            url: self.url.clone(),
            current_position: 0,
            current: None,
            cache: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bbi::{BigBedRead, BigWigRead};

    #[ignore]
    #[test]
    fn test_remote() {
        let f = RemoteFile::new("https://encode-public.s3.amazonaws.com/2020/01/17/7d2573b1-86f4-4592-a68a-ac3d5d0372d6/ENCFF592UJG.bigBed");
        let mut remote = BigBedRead::open(f).unwrap();

        let remote_intervals: Vec<_> = remote
            .get_interval("chr10", 100000000, 100010000)
            .unwrap()
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(remote_intervals.len(), 5);
    }

    #[ignore]
    #[test]
    fn test_remote2() {
        let f = RemoteFile::new("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig");
        let mut remote = BigWigRead::open(f).unwrap();

        let interval = remote.get_zoom_interval("chr17", 0, 36996442, 2048);
        let _: Vec<_> = interval.unwrap().collect();
    }

    #[ignore]
    #[test]
    fn test_remote3() {
        let f = RemoteFile::new("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig");
        let mut remote = BigWigRead::open(f).unwrap();

        let interval = remote.get_zoom_interval("chr2", 46087592, 174087320, 32768);
        let _: Vec<_> = interval.unwrap().collect();
    }

    #[ignore]
    #[test]
    fn test_remote4() {
        let f = RemoteFile::new("https://proteinpaint.stjude.org/ppdemo/hg19/bigwig/temp.bw");
        let remote = BigWigRead::open(f).unwrap();

        let _: Vec<_> = remote
            .get_interval_move("chr1", 169253475, 169257278)
            .unwrap()
            .collect();
    }
}
