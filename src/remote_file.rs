use std::fs::File;
use std::io::{self, Cursor, Read, Seek, SeekFrom, Write};

use tempfile;

use crate::seekableread::Reopen;

const READ_SIZE: usize = 20 * 1024; // 20 KB chunks

/// Remote file reads are cached to a temporary file. The size of each block is
/// `READ_SIZE + 1` bytes. The first byte of a block is `0` if the data hasn't
/// been written yet, or `1` if it has. Additionally, the data is stored by
/// tiling, so all or no data for a given block is available.

pub struct RemoteFile {
    url: String,
    last_seek: u64,
    current: Option<Cursor<Vec<u8>>>,
    cache: Option<File>,
}

impl RemoteFile {
    pub fn new(url: &str) -> RemoteFile {
        RemoteFile {
            url: url.to_string(),
            last_seek: 0,
            current: None,
            cache: None,
        }
    }
}

impl RemoteFile {
    fn block(&self) -> u64 {
        self.last_seek / READ_SIZE as u64
    }

    fn read_current_block(&mut self, recommended_size: Option<u64>) -> io::Result<u64> {
        let block = self.block();
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
            self.current = Some(Cursor::new(bytes.to_vec()));
            return Ok(READ_SIZE as u64);
        } else if status == 2 {
            let bytes_available = cache.read_u64::<byteorder::BigEndian>()?;
            let mut bytes = vec![0u8; bytes_available as usize];
            cache.read_exact(&mut bytes)?;
            self.current = Some(Cursor::new(bytes.to_vec()));
            return Ok(bytes_available);
        }

        let read_len = recommended_size
            .unwrap_or(READ_SIZE as u64)
            .max(READ_SIZE as u64);

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
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
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
        self.current = Some(Cursor::new(bytes));
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
            let bytes_available = if self.current.is_none() {
                let cursor_start = (self.last_seek / READ_SIZE as u64) * READ_SIZE as u64;
                let in_block = self.last_seek - cursor_start;
                // If we not at the start of the block, then the length that we need
                // is longer than the length of the buf itself, since we have to
                // acount for how far into the block we are.
                let bytes_available =
                    self.read_current_block(Some(in_block + remaining_buf.len() as u64))?;
                // If we are at the beginning of a block, then skip to where
                // we need to be.
                if in_block > 0 {
                    self.current
                        .as_mut()
                        .unwrap()
                        .seek(SeekFrom::Start(in_block))?;
                }
                bytes_available - in_block.min(bytes_available)
            } else {
                let bytes_in_cursor = self.current.as_ref().unwrap().get_ref().len() as u64;
                let cursor_position = self.current.as_ref().unwrap().position();
                let bytes_available = bytes_in_cursor - cursor_position;
                // We don't have enough bytes in the cursor. Let's reload this
                // block just to ensure that we have the data loaded.
                if bytes_available < remaining_buf.len() as u64 {
                    self.last_seek = ((self.last_seek / READ_SIZE as u64) * READ_SIZE as u64)
                        + self.current.as_ref().unwrap().position();
                    self.current = None;
                    let cursor_start = (self.last_seek / READ_SIZE as u64) * READ_SIZE as u64;
                    let in_block = self.last_seek - cursor_start;
                    let bytes_available =
                        self.read_current_block(Some(in_block + remaining_buf.len() as u64))?;
                    if in_block > 0 {
                        self.current
                            .as_mut()
                            .unwrap()
                            .seek(SeekFrom::Start(in_block))?;
                    }
                    bytes_available - in_block.min(bytes_available)
                } else {
                    bytes_available
                }
            };
            let read = self.current.as_mut().unwrap().read(remaining_buf)?;
            total_read += read;
            if read == 0 || read == remaining_buf.len() || read == bytes_available as usize {
                break;
            }
            let cursor_start = (self.last_seek / READ_SIZE as u64) * READ_SIZE as u64;
            let in_block = self.last_seek - cursor_start;
            let remaining_in_block = READ_SIZE - in_block as usize;
            // If we didn't read everything, we *must* have at least read until
            // the end of the block
            assert!(read >= remaining_in_block);
            // We don't want to skip into the middle of a chunk
            let skip = remaining_in_block + ((read - remaining_in_block) / READ_SIZE) * READ_SIZE;
            remaining_buf = &mut remaining_buf[skip..];
            total_read -= read - skip;
            self.current = None;
            self.last_seek += skip as u64;
        }
        Ok(total_read)
    }
}

impl Seek for RemoteFile {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        // If we seek outside the current block, we have to invalidate.
        let last_seek = self.last_seek;
        self.last_seek = match pos {
            SeekFrom::Start(s) => s,
            SeekFrom::End(_) => unimplemented!(),
            SeekFrom::Current(s) => {
                if s >= 0 {
                    last_seek + (s as u64)
                } else {
                    if last_seek < s.checked_neg().unwrap() as u64 {
                        panic!("Seeked to <0");
                    }
                    last_seek - s.checked_neg().unwrap() as u64
                }
            }
        };
        if let Some(_cursor) = self.current.as_mut() {
            /*
            // FIXME: this isn't correct
            let cursor_end = cursor.get_ref().len() as u64;
            if self.last_seek >= last_seek && self.last_seek < cursor_end {
                let new_position = self.last_seek - last_seek;
                cursor.set_position(new_position);
                return Ok(self.last_seek);
            }
            */
            self.current = None;
        }
        Ok(self.last_seek)
    }
}

impl Clone for RemoteFile {
    fn clone(&self) -> Self {
        RemoteFile {
            url: self.url.clone(),
            last_seek: 0,
            current: None,
            cache: None,
        }
    }
}

impl Reopen<RemoteFile> for RemoteFile {
    fn reopen(&self) -> io::Result<RemoteFile> {
        Ok(RemoteFile {
            url: self.url.clone(),
            last_seek: 0,
            current: None,
            cache: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bigbedread::BigBedRead;

    #[ignore]
    #[test]
    fn test_remote() {
        let f = RemoteFile::new("https://encode-public.s3.amazonaws.com/2020/01/17/7d2573b1-86f4-4592-a68a-ac3d5d0372d6/ENCFF592UJG.bigBed");
        let mut remote = BigBedRead::from(f).unwrap();

        let remote_intervals: Vec<_> = remote
            .get_interval("chr10", 100000000, 100010000)
            .unwrap()
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(remote_intervals.len(), 5);
    }
}
