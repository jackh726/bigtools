use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::io::{self, Read, Seek, SeekFrom};

pub const CACHE_SIZE: usize = 4 * 1024;

/// Wraps a `Read` and a `HashMap<u64, [u8; CACHE_SIZE]` to provide in-memory
/// caching of file contents.
///
/// An important design consideration here is that both the `Read` and the
/// `HashMap` cache are mutable refs. It is designed this way in order to
/// accommodate only caching during a particular bit of logic. For example,
/// for BigWigs, it makes sense to cache the initial header and the cir tree,
/// without caching the data section.
pub struct MemCachedRead<'a, R: Read + Seek> {
    reader: &'a mut R,
    cache: &'a mut HashMap<u64, [u8; CACHE_SIZE]>,
    current_position: Option<u64>,
}

impl<'a, R: Read + Seek> MemCachedRead<'a, R> {
    pub(crate) fn new(reader: &'a mut R, cache: &'a mut HashMap<u64, [u8; CACHE_SIZE]>) -> Self {
        MemCachedRead {
            reader,
            cache,
            current_position: None,
        }
    }
}

impl<R: Read + Seek> Read for MemCachedRead<'_, R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut remaining_buf = buf;
        let mut total_read = 0;
        loop {
            if remaining_buf.is_empty() {
                break;
            }
            let current_position = match self.current_position {
                Some(current_position) => current_position,
                None => self.seek(SeekFrom::Current(0))?,
            };
            let block_start = (current_position / CACHE_SIZE as u64) * CACHE_SIZE as u64;
            let offset = current_position - block_start;
            match self.cache.entry(block_start) {
                Entry::Occupied(entry) => {
                    let cached = entry.get();
                    let bytes = &cached
                        [offset as usize..CACHE_SIZE.min(remaining_buf.len() + offset as usize)];
                    total_read += bytes.len();
                    remaining_buf[..bytes.len()].copy_from_slice(bytes);
                    remaining_buf = &mut remaining_buf[bytes.len()..];
                    self.current_position = Some(current_position + bytes.len() as u64);
                }
                Entry::Vacant(_) => {
                    let block_end = (((current_position as f64 + remaining_buf.len() as f64)
                        / CACHE_SIZE as f64)
                        .ceil() as u64)
                        * CACHE_SIZE as u64
                        - block_start;
                    self.reader.seek(SeekFrom::Start(block_start))?;
                    let mut temp_buf = vec![0u8; block_end as usize];
                    let read = self.reader.read(&mut temp_buf)?;
                    let read_buf = &temp_buf[..read];
                    for (idx, chunk) in read_buf.chunks_exact(CACHE_SIZE).enumerate() {
                        let chunk_start = idx * CACHE_SIZE + block_start as usize;
                        let mut cached = [0u8; CACHE_SIZE];
                        cached.copy_from_slice(chunk);
                        self.cache.insert(chunk_start as u64, cached);
                    }
                    let remaining_len = remaining_buf.len();
                    if read > offset as usize {
                        let copy_start = offset as usize;
                        let copy_end = read.min(remaining_len + offset as usize);
                        let copy_to_end = (read - offset as usize).min(remaining_len);
                        self.current_position = Some(current_position + copy_to_end as u64);
                        let to_copy = &temp_buf[copy_start..copy_end];
                        remaining_buf[..copy_to_end].copy_from_slice(to_copy);
                        total_read += to_copy.len();
                    }
                    break;
                }
            }
        }
        Ok(total_read)
    }
}

impl<R: Read + Seek> Seek for MemCachedRead<'_, R> {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        match self.current_position.as_mut() {
            Some(last_position) => {
                match pos {
                    SeekFrom::Start(s) => *last_position = s,
                    SeekFrom::End(_) => *last_position = self.reader.seek(pos)?,
                    SeekFrom::Current(s) => {
                        if s >= 0 {
                            *last_position += s as u64
                        } else {
                            if *last_position < s.checked_neg().unwrap() as u64 {
                                panic!("Seeked to <0");
                            }
                            *last_position -= s.checked_neg().unwrap() as u64;
                        }
                    }
                };
                Ok(*last_position)
            }
            None => {
                let position = self.reader.seek(pos)?;
                self.current_position = Some(position);
                Ok(position)
            }
        }
    }
}

impl<R: Read + Seek> Drop for MemCachedRead<'_, R> {
    fn drop(&mut self) {
        if let Some(current_position) = self.current_position.as_ref() {
            let _ = self.reader.seek(SeekFrom::Start(*current_position));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Cursor, Write};

    #[test]
    fn test_mem_cached_file() -> io::Result<()> {
        let mut data = Cursor::new(vec![0u8; CACHE_SIZE * 16]);
        let test_data = &[1, 2, 3, 4, 5, 6, 7, 8];
        data.seek(SeekFrom::Start(150))?;
        data.write(test_data)?;
        data.seek(SeekFrom::Start(CACHE_SIZE as u64 * 2 - 4))?;
        data.write(test_data)?;
        data.seek(SeekFrom::Start(CACHE_SIZE as u64 * 7 - 4))?;
        data.write(&vec![1u8; CACHE_SIZE * 5])?;
        data.seek(SeekFrom::Start(0))?;
        let mut cache = HashMap::new();
        let mut mem_cached_file = MemCachedRead::new(&mut data, &mut cache);
        let mut test = vec![0u8; 8];

        // Read from the middle of a block
        mem_cached_file.seek(SeekFrom::Start(150))?;
        mem_cached_file.read_exact(&mut test)?;
        assert_eq!(&test, test_data);
        // Read across a block boundary
        mem_cached_file.seek(SeekFrom::Start(CACHE_SIZE as u64 * 2 - 4))?;
        mem_cached_file.read_exact(&mut test)?;
        assert_eq!(&test, test_data);
        // Read first into cache, then second from cache + across multipe blocks (cached + not cached)
        mem_cached_file.seek(SeekFrom::Start(CACHE_SIZE as u64 * 7 - 4))?;
        mem_cached_file.read_exact(&mut vec![0u8; CACHE_SIZE + 1000])?;
        let mut test = vec![0u8; CACHE_SIZE * 5];
        mem_cached_file.seek(SeekFrom::Start(CACHE_SIZE as u64 * 7 - 4))?;
        mem_cached_file.read_exact(&mut test)?;
        let sum: usize = test.iter().map(|i| *i as usize).sum();
        assert_eq!(sum, CACHE_SIZE * 5);
        // FIXME: test repeated reads

        Ok(())
    }
}
