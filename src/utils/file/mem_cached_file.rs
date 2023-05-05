use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::io::{self, Read, Seek, SeekFrom};

pub const CACHE_SIZE: usize = 64;

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
    cache: &'a mut HashMap<usize, [u8; CACHE_SIZE]>,
    current_buffer: Option<(usize, [u8; CACHE_SIZE])>,
    current_position: Option<u64>,
}

impl<'a, R: Read + Seek> MemCachedRead<'a, R> {
    pub(crate) fn new(reader: &'a mut R, cache: &'a mut HashMap<usize, [u8; CACHE_SIZE]>) -> Self {
        MemCachedRead {
            reader,
            cache,
            current_buffer: None,
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
            let block_start = (current_position as usize / CACHE_SIZE) * CACHE_SIZE;
            let block_offset = current_position as usize - block_start;
            match self.current_buffer.as_ref() {
                // If there is a buffer stored and it's the right block, use it
                Some((buffer_block, buffer)) if block_start == *buffer_block => {
                    let bytes =
                        &buffer[block_offset..CACHE_SIZE.min(remaining_buf.len() + block_offset)];
                    total_read += bytes.len();
                    remaining_buf[..bytes.len()].copy_from_slice(bytes);
                    remaining_buf = &mut remaining_buf[bytes.len()..];
                    self.current_position = Some(current_position + bytes.len() as u64);
                    continue;
                }
                // Otherwise, just fall through
                _ => {}
            }

            match self.cache.entry(block_start) {
                Entry::Occupied(entry) => {
                    let cached = entry.get();
                    let bytes =
                        &cached[block_offset..CACHE_SIZE.min(remaining_buf.len() + block_offset)];
                    total_read += bytes.len();
                    remaining_buf[..bytes.len()].copy_from_slice(bytes);
                    remaining_buf = &mut remaining_buf[bytes.len()..];
                    self.current_position = Some(current_position + bytes.len() as u64);
                    self.current_buffer = Some((block_start, *cached));
                }
                Entry::Vacant(_) => {
                    // We want read enough to cover all the cache blocks that
                    // aren't available to cover the remaining buffer.

                    // The last position in the file that we want to read
                    let position_end = current_position + remaining_buf.len() as u64;
                    // The cache block with that position
                    let last_block = (position_end as f64 / CACHE_SIZE as f64) as usize;

                    // The actual end block is the last one where the cache
                    // doesn't contain an entry for the next block
                    let first_block = current_position as usize / CACHE_SIZE;
                    let last_block = (first_block + 1..=last_block)
                        .take_while(|block| !self.cache.contains_key(&(*block * CACHE_SIZE)))
                        .last()
                        .unwrap_or(first_block);

                    // The end block with the last position we want to read
                    // FIXME: Why is this +1
                    let block_end = last_block * CACHE_SIZE + CACHE_SIZE + 1;

                    // The read length is from the beginning of the start block
                    // to the end of the end block
                    let read_length = (block_end - block_start) as usize;

                    // This *might* not be strictly necessary is some cases, but oh well.
                    // I.e., we might not to re-seek, if we've already seeked previously.
                    self.reader.seek(SeekFrom::Start(block_start as u64))?;

                    // Read all the new cache blocks
                    let mut temp_buf = vec![0u8; read_length];
                    let read = try_read_exact(self.reader, &mut temp_buf)?;
                    let read_buf = &temp_buf[..read];

                    // For each *full* chunk that was read, insert that data into
                    // the cache.
                    // FIXME: Should use `slice::array_chunks` when that's stable,
                    // in order to avoid creating the zero-initialized array.
                    // NOTE: We can end up doing extra work if there are cached
                    // blocks after the current block and we are trying to read
                    // into those. The alternative is
                    for (idx, chunk) in read_buf.chunks_exact(CACHE_SIZE).enumerate() {
                        let chunk_start = idx * CACHE_SIZE + block_start;
                        let mut cached = [0u8; CACHE_SIZE];
                        cached.copy_from_slice(chunk);
                        self.cache.insert(chunk_start, cached);
                    }

                    // Now we need to figure out what to copy from `read_buf` to
                    // `remaining_buf`.

                    // If we didn't even read enough to reach the
                    // `current_position`, then there is nothing to do. We will
                    // never be able to read enough.
                    if read <= block_offset {
                        break;
                    }

                    let remaining_len = remaining_buf.len();
                    // We want to copy, starting at the block offset...
                    let copy_start = block_offset;
                    // ...and going to either the end of the read bytes or
                    // to the number of bytes we need:
                    // `(block_offset + remaining_len) - block_offset`
                    let copy_end = read.min(remaining_len + block_offset);

                    let bytes = &read_buf[copy_start..copy_end];
                    total_read += bytes.len();
                    remaining_buf[..bytes.len()].copy_from_slice(bytes);
                    remaining_buf = &mut remaining_buf[bytes.len()..];
                    self.current_position = Some(current_position + bytes.len() as u64);
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

/// Repeatedly calls `read` until either the buffer is filled or `read` returns
/// `0`, ending EOF. This roughly follows the default implementation of
/// `read_exact`, except it does not error on EOF.
fn try_read_exact<R: Read + ?Sized>(this: &mut R, mut buf: &mut [u8]) -> io::Result<usize> {
    let mut total_read = 0;
    while !buf.is_empty() {
        match this.read(buf) {
            Ok(0) => break,
            Ok(n) => {
                total_read += n;
                let tmp = buf;
                buf = &mut tmp[n..];
            }
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(total_read)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Cursor, Write};

    #[test]
    fn test_mem_cached_file() -> io::Result<()> {
        // This test is specifically tailored to the cache size; if it changes,
        // make sure to update the test
        assert_eq!(CACHE_SIZE, 64);

        let mut data = Cursor::new(vec![0u8; CACHE_SIZE * 16]);
        let test_data = &[1, 2, 3, 4, 5, 6, 7, 8];
        // First write to middle of a block
        data.seek(SeekFrom::Start(10)).unwrap();
        data.write(test_data).unwrap();
        // Then write across a block boundary
        data.seek(SeekFrom::Start(CACHE_SIZE as u64 * 2 - 4))
            .unwrap();
        data.write(test_data).unwrap();
        data.seek(SeekFrom::Start(CACHE_SIZE as u64 * 7 - 4))
            .unwrap();
        data.write(&vec![1u8; CACHE_SIZE * 5]).unwrap();
        data.seek(SeekFrom::Start(0)).unwrap();
        let mut cache = HashMap::new();
        let mut mem_cached_file = MemCachedRead::new(&mut data, &mut cache);
        let mut test = vec![0u8; 8];

        // Read from the middle of a block
        mem_cached_file.seek(SeekFrom::Start(10)).unwrap();
        mem_cached_file.read_exact(&mut test).unwrap();
        assert_eq!(&test, test_data);
        // Read across a block boundary
        mem_cached_file
            .seek(SeekFrom::Start(CACHE_SIZE as u64 * 2 - 4))
            .unwrap();
        mem_cached_file.read_exact(&mut test).unwrap();
        assert_eq!(&test, test_data);
        // Read first into cache, then second from cache + across multipe blocks (cached + not cached)
        mem_cached_file
            .seek(SeekFrom::Start(CACHE_SIZE as u64 * 7 - 4))
            .unwrap();
        mem_cached_file
            .read_exact(&mut vec![0u8; CACHE_SIZE * 3])
            .unwrap();
        let mut test = vec![0u8; CACHE_SIZE * 5];
        mem_cached_file
            .seek(SeekFrom::Start(CACHE_SIZE as u64 * 7 - 4))
            .unwrap();
        mem_cached_file.read_exact(&mut test).unwrap();
        let sum: usize = test.iter().map(|i| *i as usize).sum();
        assert_eq!(sum, CACHE_SIZE * 5);
        // FIXME: test repeated reads

        Ok(())
    }
}
