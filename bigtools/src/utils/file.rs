pub mod file_view;
pub mod reopen;
pub mod streaming_linereader;
pub mod tell;
#[cfg(feature = "write")]
pub mod tempfilebuffer;

#[cfg(feature = "remote")]
pub mod remote_file;

use std::fs::File;
use std::io::{self, BufRead, Seek};

/// Split a file into n `chunks` of roughly equal sizes. The file is split into
/// roughly equal sizes (by bytes), then each chunk is extended to the end of
/// the last line in that chunk. Importantly, this doesn't split the file into
/// equal number of lines in each chunk. Also, the number of chunks may be lower
/// than the provided number (which, is most obvious when each chunk is only a single line).
pub fn split_file_into_chunks_by_size(f: File, chunks: u64) -> io::Result<Vec<(u64, u64)>> {
    let file_size = f.metadata()?.len();
    let mut file_reader = io::BufReader::new(f);
    let chunk_size = file_size / chunks;
    let mut chunk_vec = Vec::with_capacity(chunks as usize);
    let mut chunk_start = 0;
    let mut chunk_end = chunk_size;
    loop {
        file_reader.seek(io::SeekFrom::Start(chunk_end))?;
        file_reader.read_line(&mut String::new())?;
        let line_end = file_reader.seek(io::SeekFrom::Current(0))?;
        chunk_end = line_end;
        chunk_vec.push((chunk_start, chunk_end));
        (chunk_start, chunk_end) = (
            chunk_end,
            chunk_end.max(chunk_start + chunk_size + chunk_size),
        );
        chunk_end = chunk_end.min(file_size);

        if chunk_start >= file_size {
            break;
        }
    }

    Ok(chunk_vec)
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use super::*;

    #[test]
    fn test_split_file_into_chunks_by_size() {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        dir.push("small2.bed");

        let f = File::open(dir.clone()).unwrap();
        let chunks = split_file_into_chunks_by_size(f, 1).unwrap();
        assert_eq!(chunks, vec![(0, 3103)]);

        let f = File::open(dir.clone()).unwrap();
        let chunks = split_file_into_chunks_by_size(f, 2).unwrap();
        assert_eq!(chunks, vec![(0, 1614), (1614, 3103)]);

        let f = File::open(dir.clone()).unwrap();
        let chunks = split_file_into_chunks_by_size(f, 3).unwrap();
        assert_eq!(chunks, vec![(0, 1056), (1056, 2112), (2112, 3103)]);

        let f = File::open(dir.clone()).unwrap();
        let chunks = split_file_into_chunks_by_size(f, 200).unwrap();
        assert_eq!(chunks.len(), 50);
    }
}
