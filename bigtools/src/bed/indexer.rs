/*!
Utilities for indexing a bed file that is start-sorted (chromosomes may be out of order).
*/
use std::fs::File;
use std::io::{self, BufRead, BufReader, Seek, SeekFrom};

use index_list::{IndexList, ListIndex};

use crate::utils::tell::Tell;

/// Returns a Vec of offsets into a bed file, and the chromosome starting at each offset.
pub fn index_chroms(file: File) -> io::Result<Option<Vec<(u64, String)>>> {
    let mut file = BufReader::new(file);

    let mut line = String::new();

    file.read_line(&mut line)?;

    if line.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Empty file".to_string(),
        ));
    }

    let mut chroms = IndexList::new();

    fn parse_line(s: &str) -> Result<Option<String>, io::Error> {
        if s.is_empty() {
            return Ok(None);
        }
        let mut split = s.trim_end().splitn(4, '\t');
        let chrom = split
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Empty file".to_string()))?;
        let s = split.next().ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Missing start: {:}", s))
        })?;
        let _start = s.parse::<u32>().map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Invalid start: {:}", s))
        })?;
        let s = split.next().ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Missing end: {:}", s))
        })?;
        let _end = s.parse::<u32>().map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Invalid end: {:}", s))
        })?;

        Ok(Some(chrom.to_string()))
    }

    let chrom = parse_line(&line)?.unwrap();
    let first = chroms.insert_first((0, chrom));
    let file_size = file.seek(SeekFrom::End(0))?;

    fn do_index(
        file_size: u64,
        file: &mut BufReader<File>,
        chroms: &mut IndexList<(u64, String)>,
        line: &mut String,
        prev: ListIndex,
        next: Option<ListIndex>,
        limit: usize,
    ) -> Result<(), io::Error> {
        if limit == 0 {
            panic!("Recursive depth limit reached");
        }

        let next_tell = next
            .map(|next| chroms.get(next).unwrap().0)
            .unwrap_or(file_size);
        let mid = (next_tell + chroms.get(prev).unwrap().0) / 2;
        file.seek(SeekFrom::Start(mid))?;
        file.read_line(line)?;
        line.clear();
        let tell = file.tell()?;
        file.read_line(line)?;
        let chrom = parse_line(&*line)?;
        let chrom = match chrom {
            Some(chrom) => chrom,
            None => return Ok(()),
        };

        // There are three options:
        // 1) The chrom is the same as the previous one. We need to index
        //    between the current and next, since they must be different.
        // 2) The chrom is the same as the next one. This means that we need to
        //    continuing indexing between the previous index and the current,
        //    but not between the current and next index. There is one exceptional
        //    case to think about, when we're *at or passed* the "next" index:
        //    |1----|2----|2----|
        //    ^           ^      Indexed
        //            ^          Mid
        //                ^      New
        //    Here, the "new" index is the same as the previous last. It's
        //    hopefully clear that we one more cycle of indexing between the
        //    start and mid will mark the second line (first line of 2) as the
        //    start correctly.
        // 3) The chrom is different from both the previous and next. We need
        //    to continue to index between the previous and current as well as
        //    between the current and next.

        let curr = chroms.insert_after(prev, (tell, chrom));

        let left = chroms.get(curr).unwrap().1 != chroms.get(prev).unwrap().1 && tell < next_tell;
        let right = next
            .map(|next| {
                chroms.get(curr).unwrap().1 != chroms.get(next).unwrap().1
                    && tell < chroms.get(next).unwrap().0
            })
            .unwrap_or(true);

        if left {
            do_index(file_size, file, chroms, line, prev, Some(curr), limit - 1)?;
        }

        if right {
            do_index(file_size, file, chroms, line, curr, next, limit - 1)?;
        }

        if chroms.get(curr).unwrap().1 != chroms.get(prev).unwrap().1 && tell == next_tell {
            file.seek(SeekFrom::Start(chroms.get(prev).unwrap().0))?;
            line.clear();
            file.read_line(line)?;
            line.clear();
            let tell = file.tell()?;
            file.read_line(line)?;
            let chrom = parse_line(&*line)?
                .expect("Bad logic. Must at least find last entry for chromosome.");
            chroms.insert_after(prev, (tell, chrom));
        }
        Ok(())
    }

    do_index(
        file_size,
        &mut file,
        &mut chroms,
        &mut line,
        first,
        None,
        100,
    )?;

    let mut chroms: Vec<_> = chroms.drain_iter().collect();
    chroms.dedup_by_key(|index| index.1.clone());
    let mut deduped_chroms = chroms.clone();
    deduped_chroms.sort();
    deduped_chroms.dedup_by_key(|index| index.1.clone());
    if chroms.len() != deduped_chroms.len() {
        return Ok(None);
    }

    Ok(Some(chroms))
}

#[allow(unused)]
fn linear_index(reader: &mut BufReader<File>) -> io::Result<Vec<(u64, String)>> {
    reader.seek(SeekFrom::Start(0))?;
    let mut chroms: Vec<(u64, String)> = vec![];
    loop {
        let mut line = String::new();
        let tell = reader.seek(SeekFrom::Current(0))?;
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            break;
        }

        let chrom = line
            .split_ascii_whitespace()
            .next()
            .ok_or(io::Error::from(io::ErrorKind::InvalidData))?;
        if chroms.last().map(|c| &c.1 != chrom).unwrap_or(true) {
            chroms.push((tell, chrom.to_string()));
        }
        line.clear();
    }
    Ok(chroms)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io;
    use std::path::PathBuf;

    #[test]
    fn test_index() -> io::Result<()> {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        //dir.push("small.bed");
        dir.push("multi_chrom.bedGraph");

        let f = File::open(dir.clone())?;
        let mut f = BufReader::new(f);
        let mut line = String::new();
        let mut chroms = vec![];
        loop {
            let line_start = f.tell()?;
            match f.read_line(&mut line) {
                Ok(size) if size == 0 => break,
                Ok(_) => {
                    let mut split = line.splitn(4, '\t');
                    let chrom = split.next().unwrap().to_string();
                    match chroms.last() {
                        None => chroms.push((line_start, chrom.clone())),
                        Some(last_chrom) if last_chrom.1 != chrom => {
                            chroms.push((line_start, chrom.clone()));
                        }
                        _ => {}
                    }
                }
                Err(e) => return Err(e),
            }
            line.clear();
        }

        let f = File::open(dir)?;
        let indexed_chroms = index_chroms(f)?.unwrap();
        assert_eq!(chroms, indexed_chroms);

        Ok(())
    }
}
