use std::io::{self, BufRead};

#[derive(Debug)]
pub struct StreamingLineReader<B> {
    current_line: String,
    buf_read: B,
}

impl<B: BufRead> StreamingLineReader<B> {
    pub fn new(bf: B) -> StreamingLineReader<B> {
        StreamingLineReader {
            current_line: String::new(),
            buf_read: bf,
        }
    }

    pub fn read(&mut self) -> Option<io::Result<&'_ str>> {
        self.current_line.clear();
        match self.buf_read.read_line(&mut self.current_line) {
            Ok(size) if size == 0 => None,
            Ok(_) => Some(Ok(self.current_line.trim_end())),
            Err(e) => Some(Err(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::{self, BufReader};
    use std::path::PathBuf;

    #[test]
    fn test_works() -> io::Result<()> {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        dir.push("small.bedGraph");
        let f = File::open(dir)?;
        let bf = BufReader::new(f);
        let mut slr = StreamingLineReader::new(bf);
        assert_eq!("chr17\t1\t100\t0.5", slr.read().unwrap().unwrap());
        assert_eq!("chr17\t101\t200\t0.5", slr.read().unwrap().unwrap());
        assert_eq!("chr17\t201\t300\t0.5", slr.read().unwrap().unwrap());
        assert_eq!("chr18\t1\t100\t0.5", slr.read().unwrap().unwrap());
        assert_eq!("chr18\t101\t200\t0.5", slr.read().unwrap().unwrap());
        assert_eq!("chr19\t1\t100\t0.5", slr.read().unwrap().unwrap());
        assert!(slr.read().is_none());
        Ok(())
    }
}
