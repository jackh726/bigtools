use std::io::{self, BufRead};

#[derive(Debug)]
pub struct StreamingLineReader<B: BufRead + std::fmt::Debug> {
    current_line: String,
    buf_read: B,
}

impl<B: BufRead + std::fmt::Debug> StreamingLineReader<B> {
    pub fn new(bf: B) -> StreamingLineReader<B> {
        return StreamingLineReader {
            current_line: String::new(),
            buf_read: bf,
        };
    }

    pub fn read<'a>(&'a mut self) -> io::Result<Option<&'a str>> {
        self.current_line.clear();
        let size = self.buf_read.read_line(&mut self.current_line)?;
        if size == 0 {
            return Ok(None);
        }
        return Ok(Some(&self.current_line));
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::{self, BufReader};
    extern crate test;

    #[test]
    fn test_works() -> io::Result<()> {
        let f = File::open("/home/hueyj/temp/final.min.chr17.bedGraph")?;
        let bf = BufReader::new(f);
        let mut slr = StreamingLineReader::new(bf);
        for _i in 0..=1 {
            let line = slr.read()?;
            if let Some(l) = line {
                println!("Line: {:?}", l);
                let mut split = l.split_whitespace();
                let chrom = split.next().expect("Missing chrom").to_owned();
                let start = split.next().expect("Missing start").parse::<u32>().unwrap();
                let end = split.next().expect("Missing end").parse::<u32>().unwrap();
                let value = split.next().expect("Missing value").parse::<f32>().unwrap();
                println!("{:?} {:?} {:?} {:?}", chrom, start, end, value);
            }
        }
        Ok(())
    }

}
