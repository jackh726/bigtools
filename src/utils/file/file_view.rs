use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};

pub struct FileView {
    file: File,
    start: u64,
    end: u64,
    current: Option<u64>,
}

impl FileView {
    pub fn new(mut file: File, start: u64, end: u64) -> io::Result<FileView> {
        let file_end = file.seek(io::SeekFrom::End(0))?;
        file.seek(io::SeekFrom::Start(start))?;
        let end = end.min(file_end);
        Ok(FileView {
            file,
            start,
            end,
            current: Some(start),
        })
    }
}

impl Read for FileView {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let current = match self.current {
            Some(current) => current,
            None => {
                let current = self.seek(io::SeekFrom::Current(0))? + self.start;
                self.current = Some(current);
                current
            }
        };
        let to_read = buf.len().min((self.end - current) as usize);
        let buf = &mut buf[..to_read];
        match self.file.read(buf) {
            Ok(read) => {
                self.current = Some(current + read as u64);
                Ok(read)
            }
            Err(e) => {
                self.current = None;
                Err(e)
            }
        }
    }
}

impl Seek for FileView {
    fn seek(&mut self, pos: io::SeekFrom) -> io::Result<u64> {
        match pos {
            SeekFrom::Start(start) => {
                let seek_from = io::SeekFrom::Start(self.end.min(self.start + start));
                match self.file.seek(seek_from) {
                    Ok(new_pos) => {
                        assert!(new_pos >= self.start && new_pos <= self.end);
                        self.current = Some(new_pos);
                        let new_pos = new_pos - self.start;
                        Ok(new_pos)
                    }
                    Err(e) => {
                        self.current = None;
                        Err(e)
                    }
                }
            }
            SeekFrom::End(end) => {
                let end = end.min(0);
                let new_pos = (self.end as i64) + end;
                match self.file.seek(io::SeekFrom::Start(new_pos.max(0) as u64)) {
                    Ok(new_pos) => {
                        assert!(new_pos >= self.start && new_pos <= self.end);
                        self.current = Some(new_pos);
                        let new_pos = new_pos - self.start;
                        Ok(new_pos)
                    }
                    Err(e) => {
                        self.current = None;
                        Err(e)
                    }
                }
            }
            SeekFrom::Current(offset) => {
                let current = match self.current {
                    Some(current) => current,
                    None => {
                        let current = self.seek(io::SeekFrom::Current(0))? + self.start;
                        self.current = Some(current);
                        current
                    }
                };

                let new_pos = (current as i64) + offset;
                let new_pos = new_pos.min(self.end as i64).max(self.start as i64);
                match self.file.seek(io::SeekFrom::Start(new_pos as u64)) {
                    Ok(new_pos) => {
                        assert!(new_pos >= self.start && new_pos <= self.end);
                        self.current = Some(new_pos);
                        let new_pos = new_pos - self.start;
                        Ok(new_pos)
                    }
                    Err(e) => {
                        self.current = None;
                        Err(e)
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use std::io::{BufRead, BufReader};
    use std::path::PathBuf;

    use super::*;

    #[test]
    fn test_file_view_lines() {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        dir.push("small2.bed");

        let f = File::open(dir.clone()).unwrap();
        let file_len = f.metadata().unwrap().len();
        let view = FileView::new(f, 0, file_len).unwrap();
        let mut buf_read = BufReader::new(view);
        let mut str = String::new();

        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(
            &str,
            "chr5\t65724845\t65725305\t.\t689\t.\t11.08521\t-1.00000\t0.50517\t230\n"
        );

        buf_read.seek(io::SeekFrom::Start(1117)).unwrap();
        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(
            &str,
            "chr5\t76372720\t76373180\t.\t870\t.\t11.19471\t-1.00000\t0.51655\t230\n"
        );

        buf_read.seek(io::SeekFrom::End(-1366)).unwrap();
        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(
            &str,
            "chr1\t113141313\t113141773\t.\t562\t.\t11.24802\t-1.00000\t0.52177\t230\n"
        );

        buf_read.seek(io::SeekFrom::Start(1737)).unwrap();
        buf_read.seek(io::SeekFrom::Current(-620)).unwrap();
        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(
            &str,
            "chr5\t76372720\t76373180\t.\t870\t.\t11.19471\t-1.00000\t0.51655\t230\n"
        );

        buf_read.seek(io::SeekFrom::End(50)).unwrap();
        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(&str, "");

        let f = File::open(dir.clone()).unwrap();
        let view = FileView::new(f, 0, 1366).unwrap();
        let mut buf_read = BufReader::new(view);
        let mut str = String::new();

        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(
            &str,
            "chr5\t65724845\t65725305\t.\t689\t.\t11.08521\t-1.00000\t0.50517\t230\n"
        );

        buf_read.seek(io::SeekFrom::Start(1117)).unwrap();
        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(
            &str,
            "chr5\t76372720\t76373180\t.\t870\t.\t11.19471\t-1.00000\t0.51655\t230\n"
        );

        buf_read.seek(io::SeekFrom::Start(1366)).unwrap();
        buf_read.seek(io::SeekFrom::Current(-249)).unwrap();
        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(
            &str,
            "chr5\t76372720\t76373180\t.\t870\t.\t11.19471\t-1.00000\t0.51655\t230\n"
        );

        buf_read.seek(io::SeekFrom::End(-1500)).unwrap();
        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(
            &str,
            "chr5\t65724845\t65725305\t.\t689\t.\t11.08521\t-1.00000\t0.50517\t230\n"
        );

        buf_read.seek(io::SeekFrom::End(50)).unwrap();
        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(&str, "");

        let f = File::open(dir.clone()).unwrap();
        let view = FileView::new(f, 1117, file_len).unwrap();
        let mut buf_read = BufReader::new(view);
        let mut str = String::new();

        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(
            &str,
            "chr5\t76372720\t76373180\t.\t870\t.\t11.19471\t-1.00000\t0.51655\t230\n"
        );

        buf_read.seek(io::SeekFrom::Start(620)).unwrap();
        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(
            &str,
            "chr1\t113141313\t113141773\t.\t562\t.\t11.24802\t-1.00000\t0.52177\t230\n"
        );

        buf_read.seek(io::SeekFrom::End(-1366)).unwrap();
        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(
            &str,
            "chr1\t113141313\t113141773\t.\t562\t.\t11.24802\t-1.00000\t0.52177\t230\n"
        );

        buf_read.seek(io::SeekFrom::End(50)).unwrap();
        str.clear();
        let _ = buf_read.read_line(&mut str).unwrap();
        assert_eq!(&str, "");
    }

    #[test]
    fn test_file_view_bufs() {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        dir.push("small2.bed");

        let f = File::open(dir.clone()).unwrap();
        let file_len = f.metadata().unwrap().len();

        let mut view = FileView::new(f, 0, file_len).unwrap();

        assert_eq!(
            view.read(&mut vec![0; file_len as usize]).unwrap(),
            file_len as usize
        );

        view.seek(io::SeekFrom::Start(0)).unwrap();
        assert_eq!(
            view.read(&mut vec![0; file_len as usize + 500]).unwrap(),
            file_len as usize
        );

        view.seek(io::SeekFrom::Start(0)).unwrap();
        assert_eq!(view.read(&mut vec![0; 500]).unwrap(), 500);

        let f = File::open(dir.clone()).unwrap();
        let mut view = FileView::new(f, 500, file_len).unwrap();

        assert_eq!(
            view.read(&mut vec![0; file_len as usize]).unwrap(),
            file_len as usize - 500
        );

        view.seek(io::SeekFrom::Start(0)).unwrap();
        assert_eq!(
            view.read(&mut vec![0; file_len as usize + 500]).unwrap(),
            file_len as usize - 500
        );

        view.seek(io::SeekFrom::Start(0)).unwrap();
        assert_eq!(view.read(&mut vec![0; 500]).unwrap(), 500);

        let f = File::open(dir.clone()).unwrap();
        let mut view = FileView::new(f, 500, 1000).unwrap();

        assert_eq!(view.read(&mut vec![0; file_len as usize]).unwrap(), 500);

        view.seek(io::SeekFrom::Start(0)).unwrap();
        assert_eq!(view.read(&mut vec![0; 1000]).unwrap(), 500);

        view.seek(io::SeekFrom::Start(0)).unwrap();
        assert_eq!(view.read(&mut vec![0; 250]).unwrap(), 250);
    }
}
