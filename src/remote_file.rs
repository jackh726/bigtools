use std::io::{self, Cursor, Read, Seek, SeekFrom};

use reqwest;

use crate::seekableread::Reopen;


const READ_SIZE: usize = 8 * 1_048_576; // 8 MB chunks

pub struct RemoteFile {
    url: reqwest::Url,
    last_seek: u64,
    current: Option<Cursor<Vec<u8>>>,
}

impl RemoteFile {
    pub fn new(url: &str) -> RemoteFile {
        RemoteFile {
            url: reqwest::Url::parse(url).unwrap(),
            last_seek: 0,
            current: None,
        }
    }
}

impl Read for RemoteFile {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.current.is_none() {
            self.seek(SeekFrom::Start(self.last_seek))?;
        }
        let initial_read = self.current.as_mut().unwrap().read(buf)?;
        let parts_left: usize = (((buf.len() - initial_read) as f32) / (READ_SIZE as f32)).ceil() as usize;
        let initial_pos = self.last_seek;
        let mut total_read = initial_read;
        for part in 0..parts_left {
            let start = initial_read + part * READ_SIZE;
            self.seek(SeekFrom::Start(initial_pos + start as u64))?;
            let read = self.current.as_mut().unwrap().read(&mut buf[start..])?;
            total_read += read;
            if read <= READ_SIZE as usize {
                break
            }
        }
        Ok(total_read)
    }
}

impl Seek for RemoteFile {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.last_seek = match pos {
            SeekFrom::Start(s) => s,
            SeekFrom::End(_) => unimplemented!(),
            SeekFrom::Current(s) => {
                if s == 0 {
                    self.last_seek
                } else if s > 0 {
                    self.last_seek + (s as u64)
                } else {
                    if self.last_seek < s.checked_neg().unwrap() as u64 {
                        panic!("Seeked to <0");
                    }
                    self.last_seek - s.checked_neg().unwrap() as u64
                }
            }
        };
        let client = reqwest::Client::new();
        let mut r = client
            .get(self.url.clone())
            .header(reqwest::header::RANGE, format!("bytes={}-{}", self.last_seek, self.last_seek + READ_SIZE as u64))
            .send()
            .map_err(|_| io::Error::from(io::ErrorKind::Other))?;
        let mut buf = vec![0u8; READ_SIZE as usize];
        r.read(&mut buf).unwrap();
        self.current = Some(Cursor::new(buf));
        Ok(self.last_seek)
    }
}

impl Clone for RemoteFile {
    fn clone(&self) -> Self {
        RemoteFile {
            url: self.url.clone(),
            last_seek: 0,
            current: None,
        }
    }
}

impl Reopen<RemoteFile> for RemoteFile {
    fn reopen(&self) -> io::Result<RemoteFile> {
        Ok(RemoteFile {
            url: self.url.clone(),
            last_seek: 0,
            current: None,
        })
    }
}