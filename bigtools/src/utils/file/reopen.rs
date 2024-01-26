use std::fs::File;
use std::io::{self, Read, Seek};

/// A helper trait that for things that implement `Read`, `Seek`, and `Send`
pub trait SeekableRead: Seek + Read {}
impl<T> SeekableRead for T where T: Seek + Read {}

/// Indicates something that can be *reopened*. Importantly, reopening should be independent
/// with respect to seeks and reads from the original object.
pub trait Reopen: Sized {
    fn reopen(&self) -> io::Result<Self>;
}

pub struct ReopenableFile {
    pub path: String,
    pub file: File,
}

impl Reopen for ReopenableFile {
    fn reopen(&self) -> io::Result<Self> {
        Ok(ReopenableFile {
            path: self.path.clone(),
            file: File::open(&self.path)?,
        })
    }
}

impl Seek for ReopenableFile {
    fn seek(&mut self, pos: io::SeekFrom) -> io::Result<u64> {
        self.file.seek(pos)
    }
}

impl Read for ReopenableFile {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.file.read(buf)
    }

    fn read_vectored(&mut self, bufs: &mut [io::IoSliceMut<'_>]) -> io::Result<usize> {
        self.file.read_vectored(bufs)
    }

    fn read_to_end(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        self.file.read_to_end(buf)
    }

    fn read_to_string(&mut self, buf: &mut String) -> io::Result<usize> {
        self.file.read_to_string(buf)
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        self.file.read_exact(buf)
    }
}
