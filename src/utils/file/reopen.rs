use std::fs::File;
use std::io::{self, Read, Seek};
use std::marker::Send;

/// A helper trait that for things that implement `Read`, `Seek`, and `Send`
pub trait SeekableRead: Seek + Read + Send {}
impl<T> SeekableRead for T where T: Seek + Read + Send {}

/// Indicates something that can be *reopened*. Importantly, reopening should be independent
/// with respect to seeks and reads from the original object.
pub trait Reopen<S>: Clone + Send
where
    S: SeekableRead,
{
    fn reopen(&self) -> io::Result<S>;
}

#[derive(Clone)]
pub struct ReopenableFile {
    pub path: String,
}

impl Reopen<File> for ReopenableFile {
    fn reopen(&self) -> io::Result<File> {
        File::open(self.path.clone())
    }
}
