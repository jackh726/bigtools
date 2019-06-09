use std::fs::File;
use std::io::{self, Read, Write};
use std::sync::{Arc};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use parking_lot::{Condvar, Mutex};

// This struct provides a way to buffer data in a file, instead of memory.
#[derive(Debug)]
pub struct TempFileWrite {
    closed: bool,
}

#[derive(Debug)]
pub struct TempFileWriteWriter {
    statemutex: Arc<(Mutex<TempFileWrite>, Condvar)>,
    file_bytes: Arc<AtomicUsize>,
    iswaiting: Arc<AtomicBool>,
    boundarylocked: Arc<AtomicBool>,
    writefile: File,
}

pub struct TempFileWriteReader {
    statemutex: Arc<(Mutex<TempFileWrite>, Condvar)>,
    file_bytes: Arc<AtomicUsize>,
    iswaiting: Arc<AtomicBool>,
    boundarylocked: Arc<AtomicBool>,
    readfile: File,
}

impl TempFileWrite {
    pub fn new() -> io::Result<(TempFileWriteWriter, TempFileWriteReader)> {
        let temp = tempfile::NamedTempFile::new()?;

        let write = TempFileWrite {
            closed: false,
        };
        let pair = Arc::new((Mutex::new(write), Condvar::new()));
        let pair2 = pair.clone();
        let file_bytes = Arc::new(AtomicUsize::new(0));
        let iswaiting = Arc::new(AtomicBool::new(false));
        let boundarylocked = Arc::new(AtomicBool::new(false));
        Ok((
            TempFileWriteWriter {
                statemutex: pair,
                file_bytes: file_bytes.clone(),
                iswaiting: iswaiting.clone(),
                boundarylocked: boundarylocked.clone(),
                writefile: temp.reopen()?,
            },
            TempFileWriteReader {
                statemutex: pair2,
                file_bytes: file_bytes,
                iswaiting: iswaiting,
                boundarylocked: boundarylocked,
                readfile: temp.into_file(),
            }
        ))
    }
}

impl Drop for TempFileWriteWriter {
    fn drop(&mut self) {
        let &(ref lock, ref cvar) = &*self.statemutex;
        let mut state = lock.lock();
        state.closed = true;
        cvar.notify_one();
    }
}

impl Write for TempFileWriteWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        // We use this atomicbool as a sort of mutex. The idea is that in most cases, we aren't reading at the write/read boundary, so we don't need to "lock".
        while self.boundarylocked.compare_and_swap(false, true, Ordering::SeqCst) {
            // This means that the reader wants to read up to the write/read boundary
            // We can't write currently, otherwise we might get a race condition
            let &(ref lock, ref cvar) = &*self.statemutex;
            let mut state = lock.lock();
            self.iswaiting.store(true, Ordering::SeqCst);
            cvar.wait(&mut state);
            self.iswaiting.store(false, Ordering::SeqCst);
        }

        let writtenbytes = self.writefile.write(buf)?;
        self.file_bytes.fetch_add(writtenbytes, Ordering::SeqCst);

        if self.iswaiting.load(Ordering::SeqCst) {
            let &(_, ref cvar) = &*self.statemutex;
            self.boundarylocked.store(false, Ordering::SeqCst);
            cvar.notify_one();
        } else {
            self.boundarylocked.store(false, Ordering::SeqCst);
        }
        Ok(writtenbytes)
    }
    fn flush(&mut self) -> io::Result<()> { self.writefile.flush() }
}

impl Read for TempFileWriteReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let bytes_to_read = buf.len();
        if bytes_to_read > self.file_bytes.load(Ordering::SeqCst) {
            while self.boundarylocked.compare_and_swap(false, true, Ordering::SeqCst) {
                // We could go over, and we are currently writing
                let &(ref lock, ref cvar) = &*self.statemutex;
                let mut state = lock.lock();
                self.iswaiting.store(true, Ordering::SeqCst);
                println!("Waiting for boundary, in order to read.");
                cvar.wait(&mut state);
                self.iswaiting.store(false, Ordering::SeqCst);
            }
            // At this point, we know we have the boundarylock
            let bytes_read = self.readfile.read(buf)?;
            debug_assert!(self.file_bytes.load(Ordering::SeqCst) >= bytes_read, "Read more bytes than expected. Expected: {:?}. Actual: {:?}", self.file_bytes.load(Ordering::SeqCst), bytes_read);
            self.file_bytes.fetch_sub(bytes_read, Ordering::SeqCst);

            if self.iswaiting.load(Ordering::SeqCst) {
                println!("Iswaiting, unlock.");
                let &(_, ref cvar) = &*self.statemutex;
                self.boundarylocked.store(false, Ordering::SeqCst);
                cvar.notify_one();
            } else {
                self.boundarylocked.store(false, Ordering::SeqCst);
            }
            return Ok(bytes_read);
        } else {
            let bytes_read = self.readfile.read(buf)?;
            self.file_bytes.fetch_sub(bytes_read, Ordering::SeqCst);
            return Ok(bytes_read);
        }
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        let bytes_to_read = buf.len();
        if bytes_to_read > self.file_bytes.load(Ordering::SeqCst) {
            println!("read_exact lock");
            let &(ref lock, ref cvar) = &*self.statemutex;
            let mut state = lock.lock();

            self.iswaiting.store(true, Ordering::SeqCst);
            while !state.closed && self.file_bytes.load(Ordering::SeqCst) < bytes_to_read {
                cvar.wait(&mut state);
            }
            self.iswaiting.store(false, Ordering::SeqCst);

            self.readfile.read_exact(buf)?;

            let bytes_read = bytes_to_read;
            self.file_bytes.fetch_sub(bytes_read, Ordering::SeqCst);

            return Ok(());
        } else {
            self.readfile.read_exact(buf)?;
            self.file_bytes.fetch_sub(bytes_to_read, Ordering::SeqCst);

            return Ok(());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate test;

    #[test]
    fn test_works() -> io::Result<()> {
        let (mut writer, mut reader) = TempFileWrite::new()?;
        let test = &[1u8; 4];
        writer.write(test)?;
        let buf = &mut [0u8; 4];
        reader.read_exact(buf)?;
        assert!(buf == &[1u8; 4]);
        Ok(())
    }

    #[test]
    fn test_no_overread() -> io::Result<()> {
        let (mut writer, mut reader) = TempFileWrite::new()?;
        writer.write(&[1u8; 4])?;
        std::thread::spawn(move || {
            std::thread::sleep(std::time::Duration::from_millis(500));
            writer.write(&[1u8; 1]).unwrap();
        });
        let buf = &mut [0u8; 5];
        reader.read_exact(buf)?;
        assert!(buf.len() == 5);
        Ok(())
    }

    #[test]
    fn test_fails_on_writer_close() -> io::Result<()>  {
        let (mut writer, mut reader) = TempFileWrite::new()?;
        writer.write(&[1u8; 2])?;
        std::thread::spawn(move || {
            writer.write(&[1u8; 1]).unwrap();
            std::thread::sleep(std::time::Duration::from_millis(500));
            writer.write(&[1u8; 1]).unwrap();
        });
        let buf = &mut [0u8; 5];
        assert!(reader.read_exact(buf).is_err());
        Ok(())
    }

}
