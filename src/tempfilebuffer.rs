use std::fs::File;
use std::io::{self, Seek, Write};
use parking_lot::{Condvar, Mutex};
use std::sync::{Arc};

use crossbeam::atomic::AtomicCell;

#[derive(Debug)]
enum BufferState {
    NotStarted,
    Temp(Option<File>),
    Real(Option<File>),
}

/// This struct provides a way to "buffer" data in a temporary file
/// until the "real" file is ready to be written.
/// This is useful if you have parallel generation of data that needs
/// to be written to a file, but you don't know the size of each part.
pub struct TempFileBuffer {
    closed: Arc<(Mutex<bool>, Condvar)>,
    closed_file: Arc<AtomicCell<Option<File>>>,
    real_file: Arc<AtomicCell<Option<File>>>,
    has_switched: bool,
}

pub struct TempFileBufferWriter {
    closed: Arc<(Mutex<bool>, Condvar)>,
    buffer_state: BufferState,
    closed_file: Arc<AtomicCell<Option<File>>>,
    real_file: Arc<AtomicCell<Option<File>>>,
}

impl TempFileBuffer {
    pub fn new() -> io::Result<(TempFileBuffer, TempFileBufferWriter)> {
        let closed = Arc::new((Mutex::new(false), Condvar::new()));
        let buffer_state = BufferState::NotStarted;
        let closed_file = Arc::new(AtomicCell::new(None));
        let real_file = Arc::new(AtomicCell::new(None));
        Ok((
            TempFileBuffer { closed: closed.clone(), closed_file: closed_file.clone(), real_file: real_file.clone(), has_switched: false },
            TempFileBufferWriter { closed, buffer_state, closed_file, real_file },
        ))
    }

    #[allow(dead_code)]
    pub fn new_from_real(file: File) -> io::Result<(TempFileBuffer, TempFileBufferWriter)> {
        let closed = Arc::new((Mutex::new(false), Condvar::new()));
        let buffer_state = BufferState::Real(Some(file));
        let closed_file = Arc::new(AtomicCell::new(None));
        let real_file = Arc::new(AtomicCell::new(None));
        Ok((
            TempFileBuffer { closed: closed.clone(), closed_file: closed_file.clone(), real_file: real_file.clone(), has_switched: false },
            TempFileBufferWriter { closed, buffer_state, closed_file, real_file },
        ))
    }

    pub fn switch(&mut self, new_file: File) -> io::Result<()> {
        self.has_switched = true;
        self.real_file.swap(Some(new_file));
        Ok(())
    }

    pub fn await_file(self) -> File {
        let &(ref lock, ref cvar) =  &*self.closed;
        let mut closed = lock.lock();

        while !*closed {
            cvar.wait(&mut closed);
        }

        let real_file = self.real_file.swap(None);
        let closed_file = self.closed_file.swap(None);

        match real_file {
            Some(mut real_file) => {
                // Switch was previously called, but no writes have happened
                match closed_file {
                    Some(mut closed_file) => {
                        // The writer was dropped before the temp file was copied/dropped
                        closed_file.seek(io::SeekFrom::Start(0)).unwrap();
                        io::copy(&mut closed_file, &mut real_file).unwrap();
                        return real_file;
                    },
                    None => {
                        // The writer was dropped before any data was written and a temp file was created
                        return tempfile::tempfile().expect("Couldn't create tempfile");
                    }
                }
            },
            None => {
                // Either switch has not been called, the temp file has been copied/dropped, or no data has been written
                match closed_file {
                    Some(closed_file) => {
                        // Either switch has not been called (and this is the temp file), or the temp file has been copied/dropped (and this is the real file)
                        return closed_file;
                    },
                    None => {
                        // No data has been written yet
                        panic!("No data was written.");
                    }
                }
            }
        }
    }

    pub fn expect_closed_write<O>(self, mut real: &mut O) -> io::Result<()> where O: Write {
        let &(ref lock, ref cvar) =  &*self.closed;
        let mut closed = lock.lock();

        while !*closed {
            cvar.wait(&mut closed);
        }

        if self.has_switched {
            panic!("Should only be writing to real file.");
        }

        let real_file = self.real_file.swap(None);
        assert!(real_file.is_none(), "Should only be writing to real file.");

        let closed_file = self.closed_file.swap(None);
        match closed_file {
            Some(mut closed_file) => {
                // The temp file has been copied/dropped (and this is the real file)
                closed_file.seek(io::SeekFrom::Start(0))?;
                io::copy(&mut closed_file, &mut real)?;
            },
            None => {
                // No data has been written yet
                panic!("No data was written.");
            }
        }
        return Ok(());
    }
}

impl TempFileBufferWriter {
    fn update(&mut self) -> io::Result<()> {
        match &mut self.buffer_state {
            BufferState::NotStarted => {
                self.buffer_state = BufferState::Temp(Some(tempfile::tempfile()?));
            },
            BufferState::Temp(ref mut file) => {
                let real_file = self.real_file.swap(None).take();
                if let Some(mut new_file) = real_file {
                    let mut file = file.as_mut().unwrap();
                    file.seek(io::SeekFrom::Start(0))?;

                    io::copy(&mut file, &mut new_file)?;
                    self.buffer_state = BufferState::Real(Some(new_file));
                }
            },
            BufferState::Real(_) => {},
        }
        return Ok(())
    }
}

impl Drop for TempFileBufferWriter {
    fn drop(&mut self) {
        let &(ref lock, ref cvar) = &*self.closed;
        let mut closed = lock.lock();
        match &mut self.buffer_state {
            BufferState::NotStarted => {},
            BufferState::Temp(None) => panic!(),
            BufferState::Temp(Some(ref mut file)) => {
                self.closed_file.swap(Some(std::mem::replace(file, tempfile::tempfile().unwrap())));
            },
            BufferState::Real(None) => panic!(),
            BufferState::Real(Some(ref mut file)) => {
                self.closed_file.swap(Some(std::mem::replace(file, tempfile::tempfile().unwrap())));
            },
        }
        *closed = true;
        cvar.notify_one();
        drop(closed);
    }
}

impl Write for TempFileBufferWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.update()?;
        loop {
            match self.buffer_state {
                BufferState::NotStarted => unreachable!(),
                BufferState::Temp(ref mut file) => return file.as_mut().unwrap().write(buf),
                BufferState::Real(ref mut file) => return file.as_mut().unwrap().write(buf),
            }
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self.buffer_state {
            BufferState::NotStarted => Ok(()), // No data has been written, nothing to flush
            BufferState::Temp(ref mut file) => file.as_mut().unwrap().flush(),
            BufferState::Real(ref mut file) => file.as_mut().unwrap().flush(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    extern crate test;

    #[test]
    fn test_works() -> io::Result<()> {
        let (mut buf, mut writer) = TempFileBuffer::new()?;

        const NUM_BYTES: usize = 50;
        let _writethread = std::thread::spawn(move || {
            for i in 0..NUM_BYTES {
                std::thread::sleep(std::time::Duration::from_millis(50));
                let writebuf = &mut [(i % 8) as u8; 1];
                writer.write(writebuf).unwrap();
            }
        });

        std::thread::sleep(std::time::Duration::from_millis(250));

        let outfile = tempfile::tempfile()?;
        buf.switch(outfile)?;

        let mut file = buf.await_file();

        use std::io::Seek;
        file.seek(io::SeekFrom::Start(0))?;

        let mut out_bytes = vec![];
        file.read_to_end(&mut out_bytes)?;

        assert_eq!(out_bytes.len(), NUM_BYTES, "All bytes not accounted for.");
        Ok(())
    }
}