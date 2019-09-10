use std::fs::File;
use std::io::{self, Seek, Write};
use std::marker::Send;
use std::sync::{Arc};

use parking_lot::{Condvar, Mutex};

use crossbeam::atomic::AtomicCell;

// TODO: Change TempFileBuffer to memory buffer if under a certain threshold

enum BufferState<R: Write + Send + 'static> {
    NotStarted,
    Temp(Option<File>),
    Real(Option<R>),
}

pub enum ClosedFile<R: Write + Send + 'static> {
    Temp(File),
    Real(R)
}

/// This struct provides a way to "buffer" data in a temporary file
/// until the "real" file is ready to be written.
/// This is useful if you have parallel generation of data that needs
/// to be written to a file, but you don't know the size of each part.
pub struct TempFileBuffer<R: Write + Send + 'static> {
    closed: Arc<(Mutex<bool>, Condvar)>,
    closed_file: Arc<AtomicCell<Option<ClosedFile<R>>>>,
    real_file: Arc<AtomicCell<Option<R>>>,
    has_switched: bool,
}

pub struct TempFileBufferWriter<R: Write + Send + 'static> {
    closed: Arc<(Mutex<bool>, Condvar)>,
    buffer_state: BufferState<R>,
    closed_file: Arc<AtomicCell<Option<ClosedFile<R>>>>,
    real_file: Arc<AtomicCell<Option<R>>>,
}

impl<R: Write + Send + 'static> TempFileBuffer<R> {
    /// Creates a new `TempFileBuffer`/`TempFileBufferWriter` pair where the writer is writing to a temporary file
    pub fn new() -> io::Result<(TempFileBuffer<R>, TempFileBufferWriter<R>)> {
        let closed = Arc::new((Mutex::new(false), Condvar::new()));
        let buffer_state = BufferState::NotStarted;
        let closed_file = Arc::new(AtomicCell::new(None));
        let real_file = Arc::new(AtomicCell::new(None));
        Ok((
            TempFileBuffer { closed: closed.clone(), closed_file: closed_file.clone(), real_file: real_file.clone(), has_switched: false },
            TempFileBufferWriter { closed, buffer_state, closed_file, real_file },
        ))
    }

    /// Creates a new `TempFileBuffer`/`TempFileBufferWriter` pair where the writer is writing to the "real" file.
    /// This is more or less equivalent to
    /// ```no_run
    /// # use std::fs::File;
    /// # use std::io;
    /// # use bigwig2::tempfilebuffer::TempFileBuffer;
    /// let mut file = File::open("foo")?;
    /// let (mut buf, writer) = TempFileBuffer::new()?;
    /// buf.switch(file);
    /// # Ok::<(), io::Error>(())
    /// ```
    #[allow(dead_code)]
    pub fn new_from_real(file: R) -> io::Result<(TempFileBuffer<R>, TempFileBufferWriter<R>)> {
        let closed = Arc::new((Mutex::new(false), Condvar::new()));
        let buffer_state = BufferState::Real(Some(file));
        let closed_file = Arc::new(AtomicCell::new(None));
        let real_file = Arc::new(AtomicCell::new(None));
        Ok((
            TempFileBuffer { closed: closed.clone(), closed_file: closed_file.clone(), real_file: real_file.clone(), has_switched: true },
            TempFileBufferWriter { closed, buffer_state, closed_file, real_file },
        ))
    }

    /// Switches to the "real" file. This doesn't actually do any processing, only sets up everything to for next write or close.
    pub fn switch(&mut self, new_file: R) {
        if self.has_switched {
            panic!("Can only switch once.");
        }
        self.has_switched = true;
        if self.real_file.swap(Some(new_file)).is_some() {
            panic!("(Invalid state) Can only switch once.");
        }
    }

    pub fn await_file(self) -> ClosedFile<R> {
        let &(ref lock, ref cvar) =  &*self.closed;
        let mut closed = lock.lock();

        while !*closed {
            cvar.wait(&mut closed);
        }

        let real_file = self.real_file.swap(None);
        let closed_file = self.closed_file.swap(None);

        match (real_file, closed_file) {
            (Some(mut real_file), Some(ClosedFile::Temp(mut closed_file))) => {
                // Switch was called but no writes have happened
                // Writer was dropped with temp file having been written
                closed_file.seek(io::SeekFrom::Start(0)).unwrap();
                io::copy(&mut closed_file, &mut real_file).unwrap();
                ClosedFile::Real(real_file)
            },
            (Some(_), Some(ClosedFile::Real(_))) => unreachable!(),
            (Some(real_file), None) => {
                // Switch was called but no writes have happened
                // Writer was dropped with no tempfile being created (or written to)
                ClosedFile::Real(real_file)
            },
            (None, Some(closed)) => {
                // Switch was not called or called with subsequent writes
                closed
            },
            (None, None) => panic!("No data was written."),
        }
    }

    pub fn await_real_file(self) -> R {
        let &(ref lock, ref cvar) =  &*self.closed;
        let mut closed = lock.lock();

        while !*closed {
            cvar.wait(&mut closed);
        }

        if !self.has_switched {
            panic!("Should have switched already.");
        }

        let real_file = self.real_file.swap(None);
        let closed_file = self.closed_file.swap(None);

        match (real_file, closed_file) {
            (Some(mut real_file), Some(ClosedFile::Temp(mut closed_file))) => {
                // Switch was called but no writes have happened
                // Writer was dropped with temp file having been written
                closed_file.seek(io::SeekFrom::Start(0)).unwrap();
                io::copy(&mut closed_file, &mut real_file).unwrap();
                real_file
            },
            (Some(_), Some(ClosedFile::Real(_))) => unreachable!(),
            (Some(real_file), None) => {
                // Switch was called but no writes have happened
                // Writer was dropped with no tempfile being created (or written to)
                real_file
            },
            (None, Some(ClosedFile::Temp(_))) => unreachable!(), // Checked self.has_switched
            (None, Some(ClosedFile::Real(real_file))) => {
                real_file
            },
            (None, None) => panic!("No data was written."),
        }
    }

    pub fn await_temp_file(self) -> File {
        let &(ref lock, ref cvar) =  &*self.closed;
        let mut closed = lock.lock();

        while !*closed {
            cvar.wait(&mut closed);
        }

        if self.has_switched {
            panic!("Should not have switched already.");
        }

        let real_file = self.real_file.swap(None);
        assert!(real_file.is_none(), "Should not have switched already.");
        let closed_file = self.closed_file.swap(None);

        match closed_file {
            Some(ClosedFile::Real(_)) => unreachable!(),
            Some(ClosedFile::Temp(t)) => {
                t
            },
            None => panic!("No data was written."),
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
            Some(ClosedFile::Temp(mut closed_file)) => {
                closed_file.seek(io::SeekFrom::Start(0))?;
                io::copy(&mut closed_file, &mut real)?;
            },
            Some(ClosedFile::Real(_)) => unreachable!(),
            None => panic!("No data was written."),
            
        }
        Ok(())
    }
}

impl<R: Write + Send + 'static> TempFileBufferWriter<R> {
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
        Ok(())
    }
}

impl<R: Write + Send + 'static> Drop for TempFileBufferWriter<R> {
    fn drop(&mut self) {
        let &(ref lock, ref cvar) = &*self.closed;
        let mut closed = lock.lock();
        match &mut self.buffer_state {
            BufferState::NotStarted => {},
            BufferState::Temp(f) => {
                let temp = f.take();
                self.closed_file.swap(Some(ClosedFile::Temp(temp.unwrap())));
            },
            BufferState::Real(f) => {
                let real = f.take();
                self.closed_file.swap(Some(ClosedFile::Real(real.unwrap())));
            },
        }
        *closed = true;
        cvar.notify_one();
        drop(closed);
    }
}

impl<R: Write + Send + 'static> Write for TempFileBufferWriter<R> {
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
        buf.switch(outfile);

        let mut file = buf.await_real_file();

        use std::io::Seek;
        file.seek(io::SeekFrom::Start(0))?;

        let mut out_bytes = vec![];
        file.read_to_end(&mut out_bytes)?;

        assert_eq!(out_bytes.len(), NUM_BYTES, "All bytes not accounted for.");
        Ok(())
    }
}