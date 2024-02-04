use std::fs::File;
use std::io::{self, Seek, Write};
use std::marker::Send;
use std::sync::{Arc, Condvar, Mutex};

use crossbeam_utils::atomic::AtomicCell;

// In theory, we could have a `Memory` variant that temporarily buffers in memory until a given size.
enum BufferState<R> {
    NotStarted,
    InMemory(Vec<u8>),
    Temp(File),
    Real(R),
}

/// This struct provides a way to "buffer" data in a temporary file
/// until the "real" file is ready to be written.
/// This is useful if you have parallel generation of data that needs
/// to be written to a file, but you don't know the size of each part.
pub struct TempFileBuffer<R> {
    closed: Arc<(Mutex<Option<BufferState<R>>>, Condvar)>,
    real_file: Arc<AtomicCell<Option<R>>>,
}

pub struct TempFileBufferWriter<R> {
    closed: Arc<(Mutex<Option<BufferState<R>>>, Condvar)>,
    buffer_state: BufferState<R>,
    real_file: Arc<AtomicCell<Option<R>>>,
    inmemory: bool,
}

impl<R> TempFileBuffer<R> {
    /// Creates a new `TempFileBuffer`/`TempFileBufferWriter` pair where the writer is writing to a temporary file or in-memory
    pub fn new(inmemory: bool) -> (TempFileBuffer<R>, TempFileBufferWriter<R>) {
        let closed = Arc::new((Mutex::new(None), Condvar::new()));
        let buffer_state = BufferState::NotStarted;
        let real_file = Arc::new(AtomicCell::new(None));
        (
            TempFileBuffer {
                closed: closed.clone(),
                real_file: real_file.clone(),
            },
            TempFileBufferWriter {
                closed,
                buffer_state,
                real_file,
                inmemory,
            },
        )
    }
}

impl<R: Write + Send + 'static> TempFileBuffer<R> {
    /// Switches to the "real" file. This doesn't actually do any processing, only sets up everything to for next write or close.
    pub fn switch(&mut self, new_file: R) {
        if self.real_file.swap(Some(new_file)).is_some() {
            panic!("Can only switch once.");
        }
    }

    pub fn is_real_file_ready(&self) -> bool {
        let &(ref lock, _) = &*self.closed;
        let closed = lock.lock().unwrap();

        closed.is_some()
    }

    pub fn len(&self) -> io::Result<u64> {
        let &(ref lock, ref cvar) = &*self.closed;
        let mut closed = lock.lock().unwrap();

        while closed.is_none() {
            closed = cvar.wait(closed).unwrap();
        }
        let closed = closed.as_mut();

        match closed.unwrap() {
            BufferState::Real(_) => panic!("Should not have switched already."),
            BufferState::InMemory(data) => Ok(data.len() as u64),
            BufferState::Temp(ref mut t) => t.seek(io::SeekFrom::Current(0)),
            BufferState::NotStarted => Ok(0),
        }
    }

    pub fn await_real_file(self) -> R {
        let &(ref lock, ref cvar) = &*self.closed;
        let mut closed = lock.lock().unwrap();

        while closed.is_none() {
            closed = cvar.wait(closed).unwrap();
        }
        let closed = closed.take().unwrap();

        let real_file = self.real_file.swap(None);

        match (real_file, closed) {
            (Some(mut real_file), BufferState::InMemory(data)) => {
                // Switch was called but no writes have happened
                // Writer was dropped with data having been written
                real_file.write_all(&data).unwrap();
                real_file
            }
            (Some(mut real_file), BufferState::Temp(mut closed_file)) => {
                // Switch was called but no writes have happened
                // Writer was dropped with temp file having been written
                closed_file.seek(io::SeekFrom::Start(0)).unwrap();
                io::copy(&mut closed_file, &mut real_file).unwrap();
                real_file
            }
            (Some(_), BufferState::Real(_)) => unreachable!(),
            (Some(real_file), BufferState::NotStarted) => {
                // Switch was called but no writes have happened
                // Writer was dropped with no tempfile being created (or written to)
                real_file
            }
            (None, BufferState::Real(real_file)) => real_file,
            (None, BufferState::InMemory(_) | BufferState::Temp(_) | BufferState::NotStarted) => {
                panic!("Should have switched already.")
            }
        }
    }

    pub fn expect_closed_write<O>(self, mut real: &mut O) -> io::Result<()>
    where
        O: Write,
    {
        let &(ref lock, ref cvar) = &*self.closed;
        let mut closed = lock.lock().unwrap();

        while closed.is_none() {
            closed = cvar.wait(closed).unwrap();
        }
        let closed = closed.take().unwrap();

        let real_file = self.real_file.swap(None);
        assert!(real_file.is_none(), "Should only be writing to real file.");

        match closed {
            BufferState::Temp(mut closed_file) => {
                closed_file.seek(io::SeekFrom::Start(0))?;
                io::copy(&mut closed_file, &mut real)?;
            }
            BufferState::InMemory(data) => {
                real.write_all(&data)?;
            }
            BufferState::NotStarted => {}
            BufferState::Real(_) => panic!("Should only be writing to real file."),
        }
        Ok(())
    }
}

impl<R: Write + Send + 'static> TempFileBufferWriter<R> {
    fn update(&mut self) -> io::Result<()> {
        match &mut self.buffer_state {
            BufferState::NotStarted => {
                let real_file = self.real_file.swap(None).take();
                match real_file {
                    Some(new_file) => {
                        self.buffer_state = BufferState::Real(new_file);
                    }
                    None => {
                        if self.inmemory {
                            self.buffer_state =
                                BufferState::InMemory(Vec::with_capacity(10 * 1_000));
                        } else {
                            self.buffer_state = BufferState::Temp(tempfile::tempfile()?);
                        }
                    }
                }
            }
            BufferState::InMemory(data) => {
                let real_file = self.real_file.swap(None).take();
                if let Some(mut new_file) = real_file {
                    new_file.write_all(&data)?;
                    self.buffer_state = BufferState::Real(new_file);
                }
            }
            BufferState::Temp(ref mut file) => {
                let real_file = self.real_file.swap(None).take();
                if let Some(mut new_file) = real_file {
                    file.seek(io::SeekFrom::Start(0))?;

                    io::copy(file, &mut new_file)?;
                    self.buffer_state = BufferState::Real(new_file);
                }
            }
            BufferState::Real(_) => {}
        }
        Ok(())
    }
}

impl<R> Drop for TempFileBufferWriter<R> {
    fn drop(&mut self) {
        let &(ref lock, ref cvar) = &*self.closed;
        let mut closed = lock.lock().unwrap();
        let buffer_state = std::mem::replace(&mut self.buffer_state, BufferState::NotStarted);
        *closed = Some(buffer_state);
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
                BufferState::InMemory(ref mut data) => return data.write(buf),
                BufferState::Temp(ref mut file) => return file.write(buf),
                BufferState::Real(ref mut file) => return file.write(buf),
            }
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self.buffer_state {
            BufferState::NotStarted => Ok(()), // No data has been written, nothing to flush
            BufferState::InMemory(_) => Ok(()), // All data is written immediately to vec
            BufferState::Temp(ref mut file) => file.flush(),
            BufferState::Real(ref mut file) => file.flush(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn test_works() -> io::Result<()> {
        let (mut buf, mut writer) = TempFileBuffer::new(true);

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
