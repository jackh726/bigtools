use std::fs::File;
use std::io::{self, Read, Write};
use parking_lot::{Condvar, Mutex};
use std::sync::{Arc};

#[derive(Debug)]
enum BufferState {
    Temp(Option<File>),
    Real(Option<File>),
}

pub struct TempFileBuffer {
    state: Arc<Mutex<BufferState>>,
    closed: Arc<(Mutex<bool>, Condvar)>,
}

pub struct TempFileBufferWriter {
    state: Arc<Mutex<BufferState>>,
    closed: Arc<(Mutex<bool>, Condvar)>,
}

impl TempFileBuffer {
    pub fn new() -> io::Result<(TempFileBuffer, TempFileBufferWriter)> {
        let file = tempfile::tempfile()?;
        let state = Arc::new(Mutex::new(BufferState::Temp(Some(file)))); 
        let closed = Arc::new((Mutex::new(false), Condvar::new()));
        Ok((
            TempFileBuffer { state: state.clone(), closed: closed.clone() },
            TempFileBufferWriter { state, closed },
        ))
    }

    pub fn new_from_real(file: File) -> io::Result<(TempFileBuffer, TempFileBufferWriter)> {
        let state = Arc::new(Mutex::new(BufferState::Real(Some(file)))); 
        let closed = Arc::new((Mutex::new(false), Condvar::new()));
        Ok((
            TempFileBuffer { state: state.clone(), closed: closed.clone() },
            TempFileBufferWriter { state, closed },
        ))
    }

    pub fn switch(&mut self, mut new_file: File) -> io::Result<()> {
        let mut guard = self.state.lock();
        use std::ops::DerefMut;
        let mut state: &mut BufferState = guard.deref_mut();

        match state {
            BufferState::Temp(ref mut file) => {
                match file {
                    None => panic!(),
                    Some(file) => {
                        use std::io::Seek;
                        file.seek(io::SeekFrom::Start(0))?;

                        std::io::copy(file, &mut new_file)?;
                        *guard = BufferState::Real(Some(new_file));
                    }
                }
                Ok(())
            },
            BufferState::Real(ref mut file) => {
                panic!("Already switched!");
            }
        }
    }

    pub fn await(self) -> File {
        let &(ref lock, ref cvar) =  &*self.closed;
        let mut closed = lock.lock();

        while !*closed {
            cvar.wait(&mut closed);
        }

        let mut guard = self.state.lock();
        use std::ops::DerefMut;
        let state: &mut BufferState = guard.deref_mut();

        match state {
            BufferState::Temp(ref mut file) => {
                panic!("Must be called after switch.");
            },
            BufferState::Real(ref mut file) => {
                let file = file.take();
                match file {
                    None => panic!("Already switched!"),
                    Some(file) => file,
                }
            }
        }
    }

    pub fn await_raw(self) -> File {
        let &(ref lock, ref cvar) =  &*self.closed;
        let mut closed = lock.lock();

        while !*closed {
            cvar.wait(&mut closed);
        }

        let mut guard = self.state.lock();
        use std::ops::DerefMut;
        let state: &mut BufferState = guard.deref_mut();

        match state {
            BufferState::Temp(ref mut file) => {
                file.take().unwrap()
            },
            BufferState::Real(ref mut file) => {
                let file = file.take();
                match file {
                    None => panic!("Already switched!"),
                    Some(file) => file,
                }
            }
        }
    }
}

impl Drop for TempFileBufferWriter {
    fn drop(&mut self) {
        let &(ref lock, ref cvar) = &*self.closed;
        let mut closed = lock.lock();
        *closed = true;
        cvar.notify_one();
    }
}

impl Write for TempFileBufferWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut guard = self.state.lock();
        use std::ops::DerefMut;
        let state: &mut BufferState = guard.deref_mut();
        //println!("Write {:?}", state);

        match state {
            BufferState::Temp(ref mut file) => {
                match file {
                    None => panic!(),
                    Some(file) => file.write(buf),
                }
            },
            BufferState::Real(ref mut file) => {
                match file {
                    Some(file) => file.write(buf),
                    None => panic!("Should have been dropped by now."),
                }
                
            }
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        let mut guard = self.state.lock();
        use std::ops::DerefMut;
        let state: &mut BufferState = guard.deref_mut();

        match state {
            BufferState::Temp(ref mut file) => {
                match file {
                    None => panic!(),
                    Some(file) => file.flush(),
                }
            },
            BufferState::Real(ref mut file) => {
                match file {
                    Some(file) => file.flush(),
                    None => panic!("Should have been dropped by now."),
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate test;

    #[test]
    fn test_works() -> io::Result<()> {
        let (mut buf, mut writer) = TempFileBuffer::new()?;

        const NUM_BYTES: usize = 50;
        let writethread = std::thread::spawn(move || {
            for i in 0..NUM_BYTES {
                std::thread::sleep(std::time::Duration::from_millis(50));
                let writebuf = &mut [(i % 8) as u8; 1];
                writer.write(writebuf).unwrap();
            }
        });

        std::thread::sleep(std::time::Duration::from_millis(250));

        let outfile = tempfile::tempfile()?;
        buf.switch(outfile)?;

        let mut file = buf.await();

        use std::io::Seek;
        file.seek(io::SeekFrom::Start(0))?;

        let mut out_bytes = vec![];
        file.read_to_end(&mut out_bytes)?;

        assert_eq!(out_bytes.len(), NUM_BYTES, "All bytes not accounted for.");
        Ok(())
    }
}