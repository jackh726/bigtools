use std::fs::File;
use std::io::{self, Read, Write};
use std::sync::{Arc, Mutex, Condvar};


struct TempFileWrite {
    closed: bool,
    total_bytes: usize,
    file_bytes: usize,
}

struct TempFileWriteWriter {
    statemutex: Arc<(Mutex<TempFileWrite>, Condvar)>,
    writefile: File,
}

struct TempFileWriteReader {
    statemutex: Arc<(Mutex<TempFileWrite>, Condvar)>,
    readfile: File,
}

impl TempFileWrite {
    fn new() -> io::Result<(TempFileWriteWriter, TempFileWriteReader)> {
        let temp = tempfile::NamedTempFile::new()?;

        let write = TempFileWrite {
            closed: false,
            total_bytes: 0,
            file_bytes: 0,
        };
        let pair = Arc::new((Mutex::new(write), Condvar::new()));
        let pair2 = pair.clone();
        Ok((
            TempFileWriteWriter {
                statemutex: pair,
                writefile: temp.reopen()?,
            },
            TempFileWriteReader {
                statemutex: pair2,
                readfile: temp.into_file(),
            }
        ))
    }
}

impl Drop for TempFileWriteWriter {
    fn drop(&mut self) {
        let &(ref lock, ref cvar) = &*self.statemutex;
        let mut state = lock.lock().unwrap();
        state.closed = true;
        cvar.notify_one();
    }
}

impl Write for TempFileWriteWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let writtenbytes = self.writefile.write(buf)?;

        let &(ref lock, ref cvar) = &*self.statemutex;
        let mut state = lock.lock().unwrap();
        state.total_bytes += writtenbytes;
        state.file_bytes += writtenbytes;
        cvar.notify_one();
        
        Ok(writtenbytes)
    }
    fn flush(&mut self) -> io::Result<()> { self.writefile.flush() }
}

impl Read for TempFileWriteReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let bytes_read = self.readfile.read(buf)?;

        let &(ref lock, ref cvar) = &*self.statemutex;
        let mut state = lock.lock().unwrap();
        state.file_bytes -= bytes_read;
        cvar.notify_one();

        Ok(bytes_read)
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        let &(ref lock, ref cvar) =  &*self.statemutex;
        let mut state = lock.lock().unwrap();
        let bytes_to_read = buf.len();

        while !state.closed && state.file_bytes < bytes_to_read {
            state = cvar.wait(state).unwrap();
        }
        self.readfile.read_exact(buf)
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