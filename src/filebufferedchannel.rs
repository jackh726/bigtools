use std::fs::File;
use std::io::{self, Read, Write, Seek, SeekFrom};
use std::sync::Arc;

use crossbeam::channel::{bounded, Sender as ChannelSender, Receiver as ChannelReceiver};

use parking_lot::Mutex;

use bincode;

use serde::de::DeserializeOwned;
use serde::Serialize;

use crate::tell::Tell;

pub fn channel<T>(size: usize) -> (Sender<T>, Receiver<T>) where T: Serialize + DeserializeOwned {
    let state = Arc::new(Mutex::new(ChannelState::new(size)));
    let (channelsender, channelreceiver) = bounded(size);
    let sender = Sender { state: state.clone(), sender: channelsender };
    let receiver = Receiver { state: state, receiver: channelreceiver };
    (sender, receiver)
}

pub enum ChannelError {
    IOError(io::Error),
    Disconnected,
    Empty,
    InMemory,
}

impl From<io::Error> for ChannelError {
    fn from(error: io::Error) -> Self {
        ChannelError::IOError(error)
    }
}

enum ChannelStateStatus<T> {
    InMemory,
    OnDisk {
        serialize_size: Option<u64>,
        buffered: Option<T>,
        sender: ChannelSender<T>,
        buffer: ChannelReceiver<T>,
        file: File,
        readindex: u64,
        writeindex: u64,
    }
}

struct ChannelState<T> where T: Serialize + DeserializeOwned {
    maxsize: usize,
    status: ChannelStateStatus<T>,
}

type ChannelResult<R> = Result<R, ChannelError>;

impl<T> ChannelState<T> where T: Serialize + DeserializeOwned {
    fn new(maxsize: usize) -> ChannelState<T> {
        ChannelState {
            maxsize,
            status: ChannelStateStatus::InMemory,
        }
    }

    fn clearqueue(&mut self, sender: &mut ChannelSender<T>) -> io::Result<()> {
        match &mut self.status {
            ChannelStateStatus::InMemory => {
                let (channelsender, channelreceiver) = bounded(self.maxsize);
                let sender = std::mem::replace(sender, channelsender);
                self.status = ChannelStateStatus::OnDisk {
                    serialize_size: None,
                    buffered: None,
                    sender: sender,
                    buffer: channelreceiver,
                    file: tempfile::tempfile()?,
                    readindex: 0,
                    writeindex: 0,
                };
                return Ok(());
            },
            ChannelStateStatus::OnDisk { buffer, file, writeindex, serialize_size, .. } => {
                file.seek(SeekFrom::Start(*writeindex))?;
                for item in buffer.try_iter() {
                    if let None = serialize_size {
                        serialize_size.replace(bincode::serialized_size(&item).unwrap());
                    }
                    file.write(&bincode::serialize(&item).unwrap())?;
                }
                *writeindex = file.tell()?;
                debug_assert!(buffer.is_empty());
                return Ok(());
            }
        }
    }

    fn read(&mut self) -> ChannelResult<()> {
        match &mut self.status {
            ChannelStateStatus::InMemory => Err(ChannelError::InMemory),
            ChannelStateStatus::OnDisk { file, readindex, writeindex, buffer, sender, buffered, serialize_size } => {
                let bufferedelem = buffered.take();
                let mut sent = false;
                if let Some(elem) = bufferedelem {
                    if let Err(e) = sender.try_send(elem) {
                        use crossbeam::channel::TrySendError::*;
                        match e {
                            Disconnected(_) => return Err(ChannelError::Disconnected),
                            Full(t) => {
                                buffered.replace(t);
                                return Ok(())
                            },
                        }
                        
                    }
                    sent = true;
                }
                if readindex < writeindex {
                    file.seek(SeekFrom::Start(*readindex))?;
                }
                while readindex < writeindex {
                    let size = serialize_size.unwrap() as usize;
                    let buf = &mut vec![0u8; size];
                    file.read_exact(buf)?;
                    *readindex += size as u64;
                    let elem =  bincode::deserialize(buf).expect("Error while deserializing.");
                    if let Err(e) = sender.try_send(elem) {
                        use crossbeam::channel::TrySendError::*;
                        match e {
                            Disconnected(_) => return Err(ChannelError::Disconnected),
                            Full(t) => {
                                buffered.replace(t);
                                return Ok(())
                            },
                        }
                        
                    }
                    sent = true;
                }
                for elem in buffer.try_iter() {
                    if let Err(e) = sender.try_send(elem) {
                        use crossbeam::channel::TrySendError::*;
                        match e {
                            Disconnected(_) => return Err(ChannelError::Disconnected),
                            Full(t) => {
                                buffered.replace(t);
                                return Ok(())
                            },
                        }
                    }
                    sent = true;
                }
                if !sent {
                    return Err(ChannelError::Empty);
                } else {
                    return Ok(());
                }
            }
        }
    }

    fn read_wait(&mut self) -> ChannelResult<()> {
        match self.read() {
            Ok(_) => return Ok(()),
            Err(ChannelError::InMemory) => return Err(ChannelError::InMemory),
            Err(ChannelError::Disconnected) => return Err(ChannelError::Disconnected),
            Err(ChannelError::IOError(e)) => return Err(ChannelError::IOError(e)),
            Err(ChannelError::Empty) => {
                match &mut self.status {
                    ChannelStateStatus::InMemory => unreachable!(),
                    ChannelStateStatus::OnDisk { buffer, sender, buffered, .. } => {
                        debug_assert!(buffered.is_none());
                        use crossbeam::channel::RecvError;
                        match buffer.recv() {
                            Ok(elem) => {
                                if let Err(e) = sender.try_send(elem) {
                                    use crossbeam::channel::TrySendError::*;
                                    match e {
                                        Disconnected(_) => return Err(ChannelError::Disconnected),
                                        Full(t) => {
                                            buffered.replace(t);
                                            return Ok(())
                                        },
                                    }
                                }
                                return Ok(())
                            },
                            Err(RecvError) => return Err(ChannelError::Disconnected),
                        }
                    }
                }
            }
        }
    }
}

pub struct Sender<T> where T: Serialize + DeserializeOwned {
    state: Arc<Mutex<ChannelState<T>>>,
    sender: ChannelSender<T>,
}

pub struct Receiver<T> where T: Serialize + DeserializeOwned {
    state: Arc<Mutex<ChannelState<T>>>,
    receiver: ChannelReceiver<T>,
}

impl<T> Sender<T> where T: Serialize + DeserializeOwned {
    pub fn send(&mut self, t: T) -> io::Result<()> {
        if let Err(e) = self.sender.try_send(t) {
            use crossbeam::channel::TrySendError::*;
            match e {
                Full(t) => {
                    let mut state = self.state.lock();
                    state.clearqueue(&mut self.sender)?;
                    self.sender.try_send(t).unwrap();
                },
                Disconnected(_) => panic!("Disconnected."),
            }
        }
        Ok(())
    }
}

#[derive(Debug)]
pub enum RecvError {
    Disconnected,
    IOError(io::Error),
}

#[derive(Debug)]
pub enum TryRecvError {
    Empty,
    Disconnected,
    IOError(io::Error)
}

impl<T> Receiver<T> where T: Serialize + DeserializeOwned {
    pub fn recv(&mut self) -> Result<T, RecvError> {
        match self.receiver.try_recv() {
            Ok(t) => {
                return Ok(t)
            },
            Err(e) => {
                use crossbeam::channel::TryRecvError::*;
                match e {
                    // This will happen if we have stayed in memory and Sender is dropped
                    Disconnected => return Err(RecvError::Disconnected),
                    // We don't know if this is because we have taken all elements, or because some elements are on disk
                    Empty => {
                        let mut state = self.state.lock();
                        match state.read_wait() {
                            Ok(_) => return Ok(self.receiver.try_recv().expect("Internal error. read_wait() should only return Ok if items were sent to the receiver.")),
                            Err(ChannelError::InMemory) => {
                                drop(state);
                                match self.receiver.recv() {
                                    Ok(elem) => return Ok(elem),
                                    Err(crossbeam::channel::RecvError) => return Err(RecvError::Disconnected),
                                }
                            },
                            Err(ChannelError::Disconnected) => return Err(RecvError::Disconnected),
                            Err(ChannelError::IOError(e)) => return Err(RecvError::IOError(e)),
                            Err(ChannelError::Empty) => unreachable!(),
                        }
                    },
                }
            }
        }
    }

    pub fn try_recv(&mut self) -> Result<T, TryRecvError> {
        match self.receiver.try_recv() {
            Ok(t) => {
                return Ok(t)
            },
            Err(e) => {
                use crossbeam::channel::TryRecvError::*;
                match e {
                    // This will happen if we have stayed in memory and Sender is dropped
                    Disconnected => return Err(TryRecvError::Disconnected),
                    // We don't know if this is because we have taken all elements, or because some elements are on disk
                    Empty => {
                        let mut state = self.state.lock();
                        match state.read() {
                            Ok(_) => return Ok(self.receiver.try_recv().expect("Internal error. read() should only return Ok if items were sent to the receiver.")),
                            Err(ChannelError::InMemory) => return Err(TryRecvError::Empty),
                            Err(ChannelError::Empty) => return Err(TryRecvError::Empty),
                            Err(ChannelError::Disconnected) => return Err(TryRecvError::Disconnected),
                        Err(ChannelError::IOError(e)) => return Err(TryRecvError::IOError(e)),
                        }
                    },
                }
            }
        }
    }

    pub fn iter(&mut self) -> Iter<T> {
        Iter { receiver: self }
    }

    pub fn try_iter(&mut self) -> TryIter<T> {
        TryIter { receiver: self }
    }
}

pub struct TryIter<'a, T: 'a> where T: Serialize + DeserializeOwned {
    receiver: &'a mut Receiver<T>,
}

impl<'a, T> Iterator for TryIter<'a, T> where T: Serialize + DeserializeOwned {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.receiver.try_recv().ok()
    }
}

pub struct Iter<'a, T: 'a> where T: Serialize + DeserializeOwned {
    receiver: &'a mut Receiver<T>,
}

impl<'a, T> Iterator for Iter<'a, T> where T: Serialize + DeserializeOwned {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.receiver.recv().ok()
    }
}

pub struct IntoIter<T> where T: Serialize + DeserializeOwned {
    receiver: Receiver<T>,
}

impl<T> IntoIterator for Receiver<T> where T: Serialize + DeserializeOwned {
    type Item = T;
    type IntoIter = IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIter { receiver: self }
    }
}

impl<T> Iterator for IntoIter<T> where T: Serialize + DeserializeOwned {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.receiver.recv().ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate test;

    #[test]
    fn test_works() -> io::Result<()> {
        let (mut sender, mut receiver) = channel(50);
        for i in 0..175 {
            sender.send(i)?;
        }
        for i in 0..175 {
            assert_eq!(Some(i), receiver.try_recv().ok());
        }
        Ok(())
    }

    #[test]
    fn test_recv() -> io::Result<()> {
        let (mut sender, mut receiver) = channel(50);
        for i in 0..175 {
            sender.send(i)?;
        }
        for i in 0..175 {
            assert_eq!(Some(i), receiver.recv().ok());
        }
        Ok(())
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        let (mut sender, mut receiver) = channel(50);
        for i in 0..175 {
            sender.send(i)?;
        }
        let mut iter = receiver.iter();
        for i in 0..175 {
            assert_eq!(Some(i), iter.next());
        }
        Ok(())
    }

    #[test]
    fn test_inmemory() -> io::Result<()> {
        let (mut sender, mut receiver) = channel(50);
        for i in 0..175 {
            sender.send(i)?;
            assert_eq!(Some(i), receiver.try_recv().ok());
        }
        Ok(())
    }

    #[test]
    fn test_ondisk() -> io::Result<()> {
        let (mut sender, mut receiver) = channel(50);
        for i in 0..75 {
            sender.send(i)?;
        }
        let mut rcount = 0;
        for i in 75..175 {
            sender.send(i)?;
            assert_eq!(Some(rcount), receiver.try_recv().ok());
            rcount += 1;
        }
        while rcount < 175 {
            assert_eq!(Some(rcount), receiver.try_recv().ok());
            rcount += 1;
        }
        Ok(())
    }
}