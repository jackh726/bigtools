use std::fs::File;
use std::io::{self, Cursor, Read, Seek, SeekFrom};
use std::sync::Arc;

use crossbeam_channel::{bounded, Receiver as ChannelReceiver, Sender as ChannelSender};

use parking_lot::Mutex;

use bincode;

use serde::de::DeserializeOwned;
use serde::Serialize;

use crate::tell::Tell;

// TODO: check if crossbeam_channel need or if std::sync::mpsc okay
// TODO: validate maxsize = 0
// TODO: add an AtomicBoolean to indicate is_in_memory, to allow Receiver to not need to lock state (and update description)

// A general description of how a filebufferedchannel works
//
// The underlying implementation uses crossbeam_channel::bounded channel pairs to coordinate inter-struct passing of elements. On initial
// creation of the filebufferedchannel pair, a pair of crossbeam_channel::bounded pairs are created, which allows elements to be sent
// via the `crossbeam_channel::Sender`.
//
// There are two states to how a channel sender/receiver pair can be in: in memory or disk buffered. In short, the in memory representation
// of this channel is a simple proxy to an inner `crossbeam_channel::bounded` pair, while the disk buffered representation uses an
// `crossbeam_channel::bounded` pair to double-buffer elements.
//
// In the in memory state, the two halves of the channel delegate sending elements to the inner initial crossbeam channel. When the maximum
// capacity of the inner chanel is reached, then the channel must be switched to buffering in a temporary file. During the `send` method
// call when no elements can be added to the inner channel, a new channel and a temporary file is created. The current
//`crossbeam_channel::Sender` in `Sender` is replaced with the newly created `crossbeam_channel::Sender` and the old
// `crossbeam_channel::Sender`and newly crated `crossbeam_channel::Receiver` are stored in a `ChannelState`, which is
// stored behind an `Arc<Mutex<ChannelState>>` for both `Sender` and `Receiver`. Also stored in this ChannelState is the temporary file
// itself, as well as information regarding the last read and write position of that temporary file. At this point the state is whiched to
// file buffered, and the element currently being sent is send on the new, empty `crossbeam_channel::Sender` in `Sender`. During this state,
// the `Receiver` can simply check if there are new elements. For now, if there are no new elements state is locked and polled for possible
// new elements.
//
// In the file buffered state, the two halves of the channel form a double buffer for elements. When the `Sender` side of the channel fills,
// then the `Arc<Mutex<ChannelState>>` is locked, and elements are written to disk. When the `Receiver` side of the channel is empty, then
// the `Arc<Mutex<ChannelState>>` is locked and any elements are read from disk (up to `maxsize`), and then `Receiver` receives the next
// element. Importantly, the lock to the inner `ChannelState` is only locked if the `Sender` fills or the `Receiver` is empty.

/// Creates a filebufferedchannel sender/receiver pair.
pub fn channel<T>(size: usize) -> (Sender<T>, Receiver<T>)
where
    T: Serialize + DeserializeOwned,
{
    let state = Arc::new(Mutex::new(ChannelState::new(size)));
    let (channelsender, channelreceiver) = bounded(size);
    let sender = Sender {
        state: state.clone(),
        sender: channelsender,
    };
    let receiver = Receiver {
        state,
        receiver: channelreceiver,
    };
    (sender, receiver)
}

/// Creates a filebufferedchannel lazy sender/receiver pair. Unlike `channel`,
/// no data will be sent to the receiver until `recv` or `try_recv` is called.
/// They still will be queued to a maximum `size`. However, if the `Sender`
/// calls `flush`, then all queued items will be written to disk.
pub fn lazy_channel<T>(size: usize) -> io::Result<(Sender<T>, Receiver<T>)>
where
    T: Serialize + DeserializeOwned,
{
    let (channelsender, channelreceiver) = bounded(size);
    let (buffersender, bufferreceiver) = bounded(size);
    let state = Arc::new(Mutex::new(ChannelState::new_lazy(
        size,
        channelsender,
        bufferreceiver,
    )?));
    let sender = Sender {
        state: state.clone(),
        sender: buffersender,
    };
    let receiver = Receiver {
        state,
        receiver: channelreceiver,
    };
    Ok((sender, receiver))
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
        can_send: bool,
        serialize_size: Option<u64>,
        sender: ChannelSender<T>,
        buffer: ChannelReceiver<T>,
        file: File,
        readindex: u64,
        writeindex: u64,
    },
}

struct ChannelState<T>
where
    T: Serialize + DeserializeOwned,
{
    maxsize: usize,
    status: ChannelStateStatus<T>,
}

type ChannelResult<R> = Result<R, ChannelError>;

impl<T> ChannelState<T>
where
    T: Serialize + DeserializeOwned,
{
    fn new(maxsize: usize) -> ChannelState<T> {
        ChannelState {
            maxsize,
            status: ChannelStateStatus::InMemory,
        }
    }

    fn new_lazy(
        maxsize: usize,
        sender: ChannelSender<T>,
        buffer: ChannelReceiver<T>,
    ) -> io::Result<ChannelState<T>> {
        let status = ChannelStateStatus::OnDisk {
            can_send: false,
            serialize_size: None,
            sender,
            buffer,
            file: tempfile::tempfile()?,
            readindex: 0,
            writeindex: 0,
        };
        Ok(ChannelState { maxsize, status })
    }

    /// If we spawned this channel lazily, then no data has been sent to the
    /// receiver (all has been stored on disk or in the shared receiving queue).
    /// If we are in `OnDisk` state, this sets the `can_send` flag as `true`.
    fn mark_read(&mut self) {
        match &mut self.status {
            ChannelStateStatus::InMemory => {}
            ChannelStateStatus::OnDisk { can_send, .. } => {
                *can_send = true;
            }
        }
    }

    /// This will free *some* elements from `Sender`'s `ChannelSender`, either by writing them to the disk (if state is `OnDisk`),
    /// or by converting the `filebufferedchannel` from `InMemory` to `OnDisk` and, in the process, allocating a new buffer channel
    ///
    /// It's not guaranteed that `Sender`'s `ChannelSender` will be empty when this function returns (since new elements may be added)
    fn clearqueue(&mut self, sender: &mut ChannelSender<T>) -> ChannelResult<()> {
        // Since we already have the lock, let's push until the `Receiver`'s buffer is full
        match self.read() {
            Ok(_) => {}
            Err(ChannelError::InMemory) => {}
            Err(ChannelError::Disconnected) => return Err(ChannelError::Disconnected),
            Err(ChannelError::IOError(e)) => return Err(ChannelError::IOError(e)),
            Err(ChannelError::Empty) => {}
        }

        match &mut self.status {
            ChannelStateStatus::InMemory => {
                let (channelsender, channelreceiver) = bounded(self.maxsize);
                let sender = std::mem::replace(sender, channelsender);
                self.status = ChannelStateStatus::OnDisk {
                    can_send: true,
                    serialize_size: None,
                    sender,
                    buffer: channelreceiver,
                    file: tempfile::tempfile()?,
                    readindex: 0,
                    writeindex: 0,
                };
                Ok(())
            }
            ChannelStateStatus::OnDisk {
                buffer,
                file,
                writeindex,
                serialize_size,
                ..
            } => {
                // Decide up front the number of items to serialize and write to disk
                // This allows us to preallocate a vector to store the bytes, since we don't/can't write buffer the file
                let n = buffer.len();
                if n == 0 {
                    return Ok(());
                }
                let mut bytes = match &serialize_size {
                    Some(size) => Cursor::new(vec![0u8; (*size as usize) * n]),
                    None => Cursor::new(vec![]),
                };
                for item in buffer.try_iter().take(n) {
                    if serialize_size.is_none() {
                        serialize_size.replace(bincode::serialized_size(&item).unwrap());
                    }
                    bincode::serialize_into(&mut bytes, &item).unwrap();
                }
                file.seek(SeekFrom::Start(*writeindex))?;
                bytes.seek(SeekFrom::Start(0))?;
                io::copy(&mut bytes, file)?;
                *writeindex = file.tell()?;
                Ok(())
            }
        }
    }

    /// Result::Ok indicates that at least one element was read and sent to `Receiver` or that `Receiver` was full
    fn read(&mut self) -> ChannelResult<()> {
        match &mut self.status {
            ChannelStateStatus::InMemory => Err(ChannelError::InMemory),
            ChannelStateStatus::OnDisk {
                ref can_send,
                file,
                readindex,
                writeindex,
                buffer,
                sender,
                serialize_size,
            } => {
                if !can_send {
                    return Ok(());
                }
                if sender.capacity().unwrap() - sender.len() == 0 {
                    return Ok(());
                }

                // At this point there are two places elements can be:
                // 1) On disk
                // 2) In `buffer`, which has elements that were received from `Sender`, but haven't been buffered to disk yet.
                // We read each place in order until `buffer` is full or we run out of elements

                let mut sent = false;

                if *readindex < *writeindex {
                    // We need to seek to the read position
                    file.seek(SeekFrom::Start(*readindex))?;

                    let size = serialize_size.unwrap() as usize;
                    let diskn = (*writeindex - *readindex) as usize / size;
                    let n = std::cmp::min(diskn, sender.capacity().unwrap() - sender.len());
                    let mut buf = vec![0u8; n * size];
                    file.read_exact(&mut buf)?;
                    *readindex = file.tell()?;
                    let mut buf = Cursor::new(buf);
                    for _ in 0..n {
                        let elem = bincode::deserialize_from(&mut buf)
                            .expect("Error while deserializing.");
                        if let Err(e) = sender.try_send(elem) {
                            use crossbeam_channel::TrySendError::*;
                            match e {
                                Disconnected(_) => return Err(ChannelError::Disconnected),
                                Full(_) => {
                                    // We pre-checked for the capacity-len
                                    panic!("Buffer should not be full.");
                                }
                            }
                        }
                        sent = true;
                    }
                }
                let n = sender.capacity().unwrap() - sender.len();
                for elem in buffer.try_iter().take(n) {
                    if let Err(e) = sender.try_send(elem) {
                        use crossbeam_channel::TrySendError::*;
                        match e {
                            Disconnected(_) => return Err(ChannelError::Disconnected),
                            Full(_) => {
                                // We only take max number of remaining elements
                                panic!("Buffer should not be full.");
                            }
                        }
                    }
                    sent = true;
                }
                if !sent {
                    Err(ChannelError::Empty)
                } else {
                    Ok(())
                }
            }
        }
    }

    fn read_wait(&mut self) -> ChannelResult<()> {
        match self.read() {
            Ok(_) => Ok(()),
            Err(ChannelError::InMemory) => Err(ChannelError::InMemory),
            Err(ChannelError::Disconnected) => Err(ChannelError::Disconnected),
            Err(ChannelError::IOError(e)) => Err(ChannelError::IOError(e)),
            Err(ChannelError::Empty) => match &mut self.status {
                ChannelStateStatus::InMemory => unreachable!(),
                ChannelStateStatus::OnDisk { buffer, sender, .. } => {
                    use crossbeam_channel::RecvError;
                    match buffer.recv() {
                        Ok(elem) => {
                            if let Err(e) = sender.try_send(elem) {
                                use crossbeam_channel::TrySendError::*;
                                match e {
                                    Disconnected(_) => return Err(ChannelError::Disconnected),
                                    Full(_) => {
                                        panic!("Buffer should not be full.");
                                    }
                                }
                            }
                            Ok(())
                        }
                        Err(RecvError) => Err(ChannelError::Disconnected),
                    }
                }
            },
        }
    }
}

pub struct Sender<T>
where
    T: Serialize + DeserializeOwned,
{
    state: Arc<Mutex<ChannelState<T>>>,
    sender: ChannelSender<T>,
}

pub struct Receiver<T>
where
    T: Serialize + DeserializeOwned,
{
    state: Arc<Mutex<ChannelState<T>>>,
    receiver: ChannelReceiver<T>,
}

#[derive(Debug)]
pub enum SendError {
    Disconnected,
    IOError(io::Error),
}

impl From<ChannelError> for SendError {
    fn from(error: ChannelError) -> Self {
        match error {
            ChannelError::InMemory => unreachable!(),
            ChannelError::Disconnected => SendError::Disconnected,
            ChannelError::IOError(e) => SendError::IOError(e),
            ChannelError::Empty => unreachable!(),
        }
    }
}

impl<T> Sender<T>
where
    T: Serialize + DeserializeOwned,
{
    /// Sends a `T` over the channel. If state is `InMemory`, it will immediately
    /// be available to the `Receiver`. If state is `OnDisk`, then it will be
    /// sent the pending queue.
    ///
    /// If the item cannot be sent because the is full, then the channel will first be
    /// flushed.
    pub fn send(&mut self, t: T) -> Result<(), SendError> {
        if let Err(e) = self.sender.try_send(t) {
            use crossbeam_channel::TrySendError::*;
            match e {
                Full(t) => {
                    let mut state = self.state.lock();
                    state.clearqueue(&mut self.sender)?;
                    drop(state);
                    self.sender.try_send(t).unwrap();
                }
                Disconnected(_) => return Err(SendError::Disconnected),
            }
        }
        Ok(())
    }

    /// Forces the channel to be flushed. Any queued items will first be sent
    /// the `Receiver` (if possible), then written to disk. If there are no
    /// queued items, this is essentially a no-op (but it does still acquire
    /// a lock). If the channel was created a lazy channel and no read has yet
    /// occured, then all data sent over the channel will be on disk.
    pub fn flush(&mut self) -> Result<(), SendError> {
        let mut state = self.state.lock();
        Ok(state.clearqueue(&mut self.sender)?)
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
    IOError(io::Error),
}

impl<T> Receiver<T>
where
    T: Serialize + DeserializeOwned,
{
    pub fn recv(&mut self) -> Result<T, RecvError> {
        match self.receiver.try_recv() {
            Ok(t) => Ok(t),
            Err(e) => {
                use crossbeam_channel::TryRecvError::*;
                match e {
                    // This will happen if we have stayed in memory and Sender is dropped
                    Disconnected => Err(RecvError::Disconnected),
                    // We don't know if this is because we have taken all elements, or because some elements are on disk
                    Empty => {
                        let mut state = self.state.lock();
                        // If we were lazy, need to mark can_send as true
                        state.mark_read();
                        let read_wait = state.read_wait();
                        drop(state);
                        match read_wait {
                            Ok(_) => Ok(self.receiver.try_recv().expect("Internal error. read_wait() should only return Ok if items were sent to the receiver.")),
                            Err(ChannelError::InMemory) => {
                                match self.receiver.recv() {
                                    Ok(elem) => Ok(elem),
                                    Err(crossbeam_channel::RecvError) => Err(RecvError::Disconnected),
                                }
                            },
                            Err(ChannelError::Disconnected) => Err(RecvError::Disconnected),
                            Err(ChannelError::IOError(e)) => Err(RecvError::IOError(e)),
                            Err(ChannelError::Empty) => unreachable!(),
                        }
                    }
                }
            }
        }
    }

    pub fn try_recv(&mut self) -> Result<T, TryRecvError> {
        match self.receiver.try_recv() {
            Ok(t) => Ok(t),
            Err(e) => {
                use crossbeam_channel::TryRecvError::*;
                match e {
                    // This will happen if we have stayed in memory and Sender is dropped
                    Disconnected => Err(TryRecvError::Disconnected),
                    // We don't know if this is because we have taken all elements, or because some elements are on disk
                    Empty => {
                        let mut state = self.state.lock();
                        state.mark_read();
                        let read = state.read();
                        drop(state);
                        match read {
                            Ok(_) => Ok(self.receiver.try_recv().expect("Internal error. read() should only return Ok if items were sent to the receiver.")),
                            Err(ChannelError::InMemory) => Err(TryRecvError::Empty),
                            Err(ChannelError::Empty) => Err(TryRecvError::Empty),
                            Err(ChannelError::Disconnected) => Err(TryRecvError::Disconnected),
                            Err(ChannelError::IOError(e)) => Err(TryRecvError::IOError(e)),
                        }
                    }
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

pub struct TryIter<'a, T: 'a>
where
    T: Serialize + DeserializeOwned,
{
    receiver: &'a mut Receiver<T>,
}

impl<'a, T> Iterator for TryIter<'a, T>
where
    T: Serialize + DeserializeOwned,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.receiver.try_recv().ok()
    }
}

pub struct Iter<'a, T: 'a>
where
    T: Serialize + DeserializeOwned,
{
    receiver: &'a mut Receiver<T>,
}

impl<'a, T> Iterator for Iter<'a, T>
where
    T: Serialize + DeserializeOwned,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.receiver.recv().ok()
    }
}

pub struct IntoIter<T>
where
    T: Serialize + DeserializeOwned,
{
    receiver: Receiver<T>,
}

impl<T> IntoIterator for Receiver<T>
where
    T: Serialize + DeserializeOwned,
{
    type Item = T;
    type IntoIter = IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIter { receiver: self }
    }
}

impl<T> Iterator for IntoIter<T>
where
    T: Serialize + DeserializeOwned,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.receiver.recv().ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_works() -> Result<(), SendError> {
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
    fn test_recv() -> Result<(), SendError> {
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
    fn test_iter() -> Result<(), SendError> {
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
    fn test_into_iter() -> Result<(), SendError> {
        let (mut sender, receiver) = channel(50);
        let mut receiver_iter = receiver.into_iter();
        for i in 0..75 {
            sender.send(i)?;
        }
        let mut rcount = 0;
        for i in 75..175 {
            sender.send(i)?;
            assert_eq!(Some(rcount), receiver_iter.next());
            rcount += 1;
        }
        while rcount < 175 {
            assert_eq!(Some(rcount), receiver_iter.next());
            rcount += 1;
        }
        Ok(())
    }

    #[test]
    fn test_inmemory() -> Result<(), SendError> {
        let (mut sender, mut receiver) = channel(50);
        for i in 0..175 {
            sender.send(i)?;
            assert_eq!(Some(i), receiver.try_recv().ok());
        }
        Ok(())
    }

    #[test]
    fn test_ondisk() -> Result<(), SendError> {
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

    #[test]
    fn test_lazy() -> Result<(), SendError> {
        {
            let (mut sender, mut receiver) = lazy_channel(50).unwrap();
            for i in 0..45 {
                sender.send(i)?;
            }
            for i in 0..45 {
                assert_eq!(Some(i), receiver.try_recv().ok());
            }
        }
        {
            let (mut sender, mut receiver) = lazy_channel(50).unwrap();
            for i in 0..95 {
                sender.send(i)?;
            }
            for i in 0..95 {
                assert_eq!(Some(i), receiver.try_recv().ok());
            }
        }
        {
            let (mut sender, mut receiver) = lazy_channel(50).unwrap();
            for i in 0..45 {
                sender.send(i)?;
            }
            sender.flush()?;
            for i in 0..45 {
                assert_eq!(Some(i), receiver.try_recv().ok());
            }
        }
        {
            let (mut sender, _receiver) = lazy_channel(50).unwrap();
            for i in 0..45 {
                sender.send(i)?;
            }
            sender.flush()?;
            // TODO: test zero items in memory
        }
        Ok(())
    }
}
