pub trait Tell {
    /// Gets the current position
    fn tell(&mut self) -> std::io::Result<u64>;
}

impl<S: std::io::Seek> Tell for S {
    fn tell(&mut self) -> std::io::Result<u64> {
        self.seek(std::io::SeekFrom::Current(0))
    }
}
