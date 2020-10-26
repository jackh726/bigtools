use std::io;

pub trait ChromValues<V> {
    fn next(&mut self) -> Option<io::Result<V>>;
    fn peek(&mut self) -> Option<&V>;
}
