use std::io;

pub trait ChromGroups<V, C: ChromValues<V>> {
    fn next(&mut self) -> io::Result<Option<(String, C)>>;
    fn peek(&mut self) -> io::Result<Option<(String, C)>>;
}

pub trait ChromValues<V> {
    fn next(&mut self) -> Option<io::Result<V>>;
    fn peek(&mut self) -> Option<&V>;
}
