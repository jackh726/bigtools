use std::io;

pub trait ChromGroups<V, C: ChromValues<V>> {
    fn next(&mut self) -> Option<io::Result<(String, C)>>;
    fn peek(&mut self) -> Option<io::Result<(String, C)>>;
}

pub trait ChromValues<V> {
    fn next(&mut self) -> Option<io::Result<V>>;
    fn peek(&mut self) -> Option<&V>;
}
