
use std::io;
use crate::bigwig::Value;


pub trait ChromGroups<V: ChromValues> {
    fn next(&mut self) -> io::Result<Option<(String, V)>>;
}

pub trait ChromValues {
    fn next(&mut self) -> io::Result<Option<Value>>;
    fn peek(&mut self) -> Option<&Value>;
}