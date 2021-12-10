use std::io;

pub trait ChromValues {
    type V;

    fn next(&mut self) -> Option<io::Result<Self::V>>;
    fn peek(&mut self) -> Option<&Self::V>;
}
