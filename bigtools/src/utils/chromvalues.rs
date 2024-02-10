use std::error::Error;
use std::io;

pub trait ChromValues {
    type Value;
    type Error: Error + Send + From<io::Error> + 'static;

    fn next(&mut self) -> Option<Result<Self::Value, Self::Error>>;
    fn peek(&mut self) -> Option<Result<&Self::Value, &Self::Error>>;
}
