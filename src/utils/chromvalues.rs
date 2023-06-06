use std::io;

pub trait ChromValues {
    type Value;
    type Error: Send + From<io::Error>;

    fn next(&mut self) -> Option<Result<Self::Value, Self::Error>>;
    fn peek(&mut self) -> Option<Result<&Self::Value, &Self::Error>>;
}
