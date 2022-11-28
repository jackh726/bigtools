pub trait ChromValues {
    type Value;
    type Error: Send;

    fn next(&mut self) -> Option<Result<Self::Value, Self::Error>>;
    fn peek(&mut self) -> Option<&Self::Value>;
}
