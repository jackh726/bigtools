use std::io;

use crate::bbi::Value;

struct FillValues<I>
where
    I: Iterator<Item = io::Result<Value>>,
{
    iter: I,
    last_val: Option<Value>,
    expected_end: Option<u32>,
    last_end: u32,
}

impl<I> Iterator for FillValues<I>
where
    I: Iterator<Item = io::Result<Value>> + Send,
{
    type Item = io::Result<Value>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(last) = self.last_val.take() {
            self.last_end = last.end;
            return Some(Ok(last));
        }
        let next = self.iter.next();
        match next {
            Some(Ok(next)) => {
                if next.start > self.last_end {
                    let last = self.last_end;
                    self.last_end = next.start;
                    self.last_val.replace(next);
                    Some(Ok(Value {
                        start: last,
                        end: self.last_end,
                        value: 0.0,
                    }))
                } else {
                    self.last_end = next.end;
                    Some(Ok(next))
                }
            }
            Some(_) => next,
            None => match self.expected_end {
                None => None,
                Some(expected_end) => {
                    if self.last_end < expected_end {
                        let last = self.last_end;
                        self.last_end = expected_end;
                        Some(Ok(Value {
                            start: last,
                            end: expected_end,
                            value: 0.0,
                        }))
                    } else {
                        None
                    }
                }
            },
        }
    }
}

/// Fills any space between `Value`s with `0.0`s.
/// Note: Output values will not be merged if any input Values are `0.0`
pub fn fill<I>(iter: I) -> impl Iterator<Item = io::Result<Value>> + Send
where
    I: Iterator<Item = io::Result<Value>> + Send,
{
    FillValues {
        iter,
        last_val: None,
        expected_end: None,
        last_end: 0,
    }
}

/// Fills any space between `Value`s with `0.0`s. This will also pad the start and end with `0.0`s if they do not exist.
/// Note: Output values will not be merged if any input Values are `0.0`
///
/// If the start > the end of the first value, it will be ignored.
pub fn fill_start_to_end<I>(
    iter: I,
    start: u32,
    end: u32,
) -> impl Iterator<Item = io::Result<Value>> + Send
where
    I: Iterator<Item = io::Result<Value>> + Send,
{
    FillValues {
        iter,
        last_val: None,
        expected_end: Some(end),
        last_end: start,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_fill() {
        let intervals: Vec<io::Result<Value>> = vec![
            Ok(Value {
                start: 10,
                end: 15,
                value: 0.5,
            }),
            Ok(Value {
                start: 20,
                end: 30,
                value: 0.7,
            }),
            Ok(Value {
                start: 30,
                end: 35,
                value: 0.9,
            }),
            Err(io::Error::new(io::ErrorKind::Other, "Test error")),
        ];

        let mut iter = fill(intervals.into_iter());
        assert_eq!(
            iter.next().unwrap().unwrap(),
            Value {
                start: 0,
                end: 10,
                value: 0.0
            }
        );
        assert_eq!(
            iter.next().unwrap().unwrap(),
            Value {
                start: 10,
                end: 15,
                value: 0.5
            }
        );
        assert_eq!(
            iter.next().unwrap().unwrap(),
            Value {
                start: 15,
                end: 20,
                value: 0.0
            }
        );
        assert_eq!(
            iter.next().unwrap().unwrap(),
            Value {
                start: 20,
                end: 30,
                value: 0.7
            }
        );
        assert_eq!(
            iter.next().unwrap().unwrap(),
            Value {
                start: 30,
                end: 35,
                value: 0.9
            }
        );
        assert!(iter.next().unwrap().is_err());
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_fill_start_to_end() {
        let intervals: Vec<io::Result<Value>> = vec![
            Ok(Value {
                start: 10,
                end: 15,
                value: 0.5,
            }),
            Ok(Value {
                start: 20,
                end: 30,
                value: 0.7,
            }),
            Ok(Value {
                start: 30,
                end: 35,
                value: 0.9,
            }),
            Err(io::Error::new(io::ErrorKind::Other, "Test error")),
        ];

        let mut iter = fill_start_to_end(intervals.into_iter(), 5, 50);
        assert_eq!(
            iter.next().unwrap().unwrap(),
            Value {
                start: 5,
                end: 10,
                value: 0.0
            }
        );
        assert_eq!(
            iter.next().unwrap().unwrap(),
            Value {
                start: 10,
                end: 15,
                value: 0.5
            }
        );
        assert_eq!(
            iter.next().unwrap().unwrap(),
            Value {
                start: 15,
                end: 20,
                value: 0.0
            }
        );
        assert_eq!(
            iter.next().unwrap().unwrap(),
            Value {
                start: 20,
                end: 30,
                value: 0.7
            }
        );
        assert_eq!(
            iter.next().unwrap().unwrap(),
            Value {
                start: 30,
                end: 35,
                value: 0.9
            }
        );
        assert!(iter.next().unwrap().is_err());
        assert_eq!(
            iter.next().unwrap().unwrap(),
            Value {
                start: 35,
                end: 50,
                value: 0.0
            }
        );
        assert!(iter.next().is_none());
    }
}
