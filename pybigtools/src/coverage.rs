use bigtools::{BBIReadError, BedEntry, Value};
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// Converts sorted [`BedEntry`] intervals into non-overlapping coverage depth [`Value`] intervals.
///
/// This iterator consumes a stream of sorted, potentially overlapping BED entries and produces a
/// stream of non-overlapping intervals whose values represent the coverage depth (number of
/// overlapping entries) over its genomic positions.
///
/// The input BED entries **must** be sorted by start position, then end position.
///
/// # Notes
///
/// The iterator uses roughly constant memory by maintaining a small heap that only grows if many
/// entries share the same start or end position.
///
/// # Example
///
/// ```rust
/// use bigtools::{BedEntry, Value};
/// use pybigtools::coverage::CoverageIterator;
///
/// let entries = vec![
///     Ok(BedEntry { start: 10, end: 30, rest: "".to_string() }),
///     Ok(BedEntry { start: 20, end: 40, rest: "".to_string() }),
/// ];
///
/// let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
/// let intervals: Result<Vec<Value>, _> = coverage_iter.collect();
///
/// // Result: [10-20): coverage=1, [20-30): coverage=2, [30-40): coverage=1
/// ```
pub struct CoverageIterator<I> {
    iter: I,
    heap: BinaryHeap<Reverse<(u32, i32)>>,
    coverage: i32,
    prev_pos: Option<u32>,
    finished: bool,
}

impl<I> CoverageIterator<I>
where
    I: Iterator<Item = Result<BedEntry, BBIReadError>>,
{
    pub fn new(mut iter: I) -> Result<Self, BBIReadError> {
        let mut heap = BinaryHeap::new();

        // Initialize the heap with two entries (four positions)
        match iter.next() {
            Some(Ok(entry)) => {
                heap.push(Reverse((entry.start, 1)));
                heap.push(Reverse((entry.end, -1)));
            }
            Some(Err(e)) => {
                return Err(e);
            }
            None => {}
        }

        Ok(Self {
            iter,
            heap,
            coverage: 0,
            prev_pos: None,
            finished: false,
        })
    }
}

impl<I> Iterator for CoverageIterator<I>
where
    I: Iterator<Item = Result<BedEntry, BBIReadError>>,
{
    type Item = Result<Value, BBIReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        while !self.heap.is_empty() {
            let Reverse((pos, inc)) = self.heap.pop().unwrap();

            // Batch all events at the same position
            // First, collect all same-position events already in the heap
            let mut delta = inc;
            while let Some(&Reverse((next_pos, inc))) = self.heap.peek() {
                if next_pos == pos {
                    delta += inc;
                    self.heap.pop();
                } else {
                    break;
                }
            }

            // If we have more entries with the current position, load them.
            while let Some(result) = self.iter.next() {
                match result {
                    Ok(entry) => {
                        self.heap.push(Reverse((entry.start, 1)));
                        self.heap.push(Reverse((entry.end, -1)));

                        // If this entry contributes to the current position, process those events
                        let has_current_pos = entry.start == pos || entry.end == pos;
                        if has_current_pos {
                            while let Some(&Reverse((next_pos, inc))) = self.heap.peek() {
                                if next_pos == pos {
                                    delta += inc;
                                    self.heap.pop();
                                } else {
                                    break;
                                }
                            }
                        } else {
                            // No more events at current position, stop loading
                            break;
                        }
                    }
                    Err(e) => {
                        self.finished = true;
                        return Some(Err(e));
                    }
                }
            }

            // Check if we need to emit an interval
            let interval_to_emit = if let Some(start) = self.prev_pos {
                if self.coverage > 0 && start < pos {
                    Some(Value {
                        start,
                        end: pos,
                        value: self.coverage as f32,
                    })
                } else {
                    None
                }
            } else {
                None
            };

            // Update coverage and current position with all batched events
            self.coverage += delta;
            self.prev_pos = Some(pos);

            // Return the interval if we have one
            if let Some(interval) = interval_to_emit {
                return Some(Ok(interval));
            }
        }

        self.finished = true;
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_streaming_simple() {
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 20,
                rest: "entry1".to_string(),
            }),
            Ok(BedEntry {
                start: 30,
                end: 40,
                rest: "entry2".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        assert_eq!(intervals.len(), 2);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
                end: 20,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 30,
                end: 40,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_streaming_overlapping() {
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 30,
                rest: "entry1".to_string(),
            }),
            Ok(BedEntry {
                start: 20,
                end: 40,
                rest: "entry2".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        assert_eq!(intervals.len(), 3);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
                end: 20,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 20,
                end: 30,
                value: 2.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 30,
                end: 40,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_streaming_complex() {
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 25,
                rest: "entry1".to_string(),
            }),
            Ok(BedEntry {
                start: 15,
                end: 35,
                rest: "entry2".to_string(),
            }),
            Ok(BedEntry {
                start: 30,
                end: 45,
                rest: "entry3".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        assert_eq!(intervals.len(), 5);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
                end: 15,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 15,
                end: 25,
                value: 2.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 25,
                end: 30,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[3],
            Value {
                start: 30,
                end: 35,
                value: 2.0
            }
        );
        assert_eq!(
            intervals[4],
            Value {
                start: 35,
                end: 45,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_same_start_positions() {
        // Test multiple entries with the same start position (sorted by start, then by end)
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 20,
                rest: "entry1".to_string(),
            }),
            Ok(BedEntry {
                start: 10,
                end: 25,
                rest: "entry3".to_string(),
            }),
            Ok(BedEntry {
                start: 10,
                end: 30,
                rest: "entry2".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        // Expected: [10-20): coverage=3, [20-25): coverage=2, [25-30): coverage=1
        assert_eq!(intervals.len(), 3);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
                end: 20,
                value: 3.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 20,
                end: 25,
                value: 2.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 25,
                end: 30,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_same_end_positions() {
        // Test multiple entries with the same end position (already sorted by start)
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 30,
                rest: "entry1".to_string(),
            }),
            Ok(BedEntry {
                start: 15,
                end: 30,
                rest: "entry2".to_string(),
            }),
            Ok(BedEntry {
                start: 20,
                end: 30,
                rest: "entry3".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        // Expected: [10-15): coverage=1, [15-20): coverage=2, [20-30): coverage=3
        assert_eq!(intervals.len(), 3);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
                end: 15,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 15,
                end: 20,
                value: 2.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 20,
                end: 30,
                value: 3.0
            }
        );
    }

    #[test]
    fn test_identical_intervals() {
        // Test multiple identical intervals
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 20,
                rest: "entry1".to_string(),
            }),
            Ok(BedEntry {
                start: 10,
                end: 20,
                rest: "entry2".to_string(),
            }),
            Ok(BedEntry {
                start: 10,
                end: 20,
                rest: "entry3".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        // Expected: [10-20): coverage=3
        assert_eq!(intervals.len(), 1);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
                end: 20,
                value: 3.0
            }
        );
    }

    #[test]
    fn test_adjacent_intervals() {
        // Test intervals that are adjacent (end of one = start of next)
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 20,
                rest: "entry1".to_string(),
            }),
            Ok(BedEntry {
                start: 20,
                end: 30,
                rest: "entry2".to_string(),
            }),
            Ok(BedEntry {
                start: 30,
                end: 40,
                rest: "entry3".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        // Expected: [10-20): coverage=1, [20-30): coverage=1, [30-40): coverage=1
        assert_eq!(intervals.len(), 3);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
                end: 20,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 20,
                end: 30,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 30,
                end: 40,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_empty_input() {
        let entries: Vec<Result<BedEntry, BBIReadError>> = vec![];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        assert_eq!(intervals.len(), 0);
    }

    #[test]
    fn test_single_entry() {
        let entries = vec![Ok(BedEntry {
            start: 5,
            end: 15,
            rest: "single".to_string(),
        })];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        assert_eq!(intervals.len(), 1);
        assert_eq!(
            intervals[0],
            Value {
                start: 5,
                end: 15,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_nested_intervals() {
        // One interval completely inside another
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 30,
                rest: "outer".to_string(),
            }),
            Ok(BedEntry {
                start: 15,
                end: 25,
                rest: "inner".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        // Expected: [10-15): coverage=1, [15-25): coverage=2, [25-30): coverage=1
        assert_eq!(intervals.len(), 3);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
                end: 15,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 15,
                end: 25,
                value: 2.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 25,
                end: 30,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_mixed_same_positions() {
        // Some entries start where others end at the same position
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 20,
                rest: "entry1".to_string(),
            }),
            Ok(BedEntry {
                start: 15,
                end: 20,
                rest: "entry2".to_string(),
            }),
            Ok(BedEntry {
                start: 20,
                end: 30,
                rest: "entry3".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        // Expected: [10-15): coverage=1, [15-20): coverage=2, [20-30): coverage=1
        assert_eq!(intervals.len(), 3);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
                end: 15,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 15,
                end: 20,
                value: 2.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 20,
                end: 30,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_many_overlapping() {
        // Five intervals all overlapping each other
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 20,
                rest: "1".to_string(),
            }),
            Ok(BedEntry {
                start: 12,
                end: 22,
                rest: "2".to_string(),
            }),
            Ok(BedEntry {
                start: 14,
                end: 24,
                rest: "3".to_string(),
            }),
            Ok(BedEntry {
                start: 16,
                end: 26,
                rest: "4".to_string(),
            }),
            Ok(BedEntry {
                start: 18,
                end: 28,
                rest: "5".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        // Should produce multiple intervals with increasing then decreasing coverage
        assert!(intervals.len() > 5); // More intervals than input entries
        assert!(intervals.iter().any(|v| v.value >= 4.0)); // Should reach high coverage

        // Verify no gaps and proper ordering
        for i in 1..intervals.len() {
            assert_eq!(
                intervals[i - 1].end,
                intervals[i].start,
                "Gap detected between intervals"
            );
        }
    }

    #[test]
    fn test_single_base_intervals() {
        // Test intervals that are exactly 1 base pair long
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 11,
                rest: "1bp".to_string(),
            }),
            Ok(BedEntry {
                start: 12,
                end: 13,
                rest: "1bp".to_string(),
            }),
            Ok(BedEntry {
                start: 14,
                end: 15,
                rest: "1bp".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        assert_eq!(intervals.len(), 3);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
                end: 11,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 12,
                end: 13,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 14,
                end: 15,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_large_gaps() {
        // Test intervals with large gaps between them
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 20,
                rest: "first".to_string(),
            }),
            Ok(BedEntry {
                start: 1000,
                end: 1010,
                rest: "second".to_string(),
            }),
            Ok(BedEntry {
                start: 100000,
                end: 100010,
                rest: "third".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        assert_eq!(intervals.len(), 3);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
                end: 20,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 1000,
                end: 1010,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 100000,
                end: 100010,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_start_at_zero() {
        // Test intervals starting at position 0
        let entries = vec![
            Ok(BedEntry {
                start: 0,
                end: 10,
                rest: "zero_start".to_string(),
            }),
            Ok(BedEntry {
                start: 5,
                end: 15,
                rest: "overlapping".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        assert_eq!(intervals.len(), 3);
        assert_eq!(
            intervals[0],
            Value {
                start: 0,
                end: 5,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 5,
                end: 10,
                value: 2.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 10,
                end: 15,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_boundary_positions() {
        // Test intervals with position 0 and various boundary conditions
        let entries = vec![
            Ok(BedEntry {
                start: 0,
                end: 1,
                rest: "zero".to_string(),
            }),
            Ok(BedEntry {
                start: 1,
                end: 2,
                rest: "one".to_string(),
            }),
            Ok(BedEntry {
                start: 4294967290,
                end: 4294967295,
                rest: "high".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        assert_eq!(intervals.len(), 3);
        assert_eq!(
            intervals[0],
            Value {
                start: 0,
                end: 1,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 1,
                end: 2,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 4294967290,
                end: 4294967295,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_very_long_interval() {
        // Test an interval spanning a very large genomic region
        let entries = vec![
            Ok(BedEntry {
                start: 1000,
                end: 1000000,
                rest: "long".to_string(),
            }),
            Ok(BedEntry {
                start: 500000,
                end: 500010,
                rest: "short_inside".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        // Should handle large intervals correctly
        assert_eq!(intervals.len(), 3);
        assert_eq!(
            intervals[0],
            Value {
                start: 1000,
                end: 500000,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 500000,
                end: 500010,
                value: 2.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 500010,
                end: 1000000,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_zero_length_intervals() {
        // Test intervals with zero length (start == end)
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 10,
                rest: "zero_length1".to_string(),
            }),
            Ok(BedEntry {
                start: 15,
                end: 20,
                rest: "normal".to_string(),
            }),
            Ok(BedEntry {
                start: 25,
                end: 25,
                rest: "zero_length2".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        // Zero-length intervals should be ignored (produce no coverage)
        // Only the normal interval [15, 20) should produce coverage
        assert_eq!(intervals.len(), 1);
        assert_eq!(
            intervals[0],
            Value {
                start: 15,
                end: 20,
                value: 1.0
            }
        );
    }

    #[test]
    fn test_zero_length_with_overlaps() {
        // Test zero-length intervals mixed with overlapping normal intervals
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 20,
                rest: "normal1".to_string(),
            }),
            Ok(BedEntry {
                start: 15,
                end: 15,
                rest: "zero_length".to_string(),
            }),
            Ok(BedEntry {
                start: 15,
                end: 25,
                rest: "normal2".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        // Zero-length interval should not affect coverage calculation
        // Should get: [10-15): coverage=1, [15-20): coverage=2, [20-25): coverage=1
        assert_eq!(intervals.len(), 3);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
                end: 15,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
            Value {
                start: 15,
                end: 20,
                value: 2.0
            }
        );
        assert_eq!(
            intervals[2],
            Value {
                start: 20,
                end: 25,
                value: 1.0
            }
        );
    }
}
