use bigtools::{BBIReadError, BedEntry, Value};
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// Converts sorted [`BedEntry`] intervals into non-overlapping coverage depth [`Value`] intervals.
///
/// This iterator consumes a stream of sorted, potentially overlapping BED entries and produces a
/// stream of non-overlapping intervals whose values represent the coverage depth (number of
/// overlapping entries) over each range. Adjacent intervals with the same depth are merged, so
/// the output is a run-length-encoded coverage track.
///
/// The input BED entries **must** be sorted by start position, then end position.
///
/// # Notes
///
/// Implements a 1-D sweep-line over endpoint events: each entry contributes a `+1` at `start` and
/// a `-1` at `end`; the running sum at any scan position is the coverage depth at that base.
///
/// Memory scales with the number of intervals overlapping the current scan position (not with the
/// total number of input entries). The internal heap holds one pending `end` event for each such
/// interval until the scan moves past it.
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
pub struct CoverageIterator<I: Iterator> {
    iter: std::iter::Peekable<I>,
    ends: BinaryHeap<Reverse<u32>>,
    coverage: u32,
    run_start: Option<u32>,
    finished: bool,
}

impl<I> CoverageIterator<I>
where
    I: Iterator<Item = Result<BedEntry, BBIReadError>>,
{
    pub fn new(iter: I) -> Result<Self, BBIReadError> {
        Ok(Self {
            iter: iter.peekable(),
            ends: BinaryHeap::new(),
            coverage: 0,
            run_start: None,
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

        loop {
            // Determine the next event position: the minimum of the next entry's start
            // (peeked from the iterator) and the smallest pending end in the heap.
            let next_start = match self.iter.peek() {
                Some(Ok(e)) => Some(e.start),
                Some(Err(_)) => {
                    // Surface the iterator error
                    let err = match self.iter.next() {
                        Some(Err(e)) => e,
                        _ => unreachable!(),
                    };
                    self.finished = true;
                    return Some(Err(err));
                }
                None => None,
            };
            let next_end = self.ends.peek().map(|&Reverse(p)| p);
            let pos = match (next_start, next_end) {
                (Some(s), Some(e)) => s.min(e),
                (Some(s), None) => s,
                (None, Some(e)) => e,
                (None, None) => {
                    self.finished = true;
                    return None;
                }
            };

            // Process all events at `pos`. Start events come from the iterator (in sort order);
            // end events come from the heap. Zero-length entries (start == end == pos) push an
            // end into the heap during the start drain that is then consumed by the end drain
            // in the same iteration -- so they contribute net 0 to delta, as expected.
            let mut delta: i32 = 0;

            // Drain start events at `pos`
            while matches!(self.iter.peek(), Some(Ok(e)) if e.start == pos) {
                let entry = match self.iter.next() {
                    Some(Ok(e)) => e,
                    _ => unreachable!(),
                };
                delta += 1;
                self.ends.push(Reverse(entry.end));
            }

            // Drain end events at `pos` (including any just pushed by zero-length entries)
            while matches!(self.ends.peek(), Some(&Reverse(p)) if p == pos) {
                self.ends.pop();
                delta -= 1;
            }

            // Update coverage. `checked_add_signed` surfaces malformed input (e.g. an entry
            // with start > end, which would queue a -1 event before its matching +1) as a
            // panic rather than silent wrong output.
            let old_cov = self.coverage;
            let new_cov = old_cov.checked_add_signed(delta).expect(
                "CoverageIterator: coverage underflow (likely malformed input: start > end)",
            );

            // If coverage didn't change at this event position, the current run continues;
            // no emit, no state update. This collapses adjacent same-depth intervals
            // (run-length encoding), matching the conventional "coverage track" semantics
            if old_cov == new_cov {
                continue;
            }

            // Close the previous run if there was one. Skip emit if `run_start >= pos`,
            // which only happens under malformed or out-of-sort-order input.
            let interval_to_emit = match self.run_start {
                Some(start) if old_cov > 0 && start < pos => Some(Value {
                    start,
                    end: pos,
                    value: old_cov as f32,
                }),
                _ => None,
            };

            self.coverage = new_cov;
            self.run_start = if new_cov > 0 { Some(pos) } else { None };

            if let Some(interval) = interval_to_emit {
                return Some(Ok(interval));
            }
        }
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

        // Adjacent intervals with the same coverage are merged: [10-40): coverage=1
        assert_eq!(intervals.len(), 1);
        assert_eq!(
            intervals[0],
            Value {
                start: 10,
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

        // Adjacent intervals [0,1) and [1,2) with the same coverage merge into [0,2)
        assert_eq!(intervals.len(), 2);
        assert_eq!(
            intervals[0],
            Value {
                start: 0,
                end: 2,
                value: 1.0
            }
        );
        assert_eq!(
            intervals[1],
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

    #[test]
    fn test_malformed_start_greater_than_end() {
        // An entry with start > end is malformed. The current algorithm doesn't crash
        // (the end is only pushed to the heap after the start is consumed, so the +1 always
        // precedes the -1 and coverage stays non-negative). It emits no intervals because
        // the end position is less than prev_pos, failing the emit check. This test pins
        // that behavior so a future change that errors out is visible.
        let entries = vec![Ok(BedEntry {
            start: 20,
            end: 10,
            rest: "malformed".to_string(),
        })];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        assert_eq!(intervals, vec![]);
    }

    #[test]
    fn test_sort_violation_produces_incorrect_output() {
        // Out-of-order input is a documented precondition violation. The algorithm
        // doesn't detect it and produces incorrect (but non-panicking) output. This test
        // pins the current behavior so any future fix (e.g. error or debug assertion) is
        // visible as a test change.
        //
        // Correct output for sorted [(5,15), (10,20)] would be:
        //   [5-10): cov=1, [10-15): cov=2, [15-20): cov=1
        // With the entries presented out of order, we instead get a collapsed range at
        // cov=2 covering the entire [5, 15) span.
        let entries = vec![
            Ok(BedEntry {
                start: 10,
                end: 20,
                rest: "first".to_string(),
            }),
            Ok(BedEntry {
                start: 5,
                end: 15,
                rest: "out_of_order".to_string(),
            }),
        ];

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        assert_eq!(
            intervals,
            vec![
                Value {
                    start: 5,
                    end: 15,
                    value: 2.0,
                },
                Value {
                    start: 15,
                    end: 20,
                    value: 1.0,
                },
            ]
        );
    }

    #[test]
    fn test_many_same_position() {
        // 1000 entries all starting at the same position. Stresses the start-drain loop
        // which has to consume all of them in a single outer iteration. Output should be
        // a single interval at the full depth.
        let n = 1000u32;
        let entries: Vec<Result<BedEntry, BBIReadError>> = (0..n)
            .map(|_| {
                Ok(BedEntry {
                    start: 100,
                    end: 200,
                    rest: String::new(),
                })
            })
            .collect();

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        assert_eq!(
            intervals,
            vec![Value {
                start: 100,
                end: 200,
                value: n as f32,
            }]
        );
    }

    #[test]
    fn test_stress_many_overlapping() {
        // 10k entries with high pairwise overlap. Stresses the heap (peak depth ~ `width`)
        // and verifies correctness via algorithmic invariants:
        //   1. Total "base-coverage units" emitted equals the sum of input entry lengths.
        //   2. Intervals are non-overlapping, sorted, and contiguous (no zero-coverage gaps
        //      in a region covered by at least one entry).
        //   3. The full plateau region (positions `width-1 .. n`) is at coverage == `width`.
        //
        // We deliberately don't compare to a brute-force RLE reference because the iterator
        // emits at every event position even when coverage doesn't change between events
        // (which a coverage-based RLE would merge).
        let n: u32 = 10_000;
        let width: u32 = 100;

        let entries: Vec<Result<BedEntry, BBIReadError>> = (0..n)
            .map(|i| {
                Ok(BedEntry {
                    start: i,
                    end: i + width,
                    rest: String::new(),
                })
            })
            .collect();

        let coverage_iter = CoverageIterator::new(entries.into_iter()).unwrap();
        let result: Result<Vec<Value>, _> = coverage_iter.collect();
        let intervals = result.unwrap();

        // Invariant 1: total base-coverage units == sum of input entry lengths
        let total_base_cov: u64 = intervals
            .iter()
            .map(|v| (v.end - v.start) as u64 * v.value as u64)
            .sum();
        let expected_total: u64 = n as u64 * width as u64;
        assert_eq!(total_base_cov, expected_total);

        // Invariant 2: intervals are sorted, contiguous, all coverage > 0
        for v in &intervals {
            assert!(v.value > 0.0, "zero-coverage interval emitted");
            assert!(v.start < v.end, "degenerate interval emitted");
        }
        for w in intervals.windows(2) {
            assert_eq!(
                w[0].end, w[1].start,
                "non-contiguous intervals (no entry has cov=0 in this fixture)"
            );
        }

        // Invariant 3: plateau coverage equals `width`
        let mid_pos = n / 2;
        let at_mid = intervals
            .iter()
            .find(|v| v.start <= mid_pos && v.end > mid_pos)
            .expect("no interval covers mid_pos");
        assert_eq!(at_mid.value, width as f32);
    }
}
