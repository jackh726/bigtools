use std::io;

use crate::bigwig::Value;


/// Returns:
///  (val, None, None, overhang or None) when merging two does not break up one, and may or may not add an overhang (one.start == two.start)
///  (val, val, val or None, overhang or None) when merging two breaks up one, and may or may not add an overhang (one.start < two.start or one.end > two.end)
/// The overhang may equal the previous value
/// 
/// # Panics
/// Panics if the two Values do not overlap.
pub fn merge_into(one: Value, two: Value) -> (Value, Option<Value>, Option<Value>, Option<Value>) {
    if one.end <= two.start {
        panic!("No overlap.");
    }
    if one.start == two.start {
        // |--
        // |--
        if one.end == two.end {
            // |---|
            // |---|
            (
                Value {
                    start: one.start,
                    end: one.end,
                    value: one.value + two.value,
                },
                None,
                None,
                None,
            )
        } else if one.end < two.end {
            // |--|
            // |---|
            (
                Value {
                    start: one.start,
                    end: one.end,
                    value: one.value + two.value,
                },
                None,
                None,
                Some(Value {
                    start: one.end,
                    end: two.end,
                    value: two.value,
                }),
            )
        } else {
            // |---|
            // |--|
            if two.value == 0.0 {
                (
                    one,
                    None,
                    None,
                    None,
                )
            } else {
                (
                    Value {
                        start: two.start,
                        end: two.end,
                        value: one.value + two.value,
                    },
                    Some(Value {
                        start: two.end,
                        end: one.end,
                        value: one.value,
                    }),
                    None,
                    None,
                )
            }
        }
    } else if one.start < two.start {
        // |--
        //  |--
        if one.end == two.end {
            // |---|
            //  |--|
            if two.value == 0.0 {
                (
                    Value {
                        start: one.start,
                        end: one.end,
                        value: one.value,
                    },
                    None,
                    None,
                    None,
                )
            } else {
                (
                    Value {
                        start: one.start,
                        end: two.start,
                        value: one.value,
                    },
                    Some(Value {
                        start: two.start,
                        end: two.end,
                        value: one.value + two.value,
                    }),
                    None,
                    None,
                )
            }
        } else if one.end < two.end {
            // |---|
            //  |---|
            if one.value == 0.0 && two.value == 0.0 {
                let end = one.end;
                (
                    one,
                    None,
                    None,
                    Some(Value {
                        start: end,
                        end: two.end,
                        value: 0.0,
                    }),
                )
            } else if one.value == 0.0 {
                (
                    Value {
                        start: one.start,
                        end: two.start,
                        value: 0.0,
                    },
                    Some(Value {
                        start: two.start,
                        end: one.end,
                        value: two.value,
                    }),
                    None,
                    Some(Value {
                        start: one.end,
                        end: two.end,
                        value: two.value,
                    }),
                )
            } else if two.value == 0.0 {
                let end = one.end;
                (
                    one,
                    None,
                    None,
                    Some(Value {
                        start: end,
                        end: two.end,
                        value: 0.0,
                    }),
                )
            } else {
                (
                    Value {
                        start: one.start,
                        end: two.start,
                        value: one.value,
                    },
                    Some(Value {
                        start: two.start,
                        end: one.end,
                        value: one.value + two.value,
                    }),
                    None,
                    Some(Value {
                        start: one.end,
                        end: two.end,
                        value: two.value,
                    }),
                )
            }
        } else {
            // |----|
            //  |--|
            if two.value == 0.0 {
                (
                    one,
                    None,
                    None,
                    None,
                )
            } else {
                (
                    Value {
                        start: one.start,
                        end: two.start,
                        value: one.value,
                    },
                    Some(Value {
                        start: two.start,
                        end: two.end,
                        value: one.value + two.value,
                    }),
                    Some(Value {
                        start: two.end,
                        end: one.end,
                        value: one.value,
                    }),
                    None,
                )
            }
        }
    } else {
        //  |--
        // |--
        if one.end == two.end {
            //  |--|
            // |---|
            if one.value == 0.0 {
                (
                    two,
                    None,
                    None,
                    None,
                )
            } else {
                (
                    Value {
                        start: two.start,
                        end: one.start,
                        value: two.value,
                    },
                    Some(Value {
                        start: one.start,
                        end: one.end,
                        value: one.value + two.value,
                    }),
                    None,
                    None,
                )
            }
        } else if one.end < two.end {
            //  |--|
            // |----|
            if one.value == 0.0 {
                (
                    two,
                    None,
                    None,
                    None,
                )
            } else {
                (
                    Value {
                        start: two.start,
                        end: one.start,
                        value: two.value,
                    },
                    Some(Value {
                        start: one.start,
                        end: one.end,
                        value: one.value + two.value,
                    }),
                    None,
                    Some(Value {
                        start: one.end,
                        end: two.end,
                        value: two.value,
                    }),
                )
            }
        } else {
            //  |---|
            // |---|
            if one.value == 0.0 && two.value == 0.0 {
                (
                    Value {
                        start: two.start,
                        end: one.end,
                        value: 0.0,
                    },
                    None,
                    None,
                    None,
                )
            } else if one.value == 0.0 {
                let start = two.end;
                (
                    two,
                    Some(Value {
                        start,
                        end: one.end,
                        value: one.value,
                    }),
                    None,
                    None,
                )
            } else if two.value == 0.0 {
                (
                    Value {
                        start: two.start,
                        end: one.start,
                        value: 0.0,
                    },
                    Some(Value {
                        start: one.start,
                        end: one.end,
                        value: one.value,
                    }),
                    None,
                    None,
                )
            } else {
                (
                    Value {
                        start: two.start,
                        end: one.start,
                        value: two.value,
                    },
                    Some(Value {
                        start: one.start,
                        end: two.end,
                        value: one.value + two.value,
                    }),
                    Some(Value {
                        start: two.end,
                        end: one.end,
                        value: one.value,
                    }),
                    None,
                )
            }
        }
    }
}

struct ValueIter<I> where I : Iterator<Item=io::Result<Value>> + Send {
    error: io::Result<()>,
    sections: Vec<I>,
    queue: std::collections::VecDeque<Value>,
    buffer: Option<Box<Iterator<Item = Value> + Send>>,
}

impl<I> Iterator for ValueIter<I> where I : Iterator<Item=io::Result<Value>> + Send {
    type Item = Value;

    fn next(&mut self) -> Option<Value> {
        if let Some(buf) = &mut self.buffer {
            let next = buf.next();
            match next {
                None => self.buffer = None,
                Some(_) => return next,
            }
        }

        let queue = &mut self.queue;
        let mut out: Vec<Value> = Vec::new();

        loop {
            //?println!("\nQueue {:?}", queue);
            let mut earliest_start = None;
            'vals: for section in self.sections.iter_mut() {
                let val_result = section.next();
                let val = match val_result {
                    Some(Ok(x)) => Some(x),
                    Some(Err(e)) => {
                        self.error = Err(e);
                        None
                    }
                    None => None,
                };
                match val {
                    None => (),
                    Some(next_val) => {
                        earliest_start = match earliest_start {
                            None => Some(next_val.start),
                            Some(e) => Some(e.min(next_val.start)),
                        };

                        //?println!("q {:?} next_val {:?}", queue, next_val);
                        if queue.is_empty() || queue[queue.len() - 1].end <= next_val.start {
                            queue.push_back(next_val);
                        } else {
                            for (idx, queued) in queue.iter_mut().enumerate() {
                                if next_val.end <= queued.start {
                                    queue.insert(idx, next_val);
                                    continue 'vals;
                                }
                                if queued.end <= next_val.start {
                                    continue;
                                }
                                let nvq = std::mem::replace(queued, Value { start: 0, end: 0, value: 0.0 });
                                //?println!("Merging {:?} {:?}", nvq, nvo);
                                let (one, two, three, overhang) = merge_into(nvq, next_val);
                                //?println!("merge_into {:?} {:?} {:?} {:?}", one, two, three, overhang);
                                std::mem::replace(queued, one);

                                if let Some(th) = three {
                                    queue.insert(idx + 1, th);
                                }   
                                if let Some(tw) = two {
                                    queue.insert(idx + 1, tw);
                                }

                                let mut last_overhang = overhang;
                                'overhang: while let Some(o) = last_overhang.take() {
                                    //?println!("q {:?}", queue);
                                    //?println!("Overhang (inner): {:?}", o);
                                    if queue.is_empty() || queue[queue.len() - 1].end <= o.start {
                                        queue.push_back(o);
                                    } else {
                                        for (idx, queued) in queue.iter_mut().enumerate() {
                                            if o.end <= queued.start {
                                                queue.insert(idx, o);
                                                continue 'vals;
                                            }
                                            if queued.end <= o.start {
                                                continue;
                                            }
                                            let nvq = std::mem::replace(queued, Value { start: 0, end: 0, value: 0.0 });
                                            //?println!("Merging {:?} {:?}", nvq, o);
                                            let (one, two, three, overhang) = merge_into(nvq, o);
                                            //?println!("merge_into {:?} {:?} {:?} {:?}", one, two, three, overhang);
                                            std::mem::replace(queued, one);
                                            last_overhang = overhang;


                                            if let Some(th) = three {
                                                queue.insert(idx + 1, th);
                                            }   
                                            if let Some(tw) = two {
                                                queue.insert(idx + 1, tw);
                                            }

                                            continue 'overhang;
                                        }
                                    }
                                }
                                continue 'vals;
                            }
                            unreachable!();
                        }
                    }
                }
            }
            match earliest_start {
                None => {
                    out.extend(queue.drain(..));
                    break;
                },
                Some(start) => {
                    //?println!("earliest {:?}", start);
                    while !queue.is_empty() {
                        //?println!("q {:?}", queue);
                        let f = &queue[0];
                        if f.end > start {
                            break;
                        }
                        out.push(queue.pop_front().unwrap());
                    }
                    //?println!("out {:?}", out);
                    //?println!("q {:?}", queue);
                    if !out.is_empty() {
                        break;
                    }
                }
            }
        }

        self.buffer = Some(Box::new(out.into_iter()));
        self.buffer.as_mut().unwrap().next()
    }
}

pub fn merge_sections_many<I>(sections: Vec<I>) -> impl Iterator<Item=Value> + Send where I : Iterator<Item=io::Result<Value>> + Send {    
    ValueIter {
        // TODO: this isn't used right now
        error: Ok(()),
        sections: sections,
        queue: std::collections::VecDeque::new(),
        buffer: None,
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    extern crate test;

    #[test]
    fn test_merge_many() {
        let first = generate_sections_seq(50, 150, 1234);
        let second = generate_sections_seq(50, 150, 12345);
        println!("Running merge many with: \n{:?} \n{:?}", first, second);
        let merged = merge_sections_many(vec![first.into_iter(), second.into_iter()]).collect::<Vec<_>>();
        println!("\nMerged (many): {:?}\n", merged);
    }

    #[test]
    fn test_merge_into() {
        let one = Value {
            start: 10,
            end: 20,
            value: 0.3,
        };
        let two = Value {
            start: 12,
            end: 18,
            value: 0.5,
        };
        let (one, two, three, overhang) = merge_into(one, two);
        println!("merge_into: {:?} {:?} {:?} {:?}", one, two, three, overhang);
    }

    #[bench]
    fn bench_merge_many(b: &mut test::Bencher) {
        let first = generate_sections_seq(50, 555550, 1234);
        let second = generate_sections_seq(50, 555550, 12345);
        b.iter(|| {
            let merged = merge_sections_many(vec![first.clone().into_iter(), second.clone().into_iter()]);
            merged.for_each(drop);
        });
    }

    #[test]
    fn can_gen() {
        let _sections = generate_sections_seq(50, 150, 1234);
    }

    fn generate_sections_seq(start: u32, end: u32, seed: u64) -> Vec<Value> {
        use rand::prelude::*;

        let mut out = vec![];

        let mut rng: StdRng = SeedableRng::seed_from_u64(seed);

        let mut curr = start;
        loop {
            let value: f32 = rng.gen();
            let size = (rng.gen::<f32>() * 20.0).floor() as u32 + 5;
            let skip = 0.max((rng.gen::<f32>() * 10.0).floor() as i32 + -7) as u32;

            let curr_end = end.min(curr + size);
            out.push(Value {
                start: curr,
                end: curr_end,
                value,
            });
            if end <= curr_end + skip {
                break;
            } else {
                curr = curr + size + skip;
            }
        }
        out
    }    
}
