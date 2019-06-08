use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, Result, Seek, SeekFrom};

use byteordered::ByteOrdered;

use crate::bigwig::BigWigWriteOptions;
use crate::bigwig::ChromGroupReadStreamingIterator;
use crate::chromvalues::ChromValues;
use crate::bigwig::{BigWigRead, BigWigWrite};
use crate::bigwig::{Value as ValueSection, Value};
use crate::bigwig::ChromGroupRead;

use crate::idmap::IdMap;

/// Returns:
///  (val, None, None, overhang or None) when merging two does not break up one, and may or may not add an overhang (one.start == two.start)
///  (val, val, val or None, overhang or None) when merging two breaks up one, and may or may not add an overhang (one.start < two.start or one.end > two.end)
/// The overhang may equal the previous value
fn merge_into(one: ValueSection, two: ValueSection) -> (ValueSection, Option<ValueSection>, Option<ValueSection>, Option<ValueSection>) {
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
                ValueSection {
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
                ValueSection {
                    start: one.start,
                    end: one.end,
                    value: one.value + two.value,
                },
                None,
                None,
                Some(ValueSection {
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
                    ValueSection {
                        start: two.start,
                        end: two.end,
                        value: one.value + two.value,
                    },
                    Some(ValueSection {
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
                    ValueSection {
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
                    ValueSection {
                        start: one.start,
                        end: two.start,
                        value: one.value,
                    },
                    Some(ValueSection {
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
                    Some(ValueSection {
                        start: end,
                        end: two.end,
                        value: 0.0,
                    }),
                )
            } else if one.value == 0.0 {
                (
                    ValueSection {
                        start: one.start,
                        end: two.start,
                        value: 0.0,
                    },
                    Some(ValueSection {
                        start: two.start,
                        end: one.end,
                        value: two.value,
                    }),
                    None,
                    Some(ValueSection {
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
                    Some(ValueSection {
                        start: end,
                        end: two.end,
                        value: 0.0,
                    }),
                )
            } else {
                (
                    ValueSection {
                        start: one.start,
                        end: two.start,
                        value: one.value,
                    },
                    Some(ValueSection {
                        start: two.start,
                        end: one.end,
                        value: one.value + two.value,
                    }),
                    None,
                    Some(ValueSection {
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
                    ValueSection {
                        start: one.start,
                        end: two.start,
                        value: one.value,
                    },
                    Some(ValueSection {
                        start: two.start,
                        end: two.end,
                        value: one.value + two.value,
                    }),
                    Some(ValueSection {
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
                    ValueSection {
                        start: two.start,
                        end: one.start,
                        value: two.value,
                    },
                    Some(ValueSection {
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
                    ValueSection {
                        start: two.start,
                        end: one.start,
                        value: two.value,
                    },
                    Some(ValueSection {
                        start: one.start,
                        end: one.end,
                        value: one.value + two.value,
                    }),
                    None,
                    Some(ValueSection {
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
                    ValueSection {
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
                    Some(ValueSection {
                        start,
                        end: one.end,
                        value: one.value,
                    }),
                    None,
                    None,
                )
            } else if two.value == 0.0 {
                (
                    ValueSection {
                        start: two.start,
                        end: one.start,
                        value: 0.0,
                    },
                    Some(ValueSection {
                        start: one.start,
                        end: one.end,
                        value: one.value,
                    }),
                    None,
                    None,
                )
            } else {
                (
                    ValueSection {
                        start: two.start,
                        end: one.start,
                        value: two.value,
                    },
                    Some(ValueSection {
                        start: one.start,
                        end: two.end,
                        value: one.value + two.value,
                    }),
                    Some(ValueSection {
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

struct ValueSectionIter<I> where I : Iterator<Item=ValueSection> + std::marker::Send {
    sections: Vec<I>,
    queue: std::collections::VecDeque<ValueSection>,
    buffer: Option<Box<Iterator<Item = ValueSection> + std::marker::Send>>,
}

impl<I> Iterator for ValueSectionIter<I> where I : Iterator<Item=ValueSection> + std::marker::Send {
    type Item = ValueSection;

    fn next(&mut self) -> Option<ValueSection> {
        if let Some(buf) = &mut self.buffer {
            let next = buf.next();
            match next {
                None => self.buffer = None,
                Some(_) => return next,
            }
        }

        let queue = &mut self.queue;
        let mut out: Vec<ValueSection> = Vec::new();

        loop {
            //?println!("\nQueue {:?}", queue);
            let mut earliest_start = None;
            'vals: for section in self.sections.iter_mut() {
                let val = section.next();
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
                                let nvq = std::mem::replace(queued, ValueSection { start: 0, end: 0, value: 0.0 });
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
                                            let nvq = std::mem::replace(queued, ValueSection { start: 0, end: 0, value: 0.0 });
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

fn merge_sections_many<I>(sections: Vec<I>) -> impl Iterator<Item=ValueSection> + std::marker::Send where I : Iterator<Item=ValueSection> + std::marker::Send {    
    ValueSectionIter {
        sections: sections.into_iter().map(Box::new).collect(),
        queue: std::collections::VecDeque::new(),
        buffer: None,
    }
}

pub fn get_merged_values(bigwigs: Vec<BigWigRead>, options: BigWigWriteOptions) -> Result<impl ChromGroupReadStreamingIterator + std::marker::Send> {
    // Get sizes for each and check that all files (that have the chrom) agree
    // Check that all chrom sizes match for all files
    let mut chrom_sizes = BTreeMap::new();
    for chrom in bigwigs.iter().flat_map(BigWigRead::get_chroms).map(|c| c.name) {
        if chrom_sizes.get(&chrom).is_some() {
            continue;
        }
        let (sizes, bws): (Vec<_>, Vec<_>) = bigwigs.iter().map(|w| {
            let chroms = w.get_chroms();
            let res = chroms.iter().find(|v| v.name == chrom);
            match res {
                Some(s) => Some((s.length, w.clone())),
                None => None,
            }
        }).filter_map(|x| x).unzip();
        let size = sizes[0];
        if !sizes.iter().all(|s| *s == size) {
            eprintln!("Chrom '{:?}' had different sizes in the bigwig files. (Are you using the same assembly?)", chrom);
            return Err(io::Error::new(io::ErrorKind::Other, "Invalid input (nonmatching chroms)"));
        }

        chrom_sizes.insert(chrom.clone(), (size, bws));
    }

    let mut chrom_ids = IdMap::new();

    let pool = futures::executor::ThreadPoolBuilder::new().pool_size(6).create().expect("Unable to create thread pool.");

    let chrom_ids = chrom_sizes.iter().map(|(c, _)| chrom_ids.get_id(c.clone())).collect::<Vec<_>>().into_iter();

    struct MergingValues<I: Iterator<Item=Value> + std::marker::Send> {
        iter: std::iter::Peekable<I>,
    }

    impl<I: Iterator<Item=Value> + std::marker::Send> ChromValues for MergingValues<I> {
        fn next(&mut self) -> io::Result<Option<Value>> {
            Ok(self.iter.next())
        }

        fn peek(&mut self) -> Option<&Value> {
            self.iter.peek()
        }
    }

    struct ChromGroupReadStreamingIteratorImpl {
        pool: futures::executor::ThreadPool,
        options: BigWigWriteOptions,
        iter: Box<Iterator<Item=((String, (u32, Vec<BigWigRead>)), u32)> + std::marker::Send>,
    }

    impl ChromGroupReadStreamingIterator for ChromGroupReadStreamingIteratorImpl {
        fn next(&mut self) -> io::Result<Option<ChromGroupRead>> {
            let next = self.iter.next();
            match next {
                Some(((chrom, (size, bws)), chrom_id)) => {
                    let current_chrom = chrom.clone();

                    // Owned version of BigWigRead::get_interval
                    // TODO: how to please the borrow checker
                    // This doesn't work:
                    // let iters = bws.iter().map(|b| b.get_interval(&chrom, 1, size)).collect();

                    let iters: Vec<_> = bws.into_iter().map(move |b| {
                        let blocks = b.get_overlapping_blocks(&chrom, 1, size).unwrap();

                        let endianness = b.info.header.endianness;
                        let fp = File::open(b.path.clone()).unwrap();
                        let mut file = ByteOrdered::runtime(std::io::BufReader::new(fp), endianness);

                        if blocks.len() > 0 {
                            file.seek(SeekFrom::Start(blocks[0].offset)).unwrap();
                        }
                        let mut iter = blocks.into_iter().peekable();
                        
                        let block_iter = std::iter::from_fn(move || {
                            let next = iter.next();
                            let peek = iter.peek();
                            let next_offset = match peek {
                                None => None,
                                Some(peek) => Some(peek.offset),
                            };
                            match next {
                                None => None,
                                Some(next) => Some((next, next_offset))
                            }
                        });
                        let vals_iter = block_iter.flat_map(move |(block, next_offset)| {
                            // TODO: Could minimize this by chunking block reads
                            let vals = b.get_block_values(&mut file, &block).unwrap();
                            match next_offset {
                                None => (),
                                Some(next_offset) => {
                                    if next_offset != block.offset + block.size {
                                        file.seek(SeekFrom::Start(next_offset)).unwrap();
                                    }
                                }
                            }
                            vals
                        });
                        vals_iter
                    }).collect();

                    let mergingvalues = MergingValues { iter: merge_sections_many(iters).filter(|x| x.value != 0.0).peekable() };
                    Ok(Some(BigWigWrite::read_group(current_chrom, chrom_id, mergingvalues, self.pool.clone(), self.options.clone()).unwrap()))
                },
                None => {
                    return Ok(None)       
                },
            }
        }
    }

    let group_iter = ChromGroupReadStreamingIteratorImpl {
        pool: futures::executor::ThreadPoolBuilder::new().pool_size(4).create().expect("Unable to create thread pool."),
        options: options,
        iter: Box::new(chrom_sizes.into_iter().zip(chrom_ids)),
    };

    Ok(group_iter)

/*
    let all_values = chrom_sizes.into_iter().zip(chrom_ids).map(move |((chrom, (size, bws)), chrom_id)| {
        let current_chrom = chrom.clone();

        // Owned version of BigWigRead::get_interval
        // TODO: how to please the borrow checker
        // This doesn't work:
        // let iters = bws.iter().map(|b| b.get_interval(&chrom, 1, size)).collect();

        let iters: Vec<_> = bws.into_iter().map(move |b| {
            let blocks = b.get_overlapping_blocks(&chrom, 1, size).unwrap();

            let endianness = b.info.header.endianness;
            let fp = File::open(b.path.clone()).unwrap();
            let mut file = ByteOrdered::runtime(std::io::BufReader::new(fp), endianness);

            if blocks.len() > 0 {
                file.seek(SeekFrom::Start(blocks[0].offset)).unwrap();
            }
            let mut iter = blocks.into_iter().peekable();
            
            let block_iter = std::iter::from_fn(move || {
                let next = iter.next();
                let peek = iter.peek();
                let next_offset = match peek {
                    None => None,
                    Some(peek) => Some(peek.offset),
                };
                match next {
                    None => None,
                    Some(next) => Some((next, next_offset))
                }
            });
            let vals_iter = block_iter.flat_map(move |(block, next_offset)| {
                // TODO: Could minimize this by chunking block reads
                let vals = b.get_block_values(&mut file, &block).unwrap();
                match next_offset {
                    None => (),
                    Some(next_offset) => {
                        if next_offset != block.offset + block.size {
                            file.seek(SeekFrom::Start(next_offset)).unwrap();
                        }
                    }
                }
                vals
            });
            vals_iter
        }).collect();

        BigWigWrite::read_group(current_chrom, chrom_id, merge_sections_many(iters).filter(|x| x.value != 0.0), pool.clone(), options.clone()).unwrap()
    }); // Could be appended with `.collect::<Vec<_>>().into_iter();` to process all chroms in parallel
    // This has a huge overhead of file descriptors though
    Ok(all_values)
*/
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
        let one = ValueSection {
            start: 10,
            end: 20,
            value: 0.3,
        };
        let two = ValueSection {
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

    fn generate_sections_seq(start: u32, end: u32, seed: u64) -> Vec<ValueSection> {
        use rand::prelude::*;

        let mut out = vec![];

        let mut rng: StdRng = SeedableRng::seed_from_u64(seed);

        let mut curr = start;
        loop {
            let value: f32 = rng.gen();
            let size = (rng.gen::<f32>() * 20.0).floor() as u32 + 5;
            let skip = 0.max((rng.gen::<f32>() * 10.0).floor() as i32 + -7) as u32;

            let curr_end = end.min(curr + size);
            out.push(ValueSection {
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