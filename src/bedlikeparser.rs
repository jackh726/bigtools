use std::collections::HashMap;
use std::hash::BuildHasher;
use std::io;

use futures::future::Either;

use crate::bigwig::ChromGroupRead;
use crate::bigwig::ChromGroupReadStreamingIterator;
use crate::bigwig::WriteGroupsError;
use crate::idmap::IdMap;
use crate::chromvalues::{ChromGroups, ChromValues};


pub type ChromGroupReadFunction<C> = Box<dyn Fn(String, u32, u32, C) -> io::Result<ChromGroupRead> + Send>;

pub struct BedGraphParserChromGroupStreamingIterator<V, C: ChromValues<V> + Send, G: ChromGroups<V, C>, H: BuildHasher> {
    chrom_groups: G,
    callable: ChromGroupReadFunction<C>,
    last_chrom: Option<String>,
    chrom_ids: Option<IdMap<String>>,
    chrom_map: HashMap<String, u32, H>,
    _v: std::marker::PhantomData<V>,
    _s: std::marker::PhantomData<C>,
}

impl<V, C: ChromValues<V> + Send, G: ChromGroups<V, C>, H: BuildHasher> BedGraphParserChromGroupStreamingIterator<V, C, G, H> {
    pub fn new(vals: G, chrom_map: HashMap<String, u32, H>, callable: ChromGroupReadFunction<C>) -> Self{
        BedGraphParserChromGroupStreamingIterator {
            chrom_groups: vals,
            callable,
            last_chrom: None,
            chrom_ids: Some(IdMap::default()),
            chrom_map,
            _v: std::marker::PhantomData,
            _s: std::marker::PhantomData,
        }
    }
}


impl<V, C: ChromValues<V> + Send, G: ChromGroups<V, C>, H: BuildHasher> ChromGroupReadStreamingIterator for BedGraphParserChromGroupStreamingIterator<V, C, G, H> {
    fn next(&mut self) -> Result<Option<Either<ChromGroupRead, (IdMap<String>)>>, WriteGroupsError> {
        match self.chrom_groups.next()? {
            Some((chrom, group)) => {
                let chrom_ids = self.chrom_ids.as_mut().unwrap();
                let last = self.last_chrom.replace(chrom.clone());
                if let Some(c) = last {
                    // TODO: test this correctly fails
                    if c >= chrom {
                        return Err(WriteGroupsError::InvalidInput("Input bedGraph not sorted by chromosome. Sort with `sort -k1,1 -k2,2n`.".to_string()));
                    }
                }
                let length = match self.chrom_map.get(&chrom) {
                    Some(length) => *length,
                    None => return Err(WriteGroupsError::InvalidInput(format!("Input bedGraph contains chromosome that isn't in the input chrom sizes: {}", chrom))),
                };
                let chrom_id = chrom_ids.get_id(chrom.clone());
                let group = (self.callable)(chrom, chrom_id, length, group)?;
                Ok(Some(Either::Left(group)))
            },
            None => {
                match self.chrom_ids.take() {
                    Some(chrom_ids) => Ok(Some(Either::Right(chrom_ids))),
                    None => Ok(None),
                }
            }
        }
    }
}
