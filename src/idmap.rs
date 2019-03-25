use std::cmp::Eq;
use std::hash::Hash;
use std::collections::HashMap;
use std::sync::RwLock;

#[derive(Debug)]
struct State<K : Hash + Eq> {
    map: HashMap<K, u32>,
    next_id: u32,
}

#[derive(Debug)]
pub struct IdMap<K : Hash + Eq> {
    state: RwLock<State<K>>,
}

impl<K : Hash + Eq> IdMap<K> {
    pub fn new() -> IdMap<K> {
        IdMap {
            state: RwLock::new(State {
                map: std::collections::HashMap::new(),
                next_id: 0,
            })
        }
    }

    pub fn get_map(self) -> HashMap<K, u32> {
        let state = self.state.into_inner().expect("Poisoned lock"); 
        state.map
    } 

    /// If the key already exists in the map, this will simply return the id for it.
    /// Otherwise, it locks a mutex and returns a new id.
    /// This means that in the case of missing keys, there are two map hits.
    pub fn get_id(&mut self, key: K) -> u32 {
        {
            let state = self.state.read().expect("Poisoned lock");
            if let Some(id) = state.map.get(&key) {
                return *id;
            }
        }
        {
            let mut state = self.state.write().expect("Poisoned lock");
            let next_id = state.next_id;
            let chrom_id: u32 = *state.map.entry(key).or_insert(next_id);
            // Was the map updated between
            if chrom_id == next_id {
                state.next_id += 1;
            }
            return chrom_id;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let mut idmap: IdMap<String> = IdMap::new();
        let zero = idmap.get_id("zero".to_string());
        assert!(zero == 0);
        let one = idmap.get_id("one".to_string());
        assert!(one == 1);
        assert!(idmap.get_id("one".to_string()) == 1);
        assert!(idmap.get_id("zero".to_string()) == 0);

        let map = idmap.get_map();
        assert!(map.len() == 2);
    }
}