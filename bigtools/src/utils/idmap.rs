use std::collections::HashMap;

#[derive(Debug, Default)]
pub struct IdMap {
    map: HashMap<String, u32>,
    next_id: u32,
}

impl IdMap {
    pub fn get_map(self) -> HashMap<String, u32> {
        self.map
    }

    /// If the key already exists in the map, this will simply return the id for it.
    /// Otherwise, it locks a mutex and returns a new id.
    /// This means that in the case of missing keys, there are two map hits.
    pub fn get_id(&mut self, key: &str) -> u32 {
        if let Some(id) = self.map.get(key) {
            return *id;
        }
        let next_id = self.next_id;
        self.next_id += 1;
        let chrom_id: u32 = *self.map.entry(key.to_string()).or_insert(next_id);
        chrom_id
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let mut idmap: IdMap = IdMap::default();
        let zero = idmap.get_id("zero");
        assert!(zero == 0);
        let one = idmap.get_id("one");
        assert!(one == 1);
        assert!(idmap.get_id("one") == 1);
        assert!(idmap.get_id("zero") == 0);

        let map = idmap.get_map();
        assert!(map.len() == 2);
    }
}
