use std::io::{Read, Seek, Cursor};

use bigtools::BigWigRead;
use bigtools::utils::reopen::SeekableRead;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern {
    pub fn alert(s: &str);
}

#[wasm_bindgen]
pub fn read_bigwig(data: Vec<u8>) {
    let data = Cursor::new(data);
    run(data);
}

/*
fn run() -> () {
    use bigtools::{BigWigWrite, Value};
    use bigtools::bedchromdata::BedParserStreamingIterator;
    use tokio::runtime;

    let runtime = runtime::Builder::new_current_thread()
        .build()
        .expect("Unable to create thread pool.");

    let bigwig = BigWigWrite::create_file("test.bigWig".to_owned());
    let chrom_map = [("chr1".to_owned(), 100)]
        .into_iter()
        .collect::<std::collections::HashMap<String, u32>>();
    let data: [(String, Value); 0] = [];
    let chsi = BedParserStreamingIterator::wrap_infallible_iter(data.into_iter(), true);
    match bigwig.write(chrom_map, chsi, runtime) {
        Err(e) => println!("{}", e),
        Ok(_) => {}
    }
}
*/

fn run<R: SeekableRead>(read: R) -> () {
    let mut bigwig = BigWigRead::open(read).unwrap();
    let chroms = bigwig.chroms().to_vec();
    let data = bigwig.get_interval("chr1", 0, u32::MAX).unwrap();
    let data: Vec<_> = data.collect();
    alert(&format!("{chroms:?}"));
    alert(&format!("{data:?}"));
}
