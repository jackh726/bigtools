pub mod bbiread;
pub mod bbiwrite;
pub mod bedchromdata;
pub mod bigbedread;
pub mod bigbedwrite;
pub mod bigwigread;
pub mod bigwigwrite;

use serde::{Deserialize, Serialize};

pub(crate) const BIGWIG_MAGIC: u32 = 0x888F_FC26;
pub(crate) const BIGBED_MAGIC: u32 = 0x8789_F2EB;

pub(crate) const CIR_TREE_MAGIC: u32 = 0x2468_ACE0;
pub(crate) const CHROM_TREE_MAGIC: u32 = 0x78CA_8C91;

#[derive(Copy, Clone, Debug)]
pub struct ZoomHeader {
    pub reduction_level: u32,
    pub(crate) data_offset: u64,
    pub(crate) index_offset: u64,
}

#[derive(Copy, Clone, Debug)]
pub struct ZoomRecord {
    pub chrom: u32,
    pub start: u32,
    pub end: u32,
    pub summary: Summary,
}

#[derive(Copy, Clone, Debug)]
pub struct Summary {
    pub total_items: u64,
    pub bases_covered: u64,
    pub min_val: f64,
    pub max_val: f64,
    pub sum: f64,
    pub sum_squares: f64,
}

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Value {
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BedEntry {
    pub start: u32,
    pub end: u32,
    pub rest: String,
}

#[derive(Copy, Clone, Debug)]
pub enum BBIFile {
    BigWig,
    BigBed,
}

pub use bbiread::*;
pub use bbiwrite::*;
pub use bigbedread::*;
pub use bigbedwrite::*;
pub use bigwigread::*;
pub use bigwigwrite::*;
