#[cfg(feature = "read")]
pub(crate) mod bbiread;
#[cfg(feature = "write")]
pub(crate) mod bbiwrite;
#[cfg(feature = "write")]
pub mod bedchromdata;
#[cfg(feature = "read")]
pub(crate) mod bigbedread;
#[cfg(feature = "write")]
pub(crate) mod bigbedwrite;
#[cfg(feature = "read")]
pub(crate) mod bigwigread;
#[cfg(feature = "write")]
pub(crate) mod bigwigwrite;

#[cfg(feature = "write")]
use serde::{Deserialize, Serialize};

pub(crate) const BIGWIG_MAGIC: u32 = 0x888F_FC26;
pub(crate) const BIGBED_MAGIC: u32 = 0x8789_F2EB;

pub(crate) const CIR_TREE_MAGIC: u32 = 0x2468_ACE0;
pub(crate) const CHROM_TREE_MAGIC: u32 = 0x78CA_8C91;

/// Info on a specific zoom level in a bbi file
#[derive(Copy, Clone, Debug)]
pub struct ZoomHeader {
    pub reduction_level: u32,
    pub(crate) data_offset: u64,
    pub(crate) index_offset: u64,
    pub(crate) index_tree_offset: Option<u64>,
}

/// A single zoom item
#[derive(Copy, Clone, Debug)]
pub struct ZoomRecord {
    pub(crate) chrom: u32,
    pub start: u32,
    pub end: u32,
    pub summary: Summary,
}

/// A summary of a section of data (may be an entire file)
#[derive(Copy, Clone, Debug)]
pub struct Summary {
    pub total_items: u64,
    pub bases_covered: u64,
    pub min_val: f64,
    pub max_val: f64,
    pub sum: f64,
    pub sum_squares: f64,
}

/// Represents a single value in a bigWig file
#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "write", derive(Serialize, Deserialize))]
pub struct Value {
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

/// Represents a single entry in a bigBed file
#[derive(Clone, Debug, PartialEq)]
pub struct BedEntry {
    pub start: u32,
    pub end: u32,
    pub rest: String,
}

/// The type of bbi file
#[derive(Copy, Clone, Debug)]
pub enum BBIFile {
    BigWig,
    BigBed,
}

#[cfg(feature = "read")]
pub use bbiread::*;
#[cfg(feature = "write")]
pub use bbiwrite::*;
#[cfg(feature = "read")]
pub use bigbedread::*;
#[cfg(feature = "write")]
pub use bigbedwrite::*;
#[cfg(feature = "read")]
pub use bigwigread::*;
#[cfg(feature = "write")]
pub use bigwigwrite::*;
