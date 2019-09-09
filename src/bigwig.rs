pub(crate) const BIGWIG_MAGIC: u32 = 0x888F_FC26;
pub(crate) const BIGBED_MAGIC: u32 = 0x8789_F2EB;

pub(crate) const CIR_TREE_MAGIC: u32 = 0x2468_ACE0;
pub(crate) const CHROM_TREE_MAGIC: u32 = 0x78CA_8C91;

#[derive(Clone, Debug)]
pub(crate) struct ZoomHeader {
    pub(crate) reduction_level: u32,
    pub(crate) data_offset: u64,
    pub(crate) index_offset: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Value {
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BedEntry {
    pub start: u32,
    pub end: u32,
    pub rest: String,
}

#[derive(Clone, Debug)]
pub enum BBIFile {
    BigWig,
    BigBed,
}

pub use crate::bbiread::*;
pub use crate::bbiwrite::*;

pub use crate::bigwigread::*;
pub use crate::bigwigwrite::*;

pub use crate::bigbedread::*;
pub use crate::bigbedwrite::*;
