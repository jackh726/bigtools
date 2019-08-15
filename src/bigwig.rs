pub(crate) const BIGWIG_MAGIC_LTH: u32 = 0x888F_FC26;
pub(crate) const BIGWIG_MAGIC_HTL: u32 = 0x26FC_8F88;
#[allow(dead_code)]
pub(crate) const BIGBED_MAGIC_LTH: u32 = 0x8789_F2EB;
#[allow(dead_code)]
pub(crate) const BIGBED_MAGIC_HTL: u32 = 0xEBF2_8987;

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

pub use crate::bigwigread::*;
pub use crate::bigwigwrite::*;
