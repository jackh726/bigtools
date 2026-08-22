use std::borrow::BorrowMut;
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::vec::Vec;

use byteordered::Endianness;
use bytes::{Buf, BytesMut};
use itertools::Either;
use libdeflater::Decompressor;
use smallvec::{smallvec, SmallVec};
use thiserror::Error;

use crate::bbi::{
    BBIFile, Summary, ZoomHeader, ZoomRecord, BIGBED_MAGIC, BIGWIG_MAGIC, CHROM_TREE_MAGIC,
    CIR_TREE_MAGIC,
};
use crate::bed::bedparser::BedValueError;
use crate::utils::reopen::{Reopen, ReopenableFile, SeekableRead};
use crate::{BigBedRead, BigWigRead};

use self::internal::BBIReadInternal;

const MAX_CIR_TREE_LAYOUT_NODES: u64 = 1_000_000;
const MAX_CIR_TREE_LAYOUT_BLOCKS: u64 = 10_000_000;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct Block {
    pub(crate) offset: u64,
    pub(crate) size: u64,
}

/// Coordinate and on-disk size metadata for one primary BBI data block.
///
/// The chromosome identifiers correspond to [`ChromInfo::id`]. A block can
/// span a chromosome boundary, so both its start and end chromosome are
/// reported. Reading this layout only traverses the cir-tree index; it does not
/// read or decompress the primary data blocks.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct BBIDataBlock {
    pub start_chrom_id: u32,
    pub start_base: u32,
    pub end_chrom_id: u32,
    pub end_base: u32,
    /// Byte offset of the encoded data block in the BBI file.
    pub offset: u64,
    /// Number of bytes occupied by the encoded data block in the BBI file.
    pub data_size: u64,
}

/// Safety limits for a complete primary-data cir-tree traversal.
///
/// The defaults allow 1,000,000 index nodes and 10,000,000 primary data
/// blocks. Callers reading larger trusted files can raise either limit through
/// [`BigWigRead::data_blocks_with_limits`] or
/// [`BigBedRead::data_blocks_with_limits`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct BBIDataBlockLimits {
    /// Maximum number of cir-tree nodes visited by the traversal.
    pub max_nodes: u64,
    /// Maximum number of primary data blocks returned by the traversal.
    pub max_blocks: u64,
}

impl Default for BBIDataBlockLimits {
    fn default() -> Self {
        Self {
            max_nodes: MAX_CIR_TREE_LAYOUT_NODES,
            max_blocks: MAX_CIR_TREE_LAYOUT_BLOCKS,
        }
    }
}

impl Block {
    pub fn size(&self) -> u64 {
        self.size
    }
}

/// Header info for a bbi file
///
/// Note that info on internal properties like file offsets are not public.
/// Reading data is available through higher-level functions.
#[derive(Copy, Clone, Debug)]
pub struct BBIHeader {
    pub endianness: Endianness,
    pub version: u16,
    pub field_count: u16,
    pub defined_field_count: u16,

    pub(crate) zoom_levels: u16,
    pub(crate) chromosome_tree_offset: u64,
    pub(crate) full_data_offset: u64,
    pub(crate) full_index_offset: u64,
    pub(crate) full_index_tree_offset: Option<u64>,
    pub(crate) full_index_item_count: Option<u64>,
    pub(crate) auto_sql_offset: u64,
    pub(crate) total_summary_offset: u64,
    pub(crate) uncompress_buf_size: u32,
}

impl BBIHeader {
    pub fn is_compressed(&self) -> bool {
        self.uncompress_buf_size > 0
    }

    pub fn primary_data_size(&self) -> u64 {
        self.full_index_offset - self.full_data_offset
    }

    pub(crate) fn primary_data_bounds(
        &self,
        filetype: BBIFile,
    ) -> Result<(u64, u64), BBIReadError> {
        // BigWig stores a 32-bit section count; BigBed stores a 64-bit item
        // count. Primary blocks begin after that format-specific prefix.
        let item_count_size = match filetype {
            BBIFile::BigWig => 4,
            BBIFile::BigBed => 8,
        };
        let start = self
            .full_data_offset
            .checked_add(item_count_size)
            .ok_or_else(|| BBIReadError::InvalidFile("primary data start overflows".into()))?;
        // Kent-generated files place the chromosome tree before primary data,
        // while BigTools can place it between primary data and its cir-tree.
        // In the latter layout it is the tighter upper bound for data blocks.
        let end = if self.chromosome_tree_offset > self.full_data_offset
            && self.chromosome_tree_offset < self.full_index_offset
        {
            self.chromosome_tree_offset
        } else {
            self.full_index_offset
        };
        if start > end {
            return Err(BBIReadError::InvalidFile(
                "primary data region is invalid".into(),
            ));
        }
        Ok((start, end))
    }
}

/// Information on a chromosome in a bbi file
#[derive(Clone, Debug)]
pub struct ChromInfo {
    pub name: String,
    pub length: u32,
    pub(crate) id: u32,
}

impl ChromInfo {
    /// Return the chromosome identifier used by the BBI cir-tree index.
    pub fn id(&self) -> u32 {
        self.id
    }
}

impl PartialEq for ChromInfo {
    fn eq(&self, other: &ChromInfo) -> bool {
        self.name == other.name
    }
}

/// Info on a bbi file
#[derive(Clone, Debug)]
pub struct BBIFileInfo {
    /// The type of the bbi file - either a bigBed or a bigWig
    pub filetype: BBIFile,
    /// Header info
    pub header: BBIHeader,
    /// Info on zooms in the bbi file
    pub zoom_headers: Vec<ZoomHeader>,
    /// The chromosome info the bbi file is based on
    pub chrom_info: Vec<ChromInfo>,
}

pub(crate) struct ChromIdNotFound(pub(crate) String);

impl From<ChromIdNotFound> for BBIReadError {
    fn from(e: ChromIdNotFound) -> Self {
        BBIReadError::InvalidChromosome(e.0)
    }
}

impl BBIFileInfo {
    pub(crate) fn chrom_id(&self, chrom_name: &str) -> Result<u32, ChromIdNotFound> {
        let chrom_info = &self.chrom_info;
        let chrom = chrom_info.iter().find(|&x| x.name == chrom_name);
        match chrom {
            Some(c) => Ok(c.id),
            None => Err(ChromIdNotFound(chrom_name.to_owned())),
        }
    }
}

#[derive(Error, Debug)]
pub(crate) enum BBIFileReadInfoError {
    #[error("Invalid magic (likely not a BigWig or BigBed file)")]
    UnknownMagic,
    #[error("Invalid chromosomes section")]
    InvalidChroms,
    #[error("Error occurred: {}", .0)]
    IoError(#[from] io::Error),
}

#[derive(Error, Debug)]
pub enum CirTreeSearchError {
    #[error("The passed chromosome ({}) was incorrect.", .0)]
    InvalidChromosome(String),
    #[error("Error occurred: {}", .0)]
    IoError(#[from] io::Error),
}

/// Possible errors encountered when reading a bbi file
#[derive(Error, Debug)]
pub enum BBIReadError {
    #[error("The passed chromosome ({}) was incorrect.", .0)]
    InvalidChromosome(String),
    #[error("Invalid magic (likely a bug).")]
    UnknownMagic,
    #[error("The file was invalid: {}", .0)]
    InvalidFile(String),
    #[error("Error parsing bed-like data.")]
    BedValueError(#[from] BedValueError),
    #[error("Error occurred: {}", .0)]
    IoError(#[from] io::Error),
}

impl BBIReadError {
    /// Whether a complete primary-data layout traversal stopped at a caller
    /// safety limit rather than because the file was malformed.
    pub fn is_data_block_traversal_limit_exceeded(&self) -> bool {
        matches!(
            self,
            Self::InvalidFile(message)
                if message == "cir-tree node limit exceeded"
                    || message == "cir-tree data block limit exceeded"
        )
    }
}

impl From<CirTreeSearchError> for BBIReadError {
    fn from(value: CirTreeSearchError) -> Self {
        match value {
            CirTreeSearchError::InvalidChromosome(chrom) => BBIReadError::InvalidChromosome(chrom),
            CirTreeSearchError::IoError(e) => BBIReadError::IoError(e),
        }
    }
}

impl From<internal::FullDataCirTreeError> for BBIReadError {
    fn from(value: internal::FullDataCirTreeError) -> Self {
        match value {
            internal::FullDataCirTreeError::UnknownMagic => BBIReadError::UnknownMagic,
            internal::FullDataCirTreeError::IoError(e) => BBIReadError::IoError(e),
        }
    }
}

/// Potential errors found when trying to read data from a zoom level
#[derive(Error, Debug)]
pub enum ZoomIntervalError {
    #[error("The passed reduction level was not found")]
    ReductionLevelNotFound,
    #[error("{}", .0)]
    BBIReadError(BBIReadError),
}

impl From<ChromIdNotFound> for ZoomIntervalError {
    fn from(e: ChromIdNotFound) -> Self {
        ZoomIntervalError::BBIReadError(BBIReadError::InvalidChromosome(e.0))
    }
}

impl From<CirTreeSearchError> for ZoomIntervalError {
    fn from(e: CirTreeSearchError) -> Self {
        ZoomIntervalError::BBIReadError(e.into())
    }
}

impl From<internal::ZoomDataCirTreeError> for ZoomIntervalError {
    fn from(value: internal::ZoomDataCirTreeError) -> Self {
        match value {
            internal::ZoomDataCirTreeError::UnknownMagic => {
                ZoomIntervalError::BBIReadError(BBIReadError::UnknownMagic)
            }
            internal::ZoomDataCirTreeError::ReductionLevelNotFound => {
                ZoomIntervalError::ReductionLevelNotFound
            }
            internal::ZoomDataCirTreeError::IoError(e) => {
                ZoomIntervalError::BBIReadError(BBIReadError::IoError(e))
            }
        }
    }
}

/// What kind of cir tree in a bbi file
pub enum CirTreeIndexType {
    /// The index for the full data
    FullData,
    /// The index for zoom data of a given reduction level
    Zoom(u32),
}

/// Represents a cir tree index in a bbi file. Composed of a public
/// `CirTreeIndexType`, and a private location in the bbi file.
/// This can be passed to `search_cir_tree`.
pub struct CirTreeIndex(
    pub CirTreeIndexType,
    pub(crate) u64,
    pub(crate) u64,
    pub(crate) Option<u64>,
);

pub(crate) mod internal {
    use super::*;

    pub enum FullDataCirTreeError {
        UnknownMagic,
        IoError(io::Error),
    }

    pub enum ZoomDataCirTreeError {
        UnknownMagic,
        ReductionLevelNotFound,
        IoError(io::Error),
    }

    pub trait BBIReadInternal {
        type Read: BBIFileRead;

        /// Gets a reader to the underlying file
        fn reader(&mut self) -> &mut Self::Read;

        fn reader_and_info(&mut self) -> (&mut Self::Read, &mut BBIFileInfo);

        fn full_data_cir_tree(&mut self) -> Result<CirTreeIndex, FullDataCirTreeError> {
            let (reader, info) = self.reader_and_info();
            let index_offset = info.header.full_index_offset;
            if info.header.full_index_tree_offset.is_none()
                || info.header.full_index_item_count.is_none()
            {
                let endianness = info.header.endianness;

                reader
                    .raw_reader()
                    .seek(SeekFrom::Start(index_offset))
                    .map_err(|e| FullDataCirTreeError::IoError(e))?;

                let item_count =
                    read_cir_tree_header(endianness, reader.raw_reader()).map_err(|e| match e {
                        Either::Left(_) => FullDataCirTreeError::UnknownMagic,
                        Either::Right(e) => FullDataCirTreeError::IoError(e),
                    })?;

                info.header.full_index_tree_offset = Some(index_offset + 48);
                info.header.full_index_item_count = Some(item_count);
            }
            let root_offset = index_offset + 48;
            let index_end = info
                .zoom_headers
                .iter()
                .map(|header| header.data_offset)
                .filter(|offset| *offset > root_offset)
                .min();
            Ok(CirTreeIndex(
                CirTreeIndexType::FullData,
                root_offset,
                info.header.full_index_item_count.unwrap(),
                index_end,
            ))
        }

        fn zoom_cir_tree(
            &mut self,
            reduction_level: u32,
        ) -> Result<CirTreeIndex, ZoomDataCirTreeError> {
            let (reader, info) = self.reader_and_info();
            let zoom_header = match info
                .zoom_headers
                .iter_mut()
                .find(|h| h.reduction_level == reduction_level)
            {
                Some(h) => h,
                None => {
                    return Err(ZoomDataCirTreeError::ReductionLevelNotFound);
                }
            };

            if zoom_header.index_tree_offset.is_none() || zoom_header.index_item_count.is_none() {
                let endianness = info.header.endianness;

                reader
                    .raw_reader()
                    .seek(SeekFrom::Start(zoom_header.index_offset))
                    .map_err(|e| ZoomDataCirTreeError::IoError(e))?;

                let item_count =
                    read_cir_tree_header(endianness, reader.raw_reader()).map_err(|e| match e {
                        Either::Left(_) => ZoomDataCirTreeError::UnknownMagic,
                        Either::Right(e) => ZoomDataCirTreeError::IoError(e),
                    })?;

                zoom_header.index_tree_offset = Some(zoom_header.index_offset + 48);
                zoom_header.index_item_count = Some(item_count);
            }

            Ok(CirTreeIndex(
                CirTreeIndexType::Zoom(reduction_level),
                zoom_header.index_offset + 48,
                zoom_header.index_item_count.unwrap(),
                None,
            ))
        }
    }
}

pub struct ReductionLevelNotFound;

/// Generic methods for reading a bbi file
pub trait BBIRead: BBIReadInternal {
    /// Get basic info about the bbi file
    fn info(&self) -> &BBIFileInfo;

    fn chroms(&self) -> &[ChromInfo];
}

pub(crate) fn search_cir_tree<R: BBIFileRead>(
    info: &BBIFileInfo,
    file: &mut R,
    at: CirTreeIndex,
    chrom_name: &str,
    start: u32,
    end: u32,
) -> Result<Vec<Block>, CirTreeSearchError> {
    let chrom_ix = {
        let chrom = info.chrom_info.iter().find(|&x| x.name == chrom_name);
        match chrom {
            Some(c) => c.id,
            None => {
                return Err(CirTreeSearchError::InvalidChromosome(
                    chrom_name.to_string(),
                ));
            }
        }
    };

    let endianness = info.header.endianness;

    Ok(search_cir_tree_inner(
        endianness, file, at.1, chrom_ix, start, end,
    )?)
}

#[derive(Debug)]
pub(crate) struct UnknownMagic;

pub(crate) fn read_cir_tree_header<R: Read + Seek>(
    endianness: Endianness,
    file: &mut R,
) -> Result<u64, Either<UnknownMagic, io::Error>> {
    let mut header_data = BytesMut::zeroed(48);
    file.read_exact(&mut header_data)
        .map_err(|e| Either::Right(e))?;

    let item_count = match endianness {
        Endianness::Big => {
            let magic = header_data.get_u32();
            if magic != CIR_TREE_MAGIC {
                return Err(Either::Left(UnknownMagic));
            }

            let _blocksize = header_data.get_u32();
            let item_count = header_data.get_u64();
            let _start_chrom_idx = header_data.get_u32();
            let _start_base = header_data.get_u32();
            let _end_chrom_idx = header_data.get_u32();
            let _end_base = header_data.get_u32();
            let _end_file_offset = header_data.get_u64();
            let _item_per_slot = header_data.get_u32();
            let _reserved = header_data.get_u32();
            item_count
        }
        Endianness::Little => {
            let magic = header_data.get_u32_le();
            if magic != CIR_TREE_MAGIC {
                return Err(Either::Left(UnknownMagic));
            }

            let _blocksize = header_data.get_u32_le();
            let item_count = header_data.get_u64_le();
            let _start_chrom_idx = header_data.get_u32_le();
            let _start_base = header_data.get_u32_le();
            let _end_chrom_idx = header_data.get_u32_le();
            let _end_base = header_data.get_u32_le();
            let _end_file_offset = header_data.get_u64_le();
            let _item_per_slot = header_data.get_u32_le();
            let _reserved = header_data.get_u32_le();
            item_count
        }
    };
    Ok(item_count)
}

pub(crate) fn search_cir_tree_inner<R: BBIFileRead>(
    endianness: Endianness,
    file: &mut R,
    at: u64,
    chrom_ix: u32,
    start: u32,
    end: u32,
) -> io::Result<Vec<Block>> {
    // We currently don't check that the passed interval overlaps with *any* data.
    // We could, but would have to store this data when we check the header.
    let mut blocks = vec![];

    let mut remaining_childblocks = VecDeque::with_capacity(2048);
    remaining_childblocks.push_front(at);
    let iter = CirTreeBlockSearchIter {
        remaining_childblocks,
        file,
        endianness,
        chrom_ix,
        start,
        end,
    };

    for i in iter {
        let i = i?;
        blocks.extend(i);
    }

    Ok(blocks)
}

/// Traverse a cir-tree and return every primary data-block leaf in index order
/// (coordinate order for a well-formed cir-tree) without reading or
/// decompressing block contents.
pub(crate) fn cir_tree_data_blocks<R: BBIFileRead>(
    endianness: Endianness,
    file: &mut R,
    at: CirTreeIndex,
    limits: BBIDataBlockLimits,
    primary_data_start: u64,
    primary_data_end: u64,
) -> Result<Vec<BBIDataBlock>, BBIReadError> {
    // The header counts data records/sections, not cir-tree leaf entries. It is
    // therefore an upper bound for blocks, not a count that the layout must
    // equal. It does not bound internal nodes, whose count also depends on the
    // cir-tree fanout, so nodes use only the independent safety limit and the
    // structural index-span bound below.
    let item_count = at.2;
    let mut structural_max_nodes = u64::MAX;
    let mut structural_max_blocks = item_count;
    if let Some(index_end) = at.3 {
        let index_span = index_end.checked_sub(at.1).ok_or_else(|| {
            BBIReadError::InvalidFile("cir-tree index end precedes its root".into())
        })?;
        structural_max_nodes = index_span / 4;
        structural_max_blocks = structural_max_blocks.min(index_span / 32);
    }
    let max_nodes = limits.max_nodes.min(structural_max_nodes);
    let max_blocks = limits.max_blocks.min(structural_max_blocks);
    if item_count == 0 {
        return Ok(Vec::new());
    }

    let mut blocks = Vec::new();
    let mut remaining_nodes = VecDeque::with_capacity(2048);
    let mut visited_nodes = HashSet::new();
    let mut children = Vec::with_capacity(256);
    remaining_nodes.push_front(at.1);
    visited_nodes.insert(at.1);
    if max_nodes == 0 {
        let message = if structural_max_nodes == 0 {
            "cir-tree node count exceeds structural index bound"
        } else {
            "cir-tree node limit exceeded"
        };
        return Err(BBIReadError::InvalidFile(message.into()));
    }

    while let Some(node_offset) = remaining_nodes.pop_front() {
        children.clear();
        let first_new_block = blocks.len();
        file.data_blocks_for_cir_tree_node(
            endianness,
            node_offset,
            at.3,
            &mut blocks,
            &mut children,
        )
        .map_err(|error| {
            if matches!(
                error.kind(),
                io::ErrorKind::InvalidData | io::ErrorKind::UnexpectedEof
            ) {
                BBIReadError::InvalidFile(error.to_string())
            } else {
                BBIReadError::IoError(error)
            }
        })?;
        if blocks.len() == first_new_block && children.is_empty() {
            return Err(BBIReadError::InvalidFile(
                "cir-tree contains an empty node for a nonempty index".into(),
            ));
        }
        for block in &blocks[first_new_block..] {
            let data_end = block.offset.checked_add(block.data_size).ok_or_else(|| {
                BBIReadError::InvalidFile("cir-tree data block range overflows".into())
            })?;
            if block.data_size == 0
                || block.offset < primary_data_start
                || data_end > primary_data_end
            {
                return Err(BBIReadError::InvalidFile(format!(
                    "cir-tree data block [{}, {}) lies outside primary data [{}, {})",
                    block.offset, data_end, primary_data_start, primary_data_end
                )));
            }
        }
        let block_count = u64::try_from(blocks.len()).unwrap_or(u64::MAX);
        if block_count > max_blocks {
            let message = if block_count > structural_max_blocks {
                "cir-tree data block count exceeds structural bound"
            } else {
                "cir-tree data block limit exceeded"
            };
            return Err(BBIReadError::InvalidFile(message.into()));
        }
        // Validate every sibling before enforcing the configurable node limit so
        // malformed duplicate offsets cannot be misclassified as caller limits.
        for &child in children.iter().rev() {
            let child_header_in_range = child
                .checked_add(4)
                .is_some_and(|child_end| at.3.map_or(true, |index_end| child_end <= index_end));
            if child < at.1 || !child_header_in_range {
                return Err(BBIReadError::InvalidFile(
                    "cir-tree child node offset out of range".into(),
                ));
            }
            if !visited_nodes.insert(child) {
                return Err(BBIReadError::InvalidFile(
                    "cir-tree contains a repeated child node offset".into(),
                ));
            }
        }
        let node_count = u64::try_from(visited_nodes.len()).unwrap_or(u64::MAX);
        if node_count > max_nodes {
            let message = if node_count > structural_max_nodes {
                "cir-tree node count exceeds structural index bound"
            } else {
                "cir-tree node limit exceeded"
            };
            return Err(BBIReadError::InvalidFile(message.into()));
        }
        for child in children.drain(..).rev() {
            remaining_nodes.push_front(child);
        }
    }
    Ok(blocks)
}

fn append_node_data_blocks<L, N>(
    iter: CirTreeNodeIterator<L, N>,
    blocks: &mut Vec<BBIDataBlock>,
    children: &mut Vec<u64>,
) where
    L: Iterator<Item = CirTreeNodeLeaf>,
    N: Iterator<Item = CirTreeNodeNonLeaf>,
{
    match iter {
        CirTreeNodeIterator::Leaf(items) => {
            blocks.extend(items.map(|leaf| data_block_from_leaf(&leaf)));
        }
        CirTreeNodeIterator::NonLeaf(items) => {
            children.extend(items.map(|child| child.node_offset));
        }
    }
}

fn append_cached_node_data_blocks(
    node: &Either<Vec<CirTreeNodeLeaf>, Vec<CirTreeNodeNonLeaf>>,
    blocks: &mut Vec<BBIDataBlock>,
    children: &mut Vec<u64>,
) {
    match node {
        Either::Left(items) => {
            blocks.extend(items.iter().map(data_block_from_leaf));
        }
        Either::Right(items) => {
            children.extend(items.iter().map(|child| child.node_offset));
        }
    }
}

fn data_block_from_leaf(leaf: &CirTreeNodeLeaf) -> BBIDataBlock {
    BBIDataBlock {
        start_chrom_id: leaf.start_chrom_ix,
        start_base: leaf.start_base,
        end_chrom_id: leaf.end_chrom_ix,
        end_base: leaf.end_base,
        offset: leaf.data_offset,
        data_size: leaf.data_size,
    }
}

pub enum GenericBBIRead<R> {
    BigWig(BigWigRead<R>),
    BigBed(BigBedRead<R>),
}

impl<R: SeekableRead> BBIRead for GenericBBIRead<R> {
    fn info(&self) -> &BBIFileInfo {
        match self {
            GenericBBIRead::BigWig(b) => b.info(),
            GenericBBIRead::BigBed(b) => b.info(),
        }
    }

    fn chroms(&self) -> &[ChromInfo] {
        match self {
            GenericBBIRead::BigWig(b) => b.chroms(),
            GenericBBIRead::BigBed(b) => b.chroms(),
        }
    }
}

impl<R: SeekableRead> BBIReadInternal for GenericBBIRead<R> {
    type Read = R;

    fn reader(&mut self) -> &mut Self::Read {
        match self {
            GenericBBIRead::BigWig(b) => b.reader(),
            GenericBBIRead::BigBed(b) => b.reader(),
        }
    }

    fn reader_and_info(&mut self) -> (&mut Self::Read, &mut BBIFileInfo) {
        match self {
            GenericBBIRead::BigWig(b) => b.reader_and_info(),
            GenericBBIRead::BigBed(b) => b.reader_and_info(),
        }
    }
}
/// Possible errors encountered when opening a bigBed file to read
#[derive(Error, Debug)]
pub enum GenericBBIFileOpenError {
    #[error("File is not a bigWig or bigBed.")]
    NotABBIFile,
    #[error("The chromosomes are invalid.")]
    InvalidChroms,
    #[error("An error occurred: {}", .0)]
    IoError(#[from] io::Error),
}

impl From<BBIFileReadInfoError> for GenericBBIFileOpenError {
    fn from(error: BBIFileReadInfoError) -> Self {
        match error {
            BBIFileReadInfoError::UnknownMagic => GenericBBIFileOpenError::NotABBIFile,
            BBIFileReadInfoError::InvalidChroms => GenericBBIFileOpenError::InvalidChroms,
            BBIFileReadInfoError::IoError(e) => GenericBBIFileOpenError::IoError(e),
        }
    }
}

impl<R> GenericBBIRead<R> {
    pub fn bigwig(self) -> Option<BigWigRead<R>> {
        match self {
            GenericBBIRead::BigWig(b) => Some(b),
            GenericBBIRead::BigBed(_) => None,
        }
    }

    pub fn bigbed(self) -> Option<BigBedRead<R>> {
        match self {
            GenericBBIRead::BigBed(b) => Some(b),
            GenericBBIRead::BigWig(_) => None,
        }
    }
}

impl<R: BBIFileRead> GenericBBIRead<R> {
    /// Opens a generic bbi file for a given type that implements both `Read` and `Seek`
    pub fn open(mut read: R) -> Result<Self, GenericBBIFileOpenError> {
        let info = read_info(&mut read)?;
        match info.filetype {
            BBIFile::BigWig => Ok(GenericBBIRead::BigWig(BigWigRead { info, read })),
            BBIFile::BigBed => Ok(GenericBBIRead::BigBed(BigBedRead { info, read })),
        }
    }
}

impl GenericBBIRead<ReopenableFile> {
    /// Opens a generic bbi file
    pub fn open_file(path: &str) -> Result<Self, GenericBBIFileOpenError> {
        let reopen = ReopenableFile {
            file: File::open(&path)?,
            path: path.into(),
        };
        let b = GenericBBIRead::open(reopen);
        if b.is_err() {
            eprintln!("Error when opening: {}", path);
        }
        b
    }
}

pub trait BBIFileRead {
    type Reader: Read + Seek;

    fn get_block_data(&mut self, info: &BBIFileInfo, block: &Block) -> io::Result<Vec<u8>>;

    fn blocks_for_cir_tree_node(
        &mut self,
        endianness: Endianness,
        node_offset: u64,
        chrom_ix: u32,
        start: u32,
        end: u32,
    ) -> io::Result<(SmallVec<[u64; 4]>, SmallVec<[Block; 4]>)>;

    /// Append a cir-tree node's primary data blocks and child node offsets.
    ///
    /// Full-index traversals use this method. Implementations that cache or
    /// prefetch parsed cir-tree nodes can override the default reader behavior.
    fn data_blocks_for_cir_tree_node(
        &mut self,
        endianness: Endianness,
        node_offset: u64,
        index_end: Option<u64>,
        blocks: &mut Vec<BBIDataBlock>,
        children: &mut Vec<u64>,
    ) -> io::Result<()> {
        let iter = read_node_bounded(self.raw_reader(), node_offset, endianness, index_end)?;
        append_node_data_blocks(iter, blocks, children);
        Ok(())
    }

    fn raw_reader(&mut self) -> &mut Self::Reader;
}

impl<S: SeekableRead> BBIFileRead for S {
    type Reader = Self;

    fn get_block_data(&mut self, info: &BBIFileInfo, block: &Block) -> io::Result<Vec<u8>> {
        read_block_data(info, self, block)
    }

    fn blocks_for_cir_tree_node(
        &mut self,
        endianness: Endianness,
        node_offset: u64,
        chrom_ix: u32,
        start: u32,
        end: u32,
    ) -> io::Result<(SmallVec<[u64; 4]>, SmallVec<[Block; 4]>)> {
        let iter = match read_node(self, node_offset, endianness) {
            Ok(d) => d,
            Err(e) => return Err(e),
        };

        Ok(nodes_overlapping(iter, chrom_ix, start, end))
    }

    fn raw_reader(&mut self) -> &mut Self::Reader {
        self
    }
}

pub struct CachedBBIFileRead<S> {
    read: S,
    cir_tree_node_map: HashMap<u64, Either<Vec<CirTreeNodeLeaf>, Vec<CirTreeNodeNonLeaf>>>,
    block_data: HashMap<Block, Vec<u8>>,
}

impl<S> CachedBBIFileRead<S> {
    pub fn inner_read(&self) -> &S {
        &self.read
    }
}

impl<S: SeekableRead> CachedBBIFileRead<S> {
    pub fn new(read: S) -> Self {
        CachedBBIFileRead {
            read,
            cir_tree_node_map: HashMap::new(),
            block_data: HashMap::new(),
        }
    }
}

impl<S: SeekableRead> BBIFileRead for CachedBBIFileRead<S> {
    type Reader = S;

    fn get_block_data(&mut self, info: &BBIFileInfo, block: &Block) -> io::Result<Vec<u8>> {
        if let Some(data) = self.block_data.get(block) {
            return Ok(data.clone());
        }
        if self.block_data.len() >= 5000 {
            self.block_data.clear();
        }
        let data = read_block_data(info, &mut self.read, block)?;
        self.block_data.insert(*block, data.clone());
        Ok(data)
    }

    fn blocks_for_cir_tree_node(
        &mut self,
        endianness: Endianness,
        node_offset: u64,
        chrom_ix: u32,
        start: u32,
        end: u32,
    ) -> io::Result<(SmallVec<[u64; 4]>, SmallVec<[Block; 4]>)> {
        match self.cir_tree_node_map.entry(node_offset) {
            Entry::Occupied(node) => {
                let iter = match node.get() {
                    Either::Left(v) => CirTreeNodeIterator::Leaf(v.clone().into_iter()),
                    Either::Right(v) => CirTreeNodeIterator::NonLeaf(v.clone().into_iter()),
                };
                Ok(nodes_overlapping(iter, chrom_ix, start, end))
            }
            Entry::Vacant(e) => {
                let iter = match read_node(&mut self.read, node_offset, endianness) {
                    Ok(d) => d,
                    Err(e) => return Err(e),
                };
                let iter = match iter {
                    CirTreeNodeIterator::Leaf(v) => {
                        let v: Vec<_> = v.collect();
                        e.insert(Either::Left(v.clone()));
                        CirTreeNodeIterator::Leaf(v.into_iter())
                    }
                    CirTreeNodeIterator::NonLeaf(v) => {
                        let v: Vec<_> = v.collect();
                        e.insert(Either::Right(v.clone()));
                        CirTreeNodeIterator::NonLeaf(v.into_iter())
                    }
                };

                Ok(nodes_overlapping(iter, chrom_ix, start, end))
            }
        }
    }

    fn data_blocks_for_cir_tree_node(
        &mut self,
        endianness: Endianness,
        node_offset: u64,
        index_end: Option<u64>,
        blocks: &mut Vec<BBIDataBlock>,
        children: &mut Vec<u64>,
    ) -> io::Result<()> {
        if let Some(node) = self.cir_tree_node_map.get(&node_offset) {
            let (count, record_size) = match node {
                Either::Left(items) => (items.len(), 32),
                Either::Right(items) => (items.len(), 24),
            };
            validate_cir_tree_node_span(node_offset, count, record_size, index_end)?;
            append_cached_node_data_blocks(node, blocks, children);
            return Ok(());
        }

        let iter = read_node_bounded(&mut self.read, node_offset, endianness, index_end)?;
        append_node_data_blocks(iter, blocks, children);
        Ok(())
    }

    fn raw_reader(&mut self) -> &mut Self::Reader {
        &mut self.read
    }
}

impl<R: Reopen + SeekableRead> Reopen for CachedBBIFileRead<R> {
    fn reopen(&self) -> io::Result<Self> {
        Ok(Self {
            read: self.read.reopen()?,
            cir_tree_node_map: self.cir_tree_node_map.clone(),
            block_data: self.block_data.clone(),
        })
    }
}

#[cfg(test)]
mod data_block_tests {
    use std::io::Cursor;
    use std::path::PathBuf;

    use super::*;

    fn little_endian_leaf_node() -> Vec<u8> {
        let mut bytes = vec![1, 0, 1, 0];
        bytes.extend_from_slice(&2_u32.to_le_bytes());
        bytes.extend_from_slice(&10_u32.to_le_bytes());
        bytes.extend_from_slice(&2_u32.to_le_bytes());
        bytes.extend_from_slice(&20_u32.to_le_bytes());
        bytes.extend_from_slice(&128_u64.to_le_bytes());
        bytes.extend_from_slice(&64_u64.to_le_bytes());
        bytes
    }

    fn little_endian_cir_tree_header(item_count: u64) -> Vec<u8> {
        let mut bytes = vec![0; 48];
        bytes[0..4].copy_from_slice(&CIR_TREE_MAGIC.to_le_bytes());
        bytes[8..16].copy_from_slice(&item_count.to_le_bytes());
        bytes
    }

    fn little_endian_non_leaf_node(child_offset: u64) -> Vec<u8> {
        let mut bytes = little_endian_non_leaf_node_children(&[child_offset]);
        bytes.extend_from_slice(&[0; 8]);
        bytes
    }

    fn little_endian_non_leaf_node_children(child_offsets: &[u64]) -> Vec<u8> {
        let mut bytes = vec![0, 0];
        bytes.extend_from_slice(&(child_offsets.len() as u16).to_le_bytes());
        for &child_offset in child_offsets {
            bytes.extend_from_slice(&0_u32.to_le_bytes());
            bytes.extend_from_slice(&0_u32.to_le_bytes());
            bytes.extend_from_slice(&0_u32.to_le_bytes());
            bytes.extend_from_slice(&1_u32.to_le_bytes());
            bytes.extend_from_slice(&child_offset.to_le_bytes());
        }
        bytes
    }

    fn invalid_file_message(error: BBIReadError) -> String {
        match error {
            BBIReadError::InvalidFile(message) => message,
            other => panic!("expected invalid-file error, got {other}"),
        }
    }

    #[test]
    fn full_traversal_does_not_pollute_query_cache() {
        let mut reader = CachedBBIFileRead::new(Cursor::new(little_endian_leaf_node()));
        reader
            .cir_tree_node_map
            .insert(64, Either::Left(Vec::new()));
        let mut blocks = Vec::new();
        let mut children = Vec::new();

        reader
            .data_blocks_for_cir_tree_node(Endianness::Little, 0, None, &mut blocks, &mut children)
            .unwrap();

        assert!(children.is_empty());
        assert_eq!(
            blocks,
            vec![BBIDataBlock {
                start_chrom_id: 2,
                start_base: 10,
                end_chrom_id: 2,
                end_base: 20,
                offset: 128,
                data_size: 64,
            }]
        );
        assert_eq!(reader.cir_tree_node_map.len(), 1);
        assert!(reader.cir_tree_node_map.contains_key(&64));
        assert!(!reader.cir_tree_node_map.contains_key(&0));
    }

    #[test]
    fn full_traversal_reuses_complete_query_cache_nodes() {
        let mut directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        directory.push("resources/test");

        let bigwig_path = directory.join("valid.bigWig");
        let mut uncached_bigwig = BigWigRead::open_file(&bigwig_path).unwrap();
        let expected_bigwig = uncached_bigwig.data_blocks().unwrap();
        let mut cached_bigwig = BigWigRead::open_file(bigwig_path).unwrap().cached();
        let bigwig_chromosomes: Vec<_> = cached_bigwig
            .chroms()
            .iter()
            .map(|chrom| (chrom.name.clone(), chrom.length))
            .collect();
        for (name, length) in bigwig_chromosomes {
            drop(cached_bigwig.get_interval(&name, 0, length).unwrap());
        }
        assert!(!cached_bigwig.read.cir_tree_node_map.is_empty());
        assert_eq!(cached_bigwig.data_blocks().unwrap(), expected_bigwig);

        let bigbed_path = directory.join("bigGenePred.bb");
        let mut uncached_bigbed = BigBedRead::open_file(&bigbed_path).unwrap();
        let expected_bigbed = uncached_bigbed.data_blocks().unwrap();

        let mut narrow_cached_bigbed = BigBedRead::open_file(&bigbed_path).unwrap().cached();
        let first_chromosome = narrow_cached_bigbed.chroms()[0].name.clone();
        drop(
            narrow_cached_bigbed
                .get_interval(&first_chromosome, 1_000_000, 1_001_000)
                .unwrap(),
        );
        assert!(narrow_cached_bigbed.read.cir_tree_node_map.len() >= 2);
        assert_eq!(narrow_cached_bigbed.data_blocks().unwrap(), expected_bigbed);

        let mut cached_bigbed = BigBedRead::open_file(bigbed_path).unwrap().cached();
        let bigbed_chromosomes: Vec<_> = cached_bigbed
            .chroms()
            .iter()
            .map(|chrom| (chrom.name.clone(), chrom.length))
            .collect();
        for (name, length) in bigbed_chromosomes {
            drop(cached_bigbed.get_interval(&name, 0, length).unwrap());
        }
        assert!(cached_bigbed.read.cir_tree_node_map.len() >= 2);
        assert_eq!(cached_bigbed.data_blocks().unwrap(), expected_bigbed);
    }

    #[test]
    fn full_traversal_rejects_cyclic_index() {
        let mut bytes = little_endian_cir_tree_header(2);
        bytes.extend_from_slice(&little_endian_non_leaf_node(84));
        bytes.extend_from_slice(&little_endian_non_leaf_node(48));
        let mut reader = Cursor::new(bytes);

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 2, Some(120)),
            BBIDataBlockLimits::default(),
            128,
            192,
        )
        .unwrap_err();

        assert_eq!(
            invalid_file_message(error),
            "cir-tree contains a repeated child node offset"
        );
    }

    #[test]
    fn full_traversal_detects_duplicate_siblings_before_caller_limit() {
        let mut bytes = little_endian_cir_tree_header(1);
        bytes.extend_from_slice(&little_endian_non_leaf_node_children(&[100, 100]));
        bytes.extend_from_slice(&[1, 0, 0, 0]);
        let mut reader = Cursor::new(bytes);

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 1, Some(104)),
            BBIDataBlockLimits {
                max_nodes: 1,
                ..BBIDataBlockLimits::default()
            },
            128,
            192,
        )
        .unwrap_err();

        assert!(!error.is_data_block_traversal_limit_exceeded());
        assert_eq!(
            invalid_file_message(error),
            "cir-tree contains a repeated child node offset"
        );
    }

    #[test]
    fn full_traversal_rejects_empty_nodes_in_nonempty_index() {
        for is_leaf in [false, true] {
            let mut bytes = little_endian_cir_tree_header(1);
            bytes.extend_from_slice(&[u8::from(is_leaf), 0, 0, 0]);
            let mut reader = Cursor::new(bytes);

            let error = cir_tree_data_blocks(
                Endianness::Little,
                &mut reader,
                CirTreeIndex(CirTreeIndexType::FullData, 48, 1, Some(52)),
                BBIDataBlockLimits::default(),
                128,
                192,
            )
            .unwrap_err();

            assert!(!error.is_data_block_traversal_limit_exceeded());
            assert_eq!(
                invalid_file_message(error),
                "cir-tree contains an empty node for a nonempty index"
            );
        }
    }

    #[test]
    fn full_traversal_enforces_declared_item_count() {
        let mut bytes = little_endian_cir_tree_header(1);
        let leaf = little_endian_leaf_node();
        bytes.extend_from_slice(&[1, 0, 2, 0]);
        bytes.extend_from_slice(&leaf[4..]);
        bytes.extend_from_slice(&leaf[4..]);
        let mut reader = Cursor::new(bytes);

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 1, Some(116)),
            BBIDataBlockLimits::default(),
            128,
            192,
        )
        .unwrap_err();

        assert!(!error.is_data_block_traversal_limit_exceeded());
        assert_eq!(
            invalid_file_message(error),
            "cir-tree data block count exceeds structural bound"
        );
    }

    #[test]
    fn full_traversal_enforces_node_limit() {
        let mut bytes = little_endian_cir_tree_header(1);
        bytes.extend_from_slice(&little_endian_non_leaf_node(84));
        bytes.extend_from_slice(&little_endian_non_leaf_node(120));
        bytes.extend_from_slice(&little_endian_non_leaf_node(156));
        let mut reader = Cursor::new(bytes);

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 1, None),
            BBIDataBlockLimits {
                max_nodes: 2,
                ..BBIDataBlockLimits::default()
            },
            128,
            192,
        )
        .unwrap_err();

        assert_eq!(invalid_file_message(error), "cir-tree node limit exceeded");
    }

    #[test]
    fn full_traversal_does_not_classify_structural_node_bound_as_limit() {
        let mut reader = Cursor::new(little_endian_cir_tree_header(1));

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 1, Some(51)),
            BBIDataBlockLimits::default(),
            128,
            192,
        )
        .unwrap_err();

        assert!(!error.is_data_block_traversal_limit_exceeded());
        assert_eq!(
            invalid_file_message(error),
            "cir-tree node count exceeds structural index bound"
        );
    }

    #[test]
    fn full_traversal_does_not_derive_node_limit_from_item_count() {
        let mut bytes = little_endian_cir_tree_header(1);
        bytes.extend_from_slice(&little_endian_non_leaf_node(84));
        bytes.extend_from_slice(&little_endian_non_leaf_node(120));
        bytes.extend_from_slice(&little_endian_leaf_node());
        let mut reader = Cursor::new(bytes);

        let blocks = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 1, Some(156)),
            BBIDataBlockLimits::default(),
            128,
            192,
        )
        .unwrap();

        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].offset, 128);
    }

    #[test]
    fn full_traversal_rejects_child_outside_index() {
        let mut bytes = little_endian_cir_tree_header(2);
        bytes.extend_from_slice(&little_endian_non_leaf_node(120));
        let mut reader = Cursor::new(bytes);

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 2, Some(84)),
            BBIDataBlockLimits::default(),
            128,
            192,
        )
        .unwrap_err();

        assert_eq!(
            invalid_file_message(error),
            "cir-tree child node offset out of range"
        );
    }

    #[test]
    fn full_traversal_rejects_leaf_records_outside_index() {
        let mut bytes = little_endian_cir_tree_header(1);
        bytes.extend_from_slice(&little_endian_leaf_node());
        let mut reader = Cursor::new(bytes);

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 1, Some(83)),
            BBIDataBlockLimits::default(),
            128,
            192,
        )
        .unwrap_err();

        assert_eq!(
            invalid_file_message(error),
            "cir-tree node records extend outside index"
        );
    }

    #[test]
    fn full_traversal_validates_cached_node_span() {
        let mut reader = CachedBBIFileRead::new(Cursor::new(Vec::new()));
        reader.cir_tree_node_map.insert(
            48,
            Either::Left(vec![CirTreeNodeLeaf {
                start_chrom_ix: 2,
                start_base: 10,
                end_chrom_ix: 2,
                end_base: 20,
                data_offset: 128,
                data_size: 64,
            }]),
        );

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 1, Some(83)),
            BBIDataBlockLimits::default(),
            128,
            192,
        )
        .unwrap_err();

        assert_eq!(
            invalid_file_message(error),
            "cir-tree node records extend outside index"
        );
    }

    #[test]
    fn full_traversal_returns_empty_layout_for_empty_index() {
        let mut reader = Cursor::new(little_endian_cir_tree_header(0));

        let blocks = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 0, Some(48)),
            BBIDataBlockLimits::default(),
            8,
            48,
        )
        .unwrap();

        assert!(blocks.is_empty());
    }

    #[test]
    fn full_traversal_rejects_inverted_empty_index_span() {
        let mut reader = Cursor::new(little_endian_cir_tree_header(0));

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 0, Some(47)),
            BBIDataBlockLimits::default(),
            8,
            48,
        )
        .unwrap_err();

        assert_eq!(
            invalid_file_message(error),
            "cir-tree index end precedes its root"
        );
    }

    #[test]
    fn full_traversal_classifies_truncated_node_as_invalid_file() {
        let mut reader = Cursor::new(little_endian_cir_tree_header(1));

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 1, Some(52)),
            BBIDataBlockLimits::default(),
            128,
            192,
        )
        .unwrap_err();

        assert!(invalid_file_message(error).contains("failed to fill whole buffer"));
    }

    #[test]
    fn full_traversal_rejects_invalid_node_tag_without_panicking() {
        let mut bytes = little_endian_cir_tree_header(1);
        bytes.extend_from_slice(&[2, 0, 0, 0]);
        let mut reader = Cursor::new(bytes);

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 1, Some(52)),
            BBIDataBlockLimits::default(),
            128,
            192,
        )
        .unwrap_err();

        assert_eq!(
            invalid_file_message(error),
            "cir-tree node has invalid leaf flag 2"
        );
    }

    #[test]
    fn full_traversal_rejects_data_block_outside_primary_data() {
        let mut bytes = little_endian_cir_tree_header(1);
        bytes.extend_from_slice(&little_endian_leaf_node());
        let mut reader = Cursor::new(bytes);

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 1, Some(84)),
            BBIDataBlockLimits::default(),
            192,
            256,
        )
        .unwrap_err();

        assert_eq!(
            invalid_file_message(error),
            "cir-tree data block [128, 192) lies outside primary data [192, 256)"
        );
    }

    #[test]
    fn full_traversal_rejects_overflowing_data_block_range() {
        let mut leaf = little_endian_leaf_node();
        leaf[20..28].copy_from_slice(&u64::MAX.to_le_bytes());
        leaf[28..36].copy_from_slice(&2_u64.to_le_bytes());
        let mut bytes = little_endian_cir_tree_header(1);
        bytes.extend_from_slice(&leaf);
        let mut reader = Cursor::new(bytes);

        let error = cir_tree_data_blocks(
            Endianness::Little,
            &mut reader,
            CirTreeIndex(CirTreeIndexType::FullData, 48, 1, Some(84)),
            BBIDataBlockLimits::default(),
            128,
            256,
        )
        .unwrap_err();

        assert_eq!(
            invalid_file_message(error),
            "cir-tree data block range overflows"
        );
    }

    #[test]
    fn layout_matches_blocks_used_by_full_range_queries() {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("resources/test/valid.bigWig");
        let mut reader = BigWigRead::open_file(path).unwrap();
        let layout: HashSet<_> = reader
            .data_blocks()
            .unwrap()
            .into_iter()
            .map(|block| (block.offset, block.data_size))
            .collect();
        let chromosomes: Vec<_> = reader
            .chroms()
            .iter()
            .map(|chrom| (chrom.name.clone(), chrom.length))
            .collect();
        let mut query_blocks = HashSet::new();

        for (name, length) in chromosomes {
            let tree = match reader.full_data_cir_tree() {
                Ok(tree) => tree,
                Err(_) => panic!("fixture should contain a valid cir-tree"),
            };
            let blocks =
                search_cir_tree(&reader.info, &mut reader.read, tree, &name, 0, length).unwrap();
            query_blocks.extend(blocks.into_iter().map(|block| (block.offset, block.size)));
        }

        assert_eq!(layout, query_blocks);
    }

    #[test]
    fn primary_data_bounds_follow_supported_fixture_layouts() {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("resources/test");

        let bigwig = BigWigRead::open_file(path.join("valid.bigWig")).unwrap();
        assert_eq!(
            bigwig
                .info
                .header
                .primary_data_bounds(BBIFile::BigWig)
                .unwrap(),
            (348, 603_600)
        );

        let bigbed = BigBedRead::open_file(path.join("bigGenePred.bb")).unwrap();
        assert_eq!(
            bigbed
                .info
                .header
                .primary_data_bounds(BBIFile::BigBed)
                .unwrap(),
            (25_592, 6_252_019)
        );
    }
}

pub(crate) fn read_info<R: BBIFileRead>(file: &mut R) -> Result<BBIFileInfo, BBIFileReadInfoError> {
    let mut file = file.raw_reader();

    let mut header_data = BytesMut::zeroed(64);
    file.read_exact(&mut header_data)?;

    let magic = header_data.get_u32();
    let (filetype, endianness) = match magic {
        _ if magic == BIGWIG_MAGIC.to_le() => (BBIFile::BigWig, Endianness::Big),
        _ if magic == BIGWIG_MAGIC.to_be() => (BBIFile::BigWig, Endianness::Little),
        _ if magic == BIGBED_MAGIC.to_le() => (BBIFile::BigBed, Endianness::Big),
        _ if magic == BIGBED_MAGIC.to_be() => (BBIFile::BigBed, Endianness::Little),
        _ => return Err(BBIFileReadInfoError::UnknownMagic),
    };

    let (
        version,
        zoom_levels,
        chromosome_tree_offset,
        full_data_offset,
        full_index_offset,
        field_count,
        defined_field_count,
        auto_sql_offset,
        total_summary_offset,
        uncompress_buf_size,
    ) = match endianness {
        Endianness::Big => {
            let version = header_data.get_u16();
            let zoom_levels = header_data.get_u16();
            let chromosome_tree_offset = header_data.get_u64();
            let full_data_offset = header_data.get_u64();
            let full_index_offset = header_data.get_u64();
            let field_count = header_data.get_u16();
            let defined_field_count = header_data.get_u16();
            let auto_sql_offset = header_data.get_u64();
            let total_summary_offset = header_data.get_u64();
            let uncompress_buf_size = header_data.get_u32();
            let _reserved = header_data.get_u64();

            (
                version,
                zoom_levels,
                chromosome_tree_offset,
                full_data_offset,
                full_index_offset,
                field_count,
                defined_field_count,
                auto_sql_offset,
                total_summary_offset,
                uncompress_buf_size,
            )
        }
        Endianness::Little => {
            let version = header_data.get_u16_le();
            let zoom_levels = header_data.get_u16_le();
            let chromosome_tree_offset = header_data.get_u64_le();
            let full_data_offset = header_data.get_u64_le();
            let full_index_offset = header_data.get_u64_le();
            let field_count = header_data.get_u16_le();
            let defined_field_count = header_data.get_u16_le();
            let auto_sql_offset = header_data.get_u64_le();
            let total_summary_offset = header_data.get_u64_le();
            let uncompress_buf_size = header_data.get_u32_le();
            let _reserved = header_data.get_u64_le();

            (
                version,
                zoom_levels,
                chromosome_tree_offset,
                full_data_offset,
                full_index_offset,
                field_count,
                defined_field_count,
                auto_sql_offset,
                total_summary_offset,
                uncompress_buf_size,
            )
        }
    };

    let header = BBIHeader {
        endianness,
        version,
        zoom_levels,
        chromosome_tree_offset,
        full_data_offset,
        full_index_offset,
        full_index_tree_offset: None,
        full_index_item_count: None,
        field_count,
        defined_field_count,
        auto_sql_offset,
        total_summary_offset,
        uncompress_buf_size,
    };

    let zoom_headers = read_zoom_headers(file, &header)?;

    // TODO: could instead store this as an Option and only read when needed
    file.seek(SeekFrom::Start(header.chromosome_tree_offset))?;

    let mut header_data = BytesMut::zeroed(32);
    file.read_exact(&mut header_data)?;

    let (key_size, val_size, item_count) = match endianness {
        Endianness::Big => {
            let magic = header_data.get_u32();
            if magic != CHROM_TREE_MAGIC {
                return Err(BBIFileReadInfoError::InvalidChroms);
            }

            let _block_size = header_data.get_u32();
            let key_size = header_data.get_u32();
            let val_size = header_data.get_u32();
            let item_count = header_data.get_u64();
            let _reserved = header_data.get_u64();

            (key_size, val_size, item_count)
        }
        Endianness::Little => {
            let magic = header_data.get_u32_le();
            if magic != CHROM_TREE_MAGIC {
                return Err(BBIFileReadInfoError::InvalidChroms);
            }

            let _block_size = header_data.get_u32_le();
            let key_size = header_data.get_u32_le();
            let val_size = header_data.get_u32_le();
            let item_count = header_data.get_u64_le();
            let _reserved = header_data.get_u64_le();

            (key_size, val_size, item_count)
        }
    };

    assert_eq!(val_size, 8u32);

    let mut chrom_info = Vec::with_capacity(item_count as usize);
    read_chrom_tree_block(&mut file, endianness, &mut chrom_info, key_size)
        .map_err(|_| BBIFileReadInfoError::InvalidChroms)?;

    let info = BBIFileInfo {
        filetype,
        header,
        zoom_headers,
        chrom_info,
    };

    Ok(info)
}

fn read_zoom_headers<R: SeekableRead>(
    file: &mut R,
    header: &BBIHeader,
) -> io::Result<Vec<ZoomHeader>> {
    let endianness = header.endianness;
    let mut header_data = BytesMut::zeroed((header.zoom_levels as usize) * 24);
    file.read_exact(&mut header_data)?;

    let mut zoom_headers = vec![];
    match endianness {
        Endianness::Big => {
            for _ in 0..header.zoom_levels {
                let reduction_level = header_data.get_u32();
                let _reserved = header_data.get_u32();
                let data_offset = header_data.get_u64();
                let index_offset = header_data.get_u64();

                zoom_headers.push(ZoomHeader {
                    reduction_level,
                    data_offset,
                    index_offset,
                    index_tree_offset: None,
                    index_item_count: None,
                });
            }
        }
        Endianness::Little => {
            for _ in 0..header.zoom_levels {
                let reduction_level = header_data.get_u32_le();
                let _reserved = header_data.get_u32_le();
                let data_offset = header_data.get_u64_le();
                let index_offset = header_data.get_u64_le();

                zoom_headers.push(ZoomHeader {
                    reduction_level,
                    data_offset,
                    index_offset,
                    index_tree_offset: None,
                    index_item_count: None,
                });
            }
        }
    };

    Ok(zoom_headers)
}

#[derive(Error, Debug)]
enum ChromTreeBlockReadError {
    #[error("{}", .0)]
    InvalidFile(String),
    #[error("Error occurred: {}", .0)]
    IoError(#[from] io::Error),
}

fn read_chrom_tree_block<R: SeekableRead>(
    f: &mut R,
    endianness: Endianness,
    chroms: &mut Vec<ChromInfo>,
    key_size: u32,
) -> Result<(), ChromTreeBlockReadError> {
    let mut header_data = BytesMut::zeroed(4);
    f.read_exact(&mut header_data)?;

    let isleaf = header_data.get_u8();
    let _reserved = header_data.get_u8();
    let count = match endianness {
        Endianness::Big => header_data.get_u16(),
        Endianness::Little => header_data.get_u16_le(),
    };

    if isleaf == 1 {
        let mut bytes = BytesMut::zeroed((key_size as usize + 8) * (count as usize));
        f.read_exact(&mut bytes)?;

        for _ in 0..count {
            let key_string = match std::str::from_utf8(&bytes.as_ref()[0..(key_size as usize)]) {
                Ok(s) => s.trim_matches(char::from(0)).to_owned(),
                Err(_) => {
                    return Err(ChromTreeBlockReadError::InvalidFile(
                        "Invalid file format: Invalid utf-8 string.".to_owned(),
                    ))
                }
            };
            bytes.advance(key_size as usize);

            let (chrom_id, chrom_size) = match endianness {
                Endianness::Big => (bytes.get_u32(), bytes.get_u32()),
                Endianness::Little => (bytes.get_u32_le(), bytes.get_u32_le()),
            };
            chroms.push(ChromInfo {
                name: key_string,
                id: chrom_id,
                length: chrom_size,
            });
        }
    } else {
        // First, go through and get child blocks
        let mut children: Vec<u64> = vec![];
        children.reserve_exact(count as usize);

        let mut bytes = BytesMut::zeroed((key_size as usize + 8) * (count as usize));
        f.read_exact(&mut bytes)?;

        for _ in 0..count {
            // We don't need this, but have to read it
            bytes.advance(key_size as usize);

            // TODO: could add specific find here by comparing key string
            let child_offset = match endianness {
                Endianness::Big => bytes.get_u64(),
                Endianness::Little => bytes.get_u64_le(),
            };
            children.push(child_offset);
        }
        // Then go through each child block
        for child in children {
            f.seek(SeekFrom::Start(child))?;
            read_chrom_tree_block(f, endianness, chroms, key_size)?;
        }
    }
    Ok(())
}

#[inline]
fn compare_position(chrom1: u32, chrom1_base: u32, chrom2: u32, chrom2_base: u32) -> i8 {
    if chrom1 < chrom2 {
        -1
    } else if chrom1 > chrom2 {
        1
    } else if chrom1_base < chrom2_base {
        -1
    } else if chrom1_base > chrom2_base {
        1
    } else {
        0
    }
}

#[inline]
fn overlaps(
    chromq: u32,
    chromq_start: u32,
    chromq_end: u32,
    chromb1: u32,
    chromb1_start: u32,
    chromb2: u32,
    chromb2_end: u32,
) -> bool {
    compare_position(chromq, chromq_start, chromb2, chromb2_end) <= 0
        && compare_position(chromq, chromq_end, chromb1, chromb1_start) >= 0
}

#[derive(Copy, Clone, Debug)]
pub(crate) struct CirTreeNodeLeaf {
    start_chrom_ix: u32,
    start_base: u32,
    end_chrom_ix: u32,
    end_base: u32,
    data_offset: u64,
    data_size: u64,
}

#[derive(Copy, Clone, Debug)]
pub(crate) struct CirTreeNodeNonLeaf {
    start_chrom_ix: u32,
    start_base: u32,
    end_chrom_ix: u32,
    end_base: u32,
    node_offset: u64,
}

pub(crate) struct CirTreeLeafItemIterator {
    endianness: Endianness,
    i: usize,
    count: usize,
    bytes: Vec<u8>,
}

impl Iterator for CirTreeLeafItemIterator {
    type Item = CirTreeNodeLeaf;

    fn next(&mut self) -> Option<Self::Item> {
        let bytes = &self.bytes;
        let i = self.i;
        if i >= self.count {
            return None;
        }
        self.i += 1;

        let istart = i * 32;
        let bytes: &[u8; 32] = &bytes[istart..istart + 32].try_into().unwrap();
        let (start_chrom_ix, start_base, end_chrom_ix, end_base, data_offset, data_size) =
            match self.endianness {
                Endianness::Big => {
                    let start_chrom_ix =
                        u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
                    let start_base = u32::from_be_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]);
                    let end_chrom_ix =
                        u32::from_be_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]);
                    let end_base = u32::from_be_bytes([bytes[12], bytes[13], bytes[14], bytes[15]]);
                    let data_offset = u64::from_be_bytes([
                        bytes[16], bytes[17], bytes[18], bytes[19], bytes[20], bytes[21],
                        bytes[22], bytes[23],
                    ]);
                    let data_size = u64::from_be_bytes([
                        bytes[24], bytes[25], bytes[26], bytes[27], bytes[28], bytes[29],
                        bytes[30], bytes[31],
                    ]);

                    (
                        start_chrom_ix,
                        start_base,
                        end_chrom_ix,
                        end_base,
                        data_offset,
                        data_size,
                    )
                }
                Endianness::Little => {
                    let start_chrom_ix =
                        u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
                    let start_base = u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]);
                    let end_chrom_ix =
                        u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]);
                    let end_base = u32::from_le_bytes([bytes[12], bytes[13], bytes[14], bytes[15]]);
                    let data_offset = u64::from_le_bytes([
                        bytes[16], bytes[17], bytes[18], bytes[19], bytes[20], bytes[21],
                        bytes[22], bytes[23],
                    ]);
                    let data_size = u64::from_le_bytes([
                        bytes[24], bytes[25], bytes[26], bytes[27], bytes[28], bytes[29],
                        bytes[30], bytes[31],
                    ]);

                    (
                        start_chrom_ix,
                        start_base,
                        end_chrom_ix,
                        end_base,
                        data_offset,
                        data_size,
                    )
                }
            };

        Some(CirTreeNodeLeaf {
            start_chrom_ix,
            start_base,
            end_chrom_ix,
            end_base,
            data_offset,
            data_size,
        })
    }
}
fn cir_tree_leaf_items<R: SeekableRead>(
    file: &mut R,
    endianness: Endianness,
    count: usize,
) -> io::Result<CirTreeLeafItemIterator> {
    let mut bytes = vec![0u8; count * 32];
    file.read_exact(&mut bytes)?;

    Ok(CirTreeLeafItemIterator {
        endianness,
        i: 0,
        count,
        bytes,
    })
}

pub(crate) struct CirTreeNonLeafItemsIterator {
    endianness: Endianness,
    i: usize,
    count: usize,
    bytes: Vec<u8>,
}

impl Iterator for CirTreeNonLeafItemsIterator {
    type Item = CirTreeNodeNonLeaf;

    fn next(&mut self) -> Option<Self::Item> {
        let bytes = &self.bytes;
        let i = self.i;
        if i >= self.count {
            return None;
        }
        self.i += 1;

        let istart = i * 24;
        let bytes: &[u8; 24] = &bytes[istart..istart + 24].try_into().unwrap();
        let (start_chrom_ix, start_base, end_chrom_ix, end_base, data_offset) = match self
            .endianness
        {
            Endianness::Big => {
                let start_chrom_ix = u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
                let start_base = u32::from_be_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]);
                let end_chrom_ix = u32::from_be_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]);
                let end_base = u32::from_be_bytes([bytes[12], bytes[13], bytes[14], bytes[15]]);
                let data_offset = u64::from_be_bytes([
                    bytes[16], bytes[17], bytes[18], bytes[19], bytes[20], bytes[21], bytes[22],
                    bytes[23],
                ]);

                (
                    start_chrom_ix,
                    start_base,
                    end_chrom_ix,
                    end_base,
                    data_offset,
                )
            }
            Endianness::Little => {
                let start_chrom_ix = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
                let start_base = u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]);
                let end_chrom_ix = u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]);
                let end_base = u32::from_le_bytes([bytes[12], bytes[13], bytes[14], bytes[15]]);
                let data_offset = u64::from_le_bytes([
                    bytes[16], bytes[17], bytes[18], bytes[19], bytes[20], bytes[21], bytes[22],
                    bytes[23],
                ]);

                (
                    start_chrom_ix,
                    start_base,
                    end_chrom_ix,
                    end_base,
                    data_offset,
                )
            }
        };

        Some(CirTreeNodeNonLeaf {
            start_chrom_ix,
            start_base,
            end_chrom_ix,
            end_base,
            node_offset: data_offset,
        })
    }
}
fn cir_tree_non_leaf_items<R: SeekableRead>(
    file: &mut R,
    endianness: Endianness,
    count: usize,
) -> io::Result<CirTreeNonLeafItemsIterator> {
    let mut bytes = vec![0u8; count * 24];
    file.read_exact(&mut bytes)?;

    Ok(CirTreeNonLeafItemsIterator {
        endianness,
        i: 0,
        count,
        bytes,
    })
}

pub(crate) struct CirTreeBlockSearchIter<'a, R: BBIFileRead> {
    remaining_childblocks: VecDeque<u64>,

    file: &'a mut R,
    endianness: Endianness,
    chrom_ix: u32,
    start: u32,
    end: u32,
}

impl<'a, R: BBIFileRead> Iterator for CirTreeBlockSearchIter<'a, R> {
    type Item = io::Result<SmallVec<[Block; 4]>>;
    fn next(&mut self) -> Option<Self::Item> {
        let file = &mut *self.file;
        let endianness = self.endianness;
        let chrom_ix = self.chrom_ix;
        let start = self.start;
        let end = self.end;

        let node_offset = self.remaining_childblocks.pop_front()?;

        let (new_childblocks, blocks) =
            match file.blocks_for_cir_tree_node(endianness, node_offset, chrom_ix, start, end) {
                Ok(d) => d,
                Err(e) => return Some(Err(e)),
            };

        for child in new_childblocks.into_iter().rev() {
            self.remaining_childblocks.push_front(child);
        }

        Some(Ok(blocks))
    }
}

pub(crate) enum CirTreeNodeIterator<
    L: Iterator<Item = CirTreeNodeLeaf> = CirTreeLeafItemIterator,
    N: Iterator<Item = CirTreeNodeNonLeaf> = CirTreeNonLeafItemsIterator,
> {
    Leaf(L),
    NonLeaf(N),
}

pub(crate) fn read_node<R: SeekableRead>(
    file: &mut R,
    node_offset: u64,
    endianness: Endianness,
) -> io::Result<CirTreeNodeIterator> {
    read_node_bounded(file, node_offset, endianness, None)
}

fn validate_cir_tree_node_span(
    node_offset: u64,
    count: usize,
    record_size: u64,
    index_end: Option<u64>,
) -> io::Result<()> {
    let body_size = u64::try_from(count)
        .ok()
        .and_then(|count| count.checked_mul(record_size))
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "cir-tree node size overflows")
        })?;
    let node_end = node_offset
        .checked_add(4)
        .and_then(|offset| offset.checked_add(body_size))
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "cir-tree node range overflows")
        })?;
    if index_end.is_some_and(|index_end| node_end > index_end) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "cir-tree node records extend outside index",
        ));
    }
    Ok(())
}

fn read_node_bounded<R: SeekableRead>(
    file: &mut R,
    node_offset: u64,
    endianness: Endianness,
    index_end: Option<u64>,
) -> io::Result<CirTreeNodeIterator> {
    match file.seek(SeekFrom::Start(node_offset)) {
        Err(e) => return Err(e),
        Ok(_) => {}
    };

    let mut header_data = BytesMut::zeroed(4);
    match file.read_exact(&mut header_data) {
        Err(e) => return Err(e),
        Ok(_) => {}
    }

    let isleaf: u8 = header_data.get_u8();
    if isleaf > 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("cir-tree node has invalid leaf flag {isleaf}"),
        ));
    }
    let _reserved = header_data.get_u8();

    let count = match endianness {
        Endianness::Big => header_data.get_u16(),
        Endianness::Little => header_data.get_u16_le(),
    };
    let record_size = if isleaf == 1 { 32 } else { 24 };
    validate_cir_tree_node_span(node_offset, usize::from(count), record_size, index_end)?;

    let iter = if isleaf == 1 {
        let iter = match cir_tree_leaf_items(file, endianness, count as usize) {
            Ok(v) => v,
            Err(e) => return Err(e),
        };
        CirTreeNodeIterator::Leaf(iter)
    } else {
        let iter = match cir_tree_non_leaf_items(file, endianness, count as usize) {
            Ok(v) => v,
            Err(e) => return Err(e),
        };
        CirTreeNodeIterator::NonLeaf(iter)
    };
    Ok(iter)
}

fn nodes_overlapping<
    L: Iterator<Item = CirTreeNodeLeaf>,
    N: Iterator<Item = CirTreeNodeNonLeaf>,
>(
    iter: CirTreeNodeIterator<L, N>,
    chrom_ix: u32,
    start: u32,
    end: u32,
) -> (SmallVec<[u64; 4]>, SmallVec<[Block; 4]>) {
    match iter {
        CirTreeNodeIterator::Leaf(iter) => {
            let mut blocks: SmallVec<[_; 4]> = smallvec![];
            for child in iter {
                let block_overlaps = overlaps(
                    chrom_ix,
                    start,
                    end,
                    child.start_chrom_ix,
                    child.start_base,
                    child.end_chrom_ix,
                    child.end_base,
                );
                if block_overlaps {
                    blocks.push(Block {
                        offset: child.data_offset,
                        size: child.data_size,
                    });
                }
            }
            (smallvec![], blocks)
        }
        CirTreeNodeIterator::NonLeaf(iter) => {
            let mut new_childblocks: SmallVec<[_; 4]> = smallvec![];
            for child in iter {
                let block_overlaps = overlaps(
                    chrom_ix,
                    start,
                    end,
                    child.start_chrom_ix,
                    child.start_base,
                    child.end_chrom_ix,
                    child.end_base,
                );
                if block_overlaps {
                    new_childblocks.push(child.node_offset);
                }
            }
            (new_childblocks, smallvec![])
        }
    }
}

/// Gets the data (uncompressed, if applicable) from a given block
fn read_block_data<R: SeekableRead>(
    info: &BBIFileInfo,
    read: &mut R,
    block: &Block,
) -> io::Result<Vec<u8>> {
    let uncompress_buf_size = info.header.uncompress_buf_size as usize;

    read.seek(SeekFrom::Start(block.offset))?;

    let mut raw_data = vec![0u8; block.size as usize];
    read.read_exact(&mut raw_data)?;
    let block_data: Vec<u8> = if uncompress_buf_size > 0 {
        let mut decompressor = Decompressor::new();
        let mut outbuf = vec![0; uncompress_buf_size];
        let decompressed = decompressor
            .zlib_decompress(&raw_data, &mut outbuf)
            .unwrap();
        outbuf.truncate(decompressed);
        outbuf
    } else {
        raw_data
    };

    Ok(block_data)
}

pub(crate) fn get_zoom_block_values<B: BBIRead>(
    bbifile: &mut B,
    block: Block,
    known_offset: &mut u64,
    chrom: u32,
    start: u32,
    end: u32,
) -> Result<std::vec::IntoIter<ZoomRecord>, BBIReadError> {
    let (read, info) = bbifile.reader_and_info();
    let data = read.get_block_data(info, &block)?;
    let mut bytes = BytesMut::with_capacity(data.len());
    bytes.extend_from_slice(&data);

    let len = bytes.len();
    assert_eq!(len % (4 * 8), 0);
    let itemcount = len / (4 * 8);
    let mut records = Vec::with_capacity(itemcount);

    let endianness = bbifile.info().header.endianness;

    match endianness {
        Endianness::Big => {
            for _ in 0..itemcount {
                let chrom_id = bytes.get_u32();
                let chrom_start = bytes.get_u32();
                let chrom_end = bytes.get_u32();
                let bases_covered = u64::from(bytes.get_u32());
                let min_val = f64::from(bytes.get_f32());
                let max_val = f64::from(bytes.get_f32());
                let sum = f64::from(bytes.get_f32());
                let sum_squares = f64::from(bytes.get_f32());
                if chrom_id == chrom && chrom_end >= start && chrom_start <= end {
                    records.push(ZoomRecord {
                        chrom: chrom_id,
                        start: chrom_start,
                        end: chrom_end,
                        summary: Summary {
                            total_items: 0,
                            bases_covered,
                            min_val,
                            max_val,
                            sum,
                            sum_squares,
                        },
                    });
                }
            }
        }
        Endianness::Little => {
            for _ in 0..itemcount {
                let chrom_id = bytes.get_u32_le();
                let chrom_start = bytes.get_u32_le();
                let chrom_end = bytes.get_u32_le();
                let bases_covered = u64::from(bytes.get_u32_le());
                let min_val = f64::from(bytes.get_f32_le());
                let max_val = f64::from(bytes.get_f32_le());
                let sum = f64::from(bytes.get_f32_le());
                let sum_squares = f64::from(bytes.get_f32_le());
                if chrom_id == chrom && chrom_end >= start && chrom_start <= end {
                    records.push(ZoomRecord {
                        chrom: chrom_id,
                        start: chrom_start,
                        end: chrom_end,
                        summary: Summary {
                            total_items: 0,
                            bases_covered,
                            min_val,
                            max_val,
                            sum,
                            sum_squares,
                        },
                    });
                }
            }
        }
    }

    *known_offset = block.offset + block.size;
    Ok(records.into_iter())
}

pub struct ZoomIntervalIter<R, B> {
    _r: std::marker::PhantomData<R>,
    bbifile: B,
    known_offset: u64,
    blocks: std::vec::IntoIter<Block>,
    vals: Option<std::vec::IntoIter<ZoomRecord>>,
    // (chrom, start, end)
    interval: Option<(u32, u32, u32)>,
}

impl<R> Into<BigWigRead<R>> for ZoomIntervalIter<BigWigRead<R>, BigWigRead<R>> {
    fn into(self) -> BigWigRead<R> {
        self.bbifile
    }
}

impl<R> Into<BigBedRead<R>> for ZoomIntervalIter<BigBedRead<R>, BigBedRead<R>> {
    fn into(self) -> BigBedRead<R> {
        self.bbifile
    }
}

impl<R, B> ZoomIntervalIter<R, B> {
    pub fn new(
        bbifile: B,
        blocks: std::vec::IntoIter<Block>,
        chrom: u32,
        start: u32,
        end: u32,
    ) -> Self {
        ZoomIntervalIter {
            _r: std::marker::PhantomData,
            bbifile,
            known_offset: 0,
            blocks,
            vals: None,
            interval: Some((chrom, start, end)),
        }
    }

    pub fn empty(bbifile: B) -> Self {
        ZoomIntervalIter {
            _r: std::marker::PhantomData,
            bbifile,
            known_offset: 0,
            blocks: vec![].into_iter(),
            vals: None,
            interval: None,
        }
    }
}

impl<R, B> Iterator for ZoomIntervalIter<R, B>
where
    R: BBIRead,
    B: BorrowMut<R>,
{
    type Item = Result<ZoomRecord, BBIReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        let Some((chrom, start, end)) = self.interval else {
            return None;
        };
        loop {
            match &mut self.vals {
                Some(vals) => match vals.next() {
                    Some(v) => {
                        return Some(Ok(v));
                    }
                    None => {
                        self.vals = None;
                    }
                },
                None => {
                    let current_block = self.blocks.next()?;
                    let file = self.bbifile.borrow_mut();
                    match get_zoom_block_values(
                        file,
                        current_block,
                        &mut self.known_offset,
                        chrom,
                        start,
                        end,
                    ) {
                        Ok(vals) => {
                            self.vals = Some(vals);
                        }
                        Err(e) => {
                            return Some(Err(e));
                        }
                    }
                }
            }
        }
    }
}

pub struct ZoomIntervalCursor<B: BBIRead>(pub(crate) ZoomIntervalIter<B, B>);

impl<B: BBIRead> ZoomIntervalCursor<B> {
    pub fn interval(
        &mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
        reduction_level: u32,
    ) -> Result<(), ZoomIntervalError> {
        let cir_tree = self.0.bbifile.zoom_cir_tree(reduction_level)?;
        let chrom = self.0.bbifile.info().chrom_id(chrom_name)?;
        let (reader, info) = self.0.bbifile.reader_and_info();
        let blocks = search_cir_tree(info, reader, cir_tree, chrom_name, start, end)?;

        self.0.known_offset = 0;
        self.0.blocks = blocks.into_iter();
        self.0.vals = None;
        self.0.interval = Some((chrom, start, end));

        Ok(())
    }

    pub fn into_inner(self) -> B {
        self.0.bbifile
    }
}
