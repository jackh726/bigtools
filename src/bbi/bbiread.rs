use std::collections::hash_map::Entry;
use std::collections::{HashMap, VecDeque};
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

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct Block {
    pub(crate) offset: u64,
    pub(crate) size: u64,
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
    pub(crate) auto_sql_offset: u64,
    pub(crate) total_summary_offset: u64,
    pub(crate) uncompress_buf_size: u32,
}

/// Information on a chromosome in a bbi file
#[derive(Clone, Debug)]
pub struct ChromInfo {
    pub name: String,
    pub length: u32,
    pub(crate) id: u32,
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
pub struct CirTreeIndex(pub CirTreeIndexType, pub(crate) u64);

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
            if info.header.full_index_tree_offset.is_none() {
                let endianness = info.header.endianness;

                reader
                    .raw_reader()
                    .seek(SeekFrom::Start(index_offset))
                    .map_err(|e| FullDataCirTreeError::IoError(e))?;

                read_cir_tree_header(endianness, reader.raw_reader()).map_err(|e| match e {
                    Either::Left(_) => FullDataCirTreeError::UnknownMagic,
                    Either::Right(e) => FullDataCirTreeError::IoError(e),
                })?;

                info.header.full_index_tree_offset = Some(index_offset + 48);
            }
            Ok(CirTreeIndex(CirTreeIndexType::FullData, index_offset + 48))
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

            if zoom_header.index_tree_offset.is_none() {
                let endianness = info.header.endianness;

                reader
                    .raw_reader()
                    .seek(SeekFrom::Start(zoom_header.index_offset))
                    .map_err(|e| ZoomDataCirTreeError::IoError(e))?;

                read_cir_tree_header(endianness, reader.raw_reader()).map_err(|e| match e {
                    Either::Left(_) => ZoomDataCirTreeError::UnknownMagic,
                    Either::Right(e) => ZoomDataCirTreeError::IoError(e),
                })?;

                zoom_header.index_tree_offset = Some(zoom_header.index_offset + 48);
            }

            Ok(CirTreeIndex(
                CirTreeIndexType::Zoom(reduction_level),
                zoom_header.index_offset + 48,
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
) -> Result<(), Either<UnknownMagic, io::Error>> {
    let mut header_data = BytesMut::zeroed(48);
    file.read_exact(&mut header_data)
        .map_err(|e| Either::Right(e))?;

    match endianness {
        Endianness::Big => {
            let magic = header_data.get_u32();
            if magic != CIR_TREE_MAGIC {
                return Err(Either::Left(UnknownMagic));
            }

            let _blocksize = header_data.get_u32();
            let _item_count = header_data.get_u64();
            let _start_chrom_idx = header_data.get_u32();
            let _start_base = header_data.get_u32();
            let _end_chrom_idx = header_data.get_u32();
            let _end_base = header_data.get_u32();
            let _end_file_offset = header_data.get_u64();
            let _item_per_slot = header_data.get_u32();
            let _reserved = header_data.get_u32();
        }
        Endianness::Little => {
            let magic = header_data.get_u32_le();
            if magic != CIR_TREE_MAGIC {
                return Err(Either::Left(UnknownMagic));
            }

            let _blocksize = header_data.get_u32_le();
            let _item_count = header_data.get_u64_le();
            let _start_chrom_idx = header_data.get_u32_le();
            let _start_base = header_data.get_u32_le();
            let _end_chrom_idx = header_data.get_u32_le();
            let _end_base = header_data.get_u32_le();
            let _end_file_offset = header_data.get_u64_le();
            let _item_per_slot = header_data.get_u32_le();
            let _reserved = header_data.get_u32_le();
        }
    };
    Ok(())
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
            path: path.to_string(),
            file: File::open(path)?,
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

pub struct CachedBBIFileRead<S: SeekableRead> {
    read: S,
    cir_tree_node_map: HashMap<u64, Either<Vec<CirTreeNodeLeaf>, Vec<CirTreeNodeNonLeaf>>>,
    block_data: HashMap<Block, Vec<u8>>,
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
    let mut bytes = vec![0u8; (count as usize) * 32];
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
    assert!(isleaf == 1 || isleaf == 0, "Unexpected isleaf: {}", isleaf);
    let _reserved = header_data.get_u8();

    let count = match endianness {
        Endianness::Big => header_data.get_u16(),
        Endianness::Little => header_data.get_u16_le(),
    };

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
) -> Result<Box<dyn Iterator<Item = ZoomRecord> + Send>, BBIReadError> {
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
    Ok(Box::new(records.into_iter()))
}

pub(crate) struct ZoomIntervalIter<'a, I, B>
where
    I: Iterator<Item = Block> + Send,
    B: BBIRead,
{
    bbifile: &'a mut B,
    known_offset: u64,
    blocks: I,
    vals: Option<Box<dyn Iterator<Item = ZoomRecord> + Send + 'a>>,
    chrom: u32,
    start: u32,
    end: u32,
}

impl<'a, I, B> ZoomIntervalIter<'a, I, B>
where
    I: Iterator<Item = Block> + Send,
    B: BBIRead,
{
    pub fn new(bbifile: &'a mut B, blocks: I, chrom: u32, start: u32, end: u32) -> Self {
        ZoomIntervalIter {
            bbifile,
            known_offset: 0,
            blocks,
            vals: None,
            chrom,
            start,
            end,
        }
    }
}

impl<'a, I, B> Iterator for ZoomIntervalIter<'a, I, B>
where
    I: Iterator<Item = Block> + Send,
    B: BBIRead,
{
    type Item = Result<ZoomRecord, BBIReadError>;

    fn next(&mut self) -> Option<Self::Item> {
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
                    match get_zoom_block_values(
                        self.bbifile,
                        current_block,
                        &mut self.known_offset,
                        self.chrom,
                        self.start,
                        self.end,
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
