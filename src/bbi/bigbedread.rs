use std::borrow::BorrowMut;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Seek, SeekFrom};
use std::vec::Vec;

use byteorder::{BigEndian, LittleEndian, ReadBytesExt};
use byteordered::ByteOrdered;
use bytes::{Buf, BytesMut};
use itertools::Itertools;
use thiserror::Error;

use crate::bbi::{BBIFile, BedEntry, ZoomRecord};
use crate::bbiread::{
    read_info, BBIFileInfo, BBIFileReadInfoError, BBIRead, BBIReadError, Block, ChromInfo,
    ZoomIntervalIter,
};
use crate::internal::BBIReadInternal;
use crate::utils::reopen::{Reopen, ReopenableFile, SeekableRead};
use crate::{search_cir_tree, BBIFileRead, CachedBBIFileRead, Summary, ZoomIntervalError};

struct IntervalIter<I, R, B>
where
    I: Iterator<Item = Block> + Send,
    R: BBIFileRead,
    B: BorrowMut<BigBedRead<R>>,
{
    r: std::marker::PhantomData<R>,
    bigbed: B,
    known_offset: u64,
    blocks: I,
    vals: Option<std::vec::IntoIter<BedEntry>>,
    expected_chrom: u32,
    start: u32,
    end: u32,
}

impl<I, R, B> Iterator for IntervalIter<I, R, B>
where
    I: Iterator<Item = Block> + Send,
    R: BBIFileRead,
    B: BorrowMut<BigBedRead<R>>,
{
    type Item = Result<BedEntry, BBIReadError>;

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
                    // TODO: Could minimize this by chunking block reads
                    let current_block = self.blocks.next()?;
                    match get_block_entries(
                        self.bigbed.borrow_mut(),
                        current_block,
                        &mut self.known_offset,
                        self.expected_chrom,
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

/// Possible errors encountered when opening a bigBed file to read
#[derive(Error, Debug)]
pub enum BigBedReadOpenError {
    #[error("File is not a bigBed.")]
    NotABigBed,
    #[error("The chromosomes are invalid.")]
    InvalidChroms,
    #[error("An error occurred: {}", .0)]
    IoError(#[from] io::Error),
}

impl From<BBIFileReadInfoError> for BigBedReadOpenError {
    fn from(error: BBIFileReadInfoError) -> Self {
        match error {
            BBIFileReadInfoError::UnknownMagic => BigBedReadOpenError::NotABigBed,
            BBIFileReadInfoError::InvalidChroms => BigBedReadOpenError::InvalidChroms,
            BBIFileReadInfoError::IoError(e) => BigBedReadOpenError::IoError(e),
        }
    }
}

/// The struct used to read a bigBed file
pub struct BigBedRead<R> {
    pub(super) info: BBIFileInfo,
    pub(super) read: R,
}

impl<R: Reopen> Reopen for BigBedRead<R> {
    fn reopen(&self) -> io::Result<Self> {
        Ok(BigBedRead {
            info: self.info.clone(),
            read: self.read.reopen()?,
        })
    }
}

impl<R: BBIFileRead> BBIRead for BigBedRead<R> {
    fn info(&self) -> &BBIFileInfo {
        &self.info
    }

    fn chroms(&self) -> &[ChromInfo] {
        &self.info.chrom_info
    }
}

impl<R: BBIFileRead> BBIReadInternal for BigBedRead<R> {
    type Read = R;

    fn reader(&mut self) -> &mut R {
        &mut self.read
    }

    fn reader_and_info(&mut self) -> (&mut Self::Read, &mut BBIFileInfo) {
        (&mut self.read, &mut self.info)
    }
}

impl<R> BigBedRead<R> {
    /// Get basic info about this bigBed
    pub fn info(&self) -> &BBIFileInfo {
        &self.info
    }

    /// Gets the chromosomes present in this bigBed
    pub fn chroms(&self) -> &[ChromInfo] {
        &self.info.chrom_info
    }
}

impl BigBedRead<ReopenableFile> {
    /// Opens a new `BigBedRead` from a given path as a file.
    pub fn open_file(path: &str) -> Result<Self, BigBedReadOpenError> {
        let reopen = ReopenableFile {
            path: path.to_string(),
            file: File::open(path)?,
        };
        let b = BigBedRead::open(reopen);
        if b.is_err() {
            eprintln!("Error when opening: {}", path);
        }
        b
    }
}

impl<R> BigBedRead<R>
where
    R: SeekableRead,
{
    /// Converts this `BigBedRead`` to where the `BBIFileRead` caches index
    /// access and block data
    pub fn cached(self) -> BigBedRead<CachedBBIFileRead<R>> {
        let read = CachedBBIFileRead::new(self.read);
        BigBedRead {
            read,
            info: self.info,
        }
    }
}

impl<R: BBIFileRead> BigBedRead<R> {
    /// Opens a new `BigBedRead` for a given type that implements both `Read` and `Seek`
    pub fn open(mut read: R) -> Result<Self, BigBedReadOpenError> {
        let info = read_info(&mut read.raw_reader())?;
        match info.filetype {
            BBIFile::BigBed => {}
            _ => return Err(BigBedReadOpenError::NotABigBed),
        }

        Ok(BigBedRead { info, read })
    }

    /// Reads the autosql from this bigBed
    pub fn autosql(&mut self) -> Result<String, BBIReadError> {
        let auto_sql_offset = self.info.header.auto_sql_offset;
        let reader = self.reader().raw_reader();
        let mut reader = BufReader::new(reader);
        reader.seek(SeekFrom::Start(auto_sql_offset))?;
        let mut buffer = Vec::new();
        reader.read_until(b'\0', &mut buffer)?;
        buffer.pop();
        let autosql = String::from_utf8(buffer)
            .map_err(|_| BBIReadError::InvalidFile("Invalid autosql: not UTF-8".to_owned()))?;
        Ok(autosql)
    }

    pub fn item_count(&mut self) -> Result<u64, BBIReadError> {
        let header = self.info.header;
        let reader = self.reader().raw_reader();
        let mut reader = BufReader::new(reader);
        reader.seek(SeekFrom::Start(header.full_index_offset))?;
        Ok(match header.endianness {
            byteordered::Endianness::Big => reader.read_u64::<BigEndian>()?,
            byteordered::Endianness::Little => reader.read_u64::<LittleEndian>()?,
        })
    }

    /// Returns the summary data from bigBed
    ///
    /// Note: For version 1 of bigBeds, there is no total summary. In that
    /// case, 0 is returned for all of the summary except total items. If this
    /// matters to you, you can check the version using
    /// `info().header.version > 1`.
    pub fn get_summary(&mut self) -> io::Result<Summary> {
        let endianness = self.info.header.endianness;
        let summary_offset = self.info.header.total_summary_offset;
        let data_offset = self.info.header.full_data_offset;
        let reader = self.reader().raw_reader();
        let mut reader = ByteOrdered::runtime(reader, endianness);
        let (bases_covered, min_val, max_val, sum, sum_squares) = if summary_offset != 0 {
            reader.seek(SeekFrom::Start(summary_offset))?;
            (
                reader.read_u64()?,
                reader.read_f64()?,
                reader.read_f64()?,
                reader.read_f64()?,
                reader.read_f64()?,
            )
        } else {
            (0, 0.0, 0.0, 0.0, 0.0)
        };
        reader.seek(SeekFrom::Start(data_offset))?;
        let total_items = reader.read_u64()?;
        Ok(Summary {
            total_items,
            bases_covered,
            min_val,
            max_val,
            sum,
            sum_squares,
        })
    }

    /// For a given chromosome, start, and end, returns an `Iterator` of the
    /// intersecting `BedEntry`s. The resulting iterator takes a mutable reference
    /// of this `BigBedRead`.
    pub fn get_interval<'a>(
        &'a mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
    ) -> Result<impl Iterator<Item = Result<BedEntry, BBIReadError>> + 'a, BBIReadError> {
        let cir_tree = self.full_data_cir_tree()?;
        let blocks = search_cir_tree(&self.info, &mut self.read, cir_tree, chrom_name, start, end)?;
        // TODO: this is only for asserting that the chrom is what we expect
        let chrom_ix = self
            .info()
            .chrom_info
            .iter()
            .find(|&x| x.name == chrom_name)
            .unwrap()
            .id;
        Ok(IntervalIter {
            r: std::marker::PhantomData,
            bigbed: self,
            known_offset: 0,
            blocks: blocks.into_iter(),
            vals: None,
            expected_chrom: chrom_ix,
            start,
            end,
        })
    }

    /// For a given chromosome, start, and end, returns an `Iterator` of the
    /// intersecting `BedEntry`s. The resulting iterator takes this `BigBedRead`
    /// by value.
    pub fn get_interval_move(
        mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
    ) -> Result<impl Iterator<Item = Result<BedEntry, BBIReadError>>, BBIReadError> {
        let cir_tree = self.full_data_cir_tree()?;
        let blocks = search_cir_tree(&self.info, &mut self.read, cir_tree, chrom_name, start, end)?;
        // TODO: this is only for asserting that the chrom is what we expect
        let chrom_ix = self
            .info()
            .chrom_info
            .iter()
            .find(|&x| x.name == chrom_name)
            .unwrap()
            .id;
        Ok(IntervalIter {
            r: std::marker::PhantomData,
            bigbed: self,
            known_offset: 0,
            blocks: blocks.into_iter(),
            vals: None,
            expected_chrom: chrom_ix,
            start,
            end,
        })
    }

    /// For a given chromosome, start, and end, returns an `Iterator` of the
    /// intersecting `ZoomRecord`s.
    pub fn get_zoom_interval<'a>(
        &'a mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
        reduction_level: u32,
    ) -> Result<impl Iterator<Item = Result<ZoomRecord, BBIReadError>> + 'a, ZoomIntervalError>
    {
        let cir_tree = self
            .zoom_cir_tree(reduction_level)
            .map_err(|_| ZoomIntervalError::ReductionLevelNotFound)?;

        let chrom = self.info.chrom_id(chrom_name)?;

        let blocks = search_cir_tree(&self.info, &mut self.read, cir_tree, chrom_name, start, end)?;
        Ok(ZoomIntervalIter::new(
            self,
            blocks.into_iter(),
            chrom,
            start,
            end,
        ))
    }
}

// TODO: remove expected_chrom
fn get_block_entries<R: BBIFileRead>(
    bigbed: &mut BigBedRead<R>,
    block: Block,
    known_offset: &mut u64,
    expected_chrom: u32,
    start: u32,
    end: u32,
) -> Result<std::vec::IntoIter<BedEntry>, BBIReadError> {
    let data = bigbed.read.get_block_data(&bigbed.info, &block)?;
    let mut bytes = BytesMut::with_capacity(data.len());
    bytes.extend_from_slice(&data);
    let mut entries: Vec<BedEntry> = Vec::new();

    let mut read_entry = || -> Result<Option<BedEntry>, BBIReadError> {
        if bytes.len() < 12 {
            return Ok(None);
        }
        let (chrom_id, chrom_start, chrom_end) = match bigbed.info.header.endianness {
            byteordered::Endianness::Big => (bytes.get_u32(), bytes.get_u32(), bytes.get_u32()),
            byteordered::Endianness::Little => {
                (bytes.get_u32_le(), bytes.get_u32_le(), bytes.get_u32_le())
            }
        };
        if chrom_start == 0 && chrom_end == 0 {
            return Err(BBIReadError::InvalidFile(
                "Chrom start and end both equal 0.".to_owned(),
            ));
        }
        // FIXME: should this just return empty?
        assert_eq!(
            chrom_id, expected_chrom,
            "BUG: bigBed had multiple chroms in a section"
        );
        let nul = bytes.iter().find_position(|b| **b == b'\0');
        let s = match nul {
            Some((pos, _)) => {
                let b = bytes.split_to(pos);
                bytes.get_u8();
                b.to_vec()
            }
            None => bytes.to_vec(),
        };
        let rest = String::from_utf8(s).unwrap();
        Ok(Some(BedEntry {
            start: chrom_start,
            end: chrom_end,
            rest,
        }))
    };
    while let Some(entry) = read_entry()? {
        if entry.end >= start && entry.start <= end {
            entries.push(entry);
        }
    }

    *known_offset = block.offset + block.size;
    Ok(entries.into_iter())
}
