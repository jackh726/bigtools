/*!
Bigtools provides a high-level, performant API for reading and writing bigWig and bigBed files.

The original file format specification for bigWig and bigBed files is defined in this paper: <https://doi.org/10.1093/bioinformatics/btq351>

## Reading

The entrypoint to reading bigWigs and bigBeds is [`BigWigRead::open`] and
[`BigBedRead::open`], respectively. These take any type that implements both
[`Read`][std::io::Read] and [`Seek`][std::io::Seek]. There are also
[`BigWigRead::open_file`] and [`BigBedRead::open_file`], which take a `&str` and
will open a `File`.

Once a [`BigWigRead`] or [`BigBedRead`] have been constructed, they can be read
in a number of ways. First, the info (in the form of [`BBIFileInfo`]) is available
in `info` fields. However, to access the main data, the most common method to call
is [`BigWigRead::get_interval`] or [`BigBedRead::get_interval`], which returns an
`Iterator` of [`Value`]s or [`BedEntry`]s overlapping the provided region, respectively.

## Writing

Writing new bigWigs and bigBeds is a tad more difficult. To begin, a
[`BigWigWrite`] or [`BigBedWrite`] can be created using
[`BigWigWrite::create_file`] or [`BigBedWrite::create_file`].

Generally, bigWig and bigBed writing is done per chromosome, with compression
and io being done on an async Runtime.

The source for data to be written to bigWigs and bigBeds come from the
[`ChromData`] trait. It's effectively like a powerful `Iterator` of
`ChromData::Output` values, which itself is like peekable `Iterator` of values
(either [`Value`]s or [`BedEntry`]s for bigWigs or bigBeds respectively) over a
chromosome. The [`ChromData::advance`] method takes a function that the current
chromosome and returns a `Future` that will asynchronously process chromosomal
data. This data can be returned immediately from the `advance` method, or can
be stored (with the correct implementation) to queue multiple chromosomes
simulatenously. The
[`BedParserStreamingIterator`][crate::bbi::bedchromdata::BedParserStreamingIterator]
and
[`BedParserParallelStreamingIterator`][crate::bbi::bedchromdata::BedParserParallelStreamingIterator]
types provide serial processing of a bed-like value stream (either from a
file or an iterator) or concurrent processing from a file. Generally, these
underlying details aren't necessary unless implementing a new data source.

Given some implementation of [`ChromData`] (like [`BedParserStreamingIterator`][crate::bbi::bedchromdata::BedParserStreamingIterator]),
a bigWig can be created using [`BigWigWrite::write`] or a bigBed with
[`BigBedWrite::write`]. Both take a map of chromosome sizes, the aforementioned
data, and a `Runtime` to spawn processing on.
*/

mod bbi;
pub mod bed;
pub mod utils;

pub use bbi::*;
