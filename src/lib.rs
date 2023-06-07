/*!
Bigtools provides a high-level, performant API for reading and writing bigWig and bigBed files.

The original file format specification for bigWig and bigBed files is defined in this paper: <https://doi.org/10.1093/bioinformatics/btq351>

## Reading

## Writing

*/

pub mod bbi;
pub mod bed;
pub mod utils;

pub use bbi::*;
