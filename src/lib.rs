#![feature(test)]

pub mod bbiread;
pub mod bbiwrite;
pub mod bigbedread;
pub mod bigbedwrite;
pub mod bigwigread;
pub mod bigwigwrite;
pub mod bigwig;
pub mod idmap;
pub mod tell;
pub mod tempfilebuffer;
pub mod chromvalues;
pub mod filebufferedchannel;
pub mod seekableread;

pub mod streaming_linereader;
pub mod bedparser;

pub mod utils;

#[cfg(feature = "remote")]
pub mod remote_file;
