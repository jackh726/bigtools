//#![feature(test)]

pub mod bbiread;
pub mod bbiwrite;
pub mod bigbedread;
pub mod bigbedwrite;
pub mod bigwig;
pub mod bigwigread;
pub mod bigwigwrite;
pub mod chromvalues;
pub mod filebufferedchannel;
pub mod idmap;
pub mod seekableread;
pub mod tell;
pub mod tempfilebuffer;

pub mod bedparser;
pub mod streaming_linereader;

pub mod utils;

#[cfg(feature = "remote")]
pub mod remote_file;
