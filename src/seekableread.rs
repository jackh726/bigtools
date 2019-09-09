use std::io::{Read, Seek};
use std::marker::Send;


pub trait SeekableRead: Seek + Read + Send {}

impl<T> SeekableRead for T where T: Seek + Read + Send {}