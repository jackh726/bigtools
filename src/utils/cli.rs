use crate::bbiwrite::{DEFAULT_BLOCK_SIZE, DEFAULT_ITEMS_PER_SLOT};

use clap::Parser;

#[derive(Parser)]
pub struct BBIWriteArgs {
    /// Set the number of threads to use. This tool will typically use ~225% CPU on a HDD. SDDs may be higher. (IO bound)
    #[arg(short = 't', long)]
    #[arg(default_value_t = 6)]
    pub nthreads: usize,

    /// Set the maximum of zooms to create.
    #[arg(short = 'z', long)]
    #[arg(default_value_t = 10)]
    pub nzooms: u32,

    /// Don't use compression.
    #[arg(short = 'u', long)]
    #[arg(default_value_t = false)]
    pub uncompressed: bool,

    /// Sets whether the input is sorted. Can take `all`, `start`, or `none`.
    /// `all` means that the input bedGraph is sorted by chroms and start (`sort -k1,1 -k2,2n`).
    /// `start` means that the the chroms are out of order but the starts within a chrom is sorted.
    /// `none` means that the file is not sorted at all.
    /// `all` is default. `none` currently errors but may be supported in the future.
    /// Note that using a value other than `all` will not guarantee (though likely) support for third-party tools.
    #[arg(short = 's', long)]
    #[arg(default_value = "all")]
    pub sorted: String,

    /// Number of items to bundle in r-tree.
    #[arg(long)]
    #[arg(default_value_t = DEFAULT_BLOCK_SIZE)]
    pub block_size: u32,

    /// Number of data points bundled at lowest level.
    #[arg(long)]
    #[arg(default_value_t = DEFAULT_ITEMS_PER_SLOT)]
    pub items_per_slot: u32,
}
