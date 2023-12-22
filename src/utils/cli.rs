use std::{ffi::OsString, str::FromStr};

use crate::bbiwrite::{DEFAULT_BLOCK_SIZE, DEFAULT_ITEMS_PER_SLOT};

use clap::Args;

pub mod bedgraphtobigwig;
pub mod bedtobigbed;
pub mod bigbedtobed;
pub mod bigwigaverageoverbed;
pub mod bigwiginfo;
pub mod bigwigmerge;
pub mod bigwigtobedgraph;
pub mod bigwigvaluesoverbed;

#[derive(Debug, Args)]
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

macro_rules! compat_replace_mut {
    (
        $a:expr;
        replace:
            $($find:literal,$replace:literal);*
        ignore:
            $($ignore:literal);*
        unimplemented:
            $($unimplemented:literal);*
    ) => {{
        use std::ffi::OsString;
        use std::str::FromStr;
        match $a.to_str() {
            $(
                Some(b) if b.starts_with($find) => {
                    *($a) = OsString::from_str(&b.replace($find, $replace)).unwrap()
                }
            )*
            $(
                Some(b) if b.starts_with($ignore) => *($a) = OsString::from_str("").unwrap(),
            )*
            $(
                Some(b) if b.starts_with($unimplemented) => {
                    panic!(
                        "Unimplemented compatibility option {}.",
                        $a.to_string_lossy()
                    );
                }
            )*
            _ => {}
        }
    }}
}

pub fn compat_args(
    mut args: impl ExactSizeIterator<Item = OsString>,
) -> impl Iterator<Item = OsString> {
    let first = args.next();
    let second = args.next();
    let (command, first, second) = if first
        .as_ref()
        .map(|f| f.to_string_lossy().to_lowercase().ends_with("bigtools"))
        .unwrap_or(false)
    {
        if let Some(command) = second
            .as_ref()
            .and_then(|c| c.to_str())
            .map(|c| c.to_lowercase())
        {
            (Some(command), first, second.map(|a| a.to_ascii_lowercase()))
        } else {
            (None, first, second)
        }
    } else {
        (
            first
                .as_ref()
                .and_then(|f| f.to_str())
                .map(|f| f.to_lowercase()),
            first.map(|a| a.to_ascii_lowercase()),
            second,
        )
    };
    let mut args = {
        let mut args_vec = Vec::with_capacity(2 + args.len());
        first.map(|a| args_vec.push(a));
        second.map(|a| args_vec.push(a));
        args_vec.extend(args);
        args_vec
    };
    let iter: Box<dyn Iterator<Item = OsString>> = match command.as_deref() {
        Some("bigwigmerge") => {
            let has_input = args.iter().any(|a| {
                a.to_str()
                    .map_or(false, |a| a.starts_with("-b") || a.starts_with("-l"))
            });
            let mut args = if !has_input {
                // If there are no -l or -b, then let's see if it looks like a kent bigWigMerge call
                let in_list = args
                    .iter()
                    .any(|a| a.to_str().map_or(false, |a| a == "-inList"));
                let mut old_args = args;
                let mut args = Vec::with_capacity(old_args.len());
                old_args.reverse();
                while let Some(os_arg) = old_args.pop() {
                    let arg = os_arg.to_string_lossy();
                    if arg == "-inList" {
                        continue;
                    }
                    if arg.starts_with("-") && !arg.contains("=") {
                        args.push(os_arg);
                        args.pop().map(|a| args.push(a));
                        continue;
                    }
                    let more_args = 'more: {
                        let mut args_iter = args.iter().rev().peekable();
                        while let Some(arg) = args_iter.next() {
                            let arg = arg.to_string_lossy();
                            if arg.starts_with("-") && !arg.contains("=") {
                                args_iter.next();
                            } else {
                                break 'more true;
                            }
                        }
                        false
                    };
                    if !more_args {
                        if in_list {
                            args.push(OsString::from_str("-l").unwrap());
                        } else {
                            args.push(OsString::from_str("-b").unwrap());
                        }
                    }
                    args.push(os_arg);
                }
                args
            } else {
                args
            };

            args.iter_mut().for_each(|arg| {
                compat_replace_mut!(arg;
                    replace:
                        "-threshold", "--threshold";
                        "-adjust", "--adjust";
                        "-clip", "--clip"
                    ignore:
                        "-inList"
                    unimplemented:
                        "-max";
                        "-udcDir"
                )
            });

            Box::new(args.into_iter())
        }
        Some("bedgraphtobigwig") => {
            args.iter_mut().for_each(|arg| {
                compat_replace_mut!(arg;
                    replace:
                        "-unc", "--uncompressed";
                        "-blockSize", "--block-size";
                        "-itemsPerSlot", "--items-per-slot"
                    ignore:
                    unimplemented:
                )
            });

            Box::new(args.into_iter())
        }
        Some("bedtobigbed") => {
            args.iter_mut().for_each(|arg| {
                compat_replace_mut!(arg;
                    replace:
                        "-unc", "--uncompressed";
                        "-blockSize", "--block-size";
                        "-itemsPerSlot", "--items-per-slot";
                        "-as", "-autosql"
                    ignore:
                        "-tab";
                        "-type"
                    unimplemented:
                        "-extraIndex";
                        "-sizesIs2Bit";
                        "-sizesIsChromAliasBb";
                        "-sizesIsBb";
                        "-allow1bOverlap";
                        "-extraIndex";
                        "-udcDir"
                )
            });

            Box::new(args.into_iter())
        }
        Some("bigbedtobed") => {
            args.iter_mut().for_each(|arg| {
                compat_replace_mut!(arg;
                    replace:
                        "-chrom", "--chrom";
                        "-start", "--start";
                        "-end", "--end";
                        "-bed", "-overlap-bed"
                    ignore:
                    unimplemented:
                        "-header";
                        "-udcDir"
                )
            });

            Box::new(args.into_iter())
        }
        Some("bigwigaverageoverbed") => {
            args.iter_mut().for_each(|arg| {
                compat_replace_mut!(arg;
                    replace:
                        "-chrom", "--chrom";
                        "-start", "--start";
                        "-end", "--end"
                    ignore:
                    unimplemented:
                        "-udcDir"
                )
            });

            Box::new(args.into_iter())
        }
        Some("bigwigtobedgraph") => {
            args.iter_mut().for_each(|arg| {
                compat_replace_mut!(arg;
                    replace:
                        "-chrom", "--chrom";
                        "-start", "--start";
                        "-end", "--end"
                    ignore:
                    unimplemented:
                        "-udcDir"
                )
            });

            Box::new(args.into_iter())
        }
        _ => Box::new(args.into_iter()),
    };
    let iter: Vec<_> = iter.collect();
    iter.into_iter()
}
