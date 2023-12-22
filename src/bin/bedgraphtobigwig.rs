use bigtools::utils::cli::{
    bedgraphtobigwig::{bedgraphtobigwig, BedGraphToBigWigArgs},
    compat_args,
};
use clap::Parser;
use std::error::Error;

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    args: BedGraphToBigWigArgs,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Cli::parse_from(compat_args(std::env::args_os()));

    bedgraphtobigwig(matches.args)
}

#[test]
fn verify_cli_bedgraphtobigwig() {
    use clap::CommandFactory;
    Cli::command().debug_assert();

    let args = ["bedgraphtobigwig", "a", "b", "c"];
    let cli = Cli::try_parse_from(compat_args(args.into_iter().map(|a| a.into()))).unwrap();
    assert_eq!(cli.args.bedgraph, "a");
    assert_eq!(cli.args.chromsizes, "b");
    assert_eq!(cli.args.output, "c");
    assert_eq!(cli.args.write_args.uncompressed, false);

    let args = ["bedgraphtobigwig", "a", "b", "c", "-unc"];
    let cli = Cli::try_parse_from(compat_args(args.into_iter().map(|a| a.into()))).unwrap();
    assert_eq!(cli.args.bedgraph, "a");
    assert_eq!(cli.args.chromsizes, "b");
    assert_eq!(cli.args.output, "c");
    assert_eq!(cli.args.write_args.uncompressed, true);

    let args = ["bedgraphtobigwig", "-blockSize", "50", "a", "b", "c"];
    let cli = Cli::try_parse_from(compat_args(args.into_iter().map(|a| a.into()))).unwrap();
    assert_eq!(cli.args.bedgraph, "a");
    assert_eq!(cli.args.chromsizes, "b");
    assert_eq!(cli.args.output, "c");
    assert_eq!(cli.args.write_args.block_size, 50);

    let args = ["bedgraphtobigwig", "a", "-itemsPerSlot", "5", "b", "c"];
    let cli = Cli::try_parse_from(compat_args(args.into_iter().map(|a| a.into()))).unwrap();
    assert_eq!(cli.args.bedgraph, "a");
    assert_eq!(cli.args.chromsizes, "b");
    assert_eq!(cli.args.output, "c");
    assert_eq!(cli.args.write_args.items_per_slot, 5);
}
