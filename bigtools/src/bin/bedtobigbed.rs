include!("bigtools.rs");

#[cfg(test)]
mod test {
    use bigtools::utils::cli::compat_args;
    use clap::Parser;

    use crate::{CliCommands, SubCommands};

    #[test]
    fn verify_cli_bedtobigbed() {
        use clap::CommandFactory;
        CliCommands::command().debug_assert();

        let subcommand = |args: &str| {
            let args = args.split_whitespace();
            let cli = CliCommands::try_parse_from(compat_args(args.map(|a| a.into())))
                .map_err(|e| e.print())
                .unwrap();
            match cli {
                CliCommands::SubCommands(subcommand) => subcommand,
                CliCommands::Bigtools { .. } => panic!("Expected subcommand, parsed applet."),
            }
        };
        let applet = |args: &str| {
            let args = args.split_whitespace();
            let cli = CliCommands::try_parse_from(compat_args(args.map(|a| a.into()))).unwrap();
            match cli {
                CliCommands::Bigtools { command } => command,
                CliCommands::SubCommands(..) => panic!("Expected applet, parsed subcommand."),
            }
        };

        let args = "bedToBigBed a b c";
        let cli = subcommand(args);
        let args = match cli {
            SubCommands::BedToBigBed { args } => {
                assert_eq!(args.bed, "a");
                assert_eq!(args.chromsizes, "b");
                assert_eq!(args.output, "c");

                args
            }
            _ => panic!(),
        };

        let args_orig = args;

        macro_rules! assert_args {
            (inner; $cli: expr, $args_comp:ident; $inner:block) => {
                let args_cli = match $cli {
                    SubCommands::BedToBigBed { args } => args,
                    _ => panic!(),
                };
                #[allow(unused_mut)]
                let mut $args_comp = args_orig.clone();
                $inner
                assert_eq!(args_cli, $args_comp);
            };
            ($args:expr, |$args_comp:ident| $inner:block) => {
                let cli = subcommand($args);
                assert_args!(inner; cli, $args_comp; $inner);

                let args = &format!("bigtools {}", $args);
                let cli = applet(args);
                assert_args!(inner; cli, $args_comp; $inner);
            }
        }

        let args = "bedToBigBed a b c -unc";
        assert_args!(args, |args_comp| {
            args_comp.write_args.uncompressed = true;
        });

        let args = "bedToBigBed -blockSize 50 a b c";
        assert_args!(args, |args_comp| {
            args_comp.write_args.block_size = 50;
        });

        let args = "bedToBigBed a -itemsPerSlot 5 b c";
        assert_args!(args, |args_comp| {
            args_comp.write_args.items_per_slot = 5;
        });

        /*
        let args = "bedToBigBed a b c -type=bed3+5";
        assert_args!(args, |args_comp| {});

        let args = "bedToBigBed a b c -tab";
        assert_args!(args, |args_comp| {});

        let args = "bedToBigBed a b c -extraIndex=a,b,c";
        assert_args!(args, |args_comp| {});

        let args = "bedToBigBed a b c -sizesIs2Bit";
        assert_args!(args, |args_comp| {});

        let args = "bedToBigBed a b c -sizesIsChromAliasBb";
        assert_args!(args, |args_comp| {});

        let args = "bedToBigBed a b c -sizesIsBb";
        assert_args!(args, |args_comp| {});

        let args = "bedToBigBed a b c -udcDir";
        assert_args!(args, |args_comp| {
            args_comp.write_args.items_per_slot = 5;
        });

        let args = "bedToBigBed a b c -allow1bpOverlap";
        assert_args!(args, |args_comp| {});
        */
    }
}
