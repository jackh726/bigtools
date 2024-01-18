include!("bigtools.rs");

#[cfg(test)]
mod test {
    use bigtools::utils::cli::compat_args;
    use clap::Parser;

    use crate::{CliCommands, SubCommands};

    #[test]
    fn verify_cli_bigwigtobedgraph() {
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

        let args = "bigWigToBedGraph in out";
        let cli = subcommand(args);
        let args = match cli {
            SubCommands::BigWigToBedGraph { args } => {
                assert_eq!(args.bigwig, "in");
                assert_eq!(args.bedgraph, "out");

                args
            }
            _ => panic!("Unexpected matched subcommand."),
        };

        let args_orig = args;

        macro_rules! assert_args {
            (inner; $cli: expr, $args_comp:ident; $inner:block) => {
                let args_cli = match $cli {
                    SubCommands::BigWigToBedGraph { args } => args,
                    _ => panic!("Unexpected matched subcommand."),
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

        let args = "bigWigToBedGraph in out -chrom=chr1";
        assert_args!(args, |args_comp| {
            args_comp.chrom = Some("chr1".to_string());
        });

        let args = "bigWigToBedGraph in out -start=0";
        assert_args!(args, |args_comp| {
            args_comp.start = Some(0);
        });

        let args = "bigWigToBedGraph in out -end=1000";
        assert_args!(args, |args_comp| {
            args_comp.end = Some(1000);
        });
    }
}
