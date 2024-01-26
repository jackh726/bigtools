include!("bigtools.rs");

#[cfg(test)]
mod test {
    use bigtools::utils::cli::compat_args;
    use clap::Parser;

    use crate::{CliCommands, SubCommands};

    #[test]
    fn verify_cli_bigbedinfo() {
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
        let _applet = |args: &str| {
            let args = args.split_whitespace();
            let cli = CliCommands::try_parse_from(compat_args(args.map(|a| a.into()))).unwrap();
            match cli {
                CliCommands::Bigtools { command } => command,
                CliCommands::SubCommands(..) => panic!("Expected applet, parsed subcommand."),
            }
        };

        let args = "bigBedInfo a";
        let cli = subcommand(args);
        let _args = match cli {
            SubCommands::BigBedInfo { args } => {
                assert_eq!(args.bigbed, "a");

                args
            }
            _ => panic!(),
        };

        /*
        let args_orig = args;

        macro_rules! assert_args {
            (inner; $cli: expr, $args_comp:ident; $inner:block) => {
                let args_cli = match $cli {
                    SubCommands::BigBedInfo { args } => args,
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

        let args = "bigBedInfo a -udcDir";
        assert_args!(args, |args_comp| {});

        let args = "bigBedInfo a -chroms";
        assert_args!(args, |args_comp| {});

        let args = "bigBedInfo a -zooms";
        assert_args!(args, |args_comp| {});

        let args = "bigBedInfo a -minMax";
        assert_args!(args, |args_comp| {});
        */
    }
}
