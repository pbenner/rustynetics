use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::track_generic::GenericMutableTrack;

mod common;

fn main() {
    let matches = Command::new("bigwig-quantile-normalize")
        .about("Quantile-normalize one BigWig track to a reference track")
        .arg(Arg::new("bin-size").long("bin-size").default_value("0"))
        .arg(
            Arg::new("bin-summary")
                .long("bin-summary")
                .default_value("mean")
                .value_parser([
                    "mean",
                    "max",
                    "min",
                    "discrete mean",
                    "discrete max",
                    "discrete min",
                    "variance",
                ]),
        )
        .arg(Arg::new("initial-value").long("initial-value"))
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("reference").required(true).index(1))
        .arg(Arg::new("input").required(true).index(2))
        .arg(Arg::new("output").required(true).index(3))
        .get_matches();

    let reference = matches.get_one::<String>("reference").unwrap();
    let input = matches.get_one::<String>("input").unwrap();
    let output = matches.get_one::<String>("output").unwrap();
    let bin_size: usize = matches
        .get_one::<String>("bin-size")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid bin size: {error}");
            process::exit(1);
        });
    let summary = common::parse_bin_summary(matches.get_one::<String>("bin-summary").unwrap())
        .unwrap_or_else(|error| {
            eprintln!("{error}");
            process::exit(1);
        });
    let init = common::parse_initial_value(
        matches
            .get_one::<String>("initial-value")
            .map(String::as_str),
        f64::NAN,
    )
    .unwrap_or_else(|error| {
        eprintln!("invalid initial value: {error}");
        process::exit(1);
    });

    let track_ref = common::import_simple_track(reference, "", summary, bin_size, 0, init)
        .unwrap_or_else(|error| {
            eprintln!("importing reference failed: {error}");
            process::exit(1);
        });
    let mut track_in = common::import_simple_track(input, "", summary, bin_size, 0, init)
        .unwrap_or_else(|error| {
            eprintln!("importing input failed: {error}");
            process::exit(1);
        });

    if let Err(error) = GenericMutableTrack::wrap(&mut track_in).quantile_normalize(&track_ref) {
        eprintln!("quantile normalization failed: {error}");
        process::exit(1);
    }
    if let Err(error) = common::export_simple_track(&track_in, output) {
        eprintln!("writing output failed: {error}");
        process::exit(1);
    }
}
