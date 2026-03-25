use std::process;

use clap::{Arg, ArgAction, Command};

mod common;

fn main() {
    let matches = Command::new("bigwig-statistics")
        .about("Print summary statistics for a BigWig track")
        .arg(Arg::new("bin-size").long("bin-size").default_value("0"))
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("input").required(true).index(1))
        .get_matches();

    let input = matches.get_one::<String>("input").unwrap();
    let bin_size: usize = matches
        .get_one::<String>("bin-size")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid bin size: {error}");
            process::exit(1);
        });

    let track = common::import_lazy_track(
        input,
        "",
        common::parse_bin_summary("mean").unwrap(),
        bin_size,
        0,
        f64::NAN,
    )
    .unwrap_or_else(|error| {
        eprintln!("importing BigWig failed: {error}");
        process::exit(1);
    });
    let summary = common::summarize_track(&track).unwrap_or_else(|error| {
        eprintln!("computing statistics failed: {error}");
        process::exit(1);
    });
    println!("{summary}");
}
