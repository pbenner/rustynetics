use std::collections::BTreeMap;
use std::process;

use clap::{Arg, ArgAction, Command};

mod common;

fn main() {
    let matches = Command::new("bigwig-histogram")
        .about("Print a value histogram for a BigWig track")
        .arg(
            Arg::new("bins")
                .short('b')
                .long("bins")
                .default_value("100"),
        )
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
        .arg(
            Arg::new("cumulative")
                .short('c')
                .long("cumulative")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("input").required(true).index(1))
        .get_matches();

    let input = matches.get_one::<String>("input").unwrap();
    let bins: usize = matches
        .get_one::<String>("bins")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid histogram bin count: {error}");
            process::exit(1);
        });
    let bin_size: usize = matches
        .get_one::<String>("bin-size")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid track bin size: {error}");
            process::exit(1);
        });
    let summary = common::parse_bin_summary(matches.get_one::<String>("bin-summary").unwrap())
        .unwrap_or_else(|error| {
            eprintln!("{error}");
            process::exit(1);
        });
    let cumulative = matches.get_flag("cumulative");

    let track = common::import_simple_track(input, "", summary, bin_size, 0, f64::NAN)
        .unwrap_or_else(|error| {
            eprintln!("importing BigWig failed: {error}");
            process::exit(1);
        });
    let values = common::collect_track_values(&track).unwrap_or_else(|error| {
        eprintln!("collecting values failed: {error}");
        process::exit(1);
    });

    println!("{:>15}\t{:>15}", "x", "y");
    if values.is_empty() {
        return;
    }

    if bins == 0 {
        let mut counts: BTreeMap<u64, usize> = BTreeMap::new();
        for value in &values {
            *counts.entry(value.to_bits()).or_insert(0) += 1;
        }
        let mut running = 0usize;
        for (bits, count) in counts {
            running += count;
            let y = if cumulative {
                running as f64
            } else {
                count as f64
            };
            println!("{:15e}\t{:15.6}", f64::from_bits(bits), y);
        }
        return;
    }

    let min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let width = if max == min {
        1.0
    } else {
        (max - min) / bins as f64
    };
    let mut counts = vec![0usize; bins];
    for value in &values {
        let mut idx = if max == min {
            0
        } else {
            ((value - min) / width).floor() as usize
        };
        if idx >= bins {
            idx = bins - 1;
        }
        counts[idx] += 1;
    }

    let mut running = 0usize;
    for (i, count) in counts.into_iter().enumerate() {
        running += count;
        let x = min + i as f64 * width;
        let y = if cumulative {
            running as f64
        } else {
            count as f64
        };
        println!("{:15e}\t{:15.6}", x, y);
    }
}
