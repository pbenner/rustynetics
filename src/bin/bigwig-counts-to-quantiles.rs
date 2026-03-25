use std::cmp::Ordering;
use std::collections::HashMap;
use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::track_generic::{GenericMutableTrack, GenericTrack};

mod common;

fn main() {
    let matches = Command::new("bigwig-counts-to-quantiles")
        .about("Convert BigWig counts to empirical quantiles")
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
                .action(ArgAction::Count),
        )
        .arg(Arg::new("input").required(true).index(1))
        .arg(Arg::new("output").required(true).index(2))
        .get_matches();

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
    let verbose = matches.get_count("verbose") > 0;

    if verbose {
        eprintln!("Importing BigWig `{input}`...");
    }
    let mut track = common::import_simple_track(input, "", summary, bin_size, 0, init)
        .unwrap_or_else(|error| {
            eprintln!("importing BigWig failed: {error}");
            process::exit(1);
        });

    let mut counts = HashMap::new();
    GenericTrack::wrap(&track)
        .map(|_, _, value| {
            if !value.is_nan() {
                *counts.entry(value.to_bits()).or_insert(0usize) += 1;
            }
        })
        .unwrap_or_else(|error| {
            eprintln!("counting track values failed: {error}");
            process::exit(1);
        });

    let mut dist: Vec<(f64, usize)> = counts
        .into_iter()
        .map(|(bits, count)| (f64::from_bits(bits), count))
        .collect();
    dist.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
    let mut quantiles = HashMap::new();
    let total: usize = dist.iter().map(|(_, count)| *count).sum();
    if total > 0 {
        let mut rank = 0usize;
        for (value, count) in dist {
            rank += count;
            quantiles.insert(value.to_bits(), rank as f64 / total as f64);
        }
    }

    if verbose {
        eprintln!("Converting counts to quantiles...");
    }
    GenericMutableTrack::wrap(&mut track)
        .map(|_, _, value| {
            if value.is_nan() {
                value
            } else {
                *quantiles.get(&value.to_bits()).unwrap_or(&value)
            }
        })
        .unwrap_or_else(|error| {
            eprintln!("mapping quantiles failed: {error}");
            process::exit(1);
        });

    if verbose {
        eprintln!("Writing BigWig `{output}`...");
    }
    if let Err(error) = common::export_simple_track(&track, output) {
        eprintln!("writing output failed: {error}");
        process::exit(1);
    }
}
