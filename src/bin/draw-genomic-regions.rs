use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::genome::Genome;
use rustynetics::granges::GRanges;

mod common;

fn import_bed3(path: &str) -> GRanges {
    let mut reader = common::open_reader(Some(path)).unwrap_or_else(|error| {
        eprintln!("opening BED failed: {error}");
        process::exit(1);
    });
    let mut granges = GRanges::default();
    if let Err(error) = granges.read_bed3(&mut reader) {
        eprintln!("reading BED failed: {error}");
        process::exit(1);
    }
    granges
}

fn main() {
    let matches = Command::new("draw-genomic-regions")
        .about("Draw random genomic regions")
        .arg(
            Arg::new("exclude")
                .short('e')
                .long("exclude")
                .value_name("BED1,BED2,..."),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("genome").required(true).index(1))
        .arg(Arg::new("n").required(true).index(2))
        .arg(Arg::new("length").required(true).index(3))
        .arg(Arg::new("output").index(4))
        .get_matches();

    let genome_path = matches.get_one::<String>("genome").unwrap();
    let n: usize = matches
        .get_one::<String>("n")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid number of regions: {error}");
            process::exit(1);
        });
    let length: usize = matches
        .get_one::<String>("length")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid region length: {error}");
            process::exit(1);
        });
    let output_path = matches.get_one::<String>("output").map(String::as_str);
    let exclude = matches
        .get_one::<String>("exclude")
        .map(String::as_str)
        .unwrap_or("");
    let verbose = matches.get_count("verbose");

    let mut genome = Genome::default();
    if verbose > 0 {
        eprintln!("Reading genome `{}`...", genome_path);
    }
    if let Err(error) = genome.import(genome_path) {
        eprintln!("reading genome failed: {error}");
        process::exit(1);
    }

    let exclude_sets: Vec<GRanges> = if exclude.is_empty() {
        Vec::new()
    } else {
        exclude.split(',').map(import_bed3).collect()
    };

    let mut result = GRanges::default();
    while result.num_rows() < n {
        let mut batch = GRanges::random(n - result.num_rows(), length, &genome, false);
        for excluded in &exclude_sets {
            batch = batch.remove_overlaps_with(excluded);
        }
        result = result.append(&batch).unwrap_or_else(|error| {
            eprintln!("appending random regions failed: {error}");
            process::exit(1);
        });
    }

    let mut writer = common::open_writer(output_path).unwrap_or_else(|error| {
        eprintln!("opening output failed: {error}");
        process::exit(1);
    });
    if let Err(error) = result.write_bed3(&mut writer) {
        eprintln!("writing BED failed: {error}");
        process::exit(1);
    }
}
