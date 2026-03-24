use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::genome::Genome;
use rustynetics::track::Track;
use rustynetics::track_simple::SimpleTrack;

mod common;

fn main() {
    let matches = Command::new("bigwig-extract-chroms")
        .about("Extract selected chromosomes from a BigWig file")
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("chroms").required(true).index(1))
        .arg(Arg::new("input").required(true).index(2))
        .arg(Arg::new("output").required(true).index(3))
        .get_matches();

    let chroms: Vec<_> = matches
        .get_one::<String>("chroms")
        .unwrap()
        .split(',')
        .filter(|name| !name.is_empty())
        .map(str::to_string)
        .collect();
    let input = matches.get_one::<String>("input").unwrap();
    let output = matches.get_one::<String>("output").unwrap();

    let track = common::import_simple_track(
        input,
        "",
        common::parse_bin_summary("mean").unwrap(),
        0,
        0,
        f64::NAN,
    )
    .unwrap_or_else(|error| {
        eprintln!("importing BigWig failed: {error}");
        process::exit(1);
    });

    let mut sequences = Vec::with_capacity(chroms.len());
    let mut lengths = Vec::with_capacity(chroms.len());
    for chrom in &chroms {
        let idx = track.get_genome().get_idx(chrom).unwrap_or_else(|| {
            eprintln!("chromosome `{chrom}` not found in input BigWig");
            process::exit(1);
        });
        sequences.push(track.get_sequence(chrom).unwrap().clone_as_vec());
        lengths.push(track.get_genome().lengths[idx]);
    }
    let genome = Genome::new(chroms, lengths);
    let subset = SimpleTrack::new(String::new(), sequences, genome, track.get_bin_size())
        .unwrap_or_else(|error| {
            eprintln!("building subset track failed: {error}");
            process::exit(1);
        });

    if let Err(error) = common::export_simple_track(&subset, output) {
        eprintln!("writing BigWig failed: {error}");
        process::exit(1);
    }
}
