use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::granges::GRanges;
use rustynetics::orderedstringset::OrderedStringSet;

mod common;

fn unresolved_regions(sequences: &OrderedStringSet) -> GRanges {
    let mut seqnames = Vec::new();
    let mut from = Vec::new();
    let mut to = Vec::new();

    for name in &sequences.seqnames {
        let sequence = &sequences.sequences[name];
        let mut i = 0usize;
        while i < sequence.len() {
            if sequence[i] == b'N' || sequence[i] == b'n' {
                let start = i;
                i += 1;
                while i < sequence.len() && (sequence[i] == b'N' || sequence[i] == b'n') {
                    i += 1;
                }
                seqnames.push(name.clone());
                from.push(start);
                to.push(i);
            } else {
                i += 1;
            }
        }
    }

    GRanges::new(seqnames, from, to, Vec::new())
}

fn main() {
    let matches = Command::new("fasta-unresolved-regions")
        .about("Report runs of unresolved bases in FASTA sequences")
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("fasta").required(true).index(1))
        .arg(Arg::new("output").index(2))
        .get_matches();

    let fasta_path = matches.get_one::<String>("fasta").unwrap();
    let output_path = matches.get_one::<String>("output").map(String::as_str);
    let verbose = matches.get_count("verbose") > 0;

    let mut sequences = OrderedStringSet::empty();
    if verbose {
        eprintln!("Reading FASTA `{}`...", fasta_path);
    }
    if let Err(error) = sequences.import_fasta(fasta_path) {
        eprintln!("reading FASTA failed: {error}");
        process::exit(1);
    }

    let regions = unresolved_regions(&sequences);
    let mut writer = common::open_writer(output_path).unwrap_or_else(|error| {
        eprintln!("opening output failed: {error}");
        process::exit(1);
    });
    if let Err(error) = regions.write_bed3(&mut writer) {
        eprintln!("writing BED failed: {error}");
        process::exit(1);
    }
}
