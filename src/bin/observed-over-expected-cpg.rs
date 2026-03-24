use std::collections::HashMap;
use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::granges::GRanges;
use rustynetics::meta::MetaData;
use rustynetics::orderedstringset::OrderedStringSet;

mod common;

fn main() {
    let matches = Command::new("observed-over-expected-cpg")
        .about("Compute observed/expected CpG scores for regions or whole sequences")
        .arg(
            Arg::new("output-format")
                .long("output-format")
                .value_parser(["table", "vector"])
                .default_value("table"),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("fasta").required(true).index(1))
        .arg(Arg::new("regions").index(2))
        .get_matches();

    let fasta_path = matches.get_one::<String>("fasta").unwrap();
    let regions_path = matches.get_one::<String>("regions").map(String::as_str);
    let output_format = matches.get_one::<String>("output-format").unwrap();
    let verbose = matches.get_count("verbose") > 0;

    let mut sequences = OrderedStringSet::empty();
    if verbose {
        eprintln!("Reading FASTA `{}`...", fasta_path);
    }
    if let Err(error) = sequences.import_fasta(fasta_path) {
        eprintln!("reading FASTA failed: {error}");
        process::exit(1);
    }

    let mut regions = if let Some(path) = regions_path {
        let mut reader = common::open_reader(Some(path)).unwrap_or_else(|error| {
            eprintln!("opening BED failed: {error}");
            process::exit(1);
        });
        let mut regions = GRanges::default();
        if let Err(error) = regions.read_bed3(&mut reader) {
            eprintln!("reading BED failed: {error}");
            process::exit(1);
        }
        regions
    } else {
        let from = vec![0; sequences.seqnames.len()];
        let to = sequences
            .seqnames
            .iter()
            .map(|name| sequences.sequences[name].len())
            .collect();
        GRanges::new(sequences.seqnames.clone(), from, to, Vec::new())
    };

    let sequence_map: HashMap<String, Vec<u8>> = sequences
        .seqnames
        .iter()
        .map(|name| (name.clone(), sequences.sequences[name].clone()))
        .collect();

    let scores = regions
        .observed_over_expected_cpg(&sequence_map)
        .unwrap_or_else(|error| {
            eprintln!("computing CpG scores failed: {error}");
            process::exit(1);
        });

    match output_format.as_str() {
        "vector" => {
            for value in scores {
                println!("{value}");
            }
        }
        "table" => {
            if let Err(error) = regions.meta.add("CpG", MetaData::FloatArray(scores)) {
                eprintln!("building output table failed: {error}");
                process::exit(1);
            }
            let mut writer = common::open_writer(None).unwrap_or_else(|error| {
                eprintln!("opening stdout failed: {error}");
                process::exit(1);
            });
            if let Err(error) = regions.write_table(&mut writer, &[]) {
                eprintln!("writing table failed: {error}");
                process::exit(1);
            }
        }
        _ => unreachable!(),
    }
}
