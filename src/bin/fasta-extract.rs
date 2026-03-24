use std::error::Error;
use std::io::Write;
use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::granges::GRanges;
use rustynetics::meta::MetaData;
use rustynetics::orderedstringset::OrderedStringSet;
use rustynetics::range::Range;

mod common;

fn extract_regions(
    sequences: &OrderedStringSet,
    regions: &GRanges,
) -> Result<Vec<Vec<u8>>, Box<dyn Error>> {
    let mut result = Vec::with_capacity(regions.num_rows());

    for i in 0..regions.num_rows() {
        let sequence = sequences.get_slice(
            &regions.seqnames[i],
            Range::new(regions.ranges[i].from, regions.ranges[i].to),
        )?;
        result.push(sequence.to_vec());
    }

    Ok(result)
}

fn main() {
    let matches = Command::new("fasta-extract")
        .about("Extract FASTA subsequences for BED regions")
        .arg(
            Arg::new("output-format")
                .long("output-format")
                .value_parser(["fasta", "table"])
                .default_value("fasta"),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("fasta").required(true).index(1))
        .arg(Arg::new("input").index(2))
        .arg(Arg::new("output").index(3))
        .get_matches();

    let fasta_path = matches.get_one::<String>("fasta").unwrap();
    let input_path = matches.get_one::<String>("input").map(String::as_str);
    let output_path = matches.get_one::<String>("output").map(String::as_str);
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

    let regions = if let Some(input_path) = input_path {
        let mut regions = GRanges::default();
        let mut reader = common::open_reader(Some(input_path)).unwrap_or_else(|error| {
            eprintln!("opening BED failed: {error}");
            process::exit(1);
        });
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

    let extracted = extract_regions(&sequences, &regions).unwrap_or_else(|error| {
        eprintln!("extracting sequences failed: {error}");
        process::exit(1);
    });

    let mut writer = common::open_writer(output_path).unwrap_or_else(|error| {
        eprintln!("opening output failed: {error}");
        process::exit(1);
    });

    match output_format.as_str() {
        "fasta" => {
            for (i, sequence) in extracted.iter().enumerate() {
                let name = format!(
                    "{}|{}|{}",
                    regions.seqnames[i], regions.ranges[i].from, regions.ranges[i].to
                );
                if let Err(error) = common::write_fasta_record(&mut writer, &name, sequence) {
                    eprintln!("writing FASTA failed: {error}");
                    process::exit(1);
                }
            }
        }
        "table" => {
            let mut table = regions.clone();
            let column = extracted
                .into_iter()
                .map(|sequence| String::from_utf8_lossy(&sequence).to_string())
                .collect();
            if let Err(error) = table.meta.add("sequences", MetaData::StringArray(column)) {
                eprintln!("building output table failed: {error}");
                process::exit(1);
            }
            if let Err(error) = table.write_table(&mut writer, &[]) {
                eprintln!("writing table failed: {error}");
                process::exit(1);
            }
        }
        _ => unreachable!(),
    }

    if let Err(error) = writer.flush() {
        eprintln!("flushing output failed: {error}");
        process::exit(1);
    }
}
