use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::bigwig::BigWigFile;
use rustynetics::granges::GRanges;
use rustynetics::granges_table::{OptionPrintScientific, OptionPrintStrand};
use rustynetics::track::MutableTrack;
use rustynetics::track_simple::SimpleTrack;

mod common;

fn load_regions(path: &str) -> GRanges {
    let mut reader = common::open_reader(Some(path)).unwrap_or_else(|error| {
        eprintln!("opening BED failed: {error}");
        process::exit(1);
    });
    let mut granges = GRanges::default();
    granges.read_bed3(&mut reader).unwrap_or_else(|error| {
        eprintln!("reading BED failed: {error}");
        process::exit(1);
    });
    granges
}

fn write_table_output(
    mut regions: GRanges,
    bigwig_path: &str,
    summary: rustynetics::track_statistics::BinSummaryStatistics,
    bin_size: usize,
    bin_overlap: usize,
    init: f64,
    scientific: bool,
    output_path: Option<&str>,
) {
    regions
        .import_bigwig(
            bigwig_path,
            "counts",
            summary,
            bin_size,
            bin_overlap,
            init,
            false,
        )
        .unwrap_or_else(|error| {
            eprintln!("extracting BigWig data failed: {error}");
            process::exit(1);
        });

    let mut writer = common::open_writer(output_path).unwrap_or_else(|error| {
        eprintln!("opening output failed: {error}");
        process::exit(1);
    });
    if let Err(error) = regions.write_table(
        &mut writer,
        &[&OptionPrintStrand(true), &OptionPrintScientific(scientific)],
    ) {
        eprintln!("writing table failed: {error}");
        process::exit(1);
    }
}

fn write_bigwig_output(
    regions: GRanges,
    bigwig_path: &str,
    summary: rustynetics::track_statistics::BinSummaryStatistics,
    mut bin_size: usize,
    bin_overlap: usize,
    init: f64,
    output_path: &str,
) {
    let mut reader = BigWigFile::new_reader(bigwig_path).unwrap_or_else(|error| {
        eprintln!("opening BigWig failed: {error}");
        process::exit(1);
    });
    let genome = reader.genome().clone();
    let mut track = SimpleTrack::alloc(String::new(), genome, f64::NAN, 1.max(bin_size));

    for i in 0..regions.num_rows() {
        let (values, detected_bin_size) = reader
            .query_slice(
                &regions.seqnames[i],
                regions.ranges[i].from,
                regions.ranges[i].to,
                summary,
                bin_size,
                bin_overlap,
                init,
            )
            .unwrap_or_else(|error| {
                eprintln!("querying BigWig failed: {error}");
                process::exit(1);
            });

        if bin_size == 0 {
            bin_size = detected_bin_size;
            track = SimpleTrack::alloc(String::new(), reader.genome().clone(), f64::NAN, bin_size);
        }

        let mut seq = track
            .get_sequence_mut(&regions.seqnames[i])
            .unwrap_or_else(|error| {
                eprintln!("accessing track sequence failed: {error}");
                process::exit(1);
            });
        for (j, value) in values.into_iter().enumerate() {
            seq.set(regions.ranges[i].from + j * bin_size, value);
        }
    }

    if let Err(error) = common::export_simple_track(&track, output_path) {
        eprintln!("writing BigWig failed: {error}");
        process::exit(1);
    }
}

fn main() {
    let matches = Command::new("bigwig-extract")
        .about("Extract BigWig data for BED regions as a table or BigWig")
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
            Arg::new("bin-overlap")
                .long("bin-overlap")
                .default_value("0"),
        )
        .arg(Arg::new("initial-value").long("initial-value"))
        .arg(
            Arg::new("scientific")
                .long("scientific")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("bigwig").required(true).index(1))
        .arg(Arg::new("bed").required(true).index(2))
        .arg(Arg::new("output").index(3))
        .get_matches();

    let bigwig = matches.get_one::<String>("bigwig").unwrap();
    let bed = matches.get_one::<String>("bed").unwrap();
    let output = matches.get_one::<String>("output").map(String::as_str);
    let bin_size: usize = matches
        .get_one::<String>("bin-size")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid bin size: {error}");
            process::exit(1);
        });
    let bin_overlap: usize = matches
        .get_one::<String>("bin-overlap")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid bin overlap: {error}");
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
        0.0,
    )
    .unwrap_or_else(|error| {
        eprintln!("invalid initial value: {error}");
        process::exit(1);
    });
    let scientific = matches.get_flag("scientific");

    let regions = load_regions(bed);
    match output {
        Some(path) if path.ends_with(".bw") || path.ends_with(".bigWig") => {
            write_bigwig_output(regions, bigwig, summary, bin_size, bin_overlap, init, path);
        }
        _ => write_table_output(
            regions,
            bigwig,
            summary,
            bin_size,
            bin_overlap,
            init,
            scientific,
            output,
        ),
    }
}
