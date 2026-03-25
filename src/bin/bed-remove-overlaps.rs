use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::granges::GRanges;
use rustynetics::granges_table::OptionPrintStrand;

mod common;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum InputFormat {
    Bed3,
    Bed6,
    Bed9,
    Table,
}

fn read_table_auto(path: Option<&str>) -> Result<GRanges, Box<dyn std::error::Error>> {
    let mut reader = common::open_reader(path)?;
    let mut granges = GRanges::default();
    granges.read_table_all(&mut reader)?;
    Ok(granges)
}

fn import_input(
    path: Option<&str>,
    format: InputFormat,
) -> Result<GRanges, Box<dyn std::error::Error>> {
    match format {
        InputFormat::Bed3 => {
            let mut reader = common::open_reader(path)?;
            let mut granges = GRanges::default();
            granges.read_bed3(&mut reader)?;
            Ok(granges)
        }
        InputFormat::Bed6 => {
            let mut reader = common::open_reader(path)?;
            let mut granges = GRanges::default();
            granges.read_bed6(&mut reader)?;
            Ok(granges)
        }
        InputFormat::Bed9 => {
            let mut reader = common::open_reader(path)?;
            let mut granges = GRanges::default();
            granges.read_bed9(&mut reader)?;
            Ok(granges)
        }
        InputFormat::Table => read_table_auto(path),
    }
}

fn export_output(
    granges: &GRanges,
    path: Option<&str>,
    format: InputFormat,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = common::open_writer(path)?;
    match format {
        InputFormat::Bed3 => granges.write_bed3(&mut writer)?,
        InputFormat::Bed6 => granges.write_bed6(&mut writer)?,
        InputFormat::Bed9 => granges.write_bed9(&mut writer)?,
        InputFormat::Table => granges.write_table(&mut writer, &[&OptionPrintStrand(true)])?,
    }
    Ok(())
}

fn main() {
    let matches = Command::new("bed-remove-overlaps")
        .about("Remove regions that overlap a set of inadmissible intervals")
        .arg(
            Arg::new("input-format")
                .long("input-format")
                .default_value("bed3")
                .value_parser(["bed3", "bed6", "bed9", "table"]),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("inadmissible").required(true).index(1))
        .arg(Arg::new("input").index(2))
        .arg(Arg::new("output").index(3))
        .get_matches();

    let format = match matches.get_one::<String>("input-format").unwrap().as_str() {
        "bed3" => InputFormat::Bed3,
        "bed6" => InputFormat::Bed6,
        "bed9" => InputFormat::Bed9,
        "table" => InputFormat::Table,
        _ => unreachable!(),
    };
    let inadmissible = matches.get_one::<String>("inadmissible").unwrap();
    let input = matches.get_one::<String>("input").map(String::as_str);
    let output = matches.get_one::<String>("output").map(String::as_str);
    let verbose = matches.get_count("verbose") > 0;

    if verbose {
        eprintln!("Reading inadmissible regions `{inadmissible}`...");
    }
    let remove = import_input(Some(inadmissible), InputFormat::Bed3).unwrap_or_else(|error| {
        eprintln!("reading inadmissible regions failed: {error}");
        process::exit(1);
    });
    if verbose {
        if let Some(path) = input {
            eprintln!("Reading input `{path}`...");
        } else {
            eprintln!("Reading input from stdin...");
        }
    }
    let input_ranges = import_input(input, format).unwrap_or_else(|error| {
        eprintln!("reading input failed: {error}");
        process::exit(1);
    });
    let before = input_ranges.num_rows();
    let filtered = input_ranges.remove_overlaps_with(&remove);
    if verbose {
        eprintln!("Removed {} regions", before - filtered.num_rows());
    }
    if let Err(error) = export_output(&filtered, output, format) {
        eprintln!("writing output failed: {error}");
        process::exit(1);
    }
}
