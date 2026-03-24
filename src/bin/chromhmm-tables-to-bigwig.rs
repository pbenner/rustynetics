use std::io::BufRead;
use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::genome::Genome;
use rustynetics::track::MutableTrack;
use rustynetics::track_simple::SimpleTrack;

mod common;

fn import_table(
    reader: &mut dyn BufRead,
    track: &mut SimpleTrack,
    column_name: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut line = String::new();
    if reader.read_line(&mut line)? == 0 {
        return Err("invalid table: missing sequence information".into());
    }
    let fields: Vec<_> = line.split_whitespace().collect();
    if fields.len() != 2 {
        return Err("invalid table: first line must have two columns".into());
    }
    let seqname = fields[1];
    let mut sequence = track.get_sequence_mut(seqname)?;

    line.clear();
    if reader.read_line(&mut line)? == 0 {
        return Err("invalid table: missing header".into());
    }
    let header: Vec<_> = line.split_whitespace().collect();
    let column_idx = header
        .iter()
        .position(|name| *name == column_name)
        .ok_or_else(|| format!("column `{column_name}` not found"))?;

    for i in 0..sequence.n_bins() {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let fields: Vec<_> = line.split_whitespace().collect();
        if column_idx >= fields.len() {
            return Err(format!("invalid table row at bin {i}").into());
        }
        sequence.set_bin(i, fields[column_idx].parse()?);
    }

    Ok(())
}

fn main() {
    let matches = Command::new("chromhmm-tables-to-bigwig")
        .about("Convert ChromHMM per-chromosome tables to a BigWig track")
        .arg(Arg::new("initial-value").long("initial-value"))
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("genome").required(true).index(1))
        .arg(Arg::new("bin-size").required(true).index(2))
        .arg(Arg::new("column").required(true).index(3))
        .arg(Arg::new("output").required(true).index(4))
        .arg(Arg::new("tables").required(true).index(5).num_args(1..))
        .get_matches();

    let genome_path = matches.get_one::<String>("genome").unwrap();
    let bin_size: usize = matches
        .get_one::<String>("bin-size")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid bin size: {error}");
            process::exit(1);
        });
    let column = matches.get_one::<String>("column").unwrap();
    let output = matches.get_one::<String>("output").unwrap();
    let tables: Vec<_> = matches
        .get_many::<String>("tables")
        .unwrap()
        .map(String::as_str)
        .collect();
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

    let mut genome = Genome::default();
    genome.import(genome_path).unwrap_or_else(|error| {
        eprintln!("reading genome failed: {error}");
        process::exit(1);
    });

    let mut track = SimpleTrack::alloc(String::new(), genome, init, bin_size);
    for table in tables {
        let mut reader = common::open_reader(Some(table)).unwrap_or_else(|error| {
            eprintln!("opening table `{table}` failed: {error}");
            process::exit(1);
        });
        import_table(&mut reader, &mut track, column).unwrap_or_else(|error| {
            eprintln!("importing table `{table}` failed: {error}");
            process::exit(1);
        });
    }

    if let Err(error) = common::export_simple_track(&track, output) {
        eprintln!("writing BigWig failed: {error}");
        process::exit(1);
    }
}
