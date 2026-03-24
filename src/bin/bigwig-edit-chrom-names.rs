use std::fs;
use std::path::PathBuf;
use std::process;

use clap::{Arg, ArgAction, Command};
use regex::Regex;

use rustynetics::genome::Genome;
use rustynetics::track::Track;
use rustynetics::track_simple::SimpleTrack;

mod common;

fn main() {
    let matches = Command::new("bigwig-edit-chrom-names")
        .about("Rewrite a BigWig with chromosome names transformed by a regex")
        .arg(Arg::new("bin-size").long("bin-size").default_value("0"))
        .arg(Arg::new("output").long("output").value_name("FILE"))
        .arg(
            Arg::new("dry-run")
                .long("dry-run")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("input").required(true).index(1))
        .arg(Arg::new("regex").required(true).index(2))
        .arg(Arg::new("replacement").required(true).index(3))
        .get_matches();

    let input = matches.get_one::<String>("input").unwrap();
    let output = matches.get_one::<String>("output").map(String::as_str);
    let bin_size: usize = matches
        .get_one::<String>("bin-size")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid bin size: {error}");
            process::exit(1);
        });
    let regex = Regex::new(matches.get_one::<String>("regex").unwrap()).unwrap_or_else(|error| {
        eprintln!("invalid regular expression: {error}");
        process::exit(1);
    });
    let replacement = matches.get_one::<String>("replacement").unwrap();
    let dry_run = matches.get_flag("dry-run");

    let track = common::import_simple_track(
        input,
        "",
        common::parse_bin_summary("mean").unwrap(),
        bin_size,
        0,
        f64::NAN,
    )
    .unwrap_or_else(|error| {
        eprintln!("importing BigWig failed: {error}");
        process::exit(1);
    });

    let old_seqnames = track.get_genome().seqnames.clone();
    let new_seqnames: Vec<_> = old_seqnames
        .iter()
        .map(|name| regex.replace_all(name, replacement.as_str()).to_string())
        .collect();

    if dry_run {
        for (old_name, new_name) in old_seqnames.iter().zip(new_seqnames.iter()) {
            println!("`{old_name}` -> `{new_name}`");
        }
        return;
    }

    let lengths = track.get_genome().lengths.clone();
    let sequences: Vec<_> = old_seqnames
        .iter()
        .map(|name| track.get_sequence(name).unwrap().clone_as_vec())
        .collect();
    let genome = Genome::new(new_seqnames, lengths);
    let renamed = SimpleTrack::new(String::new(), sequences, genome, track.get_bin_size())
        .unwrap_or_else(|error| {
            eprintln!("building renamed track failed: {error}");
            process::exit(1);
        });

    if let Some(output) = output {
        if let Err(error) = common::export_simple_track(&renamed, output) {
            eprintln!("writing output failed: {error}");
            process::exit(1);
        }
        return;
    }

    let mut tmp = PathBuf::from(input);
    tmp.set_extension("chromnames.tmp.bw");
    if let Err(error) = common::export_simple_track(&renamed, tmp.to_str().unwrap()) {
        eprintln!("writing temporary file failed: {error}");
        process::exit(1);
    }
    if let Err(error) = fs::rename(&tmp, input) {
        eprintln!("replacing input file failed: {error}");
        process::exit(1);
    }
}
