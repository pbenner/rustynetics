use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::granges::GRanges;
use rustynetics::granges_table::OptionPrintStrand;
use rustynetics::meta::MetaData;
use rustynetics::orderedstringset::OrderedStringSet;
use rustynetics::range::Range;
use rustynetics::tf::{TFMatrix, PWM};

mod common;

fn log_add(mut a: f64, mut b: f64) -> f64 {
    if a > b {
        std::mem::swap(&mut a, &mut b);
    }
    if a.is_infinite() && a.is_sign_negative() {
        return b;
    }
    b + (a - b).exp().ln_1p()
}

fn import_pwm(path: &str) -> PWM {
    let mut matrix = TFMatrix::empty();
    matrix.import_matrix(path).unwrap_or_else(|error| {
        eprintln!("reading PWM `{path}` failed: {error}");
        process::exit(1);
    });
    PWM::new(matrix)
}

fn import_regions(path: Option<&str>, columns: usize, sequences: &OrderedStringSet) -> GRanges {
    if let Some(path) = path {
        let mut reader = common::open_reader(Some(path)).unwrap_or_else(|error| {
            eprintln!("opening regions failed: {error}");
            process::exit(1);
        });
        let mut granges = GRanges::default();
        granges
            .read_bed(&mut reader, columns)
            .unwrap_or_else(|error| {
                eprintln!("reading BED failed: {error}");
                process::exit(1);
            });
        granges
    } else {
        let from = vec![0; sequences.seqnames.len()];
        let to = sequences
            .seqnames
            .iter()
            .map(|name| sequences.sequences[name].len())
            .collect();
        GRanges::new(sequences.seqnames.clone(), from, to, Vec::new())
    }
}

fn scan_region(sequence: &[u8], pwm: &PWM, summary: &str) -> f64 {
    match summary {
        "max" => pwm
            .max_score(sequence, false)
            .max(pwm.max_score(sequence, true)),
        "mean" => {
            log_add(
                pwm.mean_score(sequence, false),
                pwm.mean_score(sequence, true),
            ) - (2.0f64).ln()
        }
        _ => unreachable!(),
    }
}

fn main() {
    let matches = Command::new("pwm-scan-regions")
        .about("Scan genomic regions with one or more PWMs")
        .arg(Arg::new("input").long("input").value_name("BED"))
        .arg(
            Arg::new("input-columns")
                .long("input-columns")
                .default_value("3")
                .value_parser(["3", "6", "9"]),
        )
        .arg(Arg::new("output").long("output").value_name("FILE"))
        .arg(
            Arg::new("summary")
                .long("summary")
                .default_value("max")
                .value_parser(["max", "mean"]),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("fasta").required(true).index(1))
        .arg(Arg::new("pwm").required(true).index(2).num_args(1..))
        .get_matches();

    let fasta_path = matches.get_one::<String>("fasta").unwrap();
    let region_path = matches.get_one::<String>("input").map(String::as_str);
    let output_path = matches.get_one::<String>("output").map(String::as_str);
    let columns: usize = matches
        .get_one::<String>("input-columns")
        .unwrap()
        .parse()
        .unwrap();
    let summary = matches.get_one::<String>("summary").unwrap();
    let pwm_paths: Vec<_> = matches
        .get_many::<String>("pwm")
        .unwrap()
        .map(String::as_str)
        .collect();
    let verbose = matches.get_count("verbose") > 0;

    let mut sequences = OrderedStringSet::empty();
    if verbose {
        eprintln!("Reading FASTA `{}`...", fasta_path);
    }
    sequences.import_fasta(fasta_path).unwrap_or_else(|error| {
        eprintln!("reading FASTA failed: {error}");
        process::exit(1);
    });

    let pwms: Vec<_> = pwm_paths.iter().map(|path| import_pwm(path)).collect();
    let mut granges = import_regions(region_path, columns, &sequences);

    let counts: Vec<Vec<f64>> = (0..granges.num_rows())
        .map(|i| {
            let sequence = sequences
                .get_slice(
                    &granges.seqnames[i],
                    Range::new(granges.ranges[i].from, granges.ranges[i].to),
                )
                .unwrap_or_else(|error| {
                    eprintln!("extracting region sequence failed: {error}");
                    process::exit(1);
                });
            pwms.iter()
                .map(|pwm| scan_region(sequence, pwm, summary))
                .collect()
        })
        .collect();

    granges
        .meta
        .add("counts", MetaData::FloatMatrix(counts))
        .unwrap_or_else(|error| {
            eprintln!("adding output column failed: {error}");
            process::exit(1);
        });

    let mut writer = common::open_writer(output_path).unwrap_or_else(|error| {
        eprintln!("opening output failed: {error}");
        process::exit(1);
    });
    if let Err(error) = granges.write_table(&mut writer, &[&OptionPrintStrand(true)]) {
        eprintln!("writing table failed: {error}");
        process::exit(1);
    }
}
