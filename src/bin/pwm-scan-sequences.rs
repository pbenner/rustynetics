use std::process;
use std::thread;

use clap::{Arg, ArgAction, Command};

use rustynetics::genome::Genome;
use rustynetics::stringset::StringSet;
use rustynetics::tf::{TFMatrix, PWM};
use rustynetics::track::MutableTrack;
use rustynetics::track_generic::GenericTrack;
use rustynetics::track_simple::SimpleTrack;

mod common;

fn genome_from_stringset(sequences: &StringSet) -> Genome {
    let mut seqnames: Vec<_> = sequences.keys().cloned().collect();
    seqnames.sort();
    let lengths = seqnames
        .iter()
        .map(|name| sequences[name.as_str()].len())
        .collect();
    Genome::new(seqnames, lengths)
}

fn summarize(values: &[f64], summary_name: &str) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    match summary_name {
        "mean" => values.iter().sum::<f64>() / values.len() as f64,
        "max" => values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        "min" => values.iter().copied().fold(f64::INFINITY, f64::min),
        "discrete mean" => (values.iter().sum::<f64>() / values.len() as f64).floor() + 0.5,
        _ => unreachable!(),
    }
}

fn scan_sequence_bins(pwm: &PWM, sequence: &[u8], summary_name: &str, bin_size: usize) -> Vec<f64> {
    let motif_len = pwm.matrix.length();
    if motif_len == 0 || sequence.len() < motif_len {
        return Vec::new();
    }

    let n_bins = sequence.len().div_ceil(bin_size);
    let last_start = sequence.len() - motif_len + 1;
    let mut result = Vec::with_capacity(n_bins);

    for i in 0..n_bins {
        let start = i * bin_size;
        let end = ((i + 1) * bin_size).min(last_start);
        if start >= end {
            break;
        }

        let mut scores = Vec::with_capacity((end - start) * 2);
        for j in start..end {
            scores.push(
                pwm.matrix
                    .score(&sequence[j..j + motif_len], false, 0.0, |a, b| a + b)
                    .unwrap(),
            );
            scores.push(
                pwm.matrix
                    .score(&sequence[j..j + motif_len], true, 0.0, |a, b| a + b)
                    .unwrap(),
            );
        }
        result.push(summarize(&scores, summary_name));
    }

    result
}

fn main() {
    let matches = Command::new("pwm-scan-sequences")
        .about("Scan FASTA sequences with a PWM and export a BigWig track")
        .arg(
            Arg::new("bin-summary")
                .long("bin-summary")
                .default_value("max")
                .value_parser(["mean", "max", "min", "discrete mean"]),
        )
        .arg(Arg::new("bin-size").long("bin-size").default_value("10"))
        .arg(Arg::new("threads").long("threads").default_value("1"))
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("pwm").required(true).index(1))
        .arg(Arg::new("fasta").required(true).index(2))
        .arg(Arg::new("output").required(true).index(3))
        .get_matches();

    let pwm_path = matches.get_one::<String>("pwm").unwrap();
    let fasta_path = matches.get_one::<String>("fasta").unwrap();
    let output_path = matches.get_one::<String>("output").unwrap();
    let summary_name = matches.get_one::<String>("bin-summary").unwrap();
    let bin_size: usize = matches
        .get_one::<String>("bin-size")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid bin size: {error}");
            process::exit(1);
        });
    let threads: usize = matches
        .get_one::<String>("threads")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid number of threads: {error}");
            process::exit(1);
        });
    let verbose = matches.get_count("verbose") > 0;
    if threads == 0 {
        eprintln!("invalid number of threads: must be at least 1");
        process::exit(1);
    }

    let mut matrix = TFMatrix::empty();
    if verbose {
        eprintln!("Reading PWM `{}`...", pwm_path);
    }
    if let Err(error) = matrix.import_matrix(pwm_path) {
        eprintln!("reading PWM failed: {error}");
        process::exit(1);
    }
    let pwm = PWM::new(matrix);

    let mut sequences = StringSet::empty();
    if verbose {
        eprintln!("Reading FASTA `{}`...", fasta_path);
    }
    if let Err(error) = sequences.import_fasta(fasta_path) {
        eprintln!("reading FASTA failed: {error}");
        process::exit(1);
    }

    let genome = genome_from_stringset(&sequences);
    let mut track = SimpleTrack::alloc(String::new(), genome.clone(), f64::NAN, bin_size);
    let ranges = common::worker_ranges(genome.seqnames.len(), threads);
    let results = if ranges.len() <= 1 {
        genome
            .seqnames
            .iter()
            .map(|seqname| {
                if verbose {
                    eprintln!("Scanning sequence `{seqname}`...");
                }
                (
                    seqname.clone(),
                    scan_sequence_bins(&pwm, &sequences[seqname.as_str()], summary_name, bin_size),
                )
            })
            .collect::<Vec<_>>()
    } else {
        thread::scope(|scope| {
            let mut handles = Vec::new();
            for (start, end) in ranges {
                let seqnames = &genome.seqnames[start..end];
                let sequences = &sequences;
                let pwm = pwm.clone();
                let summary_name = summary_name.to_string();
                handles.push(scope.spawn(move || {
                    seqnames
                        .iter()
                        .map(|seqname| {
                            (
                                seqname.clone(),
                                scan_sequence_bins(
                                    &pwm,
                                    &sequences[seqname.as_str()],
                                    &summary_name,
                                    bin_size,
                                ),
                            )
                        })
                        .collect::<Vec<_>>()
                }));
            }

            let mut merged = Vec::with_capacity(genome.seqnames.len());
            for handle in handles {
                merged.extend(handle.join().unwrap_or_else(|_| {
                    eprintln!("pwm-scan-sequences worker thread panicked");
                    process::exit(1);
                }));
            }
            merged
        })
    };

    for (seqname, values) in results {
        let mut target = track.get_sequence_mut(&seqname).unwrap_or_else(|error| {
            eprintln!("accessing track sequence failed: {error}");
            process::exit(1);
        });
        for (i, value) in values.into_iter().enumerate() {
            target.set_bin(i, value);
        }
    }

    if let Err(error) = GenericTrack::wrap(&track).export_bigwig(output_path, vec![]) {
        eprintln!("writing BigWig failed: {error}");
        process::exit(1);
    }
}
