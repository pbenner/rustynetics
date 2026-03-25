use std::io::{self, IsTerminal, Write};
use std::process;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering as AtomicOrdering};
use std::sync::Arc;
use std::thread;
use std::time::Duration;

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

fn spawn_status_reporter(
    total: usize,
    enabled: bool,
) -> Option<(Arc<AtomicUsize>, Arc<AtomicBool>, thread::JoinHandle<()>)> {
    if !enabled || !io::stderr().is_terminal() {
        return None;
    }

    let progress = Arc::new(AtomicUsize::new(0));
    let done = Arc::new(AtomicBool::new(false));
    let progress_clone = Arc::clone(&progress);
    let done_clone = Arc::clone(&done);

    let handle = thread::spawn(move || {
        while !done_clone.load(AtomicOrdering::Relaxed) {
            let current = progress_clone.load(AtomicOrdering::Relaxed);
            let percent = if total == 0 {
                100.0
            } else {
                current as f64 * 100.0 / total as f64
            };
            eprint!("\r\x1b[2Kregions {:>6.2}% {}/{}", percent, current, total);
            let _ = io::stderr().flush();
            thread::sleep(Duration::from_millis(100));
        }
        let current = progress_clone.load(AtomicOrdering::Relaxed);
        let percent = if total == 0 {
            100.0
        } else {
            current as f64 * 100.0 / total as f64
        };
        eprintln!("\r\x1b[2Kregions {:>6.2}% {}/{}", percent, current, total);
    });

    Some((progress, done, handle))
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
        .arg(Arg::new("status").long("status").action(ArgAction::SetTrue))
        .arg(Arg::new("threads").long("threads").default_value("1"))
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
    let status = matches.get_flag("status");
    let threads: usize = matches
        .get_one::<String>("threads")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid number of threads: {error}");
            process::exit(1);
        });
    let pwm_paths: Vec<_> = matches
        .get_many::<String>("pwm")
        .unwrap()
        .map(String::as_str)
        .collect();
    let verbose = matches.get_count("verbose") > 0;
    if threads == 0 {
        eprintln!("invalid number of threads: must be at least 1");
        process::exit(1);
    }

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
    let reporter = spawn_status_reporter(granges.num_rows(), status);
    let progress = reporter
        .as_ref()
        .map(|(progress, _, _)| Arc::clone(progress));
    let ranges = common::worker_ranges(granges.num_rows(), threads);
    let counts = if ranges.len() <= 1 {
        (0..granges.num_rows())
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
                let row: Vec<f64> = pwms
                    .iter()
                    .map(|pwm| scan_region(sequence, pwm, summary))
                    .collect();
                if let Some(progress) = &progress {
                    progress.fetch_add(1, AtomicOrdering::Relaxed);
                }
                row
            })
            .collect::<Vec<_>>()
    } else {
        thread::scope(|scope| {
            let mut handles = Vec::new();
            for (start, end) in ranges {
                let granges = &granges;
                let sequences = &sequences;
                let pwms = pwms.clone();
                let summary = summary.to_string();
                let progress = progress.clone();
                handles.push(scope.spawn(move || {
                    let mut rows = Vec::with_capacity(end - start);
                    for i in start..end {
                        let sequence = sequences
                            .get_slice(
                                &granges.seqnames[i],
                                Range::new(granges.ranges[i].from, granges.ranges[i].to),
                            )
                            .unwrap_or_else(|error| {
                                eprintln!("extracting region sequence failed: {error}");
                                process::exit(1);
                            });
                        rows.push(
                            pwms.iter()
                                .map(|pwm| scan_region(sequence, pwm, &summary))
                                .collect::<Vec<_>>(),
                        );
                        if let Some(progress) = &progress {
                            progress.fetch_add(1, AtomicOrdering::Relaxed);
                        }
                    }
                    rows
                }));
            }

            let mut merged = Vec::with_capacity(granges.num_rows());
            for handle in handles {
                merged.extend(handle.join().unwrap_or_else(|_| {
                    eprintln!("pwm-scan-regions worker thread panicked");
                    process::exit(1);
                }));
            }
            merged
        })
    };

    if let Some((progress, done, handle)) = reporter {
        progress.store(granges.num_rows(), AtomicOrdering::Relaxed);
        done.store(true, AtomicOrdering::Relaxed);
        let _ = handle.join();
    }

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
