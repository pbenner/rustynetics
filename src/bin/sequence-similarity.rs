use std::error::Error;
use std::io::Write;
use std::process;
use std::thread;

use clap::{Arg, ArgAction, Command};

use rustynetics::alphabet::{
    AmbiguousNucleotideAlphabet, ComplementableAlphabet, GappedNucleotideAlphabet,
    NucleotideAlphabet,
};
use rustynetics::granges::GRanges;
use rustynetics::kmer_counter::KmerCounter;
use rustynetics::kmer_counts::KmerCounts;
use rustynetics::meta::MetaData;
use rustynetics::orderedstringset::OrderedStringSet;
use rustynetics::range::Range;

mod common;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum OutputFormat {
    Granges,
    Table,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Measure {
    Cosine,
    Tanimoto,
    DotProduct,
}

struct Config {
    binary: bool,
    max_ambiguous: Option<Vec<usize>>,
    complement: bool,
    reverse: bool,
    revcomp: bool,
    format: OutputFormat,
    measure: Measure,
    step_size: isize,
    threads: usize,
    verbose: u8,
}

fn parse_max_ambiguous(value: &str) -> Result<Option<Vec<usize>>, String> {
    if value.trim() == "-1" {
        return Ok(None);
    }
    let values: Result<Vec<_>, _> = value.split(',').map(|field| field.trim().parse()).collect();
    let values = values.map_err(|error| format!("invalid max-ambiguous value: {error}"))?;
    Ok(Some(values))
}

fn import_regions(path: Option<&str>) -> Result<Option<GRanges>, Box<dyn Error>> {
    let Some(path) = path else {
        return Ok(None);
    };
    let mut reader = common::open_reader(Some(path))?;
    let mut granges = GRanges::default();
    granges.read_bed3(&mut reader)?;
    Ok(Some(granges))
}

fn import_fasta(path: Option<&str>) -> Result<OrderedStringSet, Box<dyn Error>> {
    let mut sequences = OrderedStringSet::empty();
    let mut reader = common::open_reader(path)?;
    sequences.read_fasta(&mut reader)?;
    Ok(sequences)
}

fn collect_sequences(
    sequences: &OrderedStringSet,
    regions: Option<&GRanges>,
) -> Result<(GRanges, Vec<Vec<u8>>), Box<dyn Error>> {
    if let Some(regions) = regions {
        let mut extracted = Vec::with_capacity(regions.num_rows());
        for i in 0..regions.num_rows() {
            let slice = sequences.get_slice(
                &regions.seqnames[i],
                Range::new(regions.ranges[i].from, regions.ranges[i].to),
            )?;
            extracted.push(slice.to_vec());
        }
        Ok((regions.clone(), extracted))
    } else {
        let from = vec![0; sequences.seqnames.len()];
        let to = sequences
            .seqnames
            .iter()
            .map(|name| sequences.sequences[name].len())
            .collect();
        let granges = GRanges::new(sequences.seqnames.clone(), from, to, Vec::new());
        let data = sequences
            .seqnames
            .iter()
            .map(|name| sequences.sequences[name].clone())
            .collect();
        Ok((granges, data))
    }
}

fn compute_similarity(counts_seq: &KmerCounts, counts_ref: &KmerCounts, measure: Measure) -> f64 {
    match measure {
        Measure::Cosine => {
            let mut dot = 0.0;
            let mut a = 0.0;
            let mut b = 0.0;
            for i in 0..counts_ref.len() {
                let x = counts_ref.at(i) as f64;
                let y = counts_seq.get_count(counts_ref.get_kmer(i)) as f64;
                a += x * x;
                dot += x * y;
            }
            for i in 0..counts_seq.len() {
                let y = counts_seq.at(i) as f64;
                b += y * y;
            }
            dot / a.sqrt() / b.sqrt()
        }
        Measure::Tanimoto => {
            let mut dot = 0.0;
            let mut a = 0.0;
            let mut b = 0.0;
            for i in 0..counts_ref.len() {
                let x = counts_ref.at(i) as f64;
                let y = counts_seq.get_count(counts_ref.get_kmer(i)) as f64;
                a += x * x;
                dot += x * y;
            }
            for i in 0..counts_seq.len() {
                let y = counts_seq.at(i) as f64;
                b += y * y;
            }
            dot / (a + b - dot)
        }
        Measure::DotProduct => {
            let mut dot = 0.0;
            for i in 0..counts_seq.len() {
                let x = counts_seq.at(i) as f64;
                let y = counts_ref.get_count(counts_seq.get_kmer(i)) as f64;
                dot += x * y;
            }
            dot / counts_ref.len() as f64
        }
    }
}

fn compute_step_size(step_size: isize, reference_len: usize) -> usize {
    if step_size == -1 {
        (reference_len / 10).max(1)
    } else {
        step_size as usize
    }
}

fn scan_sequence<A: ComplementableAlphabet + Clone>(
    counter: &mut KmerCounter<A>,
    binary: bool,
    sequence: &[u8],
) -> KmerCounts {
    if binary {
        counter.identify_kmers(sequence)
    } else {
        counter.count_kmers(sequence)
    }
}

fn compute_similarities_with_alphabet<A: ComplementableAlphabet + Clone + Send + Sync>(
    alphabet: A,
    n: usize,
    m: usize,
    config: &Config,
    regions_path: Option<&str>,
    reference_path: &str,
    fasta_path: Option<&str>,
    output_path: Option<&str>,
) -> Result<(), Box<dyn Error>> {
    let mut counter = KmerCounter::new(
        n,
        m,
        config.complement,
        config.reverse,
        config.revcomp,
        config.max_ambiguous.clone(),
        alphabet,
    )?;

    let reference_set = import_fasta(Some(reference_path))?;
    let (_, reference_sequences) = collect_sequences(&reference_set, None)?;
    if reference_sequences.len() != 1 {
        return Err("reference FASTA file should contain a single sequence".into());
    }
    let reference_sequence = reference_sequences.into_iter().next().unwrap();
    let counts_ref = scan_sequence(&mut counter, config.binary, &reference_sequence);
    if config.measure == Measure::DotProduct {
        counter.freeze();
    }

    let regions = import_regions(regions_path)?;
    let sequences = import_fasta(fasta_path)?;
    let (mut granges, extracted) = collect_sequences(&sequences, regions.as_ref())?;
    let step_size = compute_step_size(config.step_size, reference_sequence.len());
    let reference_len = reference_sequence.len();
    let ranges = common::worker_ranges(extracted.len(), config.threads);

    let similarities = if ranges.len() <= 1 {
        let mut worker_counter = counter;
        extracted
            .iter()
            .map(|sequence| {
                let delta = sequence.len().saturating_sub(reference_len);
                if delta < step_size {
                    return Vec::new();
                }
                let mut values = Vec::with_capacity(delta / step_size);
                for start in (0..=delta - step_size).step_by(step_size) {
                    let counts = scan_sequence(
                        &mut worker_counter,
                        config.binary,
                        &sequence[start..start + reference_len],
                    );
                    values.push(compute_similarity(&counts, &counts_ref, config.measure));
                }
                values
            })
            .collect::<Vec<_>>()
    } else {
        thread::scope(|scope| {
            let mut handles = Vec::new();
            for (start, end) in ranges {
                let sequences = &extracted[start..end];
                let mut counter = counter.clone();
                let counts_ref = counts_ref.clone();
                let binary = config.binary;
                let measure = config.measure;
                handles.push(scope.spawn(move || {
                    let mut rows = Vec::with_capacity(sequences.len());
                    for sequence in sequences {
                        let delta = sequence.len().saturating_sub(reference_len);
                        if delta < step_size {
                            rows.push(Vec::new());
                            continue;
                        }
                        let mut values = Vec::with_capacity(delta / step_size);
                        for start in (0..=delta - step_size).step_by(step_size) {
                            let counts = scan_sequence(
                                &mut counter,
                                binary,
                                &sequence[start..start + reference_len],
                            );
                            values.push(compute_similarity(&counts, &counts_ref, measure));
                        }
                        rows.push(values);
                    }
                    rows
                }));
            }

            let mut merged = Vec::with_capacity(extracted.len());
            for handle in handles {
                merged.extend(handle.join().unwrap_or_else(|_| {
                    eprintln!("sequence-similarity worker thread panicked");
                    process::exit(1);
                }));
            }
            merged
        })
    };

    let mut writer = common::open_writer(output_path)?;
    match config.format {
        OutputFormat::Granges => {
            granges
                .meta
                .add("similarities", MetaData::FloatMatrix(similarities))?;
            granges.write_table(&mut writer, &[])?;
        }
        OutputFormat::Table => {
            for row in similarities {
                for value in row {
                    write!(writer, "{value:8.4} ")?;
                }
                writeln!(writer)?;
            }
        }
    }
    Ok(())
}

fn main() {
    let matches = Command::new("sequence-similarity")
        .about("Compute sliding-window k-mer similarity against a reference sequence")
        .arg(
            Arg::new("alphabet")
                .long("alphabet")
                .default_value("nucleotide")
                .value_parser(["nucleotide", "gapped-nucleotide", "ambiguous-nucleotide"]),
        )
        .arg(
            Arg::new("format")
                .long("format")
                .default_value("granges")
                .value_parser(["granges", "table"]),
        )
        .arg(
            Arg::new("measure")
                .long("measure")
                .default_value("cosine")
                .value_parser(["cosine", "tanimoto", "dot-product"]),
        )
        .arg(Arg::new("step-size").long("step-size").default_value("-1"))
        .arg(Arg::new("binary").long("binary").action(ArgAction::SetTrue))
        .arg(
            Arg::new("max-ambiguous")
                .long("max-ambiguous")
                .default_value("-1"),
        )
        .arg(Arg::new("regions").long("regions").value_name("BED"))
        .arg(Arg::new("threads").long("threads").default_value("1"))
        .arg(
            Arg::new("complement")
                .long("complement")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("reverse")
                .long("reverse")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("revcomp")
                .long("revcomp")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("n").required(true).index(1))
        .arg(Arg::new("m").required(true).index(2))
        .arg(Arg::new("reference").required(true).index(3))
        .arg(Arg::new("fasta").index(4))
        .arg(Arg::new("output").index(5))
        .get_matches();

    let n: usize = matches
        .get_one::<String>("n")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid minimum k-mer length: {error}");
            process::exit(1);
        });
    let m: usize = matches
        .get_one::<String>("m")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid maximum k-mer length: {error}");
            process::exit(1);
        });
    if n == 0 || m < n {
        eprintln!("invalid k-mer size range");
        process::exit(1);
    }

    let step_size: isize = matches
        .get_one::<String>("step-size")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid step size: {error}");
            process::exit(1);
        });
    if step_size < -1 || step_size == 0 {
        eprintln!("invalid step size");
        process::exit(1);
    }

    let threads: usize = matches
        .get_one::<String>("threads")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid number of threads: {error}");
            process::exit(1);
        });
    if threads == 0 {
        eprintln!("invalid number of threads: must be at least 1");
        process::exit(1);
    }

    let config = Config {
        binary: matches.get_flag("binary"),
        max_ambiguous: parse_max_ambiguous(matches.get_one::<String>("max-ambiguous").unwrap())
            .unwrap_or_else(|error| {
                eprintln!("{error}");
                process::exit(1);
            }),
        complement: matches.get_flag("complement"),
        reverse: matches.get_flag("reverse"),
        revcomp: matches.get_flag("revcomp"),
        format: match matches.get_one::<String>("format").unwrap().as_str() {
            "granges" => OutputFormat::Granges,
            "table" => OutputFormat::Table,
            _ => unreachable!(),
        },
        measure: match matches.get_one::<String>("measure").unwrap().as_str() {
            "cosine" => Measure::Cosine,
            "tanimoto" => Measure::Tanimoto,
            "dot-product" => Measure::DotProduct,
            _ => unreachable!(),
        },
        step_size,
        threads,
        verbose: matches.get_count("verbose"),
    };

    let reference_path = matches.get_one::<String>("reference").unwrap();
    let fasta_path = matches.get_one::<String>("fasta").map(String::as_str);
    let output_path = matches.get_one::<String>("output").map(String::as_str);
    let regions_path = matches.get_one::<String>("regions").map(String::as_str);

    if config.verbose > 0 {
        eprintln!("Reading reference FASTA `{reference_path}`...");
        if let Some(path) = fasta_path {
            eprintln!("Reading target FASTA `{path}`...");
        } else {
            eprintln!("Reading target FASTA from stdin...");
        }
    }

    let result = match matches.get_one::<String>("alphabet").unwrap().as_str() {
        "nucleotide" => compute_similarities_with_alphabet(
            NucleotideAlphabet,
            n,
            m,
            &config,
            regions_path,
            reference_path,
            fasta_path,
            output_path,
        ),
        "gapped-nucleotide" => compute_similarities_with_alphabet(
            GappedNucleotideAlphabet,
            n,
            m,
            &config,
            regions_path,
            reference_path,
            fasta_path,
            output_path,
        ),
        "ambiguous-nucleotide" => compute_similarities_with_alphabet(
            AmbiguousNucleotideAlphabet,
            n,
            m,
            &config,
            regions_path,
            reference_path,
            fasta_path,
            output_path,
        ),
        _ => unreachable!(),
    };

    if let Err(error) = result {
        eprintln!("{error}");
        process::exit(1);
    }
}
