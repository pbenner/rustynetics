use std::error::Error;
use std::io;
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
use rustynetics::meta::MetaData;
use rustynetics::orderedstringset::OrderedStringSet;
use rustynetics::range::Range;

mod common;

struct Config {
    binary: bool,
    max_ambiguous: Option<Vec<usize>>,
    complement: bool,
    reverse: bool,
    revcomp: bool,
    human: bool,
    sparse: bool,
    header: bool,
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

fn count_sequences_with_alphabet<A: ComplementableAlphabet + Clone + Send + Sync>(
    alphabet: A,
    n: usize,
    m: usize,
    config: &Config,
    extracted: &[Vec<u8>],
) -> Result<Vec<rustynetics::kmer_counts::KmerCounts>, Box<dyn Error>> {
    let ranges = common::worker_ranges(extracted.len(), config.threads);
    if ranges.len() <= 1 {
        let mut counter = KmerCounter::new(
            n,
            m,
            config.complement,
            config.reverse,
            config.revcomp,
            config.max_ambiguous.clone(),
            alphabet,
        )?;
        return Ok(extracted
            .iter()
            .map(|sequence| {
                if config.binary {
                    counter.identify_kmers(sequence)
                } else {
                    counter.count_kmers(sequence)
                }
            })
            .collect());
    }

    thread::scope(
        |scope| -> Result<Vec<rustynetics::kmer_counts::KmerCounts>, Box<dyn Error>> {
            let mut handles = Vec::new();

            for (start, end) in ranges {
                let alphabet = alphabet.clone();
                let sequences = &extracted[start..end];
                let max_ambiguous = config.max_ambiguous.clone();
                let binary = config.binary;
                let complement = config.complement;
                let reverse = config.reverse;
                let revcomp = config.revcomp;

                handles.push(scope.spawn(
                    move || -> Result<Vec<rustynetics::kmer_counts::KmerCounts>, String> {
                        let mut counter = KmerCounter::new(
                            n,
                            m,
                            complement,
                            reverse,
                            revcomp,
                            max_ambiguous,
                            alphabet,
                        )?;
                        Ok(sequences
                            .iter()
                            .map(|sequence| {
                                if binary {
                                    counter.identify_kmers(sequence)
                                } else {
                                    counter.count_kmers(sequence)
                                }
                            })
                            .collect())
                    },
                ));
            }

            let mut counts = Vec::with_capacity(extracted.len());
            for handle in handles {
                let local = handle
                    .join()
                    .map_err(|_| {
                        io::Error::new(io::ErrorKind::Other, "count-kmers worker thread panicked")
                    })?
                    .map_err(|error| io::Error::new(io::ErrorKind::Other, error))?;
                counts.extend(local);
            }
            Ok(counts)
        },
    )
}

fn count_with_alphabet<A: ComplementableAlphabet + Clone + Send + Sync>(
    alphabet: A,
    n: usize,
    m: usize,
    config: &Config,
    regions_path: Option<&str>,
    fasta_path: Option<&str>,
    output_path: Option<&str>,
) -> Result<(), Box<dyn Error>> {
    let regions = import_regions(regions_path)?;
    let sequences = import_fasta(fasta_path)?;
    let (mut granges, extracted) = collect_sequences(&sequences, regions.as_ref())?;
    let counts = count_sequences_with_alphabet(alphabet, n, m, config, &extracted)?;

    let counts_list = rustynetics::kmer_counts::KmerCountsList::new(counts);

    if config.human || config.sparse {
        let values: Vec<String> = (0..counts_list.len())
            .map(|i| {
                let counts = counts_list.at(i);
                let mut parts = Vec::new();
                let mut it = counts.iter();
                while it.ok() {
                    let count = it.get_count();
                    if config.sparse && count == 0 {
                        it.next();
                        continue;
                    }
                    parts.push(format!("{}={}", it.get_kmer(), count));
                    it.next();
                }
                parts.join(",")
            })
            .collect();
        granges.meta.add("k-mers", MetaData::StringArray(values))?;
    } else {
        let values: Vec<Vec<i64>> = (0..counts_list.len())
            .map(|i| {
                let counts = counts_list.at(i);
                let mut row = Vec::with_capacity(counts.len());
                let mut it = counts.iter();
                while it.ok() {
                    row.push(it.get_count() as i64);
                    it.next();
                }
                row
            })
            .collect();
        granges.meta.add("k-mers", MetaData::IntMatrix(values))?;
    }

    let mut writer = common::open_writer(output_path)?;
    if config.header {
        for kmer in &counts_list.kmers {
            writeln!(writer, "# {}", kmer)?;
        }
    }
    granges.write_table(&mut writer, &[])?;
    Ok(())
}

fn main() {
    let matches = Command::new("count-kmers")
        .about("Count k-mers in FASTA sequences or BED regions")
        .arg(
            Arg::new("alphabet")
                .long("alphabet")
                .default_value("nucleotide")
                .value_parser(["nucleotide", "gapped-nucleotide", "ambiguous-nucleotide"]),
        )
        .arg(Arg::new("binary").long("binary").action(ArgAction::SetTrue))
        .arg(
            Arg::new("max-ambiguous")
                .long("max-ambiguous")
                .default_value("-1"),
        )
        .arg(Arg::new("regions").long("regions").value_name("BED"))
        .arg(Arg::new("header").long("header").action(ArgAction::SetTrue))
        .arg(Arg::new("human").long("human").action(ArgAction::SetTrue))
        .arg(Arg::new("sparse").long("sparse").action(ArgAction::SetTrue))
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
        .arg(Arg::new("fasta").index(3))
        .arg(Arg::new("output").index(4))
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
        human: matches.get_flag("human"),
        sparse: matches.get_flag("sparse"),
        header: matches.get_flag("header"),
        threads: matches
            .get_one::<String>("threads")
            .unwrap()
            .parse()
            .unwrap_or_else(|error| {
                eprintln!("invalid number of threads: {error}");
                process::exit(1);
            }),
        verbose: matches.get_count("verbose"),
    };
    let regions_path = matches.get_one::<String>("regions").map(String::as_str);
    let fasta_path = matches.get_one::<String>("fasta").map(String::as_str);
    let output_path = matches.get_one::<String>("output").map(String::as_str);

    if config.threads == 0 {
        eprintln!("invalid number of threads: must be at least 1");
        process::exit(1);
    }
    if config.verbose > 0 {
        if let Some(path) = fasta_path {
            eprintln!("Reading FASTA `{path}`...");
        } else {
            eprintln!("Reading FASTA from stdin...");
        }
    }

    let result = match matches.get_one::<String>("alphabet").unwrap().as_str() {
        "nucleotide" => count_with_alphabet(
            NucleotideAlphabet,
            n,
            m,
            &config,
            regions_path,
            fasta_path,
            output_path,
        ),
        "gapped-nucleotide" => count_with_alphabet(
            GappedNucleotideAlphabet,
            n,
            m,
            &config,
            regions_path,
            fasta_path,
            output_path,
        ),
        "ambiguous-nucleotide" => count_with_alphabet(
            AmbiguousNucleotideAlphabet,
            n,
            m,
            &config,
            regions_path,
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
