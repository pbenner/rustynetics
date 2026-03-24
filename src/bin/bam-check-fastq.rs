// Copyright (C) 2024 Philipp Benner
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the “Software”), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::process;

use clap::{Arg, Command};
use flate2::read::MultiGzDecoder;

use rustynetics::bam::{BamReader, BamReaderOptions};
use rustynetics::progress::{CountingReader, ProgressScope};

/* -------------------------------------------------------------------------- */

struct Config {
    filename_bam: String,
    filename_fastq: String,
    keep_pair_suffixes: bool,
    report_limit: usize,
}

/* -------------------------------------------------------------------------- */

struct CheckResult {
    fastq_records: usize,
    fastq_unique_names: usize,
    bam_records: usize,
    missing_records: usize,
    missing_unique: usize,
    missing_examples: Vec<String>,
}

/* -------------------------------------------------------------------------- */

const STATUS_RECORD_UPDATE_INTERVAL: usize = 10_000;

/* -------------------------------------------------------------------------- */

fn trim_line(line: &str) -> &str {
    line.trim_end_matches(&['\r', '\n'][..])
}

/* -------------------------------------------------------------------------- */

fn normalize_read_name(name: &str, keep_pair_suffixes: bool) -> String {
    let token = name.split_whitespace().next().unwrap_or("");
    let token = token.strip_prefix('@').unwrap_or(token);

    if keep_pair_suffixes {
        token.to_string()
    } else {
        token
            .strip_suffix("/1")
            .or_else(|| token.strip_suffix("/2"))
            .unwrap_or(token)
            .to_string()
    }
}

/* -------------------------------------------------------------------------- */

fn open_fastq_reader(
    filename: &str,
    progress: Option<rustynetics::progress::ProgressHandle>,
) -> Result<Box<dyn BufRead>, Box<dyn Error>> {
    let file = File::open(filename)?;
    let reader = CountingReader::new(file, progress);
    let mut reader = BufReader::new(reader);
    let magic = reader.fill_buf()?;

    if magic.len() >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
}

/* -------------------------------------------------------------------------- */

fn read_fastq_names(
    filename: &str,
    keep_pair_suffixes: bool,
) -> Result<(HashSet<String>, usize), Box<dyn Error>> {
    let total_bytes = File::open(filename)?.metadata()?.len();
    let mut progress = ProgressScope::new("FASTQ", total_bytes);
    let mut reader = open_fastq_reader(filename, progress.handle())?;
    let mut names = HashSet::new();
    let mut line = String::new();
    let mut line_number = 0usize;
    let mut records = 0usize;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        line_number += 1;

        let header = trim_line(&line);
        if !header.starts_with('@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid FASTQ header at line {}", line_number),
            )
            .into());
        }

        let name = normalize_read_name(header, keep_pair_suffixes);
        if name.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("empty FASTQ read name at line {}", line_number),
            )
            .into());
        }
        names.insert(name);
        records += 1;
        if records % STATUS_RECORD_UPDATE_INTERVAL == 0 {
            progress.set_records(records);
        }

        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("unexpected end of FASTQ after header in record {}", records),
            )
            .into());
        }
        line_number += 1;

        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "unexpected end of FASTQ after sequence in record {}",
                    records
                ),
            )
            .into());
        }
        line_number += 1;
        if !trim_line(&line).starts_with('+') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid FASTQ separator at line {}", line_number),
            )
            .into());
        }

        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "unexpected end of FASTQ after separator in record {}",
                    records
                ),
            )
            .into());
        }
        line_number += 1;
    }

    progress.finish(records);

    Ok((names, records))
}

/* -------------------------------------------------------------------------- */

fn check_bam_fastq(config: &Config) -> Result<CheckResult, Box<dyn Error>> {
    let (fastq_names, fastq_records) =
        read_fastq_names(&config.filename_fastq, config.keep_pair_suffixes)?;

    let options = BamReaderOptions {
        read_name: true,
        read_cigar: false,
        read_sequence: false,
        read_auxiliary: false,
        read_qual: false,
    };

    let total_bytes = File::open(&config.filename_bam)?.metadata()?.len();
    let mut progress = ProgressScope::new("BAM", total_bytes);
    let file = File::open(&config.filename_bam)?;
    let reader = BufReader::new(CountingReader::new(file, progress.handle()));
    let mut bam_reader = BamReader::new(reader, Some(options))?;

    let mut bam_records = 0usize;
    let mut missing_records = 0usize;
    let mut missing_names = HashSet::new();
    let mut missing_examples = Vec::new();

    for result in bam_reader.read_single_end() {
        let block = result?.block;
        let name = normalize_read_name(&block.read_name, config.keep_pair_suffixes);

        bam_records += 1;
        if bam_records % STATUS_RECORD_UPDATE_INTERVAL == 0 {
            progress.set_records(bam_records);
        }

        if !fastq_names.contains(&name) {
            missing_records += 1;

            if missing_names.insert(name.clone()) && missing_examples.len() < config.report_limit {
                missing_examples.push(name);
            }
        }
    }

    progress.finish(bam_records);

    Ok(CheckResult {
        fastq_records,
        fastq_unique_names: fastq_names.len(),
        bam_records,
        missing_records,
        missing_unique: missing_names.len(),
        missing_examples,
    })
}

/* -------------------------------------------------------------------------- */

fn main() {
    let matches = Command::new("bam-check-fastq")
        .version("1.0")
        .author("Philipp Benner [https://github.com/pbenner]")
        .about("Checks whether all BAM read names are present in a FASTQ file")
        .arg(
            Arg::new("bam")
                .help("The input BAM file")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("fastq")
                .help("The FASTQ file to compare against (plain text or gzip compressed)")
                .required(true)
                .index(2),
        )
        .arg(
            Arg::new("keep-pair-suffixes")
                .long("keep-pair-suffixes")
                .action(clap::ArgAction::SetTrue)
                .help("Do not strip trailing /1 or /2 suffixes before comparing read names"),
        )
        .arg(
            Arg::new("report-limit")
                .short('n')
                .long("report-limit")
                .value_parser(clap::value_parser!(usize))
                .default_value("10")
                .help("Maximum number of missing read names to print"),
        )
        .get_matches();

    let config = Config {
        filename_bam: matches.get_one::<String>("bam").unwrap().clone(),
        filename_fastq: matches.get_one::<String>("fastq").unwrap().clone(),
        keep_pair_suffixes: matches.get_flag("keep-pair-suffixes"),
        report_limit: *matches.get_one::<usize>("report-limit").unwrap(),
    };

    match check_bam_fastq(&config) {
        Ok(result) => {
            if result.missing_unique == 0 {
                println!(
                    "All {} BAM records matched read names in {} FASTQ records ({} unique names).",
                    result.bam_records, result.fastq_records, result.fastq_unique_names
                );
            } else {
                eprintln!(
                    "Found {} missing BAM records covering {} unique read names not present in {} FASTQ records ({} unique names).",
                    result.missing_records,
                    result.missing_unique,
                    result.fastq_records,
                    result.fastq_unique_names
                );

                for name in result.missing_examples {
                    eprintln!("missing read: {}", name);
                }

                process::exit(1);
            }
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use super::normalize_read_name;

    #[test]
    fn normalize_fastq_headers() {
        assert_eq!(normalize_read_name("@read/1 comment", false), "read");
        assert_eq!(normalize_read_name("@read/2", false), "read");
        assert_eq!(normalize_read_name("read", false), "read");
    }

    #[test]
    fn keep_pair_suffixes_when_requested() {
        assert_eq!(normalize_read_name("@read/1 comment", true), "read/1");
    }
}
