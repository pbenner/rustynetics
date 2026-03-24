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

use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::process;

use clap::{Arg, Command};
use flate2::write::GzEncoder;
use flate2::Compression;

use rustynetics::bam::{BamBlock, BamReader, BamReaderOptions};
use rustynetics::fastq::{FastqRecord, FastqWriter};
use rustynetics::progress::{CountingReader, ProgressScope};

/* -------------------------------------------------------------------------- */

const STATUS_RECORD_UPDATE_INTERVAL: usize = 10_000;

/* -------------------------------------------------------------------------- */

struct Config {
    filename_bam: String,
    filename_out: Option<String>,
    filename_out1: Option<String>,
    filename_out2: Option<String>,
    filename_out_single: Option<String>,
    pair_suffixes: bool,
    include_secondary: bool,
    include_supplementary: bool,
    include_qcfail: bool,
    fill_missing_quality: Option<char>,
}

/* -------------------------------------------------------------------------- */

#[derive(Default)]
struct WriteStats {
    written: usize,
    written_read1: usize,
    written_read2: usize,
    written_single: usize,
    skipped_secondary: usize,
    skipped_supplementary: usize,
    skipped_qcfail: usize,
}

/* -------------------------------------------------------------------------- */

enum OutputMode {
    Single(FastqWriter<Box<dyn Write>>),
    Split {
        read1: FastqWriter<Box<dyn Write>>,
        read2: FastqWriter<Box<dyn Write>>,
        single: Option<FastqWriter<Box<dyn Write>>>,
    },
}

/* -------------------------------------------------------------------------- */

impl OutputMode {
    fn write_record(
        &mut self,
        block: &BamBlock,
        record: &FastqRecord,
        stats: &mut WriteStats,
    ) -> Result<(), Box<dyn Error>> {
        match self {
            OutputMode::Single(writer) => {
                writer.write_record(record)?;
                stats.written += 1;
                stats.written_single += 1;
            }
            OutputMode::Split {
                read1,
                read2,
                single,
            } => {
                if block.flag.read_paired() && block.flag.first_in_pair() {
                    read1.write_record(record)?;
                    stats.written += 1;
                    stats.written_read1 += 1;
                } else if block.flag.read_paired() && block.flag.second_in_pair() {
                    read2.write_record(record)?;
                    stats.written += 1;
                    stats.written_read2 += 1;
                } else if let Some(writer) = single {
                    writer.write_record(record)?;
                    stats.written += 1;
                    stats.written_single += 1;
                } else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "record `{}` is not routable in split-output mode; provide --output-single or use -o/--output",
                            block.read_name
                        ),
                    )
                    .into());
                }
            }
        }

        Ok(())
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            OutputMode::Single(writer) => writer.flush(),
            OutputMode::Split {
                read1,
                read2,
                single,
            } => {
                read1.flush()?;
                read2.flush()?;
                if let Some(writer) = single {
                    writer.flush()?;
                }
                Ok(())
            }
        }
    }
}

/* -------------------------------------------------------------------------- */

fn open_writer(filename: Option<&str>) -> Result<FastqWriter<Box<dyn Write>>, Box<dyn Error>> {
    let writer: Box<dyn Write> = match filename {
        None | Some("-") => Box::new(BufWriter::new(io::stdout())),
        Some(path) if path.ends_with(".gz") => {
            let file = File::create(path)?;
            Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
        }
        Some(path) => Box::new(BufWriter::new(File::create(path)?)),
    };

    Ok(FastqWriter::new(writer))
}

/* -------------------------------------------------------------------------- */

fn build_output_mode(config: &Config) -> Result<OutputMode, Box<dyn Error>> {
    let split_mode = config.filename_out1.is_some() || config.filename_out2.is_some();

    if split_mode {
        if config.filename_out.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "-o/--output cannot be combined with --output1/--output2",
            )
            .into());
        }

        let out1 = config.filename_out1.as_deref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "--output1 and --output2 must be specified together",
            )
        })?;
        let out2 = config.filename_out2.as_deref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "--output1 and --output2 must be specified together",
            )
        })?;

        Ok(OutputMode::Split {
            read1: open_writer(Some(out1))?,
            read2: open_writer(Some(out2))?,
            single: match config.filename_out_single.as_deref() {
                Some(path) => Some(open_writer(Some(path))?),
                None => None,
            },
        })
    } else {
        if config.filename_out_single.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--output-single requires --output1 and --output2",
            )
            .into());
        }

        Ok(OutputMode::Single(open_writer(
            config.filename_out.as_deref(),
        )?))
    }
}

/* -------------------------------------------------------------------------- */

fn should_skip_record(block: &BamBlock, config: &Config, stats: &mut WriteStats) -> bool {
    if block.flag.secondary_alignment() && !config.include_secondary {
        stats.skipped_secondary += 1;
        return true;
    }
    if block.flag.supplementary_alignment() && !config.include_supplementary {
        stats.skipped_supplementary += 1;
        return true;
    }
    if block.flag.not_passing_filters() && !config.include_qcfail {
        stats.skipped_qcfail += 1;
        return true;
    }

    false
}

/* -------------------------------------------------------------------------- */

fn parse_fill_missing_quality(value: &str) -> Result<char, String> {
    let mut chars = value.chars();
    let fill = chars
        .next()
        .ok_or_else(|| "missing quality fill character".to_string())?;

    if chars.next().is_some() {
        return Err("fill-missing-quality expects exactly one character".to_string());
    }
    if !fill.is_ascii() {
        return Err("fill-missing-quality must be an ASCII character".to_string());
    }

    Ok(fill)
}

/* -------------------------------------------------------------------------- */

fn format_skip_summary(stats: &WriteStats) -> String {
    let mut parts = Vec::new();

    if stats.skipped_secondary > 0 {
        parts.push(format!("{} secondary", stats.skipped_secondary));
    }
    if stats.skipped_supplementary > 0 {
        parts.push(format!("{} supplementary", stats.skipped_supplementary));
    }
    if stats.skipped_qcfail > 0 {
        parts.push(format!("{} QC-failed", stats.skipped_qcfail));
    }

    if parts.is_empty() {
        String::new()
    } else {
        format!("; skipped {}", parts.join(", "))
    }
}

/* -------------------------------------------------------------------------- */

fn bam_to_fastq(config: &Config) -> Result<WriteStats, Box<dyn Error>> {
    let total_bytes = File::open(&config.filename_bam)?.metadata()?.len();
    let mut progress = ProgressScope::new("BAM", total_bytes);

    let options = BamReaderOptions {
        read_name: true,
        read_cigar: false,
        read_sequence: true,
        read_auxiliary: false,
        read_qual: true,
    };

    let file = File::open(&config.filename_bam)?;
    let reader = BufReader::new(CountingReader::new(file, progress.handle()));
    let mut bam_reader = BamReader::new(reader, Some(options))?;
    let mut output_mode = build_output_mode(config)?;
    let mut stats = WriteStats::default();
    let mut processed = 0usize;

    for result in bam_reader.read_single_end() {
        let block = result?.block;
        processed += 1;

        if processed % STATUS_RECORD_UPDATE_INTERVAL == 0 {
            progress.set_records(processed);
        }

        if should_skip_record(&block, config, &mut stats) {
            continue;
        }

        let record =
            FastqRecord::from_bam_block(&block, config.pair_suffixes, config.fill_missing_quality)?;

        output_mode.write_record(&block, &record, &mut stats)?;
    }

    output_mode.flush()?;
    progress.finish(processed);

    Ok(stats)
}

/* -------------------------------------------------------------------------- */

fn main() {
    let matches = Command::new("bam-to-fastq")
        .version("1.0")
        .author("Philipp Benner [https://github.com/pbenner]")
        .about("Reconstruct FASTQ records from a BAM file")
        .arg(
            Arg::new("bam")
                .help("The input BAM file")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("Write all FASTQ records to a single file instead of stdout")
                .num_args(1),
        )
        .arg(
            Arg::new("output1")
                .long("output1")
                .help("Write first mates to this FASTQ file")
                .num_args(1),
        )
        .arg(
            Arg::new("output2")
                .long("output2")
                .help("Write second mates to this FASTQ file")
                .num_args(1),
        )
        .arg(
            Arg::new("output-single")
                .long("output-single")
                .help("Write single-end or unclassified records to this FASTQ file in split-output mode")
                .num_args(1),
        )
        .arg(
            Arg::new("pair-suffixes")
                .long("pair-suffixes")
                .action(clap::ArgAction::SetTrue)
                .help("Append /1 or /2 to read names for paired records"),
        )
        .arg(
            Arg::new("include-secondary")
                .long("include-secondary")
                .action(clap::ArgAction::SetTrue)
                .help("Include secondary alignments"),
        )
        .arg(
            Arg::new("include-supplementary")
                .long("include-supplementary")
                .action(clap::ArgAction::SetTrue)
                .help("Include supplementary alignments"),
        )
        .arg(
            Arg::new("include-qcfail")
                .long("include-qcfail")
                .action(clap::ArgAction::SetTrue)
                .help("Include reads marked as not passing filters"),
        )
        .arg(
            Arg::new("fill-missing-quality")
                .long("fill-missing-quality")
                .num_args(1)
                .value_parser(clap::builder::ValueParser::new(parse_fill_missing_quality))
                .help("Fill missing BAM qualities with this single ASCII character"),
        )
        .get_matches();

    let config = Config {
        filename_bam: matches.get_one::<String>("bam").unwrap().clone(),
        filename_out: matches.get_one::<String>("output").cloned(),
        filename_out1: matches.get_one::<String>("output1").cloned(),
        filename_out2: matches.get_one::<String>("output2").cloned(),
        filename_out_single: matches.get_one::<String>("output-single").cloned(),
        pair_suffixes: matches.get_flag("pair-suffixes"),
        include_secondary: matches.get_flag("include-secondary"),
        include_supplementary: matches.get_flag("include-supplementary"),
        include_qcfail: matches.get_flag("include-qcfail"),
        fill_missing_quality: matches.get_one::<char>("fill-missing-quality").copied(),
    };

    match bam_to_fastq(&config) {
        Ok(stats) => match (&config.filename_out1, &config.filename_out2) {
            (Some(_), Some(_)) => {
                eprintln!(
                    "Wrote {} FASTQ records ({} to read1, {} to read2, {} to singles){}.",
                    stats.written,
                    stats.written_read1,
                    stats.written_read2,
                    stats.written_single,
                    format_skip_summary(&stats)
                );
            }
            _ => {
                eprintln!(
                    "Wrote {} FASTQ records{}.",
                    stats.written,
                    format_skip_summary(&stats)
                );
            }
        },
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    }
}
