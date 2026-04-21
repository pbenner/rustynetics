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

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{self, BufRead, BufReader, Write};
use std::process;

use clap::{Arg, Command};
use flate2::read::MultiGzDecoder;

use rustynetics::bam::{BamFlag, BamReader, BamReaderOptions};
use rustynetics::fastq::FastqReader;
use rustynetics::progress::{CountingReader, ProgressScope};

/* -------------------------------------------------------------------------- */

struct Config {
    filename_bam: String,
    filenames_fastq: Vec<String>,
    keep_pair_suffixes: bool,
    report_limit: usize,
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
struct ReadNameFingerprint {
    lo: u64,
    hi: u64,
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Copy, Default)]
struct FastqNameStats {
    fastq_count: u32,
    bam_count: u32,
    file_mask: u64,
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Copy, Default)]
struct RecordMatchStats {
    total_records: usize,
    matched_records: usize,
    missing_records: usize,
}

impl RecordMatchStats {
    fn record(&mut self, matched: bool) {
        self.total_records += 1;
        if matched {
            self.matched_records += 1;
        } else {
            self.missing_records += 1;
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Default)]
struct BamClassSummary {
    primary: RecordMatchStats,
    secondary: RecordMatchStats,
    supplementary: RecordMatchStats,
    duplicate_flagged_records: usize,
}

/* -------------------------------------------------------------------------- */

#[derive(Default)]
struct FastqFileSummary {
    filename: String,
    records: usize,
    unique_names: usize,
    bam_overlapping_unique_names: usize,
    fastq_only_unique_names: usize,
    duplicate_name_records: usize,
    duplicated_unique_names: usize,
    warning: Option<String>,
}

/* -------------------------------------------------------------------------- */

#[derive(Default)]
struct GlobalNameSummary {
    bam_unique_names: usize,
    fastq_unique_names: usize,
    shared_unique_names: usize,
    bam_only_unique_names: usize,
    fastq_only_unique_names: usize,
    bam_duplicate_name_records: usize,
    bam_duplicated_unique_names: usize,
    fastq_duplicate_name_records: usize,
    fastq_duplicated_unique_names: usize,
}

/* -------------------------------------------------------------------------- */

struct BamRecordEvaluation {
    overall: RecordMatchStats,
    classes: BamClassSummary,
    bam_only_name_counts: HashMap<ReadNameFingerprint, u32>,
    missing_examples: Vec<String>,
}

/* -------------------------------------------------------------------------- */

struct CheckResult {
    fastq_files: usize,
    fastq_records: usize,
    fastq_unique_names: usize,
    fastq_only_unique_names: usize,
    fastq_duplicate_name_records: usize,
    fastq_duplicated_unique_names: usize,
    bam_records: usize,
    bam_unique_names: usize,
    bam_only_unique_names: usize,
    bam_duplicate_name_records: usize,
    bam_duplicated_unique_names: usize,
    shared_unique_names: usize,
    record_matches: RecordMatchStats,
    bam_classes: BamClassSummary,
    missing_examples: Vec<String>,
    fastq_file_summaries: Vec<FastqFileSummary>,
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
struct FileError {
    file_type: &'static str,
    filename: String,
    source: Box<dyn Error>,
}

impl fmt::Display for FileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "failed to read {} file `{}`: {}",
            self.file_type, self.filename, self.source
        )
    }
}

impl Error for FileError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(self.source.as_ref())
    }
}

/* -------------------------------------------------------------------------- */

const STATUS_RECORD_UPDATE_INTERVAL: usize = 10_000;
const REPORT_LABEL_WIDTH: usize = 28;
const MAX_FASTQ_FILES_FOR_FAST_MASK: usize = u64::BITS as usize;

/* -------------------------------------------------------------------------- */

enum BamRecordClass {
    Primary,
    Secondary,
    Supplementary,
}

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

fn classify_bam_record(flag: &BamFlag) -> BamRecordClass {
    if flag.supplementary_alignment() {
        BamRecordClass::Supplementary
    } else if flag.secondary_alignment() {
        BamRecordClass::Secondary
    } else {
        BamRecordClass::Primary
    }
}

/* -------------------------------------------------------------------------- */

fn fingerprint_read_name(name: &str) -> ReadNameFingerprint {
    let mut lo = std::collections::hash_map::DefaultHasher::new();
    0u8.hash(&mut lo);
    name.hash(&mut lo);

    let mut hi = std::collections::hash_map::DefaultHasher::new();
    1u8.hash(&mut hi);
    name.hash(&mut hi);

    ReadNameFingerprint {
        lo: lo.finish(),
        hi: hi.finish(),
    }
}

/* -------------------------------------------------------------------------- */

fn with_file_context<T>(
    file_type: &'static str,
    filename: &str,
    result: Result<T, Box<dyn Error>>,
) -> Result<T, Box<dyn Error>> {
    result.map_err(|source| -> Box<dyn Error> {
        Box::new(FileError {
            file_type,
            filename: filename.to_string(),
            source,
        })
    })
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

fn open_bam_reader(
    filename: &str,
    progress_label: &str,
) -> Result<(BamReader<BufReader<CountingReader<File>>>, ProgressScope), Box<dyn Error>> {
    let total_bytes = File::open(filename)?.metadata()?.len();
    let progress = ProgressScope::new(progress_label, total_bytes);
    let file = File::open(filename)?;
    let reader = BufReader::new(CountingReader::new(file, progress.handle()));
    let bam_reader = BamReader::new(
        reader,
        Some(BamReaderOptions {
            read_name: true,
            read_cigar: false,
            read_sequence: false,
            read_auxiliary: false,
            read_qual: false,
        }),
    )?;

    Ok((bam_reader, progress))
}

/* -------------------------------------------------------------------------- */

fn is_tolerated_fastq_tail_error(error: &io::Error, completed_records: usize) -> bool {
    let message = error.to_string();

    if error.kind() == io::ErrorKind::UnexpectedEof || message.contains("unexpected end of file") {
        return true;
    }

    completed_records > 0
        && error.kind() == io::ErrorKind::InvalidData
        && (message.starts_with("invalid FASTQ header at line ")
            || message.starts_with("invalid FASTQ separator at line ")
            || message.starts_with("FASTQ sequence and quality lengths differ "))
}

/* -------------------------------------------------------------------------- */

fn read_fastq_file(
    filename: &str,
    file_index: usize,
    progress_label: &str,
    keep_pair_suffixes: bool,
    fastq_index: &mut HashMap<ReadNameFingerprint, FastqNameStats>,
) -> Result<FastqFileSummary, Box<dyn Error>> {
    let total_bytes = File::open(filename)?.metadata()?.len();
    let mut progress = ProgressScope::new(progress_label, total_bytes);
    let reader = open_fastq_reader(filename, progress.handle())?;
    let mut reader = FastqReader::new(reader);
    let mut seen_names = HashSet::<ReadNameFingerprint>::new();
    let mut duplicated_names = HashSet::<ReadNameFingerprint>::new();
    let mut summary = FastqFileSummary {
        filename: filename.to_string(),
        ..FastqFileSummary::default()
    };

    loop {
        let record = match reader.read_record() {
            Ok(Some(record)) => record,
            Ok(None) => break,
            Err(error) if is_tolerated_fastq_tail_error(&error, summary.records) => {
                summary.warning = Some(format!(
                    "truncated or malformed FASTQ tail ignored after {} complete record(s): {}",
                    format_count(summary.records),
                    error
                ));
                break;
            }
            Err(error) => return Err(error.into()),
        };

        let name = normalize_read_name(record.name(), keep_pair_suffixes);
        if name.is_empty() {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "empty FASTQ read name").into());
        }
        let fingerprint = fingerprint_read_name(&name);

        summary.records += 1;

        let is_new_in_file = seen_names.insert(fingerprint);
        let entry = fastq_index.entry(fingerprint).or_default();
        entry.fastq_count = entry.fastq_count.saturating_add(1);

        if is_new_in_file {
            summary.unique_names += 1;
            entry.file_mask |= 1u64 << file_index;
        } else {
            summary.duplicate_name_records += 1;
            if duplicated_names.insert(fingerprint) {
                summary.duplicated_unique_names += 1;
            }
        }

        if summary.records % STATUS_RECORD_UPDATE_INTERVAL == 0 {
            progress.set_records(summary.records);
        }
    }

    progress.finish(summary.records);

    Ok(summary)
}

/* -------------------------------------------------------------------------- */

fn read_fastq_files(
    filenames: &[String],
    keep_pair_suffixes: bool,
    fastq_index: &mut HashMap<ReadNameFingerprint, FastqNameStats>,
) -> Result<Vec<FastqFileSummary>, Box<dyn Error>> {
    if filenames.len() > MAX_FASTQ_FILES_FOR_FAST_MASK {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "bam-check-fastq currently supports at most {} FASTQ files for detailed per-file statistics",
                MAX_FASTQ_FILES_FOR_FAST_MASK
            ),
        )
        .into());
    }

    let mut summaries = Vec::with_capacity(filenames.len());

    for (index, filename) in filenames.iter().enumerate() {
        let progress_label = if filenames.len() > 1 {
            format!("FASTQ {}/{}", index + 1, filenames.len())
        } else {
            "FASTQ".to_string()
        };
        summaries.push(with_file_context(
            "FASTQ",
            filename,
            read_fastq_file(
                filename,
                index,
                &progress_label,
                keep_pair_suffixes,
                fastq_index,
            ),
        )?);
    }

    Ok(summaries)
}

/* -------------------------------------------------------------------------- */

fn summarize_fastq_and_bam_names(
    fastq_index: &HashMap<ReadNameFingerprint, FastqNameStats>,
    bam_only_name_counts: &HashMap<ReadNameFingerprint, u32>,
    fastq_file_summaries: &mut [FastqFileSummary],
) -> GlobalNameSummary {
    let mut summary = GlobalNameSummary::default();

    for stats in fastq_index.values() {
        summary.fastq_unique_names += 1;

        if stats.fastq_count > 1 {
            summary.fastq_duplicate_name_records += (stats.fastq_count - 1) as usize;
            summary.fastq_duplicated_unique_names += 1;
        }

        if stats.bam_count > 0 {
            summary.shared_unique_names += 1;
            summary.bam_unique_names += 1;

            if stats.bam_count > 1 {
                summary.bam_duplicate_name_records += (stats.bam_count - 1) as usize;
                summary.bam_duplicated_unique_names += 1;
            }

            let mut file_mask = stats.file_mask;
            while file_mask != 0 {
                let file_index = file_mask.trailing_zeros() as usize;
                fastq_file_summaries[file_index].bam_overlapping_unique_names += 1;
                file_mask &= file_mask - 1;
            }
        } else {
            summary.fastq_only_unique_names += 1;

            let mut file_mask = stats.file_mask;
            while file_mask != 0 {
                let file_index = file_mask.trailing_zeros() as usize;
                fastq_file_summaries[file_index].fastq_only_unique_names += 1;
                file_mask &= file_mask - 1;
            }
        }
    }

    for bam_count in bam_only_name_counts.values() {
        summary.bam_unique_names += 1;
        summary.bam_only_unique_names += 1;

        if *bam_count > 1 {
            summary.bam_duplicate_name_records += (*bam_count - 1) as usize;
            summary.bam_duplicated_unique_names += 1;
        }
    }

    summary
}

/* -------------------------------------------------------------------------- */

fn evaluate_bam_records(
    filename: &str,
    keep_pair_suffixes: bool,
    report_limit: usize,
    fastq_index: &mut HashMap<ReadNameFingerprint, FastqNameStats>,
) -> Result<BamRecordEvaluation, Box<dyn Error>> {
    let (mut bam_reader, mut progress) = open_bam_reader(filename, "BAM")?;
    let mut evaluation = BamRecordEvaluation {
        overall: RecordMatchStats::default(),
        classes: BamClassSummary::default(),
        bam_only_name_counts: HashMap::new(),
        missing_examples: Vec::new(),
    };

    for result in bam_reader.read_single_end() {
        let block = result?.block;
        let name = normalize_read_name(&block.read_name, keep_pair_suffixes);

        if name.is_empty() {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "empty BAM read name").into());
        }
        let fingerprint = fingerprint_read_name(&name);

        let matched = if let Some(stats) = fastq_index.get_mut(&fingerprint) {
            stats.bam_count = stats.bam_count.saturating_add(1);
            true
        } else {
            let count = evaluation
                .bam_only_name_counts
                .entry(fingerprint)
                .or_insert(0);
            *count = count.saturating_add(1);
            if *count == 1 && evaluation.missing_examples.len() < report_limit {
                evaluation.missing_examples.push(name);
            }
            false
        };

        evaluation.overall.record(matched);
        if block.flag.duplicate() {
            evaluation.classes.duplicate_flagged_records += 1;
        }

        match classify_bam_record(&block.flag) {
            BamRecordClass::Primary => evaluation.classes.primary.record(matched),
            BamRecordClass::Secondary => evaluation.classes.secondary.record(matched),
            BamRecordClass::Supplementary => evaluation.classes.supplementary.record(matched),
        }

        if evaluation.overall.total_records % STATUS_RECORD_UPDATE_INTERVAL == 0 {
            progress.set_records(evaluation.overall.total_records);
        }
    }

    progress.finish(evaluation.overall.total_records);

    Ok(evaluation)
}

/* -------------------------------------------------------------------------- */

fn check_bam_fastq(config: &Config) -> Result<CheckResult, Box<dyn Error>> {
    let mut fastq_index = HashMap::new();

    let mut fastq_file_summaries = read_fastq_files(
        &config.filenames_fastq,
        config.keep_pair_suffixes,
        &mut fastq_index,
    )?;
    let bam_evaluation = with_file_context(
        "BAM",
        &config.filename_bam,
        evaluate_bam_records(
            &config.filename_bam,
            config.keep_pair_suffixes,
            config.report_limit,
            &mut fastq_index,
        ),
    )?;
    let global_name_summary = summarize_fastq_and_bam_names(
        &fastq_index,
        &bam_evaluation.bam_only_name_counts,
        &mut fastq_file_summaries,
    );
    let fastq_records = fastq_file_summaries
        .iter()
        .map(|summary| summary.records)
        .sum();

    Ok(CheckResult {
        fastq_files: config.filenames_fastq.len(),
        fastq_records,
        fastq_unique_names: global_name_summary.fastq_unique_names,
        fastq_only_unique_names: global_name_summary.fastq_only_unique_names,
        fastq_duplicate_name_records: global_name_summary.fastq_duplicate_name_records,
        fastq_duplicated_unique_names: global_name_summary.fastq_duplicated_unique_names,
        bam_records: bam_evaluation.overall.total_records,
        bam_unique_names: global_name_summary.bam_unique_names,
        bam_only_unique_names: global_name_summary.bam_only_unique_names,
        bam_duplicate_name_records: global_name_summary.bam_duplicate_name_records,
        bam_duplicated_unique_names: global_name_summary.bam_duplicated_unique_names,
        shared_unique_names: global_name_summary.shared_unique_names,
        record_matches: bam_evaluation.overall,
        bam_classes: bam_evaluation.classes,
        missing_examples: bam_evaluation.missing_examples,
        fastq_file_summaries,
    })
}

/* -------------------------------------------------------------------------- */

fn percentage(part: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        100.0 * part as f64 / total as f64
    }
}

/* -------------------------------------------------------------------------- */

fn format_count(value: usize) -> String {
    let digits = value.to_string();
    let mut formatted = String::with_capacity(digits.len() + digits.len() / 3);

    for (index, ch) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index) % 3 == 0 {
            formatted.push(',');
        }
        formatted.push(ch);
    }

    formatted
}

/* -------------------------------------------------------------------------- */

fn format_count_and_percentage(part: usize, total: usize) -> String {
    format!("{} ({:.2}%)", format_count(part), percentage(part, total))
}

/* -------------------------------------------------------------------------- */

fn format_fraction_and_percentage(part: usize, total: usize) -> String {
    format!(
        "{} / {} ({:.2}%)",
        format_count(part),
        format_count(total),
        percentage(part, total)
    )
}

/* -------------------------------------------------------------------------- */

fn print_section_title<W: Write>(writer: &mut W, title: &str) -> io::Result<()> {
    writeln!(writer, "{title}")?;
    writeln!(writer, "{}", "-".repeat(title.len()))
}

/* -------------------------------------------------------------------------- */

fn print_stat_line<W: Write>(
    writer: &mut W,
    label: &str,
    value: impl fmt::Display,
) -> io::Result<()> {
    writeln!(
        writer,
        "  {:<width$} {}",
        label,
        value,
        width = REPORT_LABEL_WIDTH
    )
}

/* -------------------------------------------------------------------------- */

fn print_file_stat_line<W: Write>(
    writer: &mut W,
    label: &str,
    value: impl fmt::Display,
) -> io::Result<()> {
    writeln!(
        writer,
        "    {:<width$} {}",
        label,
        value,
        width = REPORT_LABEL_WIDTH.saturating_sub(2)
    )
}

/* -------------------------------------------------------------------------- */

fn print_record_class_line<W: Write>(
    writer: &mut W,
    label: &str,
    stats: &RecordMatchStats,
) -> io::Result<()> {
    writeln!(
        writer,
        "  {:<14} {:>12} {:>22} {:>22}",
        label,
        format_count(stats.total_records),
        format_count_and_percentage(stats.matched_records, stats.total_records),
        format_count_and_percentage(stats.missing_records, stats.total_records)
    )
}

/* -------------------------------------------------------------------------- */

fn print_report<W: Write>(
    writer: &mut W,
    config: &Config,
    result: &CheckResult,
    passed: bool,
) -> io::Result<()> {
    let fastq_warning_count = result
        .fastq_file_summaries
        .iter()
        .filter(|summary| summary.warning.is_some())
        .count();
    let has_warnings = fastq_warning_count > 0;
    let coverage_status = if passed { "COMPLETE" } else { "INCOMPLETE" };
    let fastq_input_status = if has_warnings {
        format!(
            "INCOMPLETE ({} of {} files truncated)",
            format_count(fastq_warning_count),
            format_count(result.fastq_files)
        )
    } else {
        "COMPLETE".to_string()
    };
    let interpretation = match (passed, has_warnings) {
        (true, false) => "all BAM reads were found in at least one provided FASTQ file",
        (true, true) => {
            "all BAM reads were found in the readable FASTQ records; truncated FASTQ files did not hide any BAM reads in this comparison"
        }
        (false, false) => "some BAM reads were not found in any provided FASTQ file",
        (false, true) => {
            "some BAM reads were not found in the readable FASTQ records; because one or more FASTQ files were truncated, some missing BAM reads may lie in unread FASTQ tails"
        }
    };

    print_section_title(writer, "Result")?;
    print_stat_line(writer, "BAM->FASTQ coverage", coverage_status)?;
    print_stat_line(
        writer,
        "Found BAM records",
        format_fraction_and_percentage(result.record_matches.matched_records, result.bam_records),
    )?;
    print_stat_line(
        writer,
        "Missing BAM records",
        format_fraction_and_percentage(result.record_matches.missing_records, result.bam_records),
    )?;
    print_stat_line(
        writer,
        "Found BAM unique names",
        format_fraction_and_percentage(result.shared_unique_names, result.bam_unique_names),
    )?;
    print_stat_line(
        writer,
        "Missing BAM unique names",
        format_fraction_and_percentage(result.bam_only_unique_names, result.bam_unique_names),
    )?;
    print_stat_line(writer, "FASTQ input status", fastq_input_status)?;
    print_stat_line(writer, "Interpretation", interpretation)?;
    print_stat_line(writer, "BAM file", &config.filename_bam)?;
    print_stat_line(writer, "FASTQ files", format_count(result.fastq_files))?;
    print_stat_line(
        writer,
        "Read name normalization",
        if config.keep_pair_suffixes {
            "keep trailing /1 and /2 suffixes"
        } else {
            "strip trailing /1 and /2 suffixes"
        },
    )?;

    if has_warnings {
        writeln!(writer)?;
        print_section_title(writer, "FASTQ Warnings")?;
        for (index, summary) in result.fastq_file_summaries.iter().enumerate() {
            if let Some(warning) = &summary.warning {
                writeln!(
                    writer,
                    "  [{}/{}] {}",
                    index + 1,
                    result.fastq_file_summaries.len(),
                    summary.filename
                )?;
                print_file_stat_line(writer, "Status", "TRUNCATED")?;
                print_file_stat_line(writer, "Warning", warning)?;
                writeln!(writer)?;
            }
        }
    }

    if !passed && !result.missing_examples.is_empty() {
        writeln!(writer)?;
        print_section_title(writer, "Missing BAM Read Names")?;
        print_stat_line(
            writer,
            "Shown",
            format!(
                "{} of {} missing unique read names",
                format_count(result.missing_examples.len()),
                format_count(result.bam_only_unique_names)
            ),
        )?;
        for name in &result.missing_examples {
            writeln!(writer, "  {}", name)?;
        }
    }

    writeln!(writer)?;
    print_section_title(writer, "BAM Record Coverage")?;
    print_stat_line(writer, "BAM records", format_count(result.bam_records))?;
    print_stat_line(
        writer,
        "Matched BAM records",
        format_count_and_percentage(result.record_matches.matched_records, result.bam_records),
    )?;
    print_stat_line(
        writer,
        "Missing BAM records",
        format_count_and_percentage(result.record_matches.missing_records, result.bam_records),
    )?;

    writeln!(writer)?;
    print_section_title(writer, "BAM Unique Read Names")?;
    print_stat_line(
        writer,
        "BAM unique read names",
        format_count(result.bam_unique_names),
    )?;
    print_stat_line(
        writer,
        "Found in FASTQ",
        format_count_and_percentage(result.shared_unique_names, result.bam_unique_names),
    )?;
    print_stat_line(
        writer,
        "Missing from FASTQ",
        format_count_and_percentage(result.bam_only_unique_names, result.bam_unique_names),
    )?;
    print_stat_line(
        writer,
        "Duplicate-name BAM records",
        format_count(result.bam_duplicate_name_records),
    )?;
    print_stat_line(
        writer,
        "Duplicated BAM read names",
        format_count(result.bam_duplicated_unique_names),
    )?;

    writeln!(writer)?;
    print_section_title(writer, "FASTQ Unique Read Names")?;
    print_stat_line(writer, "FASTQ records", format_count(result.fastq_records))?;
    print_stat_line(
        writer,
        "FASTQ unique read names",
        format_count(result.fastq_unique_names),
    )?;
    print_stat_line(
        writer,
        "Found in BAM",
        format_count_and_percentage(result.shared_unique_names, result.fastq_unique_names),
    )?;
    print_stat_line(
        writer,
        "FASTQ-only unique read names",
        format_count_and_percentage(result.fastq_only_unique_names, result.fastq_unique_names),
    )?;
    print_stat_line(
        writer,
        "Duplicate-name FASTQ records",
        format_count(result.fastq_duplicate_name_records),
    )?;
    print_stat_line(
        writer,
        "Duplicated FASTQ read names",
        format_count(result.fastq_duplicated_unique_names),
    )?;

    writeln!(writer)?;
    print_section_title(writer, "BAM Record Classes")?;
    writeln!(
        writer,
        "  {:<14} {:>12} {:>22} {:>22}",
        "Class", "Records", "Matched", "Missing"
    )?;
    print_record_class_line(writer, "Primary", &result.bam_classes.primary)?;
    print_record_class_line(writer, "Secondary", &result.bam_classes.secondary)?;
    print_record_class_line(writer, "Supplementary", &result.bam_classes.supplementary)?;
    print_stat_line(
        writer,
        "Duplicate-flagged BAM records",
        format_count(result.bam_classes.duplicate_flagged_records),
    )?;

    writeln!(writer)?;
    print_section_title(writer, "FASTQ File Summary")?;
    for (index, summary) in result.fastq_file_summaries.iter().enumerate() {
        writeln!(
            writer,
            "  [{}/{}] {}",
            index + 1,
            result.fastq_file_summaries.len(),
            summary.filename
        )?;
        print_file_stat_line(
            writer,
            "Status",
            if summary.warning.is_some() {
                "TRUNCATED"
            } else {
                "OK"
            },
        )?;
        print_file_stat_line(writer, "Records", format_count(summary.records))?;
        print_file_stat_line(
            writer,
            "Unique read names",
            format_count(summary.unique_names),
        )?;
        print_file_stat_line(
            writer,
            "Unique names also in BAM",
            format_count_and_percentage(summary.bam_overlapping_unique_names, summary.unique_names),
        )?;
        print_file_stat_line(
            writer,
            "FASTQ-only unique names",
            format_count_and_percentage(summary.fastq_only_unique_names, summary.unique_names),
        )?;
        print_file_stat_line(
            writer,
            "Duplicate-name records",
            format_count(summary.duplicate_name_records),
        )?;
        print_file_stat_line(
            writer,
            "Duplicated read names",
            format_count(summary.duplicated_unique_names),
        )?;
        if let Some(warning) = &summary.warning {
            print_file_stat_line(writer, "Warning", warning)?;
        }
        if index + 1 < result.fastq_file_summaries.len() {
            writeln!(writer)?;
        }
    }

    Ok(())
}

/* -------------------------------------------------------------------------- */

fn main() {
    let matches = Command::new("bam-check-fastq")
        .version("1.0")
        .author("Philipp Benner [https://github.com/pbenner]")
        .about("Checks whether all BAM read names are present in one or more FASTQ files")
        .arg(
            Arg::new("bam")
                .help("The input BAM file")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("fastq")
                .help("One or more FASTQ files to compare against (plain text or gzip compressed)")
                .required(true)
                .index(2)
                .num_args(1..),
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
        filenames_fastq: matches
            .get_many::<String>("fastq")
            .unwrap()
            .cloned()
            .collect(),
        keep_pair_suffixes: matches.get_flag("keep-pair-suffixes"),
        report_limit: *matches.get_one::<usize>("report-limit").unwrap(),
    };

    match check_bam_fastq(&config) {
        Ok(result) => {
            if result.bam_only_unique_names == 0 {
                let mut stdout = io::stdout().lock();
                print_report(&mut stdout, &config, &result, true).unwrap();
            } else {
                let mut stderr = io::stderr().lock();
                print_report(&mut stderr, &config, &result, false).unwrap();
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
