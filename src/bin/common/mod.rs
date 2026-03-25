#![allow(dead_code)]

use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use rustynetics::track::Track;
use rustynetics::track_bigwig::LazyTrackFile;
use rustynetics::track_simple::SimpleTrack;
use rustynetics::track_statistics::{bin_summary_statistics_from_string, BinSummaryStatistics};

pub fn parse_bin_summary(name: &str) -> Result<BinSummaryStatistics, String> {
    bin_summary_statistics_from_string(name)
        .ok_or_else(|| format!("invalid bin summary statistic `{name}`"))
}

pub fn parse_initial_value(value: Option<&str>, default: f64) -> Result<f64, Box<dyn Error>> {
    match value {
        Some(value) => Ok(value.parse()?),
        None => Ok(default),
    }
}

pub fn open_reader(path: Option<&str>) -> Result<Box<dyn BufRead>, Box<dyn Error>> {
    match path {
        None | Some("-") => Ok(Box::new(BufReader::new(io::stdin()))),
        Some(path) if path.ends_with(".gz") => {
            let file = File::open(path)?;
            Ok(Box::new(BufReader::new(GzDecoder::new(file))))
        }
        Some(path) => Ok(Box::new(BufReader::new(File::open(path)?))),
    }
}

pub fn open_writer(path: Option<&str>) -> Result<Box<dyn Write>, Box<dyn Error>> {
    match path {
        None | Some("-") => Ok(Box::new(BufWriter::new(io::stdout()))),
        Some(path) if path.ends_with(".gz") => {
            let file = File::create(path)?;
            Ok(Box::new(BufWriter::new(GzEncoder::new(
                file,
                Compression::default(),
            ))))
        }
        Some(path) => Ok(Box::new(BufWriter::new(File::create(path)?))),
    }
}

pub fn import_simple_track(
    filename: &str,
    name: &str,
    summary: BinSummaryStatistics,
    bin_size: usize,
    bin_overlap: usize,
    init: f64,
) -> Result<SimpleTrack, Box<dyn Error>> {
    let mut track = SimpleTrack::empty(name.to_string());
    track.import_bigwig(filename, name, summary, bin_size, bin_overlap, init)?;
    Ok(track)
}

pub fn export_simple_track(track: &SimpleTrack, filename: &str) -> Result<(), Box<dyn Error>> {
    track.export_bigwig(filename, vec![])?;
    Ok(())
}

pub fn import_lazy_track(
    filename: &str,
    name: &str,
    summary: BinSummaryStatistics,
    bin_size: usize,
    bin_overlap: usize,
    init: f64,
) -> Result<LazyTrackFile, Box<dyn Error>> {
    LazyTrackFile::import_bigwig(filename, name, summary, bin_size, bin_overlap, init)
}

#[derive(Clone, Copy, Debug)]
pub struct TrackSummary {
    pub valid: u64,
    pub min: f64,
    pub max: f64,
    pub sum: f64,
    pub sum_squares: f64,
}

impl Default for TrackSummary {
    fn default() -> Self {
        Self {
            valid: 0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            sum: 0.0,
            sum_squares: 0.0,
        }
    }
}

impl std::fmt::Display for TrackSummary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "(valid={}, min={:.4}, max={:.4}, sum={:.4}, sum_squares={:.4})",
            self.valid, self.min, self.max, self.sum, self.sum_squares
        )
    }
}

impl TrackSummary {
    pub fn mean(&self) -> f64 {
        if self.valid == 0 {
            f64::NAN
        } else {
            self.sum / self.valid as f64
        }
    }

    pub fn variance(&self) -> f64 {
        if self.valid == 0 {
            f64::NAN
        } else {
            self.sum_squares / self.valid as f64 - self.mean() * self.mean()
        }
    }
}

pub fn summarize_track(track: &dyn Track) -> Result<TrackSummary, Box<dyn Error>> {
    let mut summary = TrackSummary::default();

    for seqname in track.get_seq_names() {
        let sequence = track.get_sequence(&seqname)?;
        for i in 0..sequence.n_bins() {
            let value = sequence.at_bin(i);
            if value.is_nan() {
                continue;
            }
            summary.valid += 1;
            summary.min = summary.min.min(value);
            summary.max = summary.max.max(value);
            summary.sum += value;
            summary.sum_squares += value * value;
        }
    }

    Ok(summary)
}

pub fn collect_track_values(track: &dyn Track) -> Result<Vec<f64>, Box<dyn Error>> {
    let mut values = Vec::new();

    for seqname in track.get_seq_names() {
        let sequence = track.get_sequence(&seqname)?;
        for i in 0..sequence.n_bins() {
            let value = sequence.at_bin(i);
            if !value.is_nan() {
                values.push(value);
            }
        }
    }

    Ok(values)
}

pub fn write_fasta_record<W: Write>(writer: &mut W, name: &str, sequence: &[u8]) -> io::Result<()> {
    writeln!(writer, ">{name}")?;
    for chunk in sequence.chunks(80) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

pub fn worker_ranges(len: usize, threads: usize) -> Vec<(usize, usize)> {
    if len == 0 {
        return Vec::new();
    }

    let workers = threads.max(1).min(len);
    let base = len / workers;
    let remainder = len % workers;
    let mut start = 0usize;
    let mut ranges = Vec::with_capacity(workers);

    for i in 0..workers {
        let extra = usize::from(i < remainder);
        let end = start + base + extra;
        ranges.push((start, end));
        start = end;
    }

    ranges
}
