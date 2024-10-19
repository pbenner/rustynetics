/* Copyright (C) 2024 Philipp Benner
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

use std::fmt;
use std::io;
use std::error::Error;

use core::pin::Pin;
use futures::{Stream, StreamExt};
use async_stream::stream;
use futures::executor::block_on_stream;

use crate::bam::{BamFile, bam_import_genome};
use crate::genome::Genome;
use crate::reads;
use crate::infologger::Logger;
use crate::error::ArgumentError;

use crate::track::MutableTrack;
use crate::track_generic::GenericMutableTrack;
use crate::track_simple::SimpleTrack;
use crate::track_statistics::estimate_fragment_length;

/* -------------------------------------------------------------------------- */

// Define an enum to represent all the possible option values
pub enum OptionBamCoverage {
    Logger(Logger),
    BinningMethod(String),
    BinSize(usize),
    BinOverlap(i64),
    NormalizeTrack(String),
    ShiftReads([usize; 2]),
    PairedAsSingleEnd(bool),
    PairedEndStrandSpecific(bool),
    LogScale(bool),
    Pseudocounts([f64; 2]),
    EstimateFraglen(bool),
    FraglenRange((i32, i32)),
    FraglenBinSize(usize),
    FilterChroms(Vec<String>),
    RemoveFilteredChroms(bool),
    FilterMapQ(i64),
    FilterReadLengths([usize; 2]),
    FilterDuplicates(bool),
    FilterStrand(char),
    FilterPairedEnd(bool),
    FilterSingleEnd(bool),
    SmoothenControl(bool),
    SmoothenSizes(Vec<usize>),
    SmoothenMin(f64),
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for OptionBamCoverage {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OptionBamCoverage::Logger(_) => write!(f, "Logger option"),
            OptionBamCoverage::BinningMethod(s) => write!(f, "Binning Method: {}", s),
            OptionBamCoverage::BinSize(size) => write!(f, "Bin Size: {}", size),
            OptionBamCoverage::BinOverlap(overlap) => write!(f, "Bin Overlap: {}", overlap),
            OptionBamCoverage::NormalizeTrack(s) => write!(f, "Normalize Track: {}", s),
            OptionBamCoverage::ShiftReads(arr) => write!(f, "Shift Reads: {:?}", arr),
            OptionBamCoverage::PairedAsSingleEnd(b) => write!(f, "Paired as Single End: {}", b),
            OptionBamCoverage::PairedEndStrandSpecific(b) => write!(f, "Paired End Strand Specific: {}", b),
            OptionBamCoverage::LogScale(b) => write!(f, "Log Scale: {}", b),
            OptionBamCoverage::Pseudocounts(arr) => write!(f, "Pseudocounts: {:?}", arr),
            OptionBamCoverage::EstimateFraglen(b) => write!(f, "Estimate Fraglen: {}", b),
            OptionBamCoverage::FraglenRange(arr) => write!(f, "Fraglen Range: {:?}", arr),
            OptionBamCoverage::FraglenBinSize(size) => write!(f, "Fraglen Bin Size: {}", size),
            OptionBamCoverage::FilterChroms(v) => write!(f, "Filter Chroms: {:?}", v),
            OptionBamCoverage::RemoveFilteredChroms(b) => write!(f, "Remove Filtered Chroms: {}", b),
            OptionBamCoverage::FilterMapQ(q) => write!(f, "Filter MapQ: {}", q),
            OptionBamCoverage::FilterReadLengths(arr) => write!(f, "Filter Read Lengths: {:?}", arr),
            OptionBamCoverage::FilterDuplicates(b) => write!(f, "Filter Duplicates: {}", b),
            OptionBamCoverage::FilterStrand(strand) => write!(f, "Filter Strand: {}", strand),
            OptionBamCoverage::FilterPairedEnd(b) => write!(f, "Filter Paired End: {}", b),
            OptionBamCoverage::FilterSingleEnd(b) => write!(f, "Filter Single End: {}", b),
            OptionBamCoverage::SmoothenControl(b) => write!(f, "Smoothen Control: {}", b),
            OptionBamCoverage::SmoothenSizes(v) => write!(f, "Smoothen Sizes: {:?}", v),
            OptionBamCoverage::SmoothenMin(min) => write!(f, "Smoothen Min: {}", min),
        }
    }
}

/* -------------------------------------------------------------------------- */

// Define the BamCoverageConfig struct
pub struct BamCoverageConfig {
    pub logger: Logger,
    pub binning_method: String,
    pub bin_size: usize,
    pub bin_overlap: i64,
    pub normalize_track: String,
    pub shift_reads: [usize; 2],
    pub paired_as_single_end: bool,
    pub paired_end_strand_specific: bool,
    pub log_scale: bool,
    pub pseudocounts: [f64; 2],
    pub estimate_fraglen: bool,
    pub fraglen_range: (i32, i32),
    pub fraglen_bin_size: usize,
    pub filter_chroms: Vec<String>,
    pub filter_mapq: i64,
    pub filter_read_lengths: [usize; 2],
    pub filter_duplicates: bool,
    pub filter_strand: char,
    pub filter_paired_end: bool,
    pub filter_single_end: bool,
    pub remove_filtered_chroms: bool,
    pub smoothen_control: bool,
    pub smoothen_sizes: Vec<usize>,
    pub smoothen_min: f64,
}

/* -------------------------------------------------------------------------- */

// Function to insert OptionBamCoverage into BamCoverageConfig
impl BamCoverageConfig {
    pub fn insert_option(&mut self, option: OptionBamCoverage) {
        match option {
            OptionBamCoverage::Logger(logger) => {
                self.logger = logger;
            }
            OptionBamCoverage::BinningMethod(method) => {
                self.binning_method = method;
            }
            OptionBamCoverage::BinSize(size) => {
                self.bin_size = size;
            }
            OptionBamCoverage::BinOverlap(overlap) => {
                self.bin_overlap = overlap;
            }
            OptionBamCoverage::NormalizeTrack(track) => {
                self.normalize_track = track;
            }
            OptionBamCoverage::ShiftReads(reads) => {
                self.shift_reads = reads;
            }
            OptionBamCoverage::PairedAsSingleEnd(paired) => {
                self.paired_as_single_end = paired;
            }
            OptionBamCoverage::PairedEndStrandSpecific(strand_specific) => {
                self.paired_end_strand_specific = strand_specific;
            }
            OptionBamCoverage::LogScale(log_scale) => {
                self.log_scale = log_scale;
            }
            OptionBamCoverage::Pseudocounts(pseudocounts) => {
                self.pseudocounts = pseudocounts;
            }
            OptionBamCoverage::EstimateFraglen(estimate) => {
                self.estimate_fraglen = estimate;
            }
            OptionBamCoverage::FraglenRange(range) => {
                self.fraglen_range = range;
            }
            OptionBamCoverage::FraglenBinSize(size) => {
                self.fraglen_bin_size = size;
            }
            OptionBamCoverage::FilterChroms(chroms) => {
                self.filter_chroms = chroms;
            }
            OptionBamCoverage::FilterMapQ(mapq) => {
                self.filter_mapq = mapq;
            }
            OptionBamCoverage::FilterReadLengths(read_lengths) => {
                self.filter_read_lengths = read_lengths;
            }
            OptionBamCoverage::FilterDuplicates(duplicates) => {
                self.filter_duplicates = duplicates;
            }
            OptionBamCoverage::FilterStrand(strand) => {
                self.filter_strand = strand;
            }
            OptionBamCoverage::FilterPairedEnd(paired_end) => {
                self.filter_paired_end = paired_end;
            }
            OptionBamCoverage::FilterSingleEnd(single_end) => {
                self.filter_single_end = single_end;
            }
            OptionBamCoverage::RemoveFilteredChroms(remove) => {
                self.remove_filtered_chroms = remove;
            }
            OptionBamCoverage::SmoothenControl(smoothen) => {
                self.smoothen_control = smoothen;
            }
            OptionBamCoverage::SmoothenSizes(sizes) => {
                self.smoothen_sizes = sizes;
            }
            OptionBamCoverage::SmoothenMin(min) => {
                self.smoothen_min = min;
            }
        }
    }
}

/* -------------------------------------------------------------------------- */

// Function to provide default configuration
impl BamCoverageConfig {
    pub fn default() -> Self {
        BamCoverageConfig {
            logger: Logger::new_null(), // Discard logger output
            binning_method: String::from("simple"),
            bin_size: 10,
            bin_overlap: 0,
            normalize_track: String::new(),
            shift_reads: [0, 0],
            paired_as_single_end: false,
            paired_end_strand_specific: false,
            log_scale: false,
            pseudocounts: [0.0, 0.0],
            estimate_fraglen: false,
            fraglen_range: (-1, -1),
            fraglen_bin_size: 10,
            filter_chroms: Vec::new(),
            filter_mapq: 0,
            filter_read_lengths: [0, 0],
            filter_duplicates: false,
            filter_strand: '*',
            filter_paired_end: false,
            filter_single_end: false,
            remove_filtered_chroms: false,
            smoothen_control: false,
            smoothen_sizes: Vec::new(),
            smoothen_min: 20.0,
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, Default)]
pub struct FraglenEstimate {
    pub fraglen: usize,
    pub x: Vec<i32>,
    pub y: Vec<f64>,
}

/* -------------------------------------------------------------------------- */

type ReadStream<'a> = Pin<Box<dyn Stream<Item = io::Result<reads::Read>> + 'a>>;

/* -------------------------------------------------------------------------- */

pub fn filter_paired_as_single_end<'a>(
    config: &'a BamCoverageConfig,
    mut stream_in: ReadStream<'a>,
) -> ReadStream<'a> {

    // If PairedAsSingleEnd is false, return the input stream directly
    if !config.paired_as_single_end {
        return stream_in;
    }

    let output_stream = Box::pin(stream! {
        while let Some(item) = stream_in.next().await {

            match item {
                Ok(mut r) => {
                    r.paired_end = false;
                    yield Ok(r)
                },
                Err(e) => yield Err(e),
            };
        }
    });

    Box::pin(output_stream)
}

/* -------------------------------------------------------------------------- */

// Function to filter paired reads
pub fn filter_paired_end<'a>(config: &'a BamCoverageConfig, mut stream_in: ReadStream<'a>) -> ReadStream<'a> {
    if !config.filter_paired_end {
        return stream_in;
    }

    let output_stream = async_stream::stream! {
        let mut n = 0;
        let mut m = 0;

        while let Some(item) = stream_in.next().await {
            match item {
                Ok(r) => {
                    if r.paired_end {
                        yield Ok(r);
                        m += 1;
                    }
                    n += 1;
                },
                Err(e) => yield Err(e),
            }
        }

        if n != 0 {
            log!(config.logger, "Filtered out {} unpaired reads ({:.2}%)", n - m, 100.0 * (n - m) as f64 / n as f64);
        }
    };

    Box::pin(output_stream)
}

/* -------------------------------------------------------------------------- */

// Function to filter single-end reads
pub fn filter_single_end<'a>(config: &'a BamCoverageConfig, veto: bool, mut stream_in: ReadStream<'a>) -> ReadStream<'a> {
    if !config.filter_single_end && !veto {
        return stream_in;
    }

    let output_stream = async_stream::stream! {
        let mut n = 0;
        let mut m = 0;

        while let Some(item) = stream_in.next().await {
            match item {
                Ok(r) => {
                    if !r.paired_end {
                        yield Ok(r);
                        m += 1;
                    }
                    n += 1;
                },
                Err(e) => yield Err(e),
            }
        }

        if n != 0 {
            log!(config.logger, "Filtered out {} paired reads ({:.2}%)", n - m, 100.0 * (n - m) as f64 / n as f64);
        }
    };

    Box::pin(output_stream)
}

/* -------------------------------------------------------------------------- */

// Function to filter duplicates
pub fn filter_duplicates<'a>(config: &'a BamCoverageConfig, mut stream_in: ReadStream<'a>) -> ReadStream<'a> {
    if !config.filter_duplicates {
        return stream_in;
    }

    let output_stream = async_stream::stream! {
        let mut n = 0;
        let mut m = 0;

        while let Some(item) = stream_in.next().await {
            match item {
                Ok(r) => {
                    if !r.duplicate {
                        yield Ok(r);
                        m += 1;
                    }
                    n += 1;
                },
                Err(e) => yield Err(e),
            }
        }

        if n != 0 {
            log!(config.logger, "Filtered out {} duplicates ({:.2}%)", n - m, 100.0 * (n - m) as f64 / n as f64);
        }
    };

    Box::pin(output_stream)
}

/* -------------------------------------------------------------------------- */

// Function to filter based on strand
pub fn filter_strand<'a>(config: &'a BamCoverageConfig, mut stream_in: ReadStream<'a>) -> ReadStream<'a> {
    if config.filter_strand == '*' {
        return stream_in;
    }

    let output_stream = async_stream::stream! {
        let mut n = 0;
        let mut m = 0;

        while let Some(item) = stream_in.next().await {
            match item {
                Ok(r) => {
                    if r.strand == config.filter_strand {
                        yield Ok(r);
                        m += 1;
                    }
                    n += 1;
                },
                Err(e) => yield Err(e),
            }
        }

        if n != 0 {
            log!(config.logger, "Filtered out {} reads not on strand {} ({:.2}%)", n - m, config.filter_strand, 100.0 * (n - m) as f64 / n as f64);
        }
    };

    Box::pin(output_stream)
}

/* -------------------------------------------------------------------------- */

// Function to filter based on mapping quality
pub fn filter_mapq<'a>(config: &'a BamCoverageConfig, mut stream_in: ReadStream<'a>) -> ReadStream<'a> {
    if config.filter_mapq <= 0 {
        return stream_in;
    }

    let output_stream = async_stream::stream! {
        let mut n = 0;
        let mut m = 0;

        while let Some(item) = stream_in.next().await {
            match item {
                Ok(r) => {
                    if r.mapq >= config.filter_mapq {
                        yield Ok(r);
                        m += 1;
                    }
                    n += 1;
                },
                Err(e) => yield Err(e),
            }
        }

        if n != 0 {
            log!(config.logger, "Filtered out {} reads with mapping quality lower than {} ({:.2}%)", n - m, config.filter_mapq, 100.0 * (n - m) as f64 / n as f64);
        }
    };

    Box::pin(output_stream)
}

/* -------------------------------------------------------------------------- */

// Function to filter based on read length
pub fn filter_read_length<'a>(config: &'a BamCoverageConfig, mut stream_in: ReadStream<'a>) -> ReadStream<'a> {
    if config.filter_read_lengths[0] == 0 && config.filter_read_lengths[1] == 0 {
        return stream_in;
    }

    let output_stream = async_stream::stream! {
        let mut n = 0;
        let mut m = 0;

        while let Some(item) = stream_in.next().await {
            match item {
                Ok(r) => {
                    let len = r.range.to - r.range.from;
                    if len >= config.filter_read_lengths[0] &&
                       (len <= config.filter_read_lengths[1] || config.filter_read_lengths[1] == 0) {
                        yield Ok(r);
                        m += 1;
                    }
                    n += 1;
                },
                Err(e) => yield Err(e),
            }
        }

        if n != 0 {
            log!(config.logger, "Filtered out {} reads with non-admissible length ({:.2}%)", n - m, 100.0 * (n - m) as f64 / n as f64);
        }
    };

    Box::pin(output_stream)
}

/* -------------------------------------------------------------------------- */

// Function to shift reads based on strand
pub fn shift_reads<'a>(config: &'a BamCoverageConfig, mut stream_in: ReadStream<'a>) -> ReadStream<'a> {
    if config.shift_reads[0] == 0 && config.shift_reads[1] == 0 {
        return stream_in;
    }

    let output_stream = async_stream::stream! {
        while let Some(item) = stream_in.next().await {
            match item {
                Ok(mut r) => {
                    if r.strand == '+' {
                        r.range.from += config.shift_reads[0];
                        r.range.to   += config.shift_reads[0];
                    } else if r.strand == '-' {
                        r.range.from += config.shift_reads[1];
                        r.range.to   += config.shift_reads[1];
                    }

                    r.range.to  -= r.range.from;
                    r.range.from = 0;

                    yield Ok(r);
                },
                Err(e) => yield Err(e),
            }
        }

        log!(config.logger, "Shifted reads (forward strand: {}, reverse strand: {})", config.shift_reads[0], config.shift_reads[1]);
    };

    Box::pin(output_stream)
}

/* -------------------------------------------------------------------------- */

pub fn estimate_fraglen(config: &BamCoverageConfig, filename: &str, genome: &Genome) -> Result<FraglenEstimate, Box<dyn Error>> {

    log!(config.logger, "Reading tags from `{}`", filename);

    let bam_result = BamFile::open(filename, None);
    let mut bam = match bam_result {
        Ok (b)   => b,
        Err(err) => return Err(err),
    };

    // Read the reads
    let reads = Box::pin(bam.reader.read_simple_stream(false, false));

    // First round of filtering
    let reads_1 = filter_single_end(config, true, reads);
    let reads_2 = filter_read_length(config, reads_1);
    let reads_3 = filter_duplicates(config, reads_2);
    let reads_4 = filter_mapq(config, reads_3);

    let mut err_opt = None;
    // Convert stream to iterator and catch errors
    let reads_iter = block_on_stream(reads_4).map_while(|item| {
        match item {
            Ok(read) => Some(read),
            Err(err) => {
                err_opt = Some(err);
                None
            }
        }
    });

    // Estimate fragment length
    log!(config.logger, "Estimating mean fragment length");

    let r = match estimate_fragment_length(reads_iter, genome, 2000, config.fraglen_bin_size, config.fraglen_range) {
        Ok((fraglen, x, y, _n)) => {
            log!(config.logger, "Estimated mean fragment length: {}", fraglen);
            Ok(FraglenEstimate { fraglen: fraglen as usize, x, y })
        },
        Err(err) => Err(err),
    };

    // Check if reading the BAM file yielded an error
    if let Some(err) = err_opt {
        return Err(Box::new(err));
    }
    r
}

/* -------------------------------------------------------------------------- */

pub fn bam_coverage_impl(
    mut config         : BamCoverageConfig,
    filenames_treatment: &Vec<String>,
    filenames_control  : &Vec<String>,
    fraglen_treatment  : &Vec<usize>,
    fraglen_control    : &Vec<usize>,
    genome             : Genome,
) -> Result<SimpleTrack, Box<dyn Error>> {
    // Treatment data
    let mut track1 = SimpleTrack::alloc("treatment".to_string(), genome.clone(), config.bin_size);
    let mut n_treatment = 0;
    let mut n_control = 0;

    for (i, filename) in filenames_treatment.iter().enumerate() {
        let mut err_opt = None;
        let fraglen = fraglen_treatment[i];

        log!(config.logger, "Reading treatment tags from `{}`", filename);
        let mut bam = BamFile::open(filename, None)?;

        let treatment = Box::pin(bam.reader.read_simple_stream(!config.paired_as_single_end, config.paired_end_strand_specific));
        // First round of filtering
        let treatment = filter_paired_end(&config, treatment);
        let treatment = filter_single_end(&config, false, treatment);
        let treatment = filter_paired_as_single_end(&config, treatment);
        let treatment = filter_read_length(&config, treatment);
        let treatment = filter_duplicates(&config, treatment);
        let treatment = filter_mapq(&config, treatment);
        // Second round of filtering
        let treatment = filter_strand(&config, treatment);
        let treatment = shift_reads(&config, treatment);

        let treatment_iter = block_on_stream(treatment).map_while(|item| {
            match item {
                Ok(read) => Some(read),
                Err(err) => {
                    err_opt = Some(err);
                    None
                }
            }
        });
    
        n_treatment += GenericMutableTrack::wrap(&mut track1).add_reads(treatment_iter, fraglen, &config.binning_method);

        if let Some(err) = err_opt {
            return Err(Box::new(err));
        }
    }

    // Normalization for treatment
    if config.normalize_track == "rpkm" {
        log!(config.logger, "Normalizing treatment track (rpkm)");
        let c = 1_000_000.0 / (n_treatment as f64 * config.bin_size as f64);
        GenericMutableTrack::wrap(&mut track1).map(|_name, _i, x| c * x)?;
        config.pseudocounts[0] *= c;
    }

    if config.normalize_track == "cpm" {
        log!(config.logger, "Normalizing treatment track (cpm)");
        let c = 1_000_000.0 / n_treatment as f64;
        GenericMutableTrack::wrap(&mut track1).map(|_name, _i, x| c * x)?;
        config.pseudocounts[0] *= c;
    }

    if !filenames_control.is_empty() {
        // Control data
        let mut track2 = SimpleTrack::alloc("control".to_string(), genome.clone(), config.bin_size);

        for (i, filename) in filenames_control.iter().enumerate() {
            let mut err_opt = None;
            let fraglen = fraglen_control[i];

            log!(config.logger, "Reading control tags from `{}`", filename);
            let mut bam = BamFile::open(filename, None)?;
            let control = Box::pin(bam.reader.read_simple_stream(!config.paired_as_single_end, config.paired_end_strand_specific));

            // First round of filtering
            let control = filter_paired_end(&config, control);
            let control = filter_single_end(&config, false, control);
            let control = filter_paired_as_single_end(&config, control);
            let control = filter_read_length(&config, control);
            let control = filter_duplicates(&config, control);
            let control = filter_mapq(&config, control);
            // Second round of filtering
            let control = filter_strand(&config, control);
            let control = shift_reads(&config, control);

            let control_iter = block_on_stream(control).map_while(|item| {
                match item {
                    Ok(read) => Some(read),
                    Err(err) => {
                        err_opt = Some(err);
                        None
                    }
                }
            });
    
            n_control += GenericMutableTrack::wrap(&mut track2).add_reads(control_iter, fraglen, &config.binning_method);

            if let Some(err) = err_opt {
                return Err(Box::new(err));
            }
        }

        // Normalization for control
        if config.normalize_track == "rpkm" {
            log!(config.logger, "Normalizing control track (rpkm)");
            let c = 1_000_000.0 / (n_control as f64 * config.bin_size as f64);
            GenericMutableTrack::wrap(&mut track2).map(|_name, _i, x| c * x)?;
            config.pseudocounts[1] *= c;
        }

        if config.normalize_track == "cpm" {
            log!(config.logger, "Normalizing control track (cpm)");
            let c = 1_000_000.0 / n_control as f64;
            GenericMutableTrack::wrap(&mut track2).map(|_name, _i, x| c * x)?;
            config.pseudocounts[1] *= c;
        }

        if config.smoothen_control {
            GenericMutableTrack::wrap(&mut track2).smoothen(config.smoothen_min, config.smoothen_sizes.clone())?;
        }

        log!(config.logger, "Combining treatment and control tracks...");
        GenericMutableTrack::wrap(&mut track1).normalize(&track2, config.pseudocounts[0], config.pseudocounts[1], config.log_scale)?;
    } else {
        // No control data
        if config.pseudocounts[0] != 0.0 {
            log!(config.logger, "Adding pseudocount `{}`", config.pseudocounts[0]);
            GenericMutableTrack::wrap(&mut track1).map(|_name, _i, x| x + config.pseudocounts[0])?;
        }
        if config.log_scale {
            log!(config.logger, "Log-transforming data");
            GenericMutableTrack::wrap(&mut track1).map(|_name, _i, x| x.ln())?;
        }
    }

    // Filtering chromosomes
    if config.remove_filtered_chroms {
        if !config.filter_chroms.is_empty() {
            log!(config.logger, "Removing chromosomes `{}`", config.filter_chroms.join(", "));
            track1.filter_genome(|name, _length| {
                !config.filter_chroms.contains(&name.to_string())
            });
        }
    } else {
        if !config.filter_chroms.is_empty() {
            log!(config.logger, "Removing all reads from `{}`", config.filter_chroms.join(", "));
            for chr in &config.filter_chroms {
                if let Ok(mut s) = track1.get_sequence_mut(chr) {
                    for i in 0..s.n_bins() {
                        s.set_bin(i, 0.0);
                    }
                }
            }
        }
    }

    Ok(track1)
}

/* -------------------------------------------------------------------------- */

pub fn bam_coverage(
    filenames_treatment: &Vec<String>,
    filenames_control  : &Vec<String>,
    fraglen_treatment  : &Vec<Option<usize>>,
    fraglen_control    : &Vec<Option<usize>>,
    options            : Vec<OptionBamCoverage>,
) -> Result<(SimpleTrack, Vec<FraglenEstimate>, Vec<FraglenEstimate>), CoverageError> {

    let mut config = BamCoverageConfig::default();

    // Parse options
    for option in options {
        config.insert_option(option);
    }

    let mut fraglen_treatment = fraglen_treatment.clone();
    let mut fraglen_control   = fraglen_control  .clone();

    // Check fraglen arguments
    if fraglen_treatment.is_empty() {
        fraglen_treatment = vec![None; filenames_treatment.len()];
    }
    if fraglen_control.is_empty() {
        fraglen_control   = vec![None; filenames_control.len()];
    }
    if fraglen_treatment.len() != filenames_treatment.len() {
        let e = ArgumentError(
            format!("Number of provided treatment fragment lengths `{}` does not match number of treatment files `{}`",
                fraglen_treatment.len(), filenames_treatment.len())
        );
        return Err(CoverageError::new_empty(Box::new(e)));
    }
    if fraglen_control.len() != filenames_control.len() {
        let e = ArgumentError(
            format!("Number of provided control fragment lengths `{}` does not match number of control files `{}`",
                fraglen_control.len(), filenames_control.len())
        );
        return Err(CoverageError::new_empty(Box::new(e)));
    }

    // Read genome
    let mut genome: Genome = Genome::default();
    for filename in filenames_treatment.iter().chain(filenames_control.iter()) {
        let g = bam_import_genome(filename).map_err(
            |e| CoverageError::new_empty(e)
        )?;

        if genome.len() == 0 {
            genome = g;
        } else if genome != g {
            let e = ArgumentError("Treatment and control tracks have different genomes".to_string());
            return Err(CoverageError::new_empty(Box::new(e)));
        }
    }

    let mut treatment_fraglen_estimates = vec![FraglenEstimate::default(); filenames_treatment.len()];
    let mut   control_fraglen_estimates = vec![FraglenEstimate::default(); filenames_control  .len()];

    // No fragment length estimation, copy fragment lengths from arguments
    if !config.estimate_fraglen {
        for (i, fraglen) in fraglen_treatment.iter().enumerate() {
            if let Some(k) = *fraglen {
                treatment_fraglen_estimates[i].fraglen = k;
            } else {
                // If no fragment length is provided and also estimation is
                // switched off, then do not extend reads
                treatment_fraglen_estimates[i].fraglen = 0;
            }
        }
        for (i, fraglen) in fraglen_control.iter().enumerate() {
            if let Some(k) = *fraglen {
                control_fraglen_estimates[i].fraglen = k;
            } else {
                // If no fragment length is provided and also estimation is
                // switched off, then do not extend reads
                control_fraglen_estimates[i].fraglen = 0;
            }
        }
    }
    // Fragment length estimation
    else {

        for (i, filename) in filenames_treatment.iter().enumerate() {

            match fraglen_treatment[i] {
                Some(k) => treatment_fraglen_estimates[i].fraglen = k,
                None    => {
                    let estimate = estimate_fraglen(&config, filename, &genome).map_err(
                        |e| CoverageError::new_empty(e)
                    )?;

                    treatment_fraglen_estimates[i] = estimate;
                }
            }
        }

        for (i, filename) in filenames_control.iter().enumerate() {

            match fraglen_control[i] {
                Some(k) => control_fraglen_estimates[i].fraglen = k,
                None    => {
                    let estimate = estimate_fraglen(&config, filename, &genome).map_err(
                        |e| CoverageError::new_empty(e)
                    )?;

                    control_fraglen_estimates[i] = estimate;
                }
            }
        }
    }
    let fraglen_treatment_arg : Vec<usize> = treatment_fraglen_estimates.iter().map(|x| x.fraglen as usize).collect();
    let fraglen_control_arg   : Vec<usize> =   control_fraglen_estimates.iter().map(|x| x.fraglen as usize).collect();

    let track = bam_coverage_impl(config, filenames_treatment, filenames_control, &fraglen_treatment_arg, &fraglen_control_arg, genome).map_err(
        |e| CoverageError::new(e, treatment_fraglen_estimates.clone(), control_fraglen_estimates.clone())
    )?;

    Ok((track, treatment_fraglen_estimates, control_fraglen_estimates))
}

/* Coverage error type
 * -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct CoverageError{
    pub error                      : Box<dyn Error>,
    pub treatment_fraglen_estimates: Vec<FraglenEstimate>,
    pub   control_fraglen_estimates: Vec<FraglenEstimate>,
}

impl CoverageError {

    pub fn new(
        error                      : Box<dyn Error>,
        treatment_fraglen_estimates: Vec<FraglenEstimate>,
        control_fraglen_estimates  : Vec<FraglenEstimate>
    ) -> CoverageError {
        CoverageError{error, treatment_fraglen_estimates, control_fraglen_estimates}
    }

    pub fn new_empty(error: Box<dyn Error>) -> CoverageError {
        CoverageError{error, treatment_fraglen_estimates: vec![], control_fraglen_estimates: vec![]}
    }

}

impl fmt::Display for CoverageError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.error)
    }
}

impl std::error::Error for CoverageError {}
