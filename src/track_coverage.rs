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

use std::fmt;
use std::error::Error;

use futures::executor::block_on_stream;

use crate::bam::{BamFile, bam_import_genome};
use crate::genome::Genome;
use crate::read_stream::ReadStream;
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
    InitialValue(f64),
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
            OptionBamCoverage::InitialValue(v) => write!(f, "Initial Value: {}", v),
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
    pub initial_value: f64,
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
            OptionBamCoverage::InitialValue(value) => {
                self.initial_value = value;
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
            initial_value: 0.0,
            log_scale: false,
            pseudocounts: [1.0, 1.0],
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
    let reads = ReadStream::filter_single_end (reads, Some(&config.logger), true);
    let reads = ReadStream::filter_read_length(reads, Some(&config.logger), &config.filter_read_lengths);
    let reads = ReadStream::filter_duplicates (reads, Some(&config.logger),  config.filter_duplicates);
    let reads = ReadStream::filter_mapq       (reads, Some(&config.logger),  config.filter_mapq);

    let mut err_opt = None;
    // Convert stream to iterator and catch errors
    let reads_iter = block_on_stream(reads).map_while(|item| {
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
    filenames_treatment: &Vec<&str>,
    filenames_control  : &Vec<&str>,
    fraglen_treatment  : &Vec<usize>,
    fraglen_control    : &Vec<usize>,
    genome             : Genome,
) -> Result<SimpleTrack, Box<dyn Error>> {
    // Treatment data
    let mut track1 = SimpleTrack::alloc("treatment".to_string(), genome.clone(), config.initial_value, config.bin_size);
    let mut n_treatment = 0;
    let mut n_control = 0;

    for (i, filename) in filenames_treatment.iter().enumerate() {
        let mut err_opt = None;
        let fraglen = fraglen_treatment[i];

        log!(config.logger, "Reading treatment tags from `{}`", filename);
        let mut bam = BamFile::open(filename, None)?;

        let treatment = Box::pin(bam.reader.read_simple_stream(!config.paired_as_single_end, config.paired_end_strand_specific));
        // First round of filtering
        let treatment = ReadStream::filter_paired_end   (treatment, Some(&config.logger),  config.filter_paired_end);
        let treatment = ReadStream::filter_single_end   (treatment, Some(&config.logger),  config.filter_single_end);
        let treatment = ReadStream::paired_as_single_end(treatment, Some(&config.logger),  config.paired_as_single_end);
        let treatment = ReadStream::filter_read_length  (treatment, Some(&config.logger), &config.filter_read_lengths);
        let treatment = ReadStream::filter_duplicates   (treatment, Some(&config.logger),  config.filter_duplicates);
        let treatment = ReadStream::filter_mapq         (treatment, Some(&config.logger),  config.filter_mapq);
        // Second round of filtering
        let treatment = ReadStream::filter_strand       (treatment, Some(&config.logger),  config.filter_strand);
        let treatment = ReadStream::shift_reads         (treatment, Some(&config.logger), &config.shift_reads);

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
        let mut track2 = SimpleTrack::alloc("control".to_string(), genome.clone(), config.initial_value, config.bin_size);

        for (i, filename) in filenames_control.iter().enumerate() {
            let mut err_opt = None;
            let fraglen = fraglen_control[i];

            log!(config.logger, "Reading control tags from `{}`", filename);
            let mut bam = BamFile::open(filename, None)?;
            let control = Box::pin(bam.reader.read_simple_stream(!config.paired_as_single_end, config.paired_end_strand_specific));

            // First round of filtering
            let control = ReadStream::filter_paired_end   (control, Some(&config.logger),  config.filter_paired_end);
            let control = ReadStream::filter_single_end   (control, Some(&config.logger),  config.filter_single_end);
            let control = ReadStream::paired_as_single_end(control, Some(&config.logger),  config.paired_as_single_end);
            let control = ReadStream::filter_read_length  (control, Some(&config.logger), &config.filter_read_lengths);
            let control = ReadStream::filter_duplicates   (control, Some(&config.logger),  config.filter_duplicates);
            let control = ReadStream::filter_mapq         (control, Some(&config.logger),  config.filter_mapq);
            // Second round of filtering
            let control = ReadStream::filter_strand       (control, Some(&config.logger),  config.filter_strand);
            let control = ReadStream::shift_reads         (control, Some(&config.logger), &config.shift_reads);

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
    filenames_treatment: &Vec<&str>,
    filenames_control  : &Vec<&str>,
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
