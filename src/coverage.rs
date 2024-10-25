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

use crate::infologger::Logger;

/* -------------------------------------------------------------------------- */

// Define an enum to represent all the possible option values
pub enum OptionCoverage {
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

impl fmt::Display for OptionCoverage {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OptionCoverage::Logger(_) => write!(f, "Logger option"),
            OptionCoverage::BinningMethod(s) => write!(f, "Binning Method: {}", s),
            OptionCoverage::BinSize(size) => write!(f, "Bin Size: {}", size),
            OptionCoverage::BinOverlap(overlap) => write!(f, "Bin Overlap: {}", overlap),
            OptionCoverage::NormalizeTrack(s) => write!(f, "Normalize Track: {}", s),
            OptionCoverage::ShiftReads(arr) => write!(f, "Shift Reads: {:?}", arr),
            OptionCoverage::PairedAsSingleEnd(b) => write!(f, "Paired as Single End: {}", b),
            OptionCoverage::PairedEndStrandSpecific(b) => write!(f, "Paired End Strand Specific: {}", b),
            OptionCoverage::InitialValue(v) => write!(f, "Initial Value: {}", v),
            OptionCoverage::LogScale(b) => write!(f, "Log Scale: {}", b),
            OptionCoverage::Pseudocounts(arr) => write!(f, "Pseudocounts: {:?}", arr),
            OptionCoverage::EstimateFraglen(b) => write!(f, "Estimate Fraglen: {}", b),
            OptionCoverage::FraglenRange(arr) => write!(f, "Fraglen Range: {:?}", arr),
            OptionCoverage::FraglenBinSize(size) => write!(f, "Fraglen Bin Size: {}", size),
            OptionCoverage::FilterChroms(v) => write!(f, "Filter Chroms: {:?}", v),
            OptionCoverage::RemoveFilteredChroms(b) => write!(f, "Remove Filtered Chroms: {}", b),
            OptionCoverage::FilterMapQ(q) => write!(f, "Filter MapQ: {}", q),
            OptionCoverage::FilterReadLengths(arr) => write!(f, "Filter Read Lengths: {:?}", arr),
            OptionCoverage::FilterDuplicates(b) => write!(f, "Filter Duplicates: {}", b),
            OptionCoverage::FilterStrand(strand) => write!(f, "Filter Strand: {}", strand),
            OptionCoverage::FilterPairedEnd(b) => write!(f, "Filter Paired End: {}", b),
            OptionCoverage::FilterSingleEnd(b) => write!(f, "Filter Single End: {}", b),
            OptionCoverage::SmoothenControl(b) => write!(f, "Smoothen Control: {}", b),
            OptionCoverage::SmoothenSizes(v) => write!(f, "Smoothen Sizes: {:?}", v),
            OptionCoverage::SmoothenMin(min) => write!(f, "Smoothen Min: {}", min),
        }
    }
}

/* -------------------------------------------------------------------------- */

// Define the CoverageConfig struct
pub struct CoverageConfig {
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

// Function to insert OptionCoverage into CoverageConfig
impl CoverageConfig {
    pub fn insert_option(&mut self, option: OptionCoverage) {
        match option {
            OptionCoverage::Logger(logger) => {
                self.logger = logger;
            }
            OptionCoverage::BinningMethod(method) => {
                self.binning_method = method;
            }
            OptionCoverage::BinSize(size) => {
                self.bin_size = size;
            }
            OptionCoverage::BinOverlap(overlap) => {
                self.bin_overlap = overlap;
            }
            OptionCoverage::InitialValue(value) => {
                self.initial_value = value;
            }
            OptionCoverage::NormalizeTrack(track) => {
                self.normalize_track = track;
            }
            OptionCoverage::ShiftReads(reads) => {
                self.shift_reads = reads;
            }
            OptionCoverage::PairedAsSingleEnd(paired) => {
                self.paired_as_single_end = paired;
            }
            OptionCoverage::PairedEndStrandSpecific(strand_specific) => {
                self.paired_end_strand_specific = strand_specific;
            }
            OptionCoverage::LogScale(log_scale) => {
                self.log_scale = log_scale;
            }
            OptionCoverage::Pseudocounts(pseudocounts) => {
                self.pseudocounts = pseudocounts;
            }
            OptionCoverage::EstimateFraglen(estimate) => {
                self.estimate_fraglen = estimate;
            }
            OptionCoverage::FraglenRange(range) => {
                self.fraglen_range = range;
            }
            OptionCoverage::FraglenBinSize(size) => {
                self.fraglen_bin_size = size;
            }
            OptionCoverage::FilterChroms(chroms) => {
                self.filter_chroms = chroms;
            }
            OptionCoverage::FilterMapQ(mapq) => {
                self.filter_mapq = mapq;
            }
            OptionCoverage::FilterReadLengths(read_lengths) => {
                self.filter_read_lengths = read_lengths;
            }
            OptionCoverage::FilterDuplicates(duplicates) => {
                self.filter_duplicates = duplicates;
            }
            OptionCoverage::FilterStrand(strand) => {
                self.filter_strand = strand;
            }
            OptionCoverage::FilterPairedEnd(paired_end) => {
                self.filter_paired_end = paired_end;
            }
            OptionCoverage::FilterSingleEnd(single_end) => {
                self.filter_single_end = single_end;
            }
            OptionCoverage::RemoveFilteredChroms(remove) => {
                self.remove_filtered_chroms = remove;
            }
            OptionCoverage::SmoothenControl(smoothen) => {
                self.smoothen_control = smoothen;
            }
            OptionCoverage::SmoothenSizes(sizes) => {
                self.smoothen_sizes = sizes;
            }
            OptionCoverage::SmoothenMin(min) => {
                self.smoothen_min = min;
            }
        }
    }
}

/* -------------------------------------------------------------------------- */

// Function to provide default configuration
impl CoverageConfig {
    pub fn default() -> Self {
        CoverageConfig {
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
 