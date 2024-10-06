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
use std::io::Write;
use std::sync::Arc;
use std::io;

use core::pin::Pin;
use futures::{Stream, StreamExt};
use async_stream::stream;

use crate::reads;

/* -------------------------------------------------------------------------- */

// Define an enum to represent all the possible option values
pub enum OptionBamCoverage {
    Logger(Arc<dyn Write>),
    BinningMethod(String),
    BinSize(i64),
    BinOverlap(i64),
    NormalizeTrack(String),
    ShiftReads([i64; 2]),
    PairedAsSingleEnd(bool),
    PairedEndStrandSpecific(bool),
    LogScale(bool),
    Pseudocounts([f64; 2]),
    EstimateFraglen(bool),
    FraglenRange([i64; 2]),
    FraglenBinSize(i64),
    FilterChroms(Vec<String>),
    RemoveFilteredChroms(bool),
    FilterMapQ(i64),
    FilterReadLengths([i64; 2]),
    FilterDuplicates(bool),
    FilterStrand(u8),
    FilterPairedEnd(bool),
    FilterSingleEnd(bool),
    SmoothenControl(bool),
    SmoothenSizes(Vec<i64>),
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
    pub logger: Arc<dyn Write>,
    pub binning_method: String,
    pub bin_size: i64,
    pub bin_overlap: i64,
    pub normalize_track: String,
    pub shift_reads: [i64; 2],
    pub paired_as_single_end: bool,
    pub paired_end_strand_specific: bool,
    pub log_scale: bool,
    pub pseudocounts: [f64; 2],
    pub estimate_fraglen: bool,
    pub fraglen_range: [i64; 2],
    pub fraglen_bin_size: i64,
    pub filter_chroms: Vec<String>,
    pub filter_map_q: i64,
    pub filter_read_lengths: [i64; 2],
    pub filter_duplicates: bool,
    pub filter_strand: u8,
    pub filter_paired_end: bool,
    pub filter_single_end: bool,
    pub remove_filtered_chroms: bool,
    pub smoothen_control: bool,
    pub smoothen_sizes: Vec<i64>,
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
            OptionBamCoverage::FilterMapQ(map_q) => {
                self.filter_map_q = map_q;
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
            logger: Arc::new(io::sink()), // Discard logger output
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
            fraglen_range: [-1, -1],
            fraglen_bin_size: 10,
            filter_chroms: Vec::new(),
            filter_map_q: 0,
            filter_read_lengths: [0, 0],
            filter_duplicates: false,
            filter_strand: b'*',
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

pub struct FraglenEstimate {
    pub fraglen: i64,
    pub x: Vec<i64>,
    pub y: Vec<f64>,
}

/* -------------------------------------------------------------------------- */

pub fn filter_paired_as_single_end(
    config: &BamCoverageConfig,
    mut stream_in: Pin<Box<dyn Stream<Item = io::Result<reads::Read>>>>,
) -> Pin<Box<dyn Stream<Item = io::Result<reads::Read>>>> {

    // If PairedAsSingleEnd is false, return the input stream directly
    if !config.paired_as_single_end {
        return stream_in;
    }

    Box::pin(stream! {
        while let Some(item) = stream_in.next().await {

            match item {
                Ok(mut r) => {
                    r.paired_end = false;
                    yield Ok(r)
                },
                Err(e) => yield Err(e),
            };
        }
    })

}
