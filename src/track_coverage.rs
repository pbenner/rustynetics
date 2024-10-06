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

/* -------------------------------------------------------------------------- */

// Define an enum to represent all the possible option values
pub enum OptionTrackCoverage {
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

impl fmt::Display for OptionTrackCoverage {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OptionTrackCoverage::Logger(_) => write!(f, "Logger option"),
            OptionTrackCoverage::BinningMethod(s) => write!(f, "Binning Method: {}", s),
            OptionTrackCoverage::BinSize(size) => write!(f, "Bin Size: {}", size),
            OptionTrackCoverage::BinOverlap(overlap) => write!(f, "Bin Overlap: {}", overlap),
            OptionTrackCoverage::NormalizeTrack(s) => write!(f, "Normalize Track: {}", s),
            OptionTrackCoverage::ShiftReads(arr) => write!(f, "Shift Reads: {:?}", arr),
            OptionTrackCoverage::PairedAsSingleEnd(b) => write!(f, "Paired as Single End: {}", b),
            OptionTrackCoverage::PairedEndStrandSpecific(b) => write!(f, "Paired End Strand Specific: {}", b),
            OptionTrackCoverage::LogScale(b) => write!(f, "Log Scale: {}", b),
            OptionTrackCoverage::Pseudocounts(arr) => write!(f, "Pseudocounts: {:?}", arr),
            OptionTrackCoverage::EstimateFraglen(b) => write!(f, "Estimate Fraglen: {}", b),
            OptionTrackCoverage::FraglenRange(arr) => write!(f, "Fraglen Range: {:?}", arr),
            OptionTrackCoverage::FraglenBinSize(size) => write!(f, "Fraglen Bin Size: {}", size),
            OptionTrackCoverage::FilterChroms(v) => write!(f, "Filter Chroms: {:?}", v),
            OptionTrackCoverage::RemoveFilteredChroms(b) => write!(f, "Remove Filtered Chroms: {}", b),
            OptionTrackCoverage::FilterMapQ(q) => write!(f, "Filter MapQ: {}", q),
            OptionTrackCoverage::FilterReadLengths(arr) => write!(f, "Filter Read Lengths: {:?}", arr),
            OptionTrackCoverage::FilterDuplicates(b) => write!(f, "Filter Duplicates: {}", b),
            OptionTrackCoverage::FilterStrand(strand) => write!(f, "Filter Strand: {}", strand),
            OptionTrackCoverage::FilterPairedEnd(b) => write!(f, "Filter Paired End: {}", b),
            OptionTrackCoverage::FilterSingleEnd(b) => write!(f, "Filter Single End: {}", b),
            OptionTrackCoverage::SmoothenControl(b) => write!(f, "Smoothen Control: {}", b),
            OptionTrackCoverage::SmoothenSizes(v) => write!(f, "Smoothen Sizes: {:?}", v),
            OptionTrackCoverage::SmoothenMin(min) => write!(f, "Smoothen Min: {}", min),
        }
    }
}
