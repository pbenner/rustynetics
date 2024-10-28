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

/* -------------------------------------------------------------------------- */

use std::fmt;

use crate::range::Range;
use crate::granges::GRanges;

/* -------------------------------------------------------------------------- */

/// Represents a single genomic range with associated sequence name, range, and strand.
///
/// `GRange` encapsulates a genomic interval defined by a chromosome or sequence name (`seqname`),
/// a `Range` struct for start and end coordinates, and a `strand` character to denote direction.
///
/// # Fields
/// - `seqname`: Name of the chromosome or sequence.
/// - `range`: Interval of the range, as a `Range` object with `from` and `to` positions.
/// - `strand`: Character representing the strand ('+', '-', or '.' for no strand).
pub struct GRange {
    pub seqname : String,
    pub range   : Range,
    pub strand  : char,
}

/* -------------------------------------------------------------------------- */

impl GRange {
    /// Creates a new `GRange` object with specified sequence name, range, and strand.
    ///
    /// # Arguments
    /// - `seqname`: The chromosome or sequence name for the range.
    /// - `from`: The starting coordinate of the range.
    /// - `to`: The ending coordinate of the range.
    /// - `strand`: Character representing the strand, which can be `+`, `-`, or `.`.
    ///
    /// # Returns
    /// A new `GRange` instance.
    pub fn new(seqname : String, from : usize, to : usize, strand : char) -> Self {
        GRange {
            seqname,
            range: Range::new(from, to),
            strand,
        }
    }
}

/* -------------------------------------------------------------------------- */

/// A row within a `GRanges` structure, providing access to a specific indexed row in `GRanges`.
///
/// `GRangesRow` allows read-only access to a particular row within the `GRanges` container, which
/// contains sequence name, range, and strand information for genomic ranges.
///
/// # Fields
/// - `granges`: Reference to the `GRanges` structure containing the data.
/// - `row`: Index of the row within the `GRanges` structure to access.
pub struct GRangesRow<'a> {
    granges: &'a GRanges,
    row    : usize,
}

/* -------------------------------------------------------------------------- */

impl<'a> GRangesRow<'a> {
    /// Creates a new `GRangesRow` referencing a specific row in a `GRanges` object.
    ///
    /// # Arguments
    /// - `granges`: Reference to the `GRanges` object to access.
    /// - `row`: Row index within the `GRanges` structure.
    ///
    /// # Returns
    /// A `GRangesRow` referencing the specified row in `GRanges`.
    pub fn new(granges: &'a GRanges, row: usize) -> Self {
        GRangesRow { granges, row }
    }

    /// Returns the sequence name (`seqname`) of the genomic range in this row.
    pub fn seqname(&self) -> &String {
        &self.granges.seqnames[self.row]
    }

    /// Returns the range coordinates (`range`) of the genomic range in this row.
    pub fn range(&self) -> &Range {
        &self.granges.ranges[self.row]
    }

    /// Returns the strand of the genomic range in this row.
    pub fn strand(&self) -> char {
        self.granges.strand[self.row]
    }
}

/* -------------------------------------------------------------------------- */

impl<'a> fmt::Display for GRangesRow<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "GRangesRow(seqname={}, range=({}, {}), strand={})",
            self.granges.seqnames[self.row],
            self.granges.ranges  [self.row].from,
            self.granges.ranges  [self.row].to,
            self.granges.strand  [self.row] as char
        )
    }
}
