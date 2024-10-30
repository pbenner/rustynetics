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

use crate::range::Range;
use crate::granges_row::GRange;

/* -------------------------------------------------------------------------- */

/// Represents a genomic read with associated metadata.
///
/// The `Read` struct holds information about a single genomic read, including
/// its chromosome (`seqname`), position range, strand orientation, mapping quality (MAPQ),
/// duplication status, and whether it is paired-end.
///
/// # Fields
///
/// - `seqname`: The name of the chromosome or sequence where the read is located.
/// - `range`: A `Range` object indicating the start and end positions of the read.
/// - `strand`: A character representing the read's strand ('+' for forward, '-' for reverse).
/// - `mapq`: Mapping quality (MAPQ) score, typically an integer value indicating the
///   confidence of the read's alignment.
/// - `duplicate`: A boolean flag indicating if the read is marked as a duplicate.
/// - `paired_end`: A boolean flag indicating if the read is part of a paired-end read.
///
/// # Examples
///
/// ```
/// use rustynetics::read::Read;
/// use rustynetics::range::Range;
///
/// let read = Read {
///     seqname: "chr1".to_string(),
///     range: Range::new(100, 150),
///     strand: '+',
///     mapq: 60,
///     duplicate: false,
///     paired_end: true,
/// };
/// println!("{}", read);
/// ```
#[derive(Clone, Debug)]
pub struct Read {
    pub seqname   : String,
    pub range     : Range,
    pub strand    : char,
    pub mapq      : i64,
    pub duplicate : bool,
    pub paired_end: bool,
}

/* -------------------------------------------------------------------------- */

impl Read {

    /// Converts the `Read` to a `GRange` object.
    ///
    /// The `GRange` object contains only the sequence name, range, and strand,
    /// omitting other fields specific to reads.
    pub fn get_grange(&self) -> GRange {
        GRange{
            seqname: self.seqname.clone(),
            range  : self.range,
            strand : self.strand
        }
    }

    /// Extends the read's range in the 3' direction if the read is not paired-end.
    ///
    /// Depending on the strand orientation, this method extends the `to` position
    /// if the strand is '+', or the `from` position if the strand is '-'. For negative strand
    /// reads, the start (`from`) is adjusted to avoid going negative.
    ///
    /// # Parameters
    ///
    /// - `d`: The extension distance in base pairs.
    ///
    /// # Errors
    ///
    /// Returns a `StrandMissingError` if the strand is not valid ('+' or '-').
    pub fn extend(&self, d: usize) -> Result<Range, Box<dyn Error>> {
        let mut from = self.range.from;
        let mut to   = self.range.to;

        if !self.paired_end && d > 0 {
            // Extend read in 3' direction
            match self.strand {
                '+' => {
                    to = from + d;
                }
                '-' => {
                    if d > to {
                        from = 0;
                    } else {
                        from = to.saturating_sub(d); // Ensure `from` doesn't go negative
                    }
                }
                _ => {
                    return Err(Box::new(StrandMissingError(self.strand)));
                }
            }
        }

        Ok(Range::new(from, to))
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for Read {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "READ: {}:[{}-{}):{}, MAPQ={}, DUPLICATE={}, PAIRED_END={}",
            self.seqname,
            self.range.from,
            self.range.to,
            self.strand,
            self.mapq,
            self.duplicate,
            self.paired_end)
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct StrandMissingError(char);

impl fmt::Display for StrandMissingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "strand information is missing for read with strand `{}`", self.0)
    }
}

impl Error for StrandMissingError {}
