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

    pub fn get_grange(&self) -> GRange {
        GRange{
            seqname: self.seqname.clone(),
            range  : self.range,
            strand : self.strand
        }
    }

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
