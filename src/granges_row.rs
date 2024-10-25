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

pub struct GRange {
    pub seqname : String,
    pub range   : Range,
    pub strand  : char,
}

/* -------------------------------------------------------------------------- */

impl GRange {
    pub fn new(seqname : String, from : usize, to : usize, strand : char) -> Self {
        GRange {
            seqname,
            range: Range::new(from, to),
            strand,
        }
    }
}

/* -------------------------------------------------------------------------- */

pub struct GRangesRow<'a> {
    granges: &'a GRanges,
    row    : usize,
}

/* -------------------------------------------------------------------------- */

impl<'a> GRangesRow<'a> {
    pub fn new(granges: &'a GRanges, row: usize) -> Self {
        GRangesRow { granges, row }
    }

    pub fn seqname(&self) -> &String {
        &self.granges.seqnames[self.row]
    }

    pub fn range(&self) -> &Range {
        &self.granges.ranges[self.row]
    }

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
