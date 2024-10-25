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
use std::cmp::Ordering;

use crate::granges::GRanges;
use crate::error::ArgumentError;

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn sorted_indices(&self, name: &str, reverse: bool) -> Result<Vec<usize>, ArgumentError> {
        let mut indices: Vec<usize> = (0..self.num_rows()).collect();
        match name {
            "seqnames" => {
                indices.sort_by(|&i, &j| {
                    let cmp = self.seqnames[i].cmp(&self.seqnames[j]);
                    if reverse {
                        cmp.reverse()
                    } else {
                        cmp
                    }
                });
                Ok(indices)
            }
            "ranges" => {
                indices.sort_by(|&i, &j| {
                    let cmp = self.ranges[i].from.cmp(&self.ranges[j].from);
                    if cmp == Ordering::Equal {
                        self.ranges[i].to.cmp(&self.ranges[j].to)
                    } else if reverse {
                        cmp.reverse()
                    } else {
                        cmp
                    }
                });
                Ok(indices)
            }
            "strand" => {
                indices.sort_by(|&i, &j| {
                    let cmp = self.strand[i].cmp(&self.strand[j]);
                    if reverse {
                        cmp.reverse()
                    } else {
                        cmp
                    }
                });
                Ok(indices)
            }
            _ => Err(ArgumentError("Invalid sort name".to_string())),
        }
    }

    pub fn sort(&self, name: &str, reverse: bool) -> Result<Self, ArgumentError> {
        if name.is_empty() {
            let mut l = GRangesSort::new(self);
            if reverse {
                l.indices.sort_by(|&i, &j| j.cmp(&i));
            } else {
                l.indices.sort();
            }
            let indices = l.indices.clone();
            Ok(self.subset(&indices))
        } else {
            let j = self.sorted_indices(name, reverse)?;
            Ok(self.subset(&j))
        }
    }

}

/* -------------------------------------------------------------------------- */

struct GRangesSort<'a> {
    granges: &'a GRanges,
    indices: Vec<usize>,
}

impl<'a> GRangesSort<'a> {
    fn new(granges: &'a GRanges) -> Self {
        let indices: Vec<usize> = (0..granges.num_rows()).collect();
        GRangesSort { granges, indices }
    }
}

impl<'a> fmt::Display for GRangesSort<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GRangesSort(indices={:?})", self.indices)
    }
}

impl<'a> PartialEq for GRangesSort<'a> {
    fn eq(&self, other: &Self) -> bool {
        self.indices == other.indices
    }
}

impl<'a> Eq for GRangesSort<'a> {}

impl<'a> PartialOrd for GRangesSort<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a> Ord for GRangesSort<'a> {
    fn cmp(&self, other: &Self) -> Ordering {
        let n = self.granges.num_rows();
        for i in 0..n {
            let cmp = self.granges.seqnames[self.indices[i]].cmp(&self.granges.seqnames[other.indices[i]]);
            if cmp != Ordering::Equal {
                return cmp;
            }
        }
        Ordering::Equal
    }
}
