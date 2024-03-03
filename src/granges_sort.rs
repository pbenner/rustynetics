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
use std::cmp::Ordering;

use crate::granges::GRanges;
use crate::error::Error;

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn sort(&self, name: &str, reverse: bool) -> Result<Self, Error> {
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
