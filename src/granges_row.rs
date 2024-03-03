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

/* -------------------------------------------------------------------------- */

use std::fmt;

use crate::granges::GRanges;

/* -------------------------------------------------------------------------- */

pub struct GRangesRow<'a> {
    granges: &'a GRanges,
    row    : usize,
}

impl<'a> GRangesRow<'a> {
    pub fn new(granges: &'a GRanges, row: usize) -> Self {
        GRangesRow { granges, row }
    }
}

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
