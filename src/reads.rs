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

    pub fn extend_read(&self, d: usize) -> Result<Range, Box<dyn Error>> {
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

#[derive(Debug)]
pub struct StrandMissingError(char);

impl fmt::Display for StrandMissingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "strand information is missing for read with strand `{}`", self.0)
    }
}

impl Error for StrandMissingError {}
