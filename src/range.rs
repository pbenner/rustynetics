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

use std::cmp::{max, min};
use std::fmt;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct Range {
    from: i32,
    to: i32,
}

impl Range {
    pub fn new(from: i32, to: i32) -> Range {
        if from > to {
            panic!("NewRange(): invalid range, i.e. from > to (from={}, to={})", from, to);
        }
        Range { from, to }
    }

    pub fn intersection(&self, other: &Range) -> Range {
        let from = max(self.from, other.from);
        let to = min(self.to, other.to);

        if to < from {
            Range::new(from, from)
        } else {
            Range::new(from, to)
        }
    }
}

impl fmt::Display for Range {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{} {})", self.from, self.to)
    }
}
