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

use std::cmp::{max, min};
use std::fmt;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Copy, Debug)]
pub struct Range {
    pub from: usize,
    pub to  : usize,
}

/* -------------------------------------------------------------------------- */

impl Range {
    pub fn new(from: usize, to: usize) -> Range {
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

/* -------------------------------------------------------------------------- */

impl PartialEq for Range {
    fn eq(&self, other: &Self) -> bool {
        self.from == other.from &&
        self.to   == other.to
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for Range {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad(&format!("[{}, {})", self.from, self.to))
    }
}
