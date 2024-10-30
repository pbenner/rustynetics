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

/// A struct representing a closed-open interval with a start (`from`) and end (`to`) position.
///
/// `Range` represents a genomic or numeric range with an inclusive `from` bound and an exclusive `to` bound.
/// This interval can be used for calculating intersections or checking overlaps with other ranges.
///
/// # Fields
///
/// - `from`: The start position of the range (inclusive).
/// - `to`: The end position of the range (exclusive).
///
/// # Examples
///
/// ```
/// use rustynetics::range::Range;
///
/// let range = Range::new(100, 200);
/// println!("{}", range); // Output: [100, 200)
/// ```
#[derive(Clone, Copy, Debug)]
pub struct Range {
    pub from: usize,
    pub to  : usize,
}

/* -------------------------------------------------------------------------- */

impl Range {
    /// Creates a new `Range` with specified `from` and `to` bounds.
    ///
    /// # Parameters
    ///
    /// - `from`: The inclusive start position of the range.
    /// - `to`: The exclusive end position of the range.
    ///
    /// # Panics
    ///
    /// This function panics if `from` is greater than `to`, as it would result in an invalid range.
    ///
    /// # Examples
    ///
    /// ```
    /// use rustynetics::range::Range;
    ///
    /// let range = Range::new(5, 10);
    /// ```
    pub fn new(from: usize, to: usize) -> Range {
        if from > to {
            panic!("NewRange(): invalid range, i.e. from > to (from={}, to={})", from, to);
        }
        Range { from, to }
    }

    /// Computes the intersection of this range with another range.
    ///
    /// The intersection of two ranges is the overlapping interval where both ranges meet.
    /// If there is no overlap, the intersection is a zero-length range at the maximum `from` position.
    ///
    /// # Parameters
    ///
    /// - `other`: A reference to another `Range` with which to calculate the intersection.
    ///
    /// # Returns
    ///
    /// A `Range` representing the intersection of the two ranges. If the ranges do not overlap,
    /// a zero-length range is returned with `from == to`.
    ///
    /// # Examples
    ///
    /// ```
    /// use rustynetics::range::Range;
    ///
    /// let range1 = Range::new(10, 20);
    /// let range2 = Range::new(15, 25);
    /// let intersection = range1.intersection(&range2);
    /// println!("{}", intersection); // Output: [15, 20)
    /// ```
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
    /// Compares two `Range` objects for equality.
    ///
    /// Two `Range` instances are considered equal if they have the same `from` and `to` values.
    ///
    /// # Examples
    ///
    /// ```
    /// use rustynetics::range::Range;
    ///
    /// let range1 = Range::new(5, 10);
    /// let range2 = Range::new(5, 10);
    /// assert_eq!(range1, range2);
    /// ```
    fn eq(&self, other: &Self) -> bool {
        self.from == other.from &&
        self.to   == other.to
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for Range {
    /// Formats the `Range` as a string for display.
    ///
    /// The format used is `[from, to)`, where `from` is inclusive and `to` is exclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use rustynetics::range::Range;
    ///
    /// let range = Range::new(5, 10);
    /// println!("{}", range); // Output: [5, 10)
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad(&format!("[{}, {})", self.from, self.to))
    }
}
