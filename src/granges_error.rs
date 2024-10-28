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

/* Generic error types
* -------------------------------------------------------------------------- */

/// A custom error type representing a missing column in the `GRanges` or `Meta` object.
///
/// `MissingColumn` is used to indicate when a required column is not found within a `GRanges`
/// or `Meta` object, allowing for more informative error handling and debugging.
///
/// # Fields
///
/// - `0`: A `String` that holds the name of the missing column.
///
/// # Example
///
/// ```rust
/// use rustynetics::granges_error::MissingColumn;
///
/// let missing_column = MissingColumn(String::from("gene_id"));
/// println!("{}", missing_column); // Outputs: "GRanges/Meta object is missing a column named `gene_id`"
/// ```
#[derive(Debug)]
pub struct MissingColumn(pub String);

impl fmt::Display for MissingColumn {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GRanges/Meta object is missing a column named `{}`", self.0)
    }
}

impl std::error::Error for MissingColumn {}
