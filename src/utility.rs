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

use std::collections::HashSet;
use std::path::Path;

use num::traits::PrimInt;

/* -------------------------------------------------------------------------- */

/// Removes duplicate integers from a slice, preserving the order of first occurrences.
///
/// This function iterates over a slice of `usize` values, and returns a `Vec<usize>` with duplicate
/// values removed. Only the first occurrence of each unique value is retained in the result, and the
/// order of appearance in the input slice is preserved.
///
/// # Parameters
/// - `s`: A slice of `usize` values from which duplicates will be removed.
///
/// # Returns
/// A `Vec<usize>` containing only the unique values from the input slice `s`, in the order they first appear.
///
/// # Examples
///
/// ```rust,ignore
/// let input = vec![1, 2, 2, 3, 4, 4, 5];
/// let result = remove_duplicates_int(&input);
/// assert_eq!(result, vec![1, 2, 3, 4, 5]);
///
/// let input = vec![10, 10, 20, 30, 20];
/// let result = remove_duplicates_int(&input);
/// assert_eq!(result, vec![10, 20, 30]);
///
/// let input = vec![1, 1, 1, 1];
/// let result = remove_duplicates_int(&input);
/// assert_eq!(result, vec![1]);
/// ```
///
/// # Complexity
/// This function has a time complexity of approximately O(n), where `n` is the length of the input slice,
/// due to the use of a `HashSet` to track unique elements.
///
/// # Note
/// The function requires the `HashSet` from the standard library to track elements seen so far.
pub fn remove_duplicates_int(s: &[usize]) -> Vec<usize> {
    let mut m: HashSet<usize> = HashSet::new();
    let mut r: Vec<usize> = Vec::new();
    for &v in s {
        if m.insert(v) {
            r.push(v);
        }
    }
    r
}

/* -------------------------------------------------------------------------- */

/// Trims trailing whitespace and removes the outermost matching quotes (either single or double) from a string if they exist.
///
/// # Parameters
/// - `input`: A string slice that may contain trailing whitespace and/or outermost quotes.
///
/// # Returns
/// A new `String` with any trailing whitespace removed and the outermost matching quotes (single or double) removed, if they exist.
/// If the outermost characters are not matching quotes, only the trailing whitespace is removed.
///
/// # Examples
///
/// ```rust,ignore
/// let input = "  'example text'   ";
/// let result = trim_and_unquote(input);
/// assert_eq!(result, "example text");
///
/// let input = " \"hello world\" ";
/// let result = trim_and_unquote(input);
/// assert_eq!(result, "hello world");
///
/// let input = "no quotes here   ";
/// let result = trim_and_unquote(input);
/// assert_eq!(result, "no quotes here");
///
/// let input = "'unmatched quotes";
/// let result = trim_and_unquote(input);
/// assert_eq!(result, "'unmatched quotes");
/// ```
pub fn trim_and_unquote(input: &str) -> String {
    // Step 1: Trim trailing whitespace
    let trimmed = input.trim_end();

    // Step 2: Remove outermost quotes if they exist
    if (trimmed.starts_with('"') && trimmed.ends_with('"')) ||
       (trimmed.starts_with('\'') && trimmed.ends_with('\'')) {
        trimmed[1..trimmed.len()-1].to_string()
    } else {
        trimmed.to_string()
    }
}

/* -------------------------------------------------------------------------- */

/// Performs integer division with rounding up.
///
/// Given two integers `a` and `b`, this function calculates `a / b` with rounding up,
/// which ensures that any remainder will result in an additional increment of the quotient.
///
/// # Type Parameters
/// - `T`: A type that implements the `PrimInt` trait, representing a primitive integer type.
///
/// # Parameters
/// - `a`: The dividend.
/// - `b`: The divisor.
///
/// # Returns
/// The result of `a / b`, rounded up to the nearest integer.
///
/// # Panics
/// Panics if `b` is zero, as division by zero is undefined.
///
/// # Examples
///
/// ```rust,ignore
/// let result = div_int_up(7, 3);
/// assert_eq!(result, 3);  // 7 / 3 rounded up is 3
///
/// let result = div_int_up(10, 2);
/// assert_eq!(result, 5);  // 10 / 2 is exactly 5
/// ```
pub fn div_int_up<T : PrimInt>(a: T, b: T) -> T {
    (a + b - T::one()) / b
}

/// Performs integer division with truncation (rounding down).
///
/// This function divides two integers `n` and `d` and rounds down, discarding any remainder,
/// which is the typical behavior of integer division.
///
/// # Type Parameters
/// - `T`: A type that implements the `PrimInt` trait, representing a primitive integer type.
///
/// # Parameters
/// - `n`: The dividend.
/// - `d`: The divisor.
///
/// # Returns
/// The result of `n / d`, rounded down to the nearest integer (truncated).
///
/// # Panics
/// Panics if `d` is zero, as division by zero is undefined.
///
/// # Examples
///
/// ```rust,ignore
/// let result = div_int_down(7, 3);
/// assert_eq!(result, 2);  // 7 / 3 rounded down is 2
///
/// let result = div_int_down(10, 2);
/// assert_eq!(result, 5);  // 10 / 2 is exactly 5
/// ```
pub fn div_int_down<T : PrimInt>(n: T, d: T) -> T {
    n / d
}

/* -------------------------------------------------------------------------- */

/// Checks if a file has a `.gz` extension, typically indicating a gzip-compressed file.
///
/// This function takes a file path and checks whether its extension is `.gz`,
/// commonly used for gzip-compressed files.
///
/// # Type Parameters
/// - `P`: A type that can be referenced as a `Path`, such as `Path` or `PathBuf`.
///
/// # Parameters
/// - `filename`: The file path to check.
///
/// # Returns
/// `true` if the file has a `.gz` extension; `false` otherwise.
///
/// # Examples
///
/// ```rust,ignore
/// let result = is_gzip("file.txt.gz");
/// assert!(result);  // file has a .gz extension
///
/// let result = is_gzip("file.txt");
/// assert!(!result);  // file does not have a .gz extension
/// ```
pub fn is_gzip<P: AsRef<Path>>(filename: P) -> bool {
    filename.as_ref().extension().map_or(false, |ext| ext == "gz")
}
