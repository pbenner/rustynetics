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
use std::collections::HashMap;
use std::error::Error;

use crate::range::Range;
use crate::genome::Genome;
use crate::granges_row::GRangesRow;
use crate::granges_find_overlaps::find_overlaps;
use crate::meta::Meta;
use crate::utility::remove_duplicates_int;

/* -------------------------------------------------------------------------- */

/// A structure representing genomic ranges with associated metadata, 
/// commonly used for handling regions of interest in genomic data processing.
///
/// The `GRanges` struct stores chromosome names (`seqnames`), genomic ranges
/// (`ranges`), strand orientation (`strand`), and additional metadata (`meta`).
///
/// # Fields
/// - `seqnames`: A vector of chromosome names where each entry corresponds
///   to a genomic range in `ranges`.
/// - `ranges`: A vector of `Range` structs representing genomic intervals,
///   each defined by start and end positions.
/// - `strand`: A vector indicating the strand orientation (`+`, `-`, or `*`)
///   for each range. A strand of `'*'` denotes an unspecified or neutral strand.
/// - `meta`: An instance of the `Meta` struct, holding additional user-defined
///   metadata, such as associated gene expression values or other annotations.
///
/// # Examples
/// ```
/// use rustynetics::granges::GRanges;
///
/// // Example of creating a new GRanges instance
/// let granges = GRanges::new(
///     vec!["chr1".to_string(), "chr2".to_string()],
///     vec![100, 200],
///     vec![150, 250],
///     vec!['+', '-']
/// );
/// ```
///
/// # Usage
/// The `GRanges` struct provides several methods for manipulating and querying
/// genomic ranges, such as `subset`, `intersection`, `filter_genome`, and more.
/// Additionally, methods for adjusting range lengths (`set_lengths`) or filtering
/// based on strand orientation are available, facilitating flexible operations 
/// on genomic data.
///
/// # Note
/// The `GRanges` struct is particularly useful for bioinformatics applications
/// where analysis of specific genomic regions or annotations is needed, and 
/// supports various metadata types via the `Meta` struct.
#[derive(Clone, Debug)]
pub struct GRanges {
    pub seqnames: Vec<String>,
    pub ranges  : Vec<Range>,
    pub strand  : Vec<char>,
    pub meta    : Meta,
}

/* -------------------------------------------------------------------------- */

impl Default for GRanges {
    fn default() -> Self {
        let seqnames = vec!["".to_string(); 0];
        let ranges   = vec![Range::new(0, 0); 0];
        let strand   = vec!['*'; 0];
        GRanges {
            seqnames,
            ranges,
            strand,
            meta: Meta::default(),
        }
    }
}

/* -------------------------------------------------------------------------- */

impl GRanges {

    /// Creates a new `GRanges` instance with specified sequence names, start and end positions, and strand information.
    ///
    /// # Arguments
    /// - `seqnames`: Vector of chromosome names.
    /// - `from`: Vector of start positions.
    /// - `to`: Vector of end positions.
    /// - `strand`: Vector of strand orientations (`+`, `-`, or `*`). If empty, defaults to `'*'`.
    ///
    /// # Panics
    /// Panics if the lengths of `seqnames`, `from`, and `to` do not match.
    pub fn new(seqnames: Vec<String>, from: Vec<usize>, to: Vec<usize>, strand: Vec<char>) -> Self {
        let n = seqnames.len();
        if from.len() != n || to.len() != n || (strand.len() != 0 && strand.len() != n) {
            panic!("NewGRanges(): invalid arguments!");
        }
        let strand = if strand.len() == 0 {
            vec!['*'; n]
        } else {
            strand
        };
        let ranges = from
            .iter()
            .zip(to.iter())
            .map(|(&f, &t)| Range::new(f, t))
            .collect();
        GRanges {
            seqnames,
            ranges,
            strand,
            meta: Meta::default(),
        }
    }

    /// Returns the number of rows in the `GRanges` instance.
    pub fn num_rows(&self) -> usize {
        self.ranges.len()
    }

    /// Returns a row of data as a `GRangesRow` at the specified index.
    ///
    /// # Arguments
    /// - `i`: The index of the row.
    ///
    /// # Panics
    /// Panics if the index is out of bounds.
    pub fn row(&self, i: usize) -> GRangesRow {
        GRangesRow::new(self, i)
    }

    /// Appends another `GRanges` instance to the current one, merging their data.
    ///
    /// # Arguments
    /// - `other`: The `GRanges` instance to append.
    ///
    /// # Returns
    /// A new `GRanges` instance containing data from both instances.
    ///
    /// # Errors
    /// Returns an error if metadata cannot be merged.
    pub fn append(&self, other: &GRanges) -> Result<Self, Box<dyn Error>> {
        let mut seqnames = self.seqnames.clone();
        seqnames.extend(other.seqnames.iter().cloned());
        let mut ranges = self.ranges.clone();
        ranges.extend(other.ranges.iter().cloned());
        let mut strand = self.strand.clone();
        strand.extend(other.strand.iter().cloned());
        let meta = self.meta.clone().append(&other.meta)?;
        Ok(GRanges {
            seqnames,
            ranges,
            strand,
            meta,
        })
    }

    /// Removes rows at the specified indices.
    ///
    /// # Arguments
    /// - `indices`: Vector of row indices to remove.
    ///
    /// # Returns
    /// A new `GRanges` instance with specified rows removed.
    pub fn remove(&self, indices: &[usize]) -> Self {
        if indices.is_empty() {
            return self.clone();
        }
        let indices = remove_duplicates_int(indices);
        let mut indices = indices.to_vec();
        indices.sort();
        let n = self.num_rows();
        let m = n - indices.len();
        let mut idx = Vec::with_capacity(m);
        let mut k = 0;
        for i in 0..n {
            while k < indices.len() - 1 && i > indices[k] {
                k += 1;
            }
            if i != indices[k] {
                idx.push(i);
            }
        }
        let result = self.subset(&idx);
        let meta = self.meta.subset(&idx);
        GRanges {
            seqnames: result.seqnames,
            ranges: result.ranges,
            strand: result.strand,
            meta,
        }
    }

    /// Removes rows that overlap with any row in the given `GRanges` instance.
    ///
    /// # Arguments
    /// - `subject`: `GRanges` instance to check for overlaps.
    ///
    /// # Returns
    /// A new `GRanges` instance with overlapping rows removed.
    pub fn remove_overlaps_with(&self, subject: &GRanges) -> Self {
        let (query_hits, _) = find_overlaps(self, subject);
        self.remove(&query_hits)
    }

    /// Keeps only rows that overlap with any row in the given `GRanges` instance.
    ///
    /// # Arguments
    /// - `subject`: `GRanges` instance to check for overlaps.
    ///
    /// # Returns
    /// A new `GRanges` instance with only overlapping rows.
    pub fn keep_overlaps_with(&self, subject: &GRanges) -> Self {
        let (query_hits, _) = find_overlaps(self, subject);
        let query_hits = remove_duplicates_int(&query_hits);
        let mut query_hits = query_hits.to_vec();
        query_hits.sort();
        self.subset(&query_hits)
    }

    /// Subsets the `GRanges` instance based on the specified row indices.
    ///
    /// # Arguments
    /// - `indices`: Vector of row indices to include in the subset.
    ///
    /// # Returns
    /// A new `GRanges` instance containing only the specified rows.
    pub fn subset(&self, indices: &[usize]) -> Self {
        let seqnames = indices.iter().map(|&i| self.seqnames[i].clone()).collect();
        let from     = indices.iter().map(|&i| self.ranges  [i].from   ).collect();
        let to       = indices.iter().map(|&i| self.ranges  [i].to     ).collect();
        let strand   = indices.iter().map(|&i| self.strand  [i]        ).collect();
        let result   = GRanges::new(seqnames, from, to, strand);
        let meta     = self.meta.subset(indices);
        GRanges {
            seqnames: result.seqnames,
            ranges  : result.ranges,
            strand  : result.strand,
            meta,
        }
    }

    /// Slices the `GRanges` instance between the specified row indices.
    ///
    /// # Arguments
    /// - `ifrom`: Start index (inclusive).
    /// - `ito`: End index (exclusive).
    ///
    /// # Returns
    /// A new `GRanges` instance containing rows within the specified range.
    pub fn slice(&self, ifrom: usize, ito: usize) -> Self {
        let ifrom = ifrom.min(self.num_rows());
        let ito   = ito  .min(self.num_rows());
        let seqnames = (ifrom..ito).map(|i| self.seqnames[i].clone()).collect();
        let from     = (ifrom..ito).map(|i| self.ranges  [i].from   ).collect();
        let to       = (ifrom..ito).map(|i| self.ranges  [i].to     ).collect();
        let strand   = (ifrom..ito).map(|i| self.strand  [i]        ).collect();
        let result   = GRanges::new(seqnames, from, to, strand);
        let meta     = self.meta.slice(ifrom, ito);
        GRanges {
            seqnames: result.seqnames,
            ranges: result.ranges,
            strand: result.strand,
            meta,
        }
    }

    /// Computes the intersection between two `GRanges` instances.
    ///
    /// # Arguments
    /// - `s`: The other `GRanges` instance to intersect with.
    ///
    /// # Returns
    /// A new `GRanges` instance containing rows that intersect between `self` and `s`.
    pub fn intersection(&self, s: &GRanges) -> Self {
        let (query_hits, subject_hits) = find_overlaps(self, s);
        let n = query_hits.len();

        let mut seqnames = Vec::with_capacity(n);
        let mut from     = Vec::with_capacity(n);
        let mut to       = Vec::with_capacity(n);
        let mut strand   = Vec::with_capacity(n);

        for i in 0..n {
            let i_q =   query_hits[i];
            let i_s = subject_hits[i];
            let gr  = self.ranges[i_q].intersection(&s.ranges[i_s]);

            seqnames.push(self.seqnames[i_q].clone());
            strand  .push(self.strand  [i_q]);
            from    .push(gr.from);
            to      .push(gr.to  );
        }
        let mut granges = GRanges::new(seqnames, from, to, strand);

        granges.meta = self.meta.subset(&query_hits);
        granges
    }

    /// Filters rows based on their inclusion in a specified genome.
    ///
    /// # Arguments
    /// - `genome`: The genome to filter by.
    ///
    /// # Returns
    /// A new `GRanges` instance containing rows included in the genome.
    pub fn filter_genome(&self, genome: &Genome) -> Self {
        let mut idx = Vec::new();
        let mut seqnames = HashMap::new();
        for i in 0..genome.len() {
            seqnames.insert(genome.seqnames[i].clone(), genome.lengths[i]);
        }
        for i in 0..self.num_rows() {
            let num_rows = seqnames.get(&self.seqnames[i]).cloned();
            if let Some(num_rows) = num_rows {
                if self.ranges[i].to <= num_rows {
                    continue;
                }
            }
            idx.push(i);
        }
        self.remove(&idx)
    }

    /// Filters rows based on a specific strand.
    ///
    /// # Arguments
    /// - `s`: The strand to filter by (`+`, `-`, or `*`).
    ///
    /// # Returns
    /// A new `GRanges` instance containing rows matching the specified strand.
    pub fn filter_strand(&self, s: char) -> Self {
        let mut idx = Vec::new();
        for i in 0..self.num_rows() {
            if self.strand[i] != s {
                idx.push(i);
            }
        }
        self.remove(&idx)
    }

    /// Adjusts the lengths of ranges in the `GRanges` instance based on strand orientation.
    ///
    /// For each row:
    /// - If the strand is `'+'`, the end position (`to`) is set to the start position (`from`) plus `n`.
    /// - If the strand is `'-'`, the start position (`from`) is set to the end position (`to`) minus `n`.
    /// - If the strand is `'*'`, no change is made to the range.
    ///
    /// # Arguments
    /// - `n`: The new length for each range.
    ///
    /// # Returns
    /// A new `GRanges` instance with modified range lengths based on strand orientation.
    ///
    /// # Example
    /// ```
    /// use rustynetics::granges::GRanges;
    ///
    /// let granges = GRanges::new(
    ///     vec!["chr1".to_string(), "chr2".to_string()],
    ///     vec![100, 200],
    ///     vec![150, 250],
    ///     vec!['+', '-']
    /// );
    /// let modified_granges = granges.set_lengths(50);
    /// ```
    /// In this example, the range on `'+'` strand becomes `[100, 150)` and on `'-'` strand becomes `[200, 250)`.
    pub fn set_lengths(&self, n: usize) -> Self {
        let mut s = self.clone();

        for i in 0..s.num_rows() {
            // forward strand
            if s.strand[i] == '+' {
                s.ranges[i].to = s.ranges[i].from + n;
            }
            // reverse strand
            if s.strand[i] == '-' {
                s.ranges[i].from = s.ranges[i].to - n;
            }
            // if strand is '*' do nothing
        }
        s
    }

}

/* -------------------------------------------------------------------------- */

impl PartialEq for GRanges {
    fn eq(&self, other: &Self) -> bool {
        self.seqnames == other.seqnames &&
        self.ranges   == other.ranges   &&
        self.strand   == other.strand   &&
        self.meta     == other.meta
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for GRanges {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad(&format!("{}", self.format_pretty(10).unwrap()))
    }
}
