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

use std::io;
use std::str::FromStr;
 
use crate::range::Range;
use crate::granges::GRanges;
 
/* -------------------------------------------------------------------------- */
 
/// A reader for processing `GRanges` data from a tabular format.
///
/// This struct reads data from a table format, interpreting it as genomic ranges. It
/// identifies the relevant columns in the header and reads subsequent lines to populate
/// the associated `GRanges` structure.
pub struct GRangesTableReader {
    col_seqname: i32,
    col_from   : i32,
    col_to     : i32,
    col_strand : i32,
    granges    : GRanges,
}

/* -------------------------------------------------------------------------- */

impl GRangesTableReader {

    /// Creates a new instance of `GRangesTableReader`.
    ///
    /// This constructor initializes the column indices to -1, indicating that they have not yet been assigned,
    /// and creates a default `GRanges` instance to store the read data.
    ///
    /// # Returns
    /// A new instance of `GRangesTableReader`.
    pub fn new() -> Self {
        GRangesTableReader{
            col_seqname: -1,
            col_from   : -1,
            col_to     : -1,
            col_strand : -1,
            granges    : GRanges::default(),
        }
    }

    /// Reads the header line of the input data to determine column positions.
    ///
    /// This method analyzes the given line to set the indices of the columns corresponding
    /// to the sequence names, start, end, and strand. It raises an error if any of the required
    /// columns are missing.
    ///
    /// # Arguments
    /// - `line`: A reference to a string containing the header line.
    ///
    /// # Returns
    /// An `io::Result<()>`, which will be `Ok(())` if the operation succeeds, or an error if a required column is missing.
    pub fn read_header(&mut self, line: &String) -> io::Result<()> {

        let fields: Vec<&str> = line.trim().split_whitespace().collect();
        for (i, field) in fields.iter().enumerate() {
            match *field {
                "seqnames" => self.col_seqname = i as i32,
                "from"     => self.col_from    = i as i32,
                "to"       => self.col_to      = i as i32,
                "strand"   => self.col_strand  = i as i32,
                "start" if self.col_from == -1 => self.col_from = i as i32,
                "end"   if self.col_to   == -1 => self.col_to   = i as i32,
                _ => (),
            }
        }

        if self.col_seqname == -1 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "is missing a seqnames column"));
        }
        if self.col_from == -1 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "is missing a from column"));
        }
        if self.col_to == -1 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "is missing a to column"));
        }
        Ok(())
    }

    /// Reads a line of data and populates the `GRanges` instance.
    ///
    /// This method processes the given line, extracts the relevant fields based on the previously
    /// determined column indices, and adds the data to the `GRanges` instance. It raises an error if the
    /// data is invalid or if parsing fails.
    ///
    /// # Arguments
    /// - `line`: A reference to a string containing the line to read.
    /// - `i`: The index of the line being processed (used for error reporting).
    ///
    /// # Returns
    /// An `io::Result<()>`, which will be `Ok(())` if the operation succeeds, or an error if parsing fails.
    pub fn read_line(&mut self, line: &String, i: i32) -> io::Result<()> {

        let fields: Vec<&str> = line.trim().split_whitespace().collect();
        if fields.len() < self.col_seqname as usize || fields.len() < self.col_from as usize || fields.len() < self.col_to as usize {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid table"));
        }

        let seqname = fields[self.col_seqname as usize].to_string();
        let from    = usize::from_str(fields[self.col_from as usize]).map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("parsing `from' column `{}` failed at line `{}`", self.col_from + 1, i + 1)))?;
        let to      = usize::from_str(fields[self.col_to   as usize]).map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("parsing `to' column `{}` failed at line `{}`"  , self.col_to   + 1, i + 1)))?;

        self.granges.seqnames.push(seqname);
        self.granges.ranges  .push(Range{ from, to });

        if self.col_strand != -1 {
            if fields.len() < self.col_strand as usize {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid table"));
            }
            let strand = fields[self.col_strand as usize].chars().next().unwrap();
            self.granges.strand.push(strand);
        } else {
            self.granges.strand.push('*');
        }
        Ok(())
    }

    /// Copies the data from the internal `GRanges` instance to another `GRanges` instance.
    ///
    /// This method takes a mutable reference to a `GRanges` instance and populates it with the
    /// sequence names, ranges, and strand information from the reader's internal `GRanges`.
    ///
    /// # Arguments
    /// - `granges`: A mutable reference to a `GRanges` instance to which the data will be copied.
    pub fn push(&mut self, granges: &mut GRanges) {
        granges.seqnames = self.granges.seqnames.clone();
        granges.ranges   = self.granges.ranges  .clone();
        granges.strand   = self.granges.strand  .clone();
    }

}
