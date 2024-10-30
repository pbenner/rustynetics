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
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

/* -------------------------------------------------------------------------- */

/// Represents a genomic structure with sequence names and corresponding lengths.
///
/// The `Genome` struct provides a way to store, manage, and query information about
/// a collection of genomic sequences. Each sequence is represented by its name (e.g.,
/// chromosome name) and its length in base pairs. The struct provides various
/// methods to manipulate and retrieve genomic information, including filtering,
/// summing lengths, and reading data from files.
///
/// # Fields
///
/// - `seqnames`: A vector of strings representing the names of the sequences.
/// - `lengths`: A vector of `usize` values representing the lengths of each sequence.
#[derive(Clone, Debug, Default)]
pub struct Genome {
    pub seqnames: Vec<String>,
    pub lengths : Vec<usize>
}

/* -------------------------------------------------------------------------- */

impl Genome {
    /// Creates a new `Genome` instance from sequence names and lengths.
    ///
    /// # Arguments
    ///
    /// - `seqnames`: A vector of sequence names.
    /// - `lengths`: A vector of sequence lengths.
    ///
    /// # Panics
    ///
    /// Panics if the lengths of `seqnames` and `lengths` are not equal.
    ///
    /// # Example
    ///
    /// ```rust
    /// use rustynetics::genome::Genome;
    ///
    /// let genome = Genome::new(vec!["chr1".to_string(), "chr2".to_string()], vec![1000, 2000]);
    /// ```
    pub fn new(seqnames: Vec<String>, lengths: Vec<usize>) -> Self {
        if seqnames.len() != lengths.len() {
            panic!("NewGenome(): invalid parameters");
        }
        Genome { seqnames, lengths }
    }

    /// Returns the number of sequences in the genome.
    pub fn len(&self) -> usize {
        self.seqnames.len()
    }

    /// Finds the index of a sequence by its name.
    ///
    /// # Arguments
    ///
    /// - `seqname`: The name of the sequence to search for.
    ///
    /// # Returns
    ///
    /// An `Option<usize>` containing the index if found, or `None` if the sequence name
    /// is not present.
    pub fn get_idx(&self, seqname: &str) -> Option<usize> {
        for (i, sn) in self.seqnames.iter().enumerate() {
            if seqname == sn {
                return Some(i);
            }
        }
        None
    }

    /// Retrieves the length of a sequence by its name.
    ///
    /// # Arguments
    ///
    /// - `seqname`: The name of the sequence.
    ///
    /// # Returns
    ///
    /// A `Result` containing the sequence length if found, or an error message if not.
    pub fn seq_length(&self, seqname: &str) -> Result<usize, String> {
        for (i, sn) in self.seqnames.iter().enumerate() {
            if seqname == sn {
                return Ok(self.lengths[i]);
            }
        }
        Err(format!("sequence `{}` not found in genome", seqname))
    }

    /// Returns the sum of all sequence lengths in the genome.
    pub fn sum_seq_length(&self) -> usize {
        self.lengths.iter().sum()
    }

    /// Adds a new sequence to the genome.
    ///
    /// # Arguments
    ///
    /// - `seqname`: The name of the new sequence.
    /// - `length`: The length of the new sequence.
    ///
    /// # Returns
    ///
    /// The index of the newly added sequence if successful, or an error if the sequence
    /// already exists.
    pub fn add_sequence(&mut self, seqname: String, length: usize) -> Result<usize, String> {
        if self.seqnames.contains(&seqname) {
            Err(format!("sequence `{}` already exists", seqname))
        } else {
            self.seqnames.push(seqname);
            self.lengths .push(length);
            Ok(self.len() - 1)
        }
    }

    /// Filters the genome based on a predicate function.
    ///
    /// # Arguments
    ///
    /// - `f`: A closure that takes a sequence name and length, returning `true` if the
    ///   sequence should be included in the filtered genome.
    ///
    /// # Returns
    ///
    /// A new `Genome` containing only the sequences that satisfy the predicate.
    pub fn filter<F>(&mut self, f: F) -> Self
    where
        F: Fn(&str, usize) -> bool,
    {
        let mut seqnames = Vec::new();
        let mut lengths  = Vec::new();

        for (seqname, length) in self.seqnames.iter().zip(self.lengths.iter()) {
            if f(seqname, *length) {
                seqnames.push(seqname.clone());
                lengths .push(*length);
            }
        }
        Genome {
            seqnames: seqnames,
            lengths : lengths
        }
    }

    /// Returns an iterator over the genome's sequences.
    ///
    /// Each item is a tuple with the sequence name and its length.
    pub fn iter(&self) -> impl Iterator<Item = (&String, &usize)> {
        self.seqnames.iter().zip(self.lengths.iter())
    }
}

/* -------------------------------------------------------------------------- */

impl Genome {

    /// Reads genome data from a reader.
    ///
    /// # Arguments
    ///
    /// - `reader`: Any object implementing the `Read` trait (e.g., a file, network stream, or other byte source).
    ///
    /// # Expected Format
    ///
    /// The input data should be a whitespace-separated text format, where each line represents a sequence.
    /// - The first item on each line is the sequence name (e.g., chromosome or contig name).
    /// - The second item is the sequence length (an unsigned integer).
    ///
    /// For example:
    /// ```text
    /// chr1 248956422
    /// chr2 242193529
    /// chrX 156040895
    /// ```
    ///
    /// Blank lines are ignored, but any line with fewer than two items will trigger an error.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on successful read or an `io::Error` if the data is malformed or if an I/O error occurs.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - A line does not contain at least two items (sequence name and length).
    /// - The length field cannot be parsed as an unsigned integer.
    pub fn read<R: Read>(&mut self, reader: R) -> io::Result<()> {
        let reader = BufReader::new(reader);
        let mut seqnames = Vec::new();
        let mut lengths  = Vec::new();

        for line in reader.lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 2 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid genome file"));
            }
            let length: usize = fields[1].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            seqnames.push(fields[0].to_string());
            lengths .push(length);
        }
        *self = Genome::new(seqnames, lengths);
        Ok(())
    }

    /// Imports genome data from a file.
    ///
    /// # Arguments
    ///
    /// - `filename`: A reference to a path that represents the file to read from.
    ///
    /// # Expected Format
    ///
    /// The file should follow a whitespace-separated text format, where each line contains:
    /// - The sequence name (e.g., chromosome or contig name).
    /// - The sequence length as an unsigned integer.
    ///
    /// For example:
    /// ```text
    /// chr1 248956422
    /// chr2 242193529
    /// chrX 156040895
    /// ```
    ///
    /// Blank lines are ignored. Any line with fewer than two items will result in an error.
    ///
    /// # Returns
    ///
    /// A result indicating success or an I/O error if the file cannot be read or is incorrectly formatted.
    pub fn import<P: AsRef<Path>>(&mut self, filename: P) -> io::Result<()> {
        let file = File::open(filename.as_ref())?;
        self.read(file).map_err(|e| io::Error::new(io::ErrorKind::Other, format!("reading genome from `{:?}` failed: {}", filename.as_ref(), e)))
    }

}

/* -------------------------------------------------------------------------- */

impl Genome {

    /// Writes the genome data to a writer.
    ///
    /// # Arguments
    ///
    /// - `writer`: Any object implementing the `Write` trait (e.g., a file, network stream, or other byte sink).
    ///
    /// # Format
    ///
    /// Each line contains:
    /// - The sequence name as a string.
    /// - The sequence length as an unsigned integer.
    ///
    /// For example:
    /// ```text
    /// chr1 248956422
    /// chr2 242193529
    /// chrX 156040895
    /// ```
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on successful write or an `io::Error` if there is an I/O issue.
    ///
    /// # Errors
    ///
    /// Returns an error if writing fails.
    pub fn write<W: Write>(&self, mut writer: W) -> io::Result<()> {
        for (seqname, length) in self.iter() {
            writeln!(writer, "{} {}", seqname, length)?;
        }
        Ok(())
    }

    /// Exports the genome data to a file.
    ///
    /// # Arguments
    ///
    /// - `filename`: A reference to a path where the genome data will be saved.
    ///
    /// # Format
    ///
    /// The output format is the same as the `write` method:
    /// - Each line contains the sequence name and sequence length separated by a space.
    ///
    /// For example:
    /// ```text
    /// chr1 248956422
    /// chr2 242193529
    /// chrX 156040895
    /// ```
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on successful export or an `io::Error` if there is an issue writing to the file.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created or written to.
    pub fn export<P: AsRef<Path>>(&self, filename: P) -> io::Result<()> {
        let file = File::create(filename.as_ref())?;
        self.write(file).map_err(|e| io::Error::new(io::ErrorKind::Other, format!("writing genome to `{:?}` failed: {}", filename.as_ref(), e)))
    }

}

/* -------------------------------------------------------------------------- */

impl PartialEq for Genome {
    /// Checks equality between two `Genome` instances.
    ///
    /// # Returns
    ///
    /// `true` if both `Genome` instances have identical sequences and lengths, `false` otherwise.
    fn eq(&self, other: &Self) -> bool {
        if self.len() != other.len() {
            return false
        }
        for (seqname, l1) in self.seqnames.iter().zip(self.lengths.iter()) {
            let l2 = match other.seq_length(seqname) {
                Ok(l2) => l2,
                Err(_) => return false,
            };
            if *l1 != l2 {
                return false
            }
        }
        true
    }
}

impl Eq for Genome {}

/* -------------------------------------------------------------------------- */

impl fmt::Display for Genome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{:<10} {:>10}", "seqnames", "lengths")?;
        for (seqname, length) in self.iter() {
            writeln!(f, "{:<10} {:>10}", seqname, length)?;
        }
        Ok(())
    }
}
