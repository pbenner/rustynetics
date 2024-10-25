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
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, Default)]
pub struct Genome {
    pub seqnames: Vec<String>,
    pub lengths : Vec<usize>
}

/* -------------------------------------------------------------------------- */

impl Genome {
    pub fn new(seqnames: Vec<String>, lengths: Vec<usize>) -> Self {
        if seqnames.len() != lengths.len() {
            panic!("NewGenome(): invalid parameters");
        }
        Genome { seqnames, lengths }
    }

    pub fn len(&self) -> usize {
        self.seqnames.len()
    }

    pub fn get_idx(&self, seqname: &str) -> Option<usize> {
        for (i, sn) in self.seqnames.iter().enumerate() {
            if seqname == sn {
                return Some(i);
            }
        }
        None
    }

    pub fn seq_length(&self, seqname: &str) -> Result<usize, String> {
        for (i, sn) in self.seqnames.iter().enumerate() {
            if seqname == sn {
                return Ok(self.lengths[i]);
            }
        }
        Err(format!("sequence `{}` not found in genome", seqname))
    }

    pub fn sum_seq_length(&self) -> usize {
        self.lengths.iter().sum()
    }

    pub fn add_sequence(&mut self, seqname: String, length: usize) -> Result<usize, String> {
        if self.seqnames.contains(&seqname) {
            Err(format!("sequence `{}` already exists", seqname))
        } else {
            self.seqnames.push(seqname);
            self.lengths .push(length);
            Ok(self.len() - 1)
        }
    }

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

    pub fn import<P: AsRef<Path>>(&mut self, filename: P) -> io::Result<()> {
        let file = File::open(filename.as_ref())?;
        self.read(file).map_err(|e| io::Error::new(io::ErrorKind::Other, format!("reading genome from `{:?}` failed: {}", filename.as_ref(), e)))
    }

    pub fn iter(&self) -> impl Iterator<Item = (&String, &usize)> {
        self.seqnames.iter().zip(self.lengths.iter())
    }
}

/* -------------------------------------------------------------------------- */

impl PartialEq for Genome {
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
