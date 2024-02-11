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

use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

/* -------------------------------------------------------------------------- */

#[derive(Clone)]
pub struct Genome {
    sequences: HashMap<String, i32>,
}

/* -------------------------------------------------------------------------- */

impl Genome {
    pub fn new(seqnames: Vec<String>, lengths: Vec<i32>) -> Self {
        if seqnames.len() != lengths.len() {
            panic!("NewGenome(): invalid parameters");
        }
        let sequences = seqnames.into_iter().zip(lengths.into_iter()).collect();
        Genome { sequences }
    }

    pub fn length(&self) -> usize {
        self.sequences.len()
    }

    pub fn seq_length(&self, seqname: &str) -> Result<i32, String> {
        self.sequences
            .get(seqname)
            .cloned()
            .ok_or_else(|| format!("sequence `{}` not found in genome", seqname))
    }

    pub fn sum_lengths(&self) -> i32 {
        self.sequences.values().sum()
    }

    pub fn add_sequence(&mut self, seqname: String, length: i32) -> Result<usize, String> {
        if self.sequences.contains_key(&seqname) {
            Err(format!("sequence `{}` already exists", seqname))
        } else {
            self.sequences.insert(seqname, length);
            Ok(self.length() - 1)
        }
    }

    pub fn filter<F>(&self, f: F) -> Self
    where
        F: Fn(&String, &i32) -> bool,
    {
        let filtered_sequences = self.sequences.iter().filter(|&(name, length)| f(name, length)).map(|(name, length)| (name.clone(), *length)).collect();
        Genome {
            sequences: filtered_sequences,
        }
    }

    pub fn equals(&self, other: &Self) -> bool {
        self.sequences == other.sequences
    }

    pub fn read<R: Read>(&mut self, reader: R) -> io::Result<()> {
        let reader = BufReader::new(reader);
        let mut seqnames = Vec::new();
        let mut lengths = Vec::new();

        for line in reader.lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 2 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid genome file"));
            }
            let length: i32 = fields[1].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            seqnames.push(fields[0].to_string());
            lengths.push(length);
        }
        *self = Genome::new(seqnames, lengths);
        Ok(())
    }

    pub fn import<P: AsRef<Path>>(&mut self, filename: P) -> io::Result<()> {
        let file = File::open(filename.as_ref())?;
        self.read(file).map_err(|e| io::Error::new(io::ErrorKind::Other, format!("reading genome from `{:?}` failed: {}", filename.as_ref(), e)))
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for Genome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{:<10} {:>10}", "seqnames", "lengths")?;
        for (seqname, length) in &self.sequences {
            writeln!(f, "{:<10} {:>10}", seqname, length)?;
        }
        Ok(())
    }
}
