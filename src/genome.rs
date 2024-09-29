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
