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

use std::fmt;
use std::error::Error;

use crate::genome::Genome;
use crate::granges_row::GRangesRow;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct TrackSequence<'a> {
    sequence: &'a Vec<f64>,
    bin_size: usize,
}

/* -------------------------------------------------------------------------- */

impl<'a> TrackSequence<'a> {

    pub fn new(sequence: &'a Vec<f64>, bin_size: usize) -> Self {
        Self {
            sequence: sequence,
            bin_size: bin_size,
        }
    }

    pub fn at(&self, i: usize) -> f64 {
        self.sequence[i / self.bin_size]
    }

    pub fn at_bin(&self, i: usize) -> f64 {
        self.sequence[i]
    }

    pub fn n_bins(&self) -> usize {
        self.sequence.len()
    }

    pub fn get_bin_size(&self) -> usize {
        self.bin_size
    }

}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct TrackMutableSequence<'a> {
    sequence: &'a mut Vec<f64>,
    bin_size: usize,
}

/* -------------------------------------------------------------------------- */

impl<'a> TrackMutableSequence<'a> {

    pub fn new(sequence: &'a mut Vec<f64>, bin_size: usize) -> Self {
        Self {
            sequence: sequence,
            bin_size: bin_size,
        }
    }

    pub fn at(&self, i: usize) -> f64 {
        self.sequence[i / self.bin_size]
    }

    pub fn at_bin(&self, i: usize) -> f64 {
        self.sequence[i]
    }

    pub fn n_bins(&self) -> usize {
        self.sequence.len()
    }

    pub fn get_bin_size(&self) -> usize {
        self.bin_size
    }

    pub fn set(&mut self, i: usize, v: f64) {
        self.sequence[i / self.bin_size] = v;
    }

    pub fn set_bin(&mut self, i: usize, v: f64) {
        self.sequence[i] = v;
    }
}

/* -------------------------------------------------------------------------- */

pub trait Track {
    fn get_name(&self) -> String;
    fn get_bin_size(&self) -> usize;
    fn get_sequence(&self, seqname: &str) -> Result<TrackSequence, Box<dyn Error>>;
    fn get_genome(&self) -> &Genome;
    fn get_seq_names(&self) -> Vec<String>;
    fn get_slice(&self, r: &GRangesRow) -> Result<Vec<f64>, Box<dyn Error>>;
}

/* -------------------------------------------------------------------------- */

pub trait MutableTrack : Track {
    fn get_sequence_mut(&mut self, seqname: &str) -> Result<TrackMutableSequence, Box<dyn Error>>;
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct SequenceNotFoundError(pub String);

impl fmt::Display for SequenceNotFoundError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Sequence `{}` was not found", self.0)
    }
}

impl Error for SequenceNotFoundError {}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct GenomeMismatchError(pub String);

impl fmt::Display for GenomeMismatchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Sequence `{}` was not found", self.0)
    }
}

impl Error for GenomeMismatchError {}
