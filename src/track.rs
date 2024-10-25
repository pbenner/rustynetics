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
use std::error::Error;
use std::{cell::RefCell, rc::Rc};

use crate::genome::Genome;
use crate::granges_row::GRangesRow;

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct TrackSequence {
    sequence: Rc<RefCell<Vec<f64>>>,
    bin_size: usize,
}

/* -------------------------------------------------------------------------- */

impl TrackSequence {

    pub fn new(sequence: Rc<RefCell<Vec<f64>>>, bin_size: usize) -> Self {
        Self {
            sequence: sequence,
            bin_size: bin_size,
        }
    }

    pub fn clone_as_vec(&self) -> Vec<f64> {
        self.sequence.borrow().clone()
    }

    pub fn at(&self, i: usize) -> f64 {
        self.sequence.borrow()[i / self.bin_size]
    }

    pub fn at_bin(&self, i: usize) -> f64 {
        self.sequence.borrow()[i]
    }

    pub fn n_bins(&self) -> usize {
        self.sequence.borrow().len()
    }

    pub fn get_bin_size(&self) -> usize {
        self.bin_size
    }

}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct TrackMutableSequence {
    sequence: Rc<RefCell<Vec<f64>>>,
    bin_size: usize,
}

/* -------------------------------------------------------------------------- */

impl TrackMutableSequence {

    pub fn new(sequence: Rc<RefCell<Vec<f64>>>, bin_size: usize) -> Self {
        Self {
            sequence: sequence,
            bin_size: bin_size,
        }
    }

    pub fn clone_as_vec(&self) -> Vec<f64> {
        self.sequence.borrow().clone()
    }

    pub fn at(&self, i: usize) -> f64 {
        self.sequence.borrow()[i / self.bin_size]
    }

    pub fn at_bin(&self, i: usize) -> f64 {
        self.sequence.borrow()[i]
    }

    pub fn n_bins(&self) -> usize {
        self.sequence.borrow().len()
    }

    pub fn get_bin_size(&self) -> usize {
        self.bin_size
    }

    pub fn set(&mut self, i: usize, v: f64) {
        self.sequence.borrow_mut()[i / self.bin_size] = v;
    }

    pub fn set_bin(&mut self, i: usize, v: f64) {
        self.sequence.borrow_mut()[i] = v;
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
