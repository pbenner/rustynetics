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

use std::error::Error;
use std::collections::HashMap;
use std::string::String;
use std::{cell::RefCell, rc::Rc};

use crate::genome::Genome;
use crate::granges_row::GRangesRow;

use crate::track::{Track, MutableTrack};
use crate::track::{TrackSequence, TrackMutableSequence, GenomeMismatchError, SequenceNotFoundError};

/* -------------------------------------------------------------------------- */

pub type TMapType = HashMap<String, Rc<RefCell<Vec<f64>>>>;

/* -------------------------------------------------------------------------- */

// A track is a container for experimental data mapped to genomic
// locations. The data is binned in order to reduce memory usage.
// The first position in a sequence is numbered 0.

#[derive(Debug)]
pub struct SimpleTrack {
    pub name    : String,
    pub genome  : Genome,
    pub data    : TMapType,
    pub bin_size: usize,
}

/* -------------------------------------------------------------------------- */

impl SimpleTrack {

    pub fn new(name: String, sequences: Vec<Vec<f64>>, genome: Genome, bin_size: usize) -> Result<Self, Box<dyn Error>> {
        if sequences.len() != genome.len() {
            return Err(Box::new(GenomeMismatchError("number of track sequences does not match genome".to_string())));
        }
        let mut data: TMapType = HashMap::new();
        for (i, sequence) in sequences.iter().enumerate() {
            if sequence.len() != genome.lengths[i] / bin_size {
                return Err(Box::new(GenomeMismatchError("genome has invalid length for the given sequence and binsize".to_string())));
            }
            data.insert(genome.seqnames[i].clone(), Rc::new(RefCell::new(sequence.clone())));
        }
        Ok(SimpleTrack { name, genome, data, bin_size })
    }

    pub fn alloc(name: String, genome: Genome, bin_size: usize) -> Self {
        let mut data: TMapType = HashMap::new();

        for i in 0..genome.len() {
            // By convention, drop the last positions if they do not fully
            // cover the last bin (i.e., round down). This is required by
            // wig-related tools.
            data.insert(
                genome.seqnames[i].clone(),
                Rc::new(RefCell::new(vec![0.0; genome.lengths[i] / bin_size])),
            );
        }
        SimpleTrack { name, genome, data, bin_size }
    }

    pub fn empty(name: String) -> Self {
        let data: TMapType = HashMap::new();
        SimpleTrack { name, genome: Genome::default(), data, bin_size: 0 }
    }

    pub fn shallow_clone(&self) -> Self {
        let name     = self.name.clone();
        let bin_size = self.bin_size;
        let data     = self.data.clone();
        let genome   = self.genome.clone();

        SimpleTrack { name, genome, data, bin_size }
    }

    pub fn index(&self, position: usize) -> usize {
        position / self.bin_size
    }

    pub fn filter_genome<F>(&mut self, f: F)
    where
        F: Fn(&str, usize) -> bool,
    {
        let retain_seqnames: Vec<String> = self
            .data
            .keys()
            .filter(|seqname| {
                let idx = self.genome.seqnames.iter().position(|x| x == *seqname).unwrap();
                f(seqname, self.genome.lengths[idx])
            })
            .cloned()
            .collect();

        // Retain the seqnames in self.data
        self.data.retain(|seqname, _| retain_seqnames.contains(seqname));

        // Perform the filtering on genome
        self.genome.filter(f);
    }
}

/* -------------------------------------------------------------------------- */

impl Clone for SimpleTrack {

    fn clone(&self) -> Self {
        let name = self.name.clone();
        let bin_size = self.bin_size;
        let mut data: TMapType = HashMap::new();
        let genome = self.genome.clone();

        for (name, sequence) in &self.data {
            let t = sequence.clone();
            data.insert(name.clone(), t);
        }
        SimpleTrack { name, genome, data, bin_size }
    }

}

/* -------------------------------------------------------------------------- */

impl Track for SimpleTrack {

    fn get_bin_size(&self) -> usize {
        self.bin_size
    }

    fn get_name(&self) -> String {
        String::from(&self.name)
    }

    fn get_seq_names(&self) -> Vec<String> {
        self.genome.seqnames.clone()
    }

    fn get_genome(&self) -> &Genome {
        &self.genome
    }

    fn get_sequence(&self, query: &str) -> Result<TrackSequence, Box<dyn Error>> {
        match self.data.get(query) {
            Some(seq) => Ok(TrackSequence::new(seq.clone(), self.bin_size)),
            None      => Err(Box::new(SequenceNotFoundError(format!("sequence `{}` not found", query)))),
        }
    }

    fn get_slice(&self, r: &GRangesRow) -> Result<Vec<f64>, Box<dyn Error>> {
        let seq = match self.data.get(r.seqname()) {
            Some(seq) => seq.borrow(),
            None => return Err(Box::new(SequenceNotFoundError(format!("GetSlice(): invalid seqname `{}`", r.seqname())))),
        };

        let from = r.range().from / self.bin_size;
        let to   = r.range().to   / self.bin_size;

        if from >= seq.len() {
            return Ok(vec![]);
        }

        let from = from.max(0);
        let to   = to  .min(seq.len());

        Ok(seq[from..to].to_vec())
    }

}

/* -------------------------------------------------------------------------- */

impl MutableTrack for SimpleTrack {

    fn get_sequence_mut(&mut self, query: &str) -> Result<TrackMutableSequence, Box<dyn Error>> {
        match self.data.get_mut(query) {
            Some(seq) => Ok(TrackMutableSequence::new(seq.clone(), self.bin_size)),
            None      => Err(Box::new(SequenceNotFoundError(format!("sequence `{}` not found", query)))),
        }
    }

}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::track::{Track, MutableTrack};
    use crate::track_simple::SimpleTrack;
    use crate::genome::Genome;

    #[test]
    fn test_simple_track_1() {

        let seq_1 = vec![1.0, 2.0, 3.0, 4.0];

        let sequences  = vec![seq_1];

        let seqnames  = vec!["test1"].into_iter().map(|x| { x.to_string() }).collect();
        let lengths   = vec![400];
        let genome    = Genome::new(seqnames, lengths);

        let mut track = SimpleTrack::new("track_name".to_string(), sequences, genome, 100).unwrap();
        let mut s1    = track.get_sequence_mut("test1").unwrap();

        s1.set(1, 100.0);

        assert_eq!(track.get_sequence("test1").unwrap().at(  0), 100.0);
        assert_eq!(track.get_sequence("test1").unwrap().at(100),   2.0);
        assert_eq!(track.get_sequence("test1").unwrap().at(200),   3.0);

    }
}
