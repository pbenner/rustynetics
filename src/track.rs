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

use crate::genome::Genome;
use crate::granges_row::GRangesRow;

/* -------------------------------------------------------------------------- */

pub struct TrackSequence {
    sequence: Vec<f64>,
    bin_size: usize,
}

impl TrackSequence {

    pub fn new(sequence: &Vec<f64>, bin_size: usize) -> TrackSequence {
        TrackSequence {
            sequence: sequence.clone(),
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

pub struct TrackMutableSequence {
    track_sequence: TrackSequence,
}

impl TrackMutableSequence {

    pub fn new(sequence: &Vec<f64>, bin_size: usize) -> TrackMutableSequence {
        TrackMutableSequence {
            track_sequence: TrackSequence::new(sequence, bin_size)
        }
    }

    pub fn set(&mut self, i: usize, v: f64) {
        self.track_sequence.sequence[i / self.track_sequence.bin_size] = v;
    }

    pub fn set_bin(&mut self, i: usize, v: f64) {
        self.track_sequence.sequence[i] = v;
    }
}

/* -------------------------------------------------------------------------- */

pub trait Track {
    fn get_name(&self) -> String;
    fn get_bin_size(&self) -> usize;
    fn get_sequence(&self, seqname: &str) -> Result<TrackSequence, String>;
    fn get_genome(&self) -> &Genome;
    fn get_seq_names(&self) -> Vec<String>;
    fn get_slice(&self, r: &GRangesRow) -> Result<Vec<f64>, String>;
}

pub trait MutableTrack: Track {
    fn get_mutable_sequence(&self, seqname: &str) -> Result<TrackMutableSequence, String>;
}

/* -------------------------------------------------------------------------- */

pub struct GenericTrack {
    track: Box<dyn Track>,
}

pub struct GenericMutableTrack {
    mutable_track: Box<dyn MutableTrack>,
}
