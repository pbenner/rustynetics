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

use crate::range::Range;
use crate::reads::Read;
use crate::track::{MutableTrack, Track};

/* -------------------------------------------------------------------------- */

pub struct GenericTrack<'a> {
    track : &'a dyn Track,
}

pub struct GenericMutableTrack<'a> {
    track : &'a mut dyn MutableTrack,
}

/* -------------------------------------------------------------------------- */

impl<'a> GenericMutableTrack<'a> {

    pub fn wrap(track : &'a mut dyn MutableTrack) -> Self {
        Self{track}
    }

    // Add a single read to the track by incrementing the value of each bin that
    // overlaps with the read. Single end reads are extended in 3' direction
    // to have a length of [d]. This is the same as the macs2 `extsize' parameter.
    // Reads are not extended if [d] is zero.
    // The function returns an error if the read's position is out of range
    pub fn add_read(&mut self, read: &Read, d: usize) -> Result<(), Box<dyn Error>> {

        let bin_size        = self.track.get_bin_size();
        let mut seq         = self.track.get_sequence_mut(&read.seqname)?;
        let Range{from, to} = read.extend(d)?;

        if from / bin_size >= seq.n_bins() {
            return Err(Box::new(ReadOutOfRangeError(read.clone())));
        }

        for j in (from / bin_size)..=(to - 1) / bin_size {
            if j >= seq.n_bins() {
                break;
            } else {
                let new_value = seq.at_bin(j) + 1.0;
                seq.set_bin(j, new_value);
            }
        }

        Ok(())
    }

    // Add a single read to the track by adding the fraction of overlap between
    // the read and each bin. Single end reads are extended in 3' direction
    // to have a length of [d]. This is the same as the macs2 `extsize' parameter.
    // Reads are not extended if [d] is zero.
    // The function returns an error if the read's position is out of range
    fn add_read_mean_overlap(&mut self, read: &Read, d: usize) -> Result<(), Box<dyn Error>> {

        let bin_size        = self.track.get_bin_size();
        let mut seq         = self.track.get_sequence_mut(&read.seqname)?;
        let Range{from, to} = read.extend(d)?;

        if from / bin_size >= seq.n_bins() {
            return Err(Box::new(ReadOutOfRangeError(read.clone())));
        }

        for j in (from / bin_size)..=((to - 1) / bin_size) {
            if j >= seq.n_bins() {
                break;
            } else {
                let jfrom = std::cmp::max(from, j * bin_size);
                let jto   = std::cmp::min(to, (j + 1) * bin_size);
                seq.set_bin(j, seq.at_bin(j) + (jto - jfrom) as f64 / bin_size as f64);
            }
        }

        Ok(())
    }

    // Add a single read to the track by adding the number of overlapping nucleotides
    // between the read and each bin. Single end reads are extended in 3' direction
    // to have a length of [d]. This is the same as the macs2 `extsize' parameter.
    // Reads are not extended if [d] is zero.
    // The function returns an error if the read's position is out of range
    fn add_read_overlap(&mut self, read: &Read, d: usize) -> Result<(), Box<dyn Error>> {

        let bin_size        = self.track.get_bin_size();
        let mut seq         = self.track.get_sequence_mut(&read.seqname)?;
        let Range{from, to} = read.extend(d)?;

        if from / bin_size >= seq.n_bins() {
            return Err(Box::new(ReadOutOfRangeError(read.clone())));
        }

        for j in (from / bin_size)..=((to - 1) / bin_size) {
            if j >= seq.n_bins() {
                break;
            } else {
                let jfrom = std::cmp::max(from, j * bin_size);
                let jto   = std::cmp::min(to, (j + 1) * bin_size);
                seq.set_bin(j, seq.at_bin(j) + (jto - jfrom) as f64);
            }
        }

        Ok(())
    }

    // Add reads to track. All single end reads are extended in 3' direction
    // to have a length of [d]. This is the same as the macs2 `extsize' parameter.
    // Reads are not extended if [d] is zero.
    // If [method] is "default", the value of each bin that overlaps the read
    // is incremented. If [method] is "overlap", each bin that overlaps the read is
    // incremented by the number of overlapping nucleotides. If [method] is "mean
    // overlap", each bin that overlaps the read is incremented by the fraction
    // of overlapping nucleotides within the bin.
    // The function returns an error if the read's position is out of range
    fn add_reads(
        &mut self,
        reads : impl Iterator<Item = Read>,
        d     : usize,
        method: &str) -> usize
    {
        let mut n = 0;

        match method {
            "" | "simple" | "default" => {
                for read in reads {
                    if self.add_read(&read, d).is_ok() {
                        n += 1;
                    }
                }
            }
            "mean overlap" => {
                for read in reads {
                    if self.add_read_mean_overlap(&read, d).is_ok() {
                        n += 1;
                    }
                }
            }
            "overlap" => {
                for read in reads {
                    if self.add_read_overlap(&read, d).is_ok() {
                        n += 1;
                    }
                }
            }
            _ => panic!("invalid binning method"),
        }

        n
    }

}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
struct ReadOutOfRangeError(Read);

impl fmt::Display for ReadOutOfRangeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "read {:?} is out of range", self.0)
    }
}

impl Error for ReadOutOfRangeError {}
