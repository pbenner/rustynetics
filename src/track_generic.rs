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

    fn extend_read(&self, read: &Read, d: usize) -> Result<(usize, usize), Box<dyn Error>> {
        let mut from = read.range.from;
        let mut to = read.range.to;

        if !read.paired_end && d > 0 {
            // Extend read in 3' direction
            match read.strand {
                '+' => {
                    to = from + d;
                }
                '-' => {
                    if d > to {
                        from = 0;
                    } else {
                        from = to.saturating_sub(d); // Ensure `from` doesn't go negative
                    }
                }
                _ => {
                    return Err(Box::new(StrandMissingError(read.strand)));
                }
            }
        }

        Ok((from, to))
    }

    // Add a single read to the track by incrementing the value of each bin that
    // overlaps with the read. Single end reads are extended in 3' direction
    // to have a length of [d]. This is the same as the macs2 `extsize' parameter.
    // Reads are not extended if [d] is zero.
    // The function returns an error if the read's position is out of range
    pub fn add_read(&mut self, read: &Read, d: usize) -> Result<(), Box<dyn Error>> {
        let mut seq = self.track.get_sequence_mut(&read.seqname)?;

        let (from, to) = self.extend_read(&read, d)?;

        let bin_size = self.track.get_bin_size();

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

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
struct StrandMissingError(char);

impl fmt::Display for StrandMissingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "strand information is missing for read with strand `{}`", self.0)
    }
}

impl Error for StrandMissingError {}
