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
use std::collections::HashMap;
use std::error::Error;

use crate::range::Range;
use crate::reads::Read;
use crate::track::{MutableTrack, Track};
use crate::utility::{div_int_up, div_int_down};
use crate::utility_cumdist::{CumDist, OrderedFloat};

/* -------------------------------------------------------------------------- */

pub struct GenericTrack<'a> {
    pub track : &'a dyn Track,
}

/* -------------------------------------------------------------------------- */

impl<'a> GenericTrack<'a> {

    pub fn wrap(track : &'a dyn Track) -> Self {
        Self{track}
    }

    pub fn reduce<F>(&self, f: F, x0: f64) -> HashMap<String, f64>
    where
        F: Fn(&str, usize, f64, f64) -> f64,
    {
        let mut result = HashMap::new();
        let bin_size = self.track.get_bin_size();

        for name in self.track.get_seq_names() {
            let sequence = match self.track.get_sequence(&name) {
                Ok(seq) => seq,
                Err(_) => continue,
            };

            if sequence.n_bins() == 0 {
                continue;
            }

            let mut tmp = f(&name, 0, x0, sequence.at_bin(0));

            for i in 1..sequence.n_bins() {
                tmp = f(&name, i * bin_size, tmp, sequence.at_bin(i));
            }

            result.insert(name, tmp);
        }

        result
    }

    pub fn map<F>(
        &self,
        mut f: F
    ) -> Result<(), Box<dyn Error>>
    where
        F: FnMut(&str, usize, f64),
    {
        let bin_size = self.track.get_bin_size();

        for name in self.track.get_seq_names() {

            let seq = self.track.get_sequence(&name)?;
            
            for i in 0..seq.n_bins() {
                // Call the function `f` with name, index and value
                f(&name, i * bin_size, seq.at_bin(i));
            }
        }

        Ok(())
    }

    pub fn window_map<F>(
        &self,
        window_size: usize,
        mut f      : F,
    ) -> Result<(), Box<dyn Error>>
    where
        F: FnMut(&str, usize, &[f64]),
    {
        if window_size == 0 {
            return Err(Box::new(InvalidWindowSizeError));
        }

        let mut v    = vec![f64::NAN; window_size];
        let bin_size = self.track.get_bin_size();

        for name in self.track.get_seq_names() {
            let seq = self.track.get_sequence(&name)?;
            for i in 0..seq.n_bins() {
                for j in 0..window_size {
                    let k = i as isize - (window_size / 2) as isize + j as isize;
                    v[j] = if k < 0 || k >= seq.n_bins() as isize {
                        f64::NAN
                    } else {
                        seq.at_bin(k as usize)
                    };
                }
                f(&name, i * bin_size, &v);
            }
        }

        Ok(())
    }

    pub fn map_list<F>(
        tracks: &[&dyn Track],
        mut f : F,
    ) -> Result<(), Box<dyn Error>>
    where
        F: FnMut(&str, usize, &[f64]) -> f64,
    {
        if tracks.is_empty() {
            return Ok(());
        }
    
        let n        = tracks.len();
        let bin_size = tracks[0].get_bin_size();
        let mut v    = vec![f64::NAN; n];
    
        // Check bin sizes
        for i in 1..n {
            if tracks[0].get_bin_size() != tracks[i].get_bin_size() {
                return Err(Box::new(BinSizeMismatchError));
            }
        }

        for name in tracks[0].get_seq_names() {
            let mut sequences = Vec::new();
            let mut nbins     = None;

            // Collect source sequences
            for (k, t) in tracks.iter().enumerate() {
                if let Ok(seq) = t.get_sequence(&name) {
                    if nbins.is_none() {
                        nbins = Some(seq.n_bins());
                    }
                    if seq.n_bins() != nbins.unwrap() {
                        return Err(Box::new(SequenceLengthMismatchError(format!(
                            "sequence `{}` in track `{}` has invalid length (`{}` instead of `{}`)",
                            name,
                            k,
                            seq.n_bins(),
                            nbins.unwrap()
                        ))));
                    }
                    sequences.push(seq);
                }
            }

            // Reduce length of v if some tracks are missing a sequence
            let v_len = sequences.len();
            v.truncate(v_len);

            // Loop over sequence
            for i in 0..nbins.unwrap() {
                // Copy values to local vector
                for (j, seq) in sequences.iter().enumerate() {
                    v[j] = seq.at_bin(i);
                }
                // Apply function
                f(&name, i * bin_size, &v);
            }
        }
    
        Ok(())
    }

    pub fn window_map_list<F>(
        tracks     : &[&dyn Track],
        window_size: usize,
        mut f      : F,
    ) -> Result<(), Box<dyn Error>>
    where
        F: FnMut(&str, usize, &[Vec<f64>]) -> f64,
    {
        if tracks.is_empty() {
            return Ok(());
        }
        if window_size == 0 {
            return Err(Box::new(InvalidWindowSizeError));
        }
    
        let n        = tracks.len();
        let bin_size = tracks[0].get_bin_size();
        let mut v: Vec<Vec<f64>> = vec![vec![f64::NAN; window_size]; n];
    
        // Check bin sizes
        for i in 1..n {
            if bin_size != tracks[i].get_bin_size() {
                return Err(Box::new(BinSizeMismatchError));
            }
        }

        for name in tracks[0].get_seq_names() {
            let mut sequences = Vec::new();
            let mut nbins     = None;

            // Collect source sequences
            for (k, t) in tracks.iter().enumerate() {
                if let Ok(seq) = t.get_sequence(&name) {
                    if nbins.is_none() {
                        nbins = Some(seq.n_bins());
                    }
                    if seq.n_bins() != nbins.unwrap() {
                        return Err(Box::new(SequenceLengthMismatchError(format!(
                            "sequence `{}` in track `{}` has invalid length (`{}` instead of `{}`)",
                            name,
                            k,
                            seq.n_bins(),
                            nbins.unwrap()
                        ))));
                    }
                    sequences.push(seq);
                }
            }

            // Reduce length of v if some tracks are missing a sequence
            let v_len = sequences.len();
            v.truncate(v_len);

            // Loop over sequence
            for i in 0..nbins.unwrap() {
                // Copy values to local vector
                for (j, seq) in sequences.iter().enumerate() {
                    for k in 0..window_size {
                        let t = i as isize - (window_size / 2) as isize + k as isize;
                        if t < 0 || t >= seq.n_bins() as isize {
                            v[j][k] = f64::NAN;
                        } else {
                            v[j][k] = seq.at_bin(t as usize);
                        }
                    }
                }
                // Apply function
                f(&name, i * bin_size, &v);
            }
        }
    
        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

pub struct GenericMutableTrack<'a> {
    pub track : &'a mut dyn MutableTrack,
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
                let mut v = seq.at_bin(j); if v.is_nan() {
                    v = 0.0;
                }
                seq.set_bin(j, v + 1.0);
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
                let mut v = seq.at_bin(j); if v.is_nan() {
                    v = 0.0;
                }
                let jfrom = std::cmp::max(from, j * bin_size);
                let jto   = std::cmp::min(to, (j + 1) * bin_size);

                seq.set_bin(j, v + (jto - jfrom) as f64 / bin_size as f64);
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
                let mut v = seq.at_bin(j); if v.is_nan() {
                    v = 0.0;
                }
                let jfrom = std::cmp::max(from, j * bin_size);
                let jto   = std::cmp::min(to, (j + 1) * bin_size);
                seq.set_bin(j, v + (jto - jfrom) as f64);
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
    pub fn add_reads(
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

    // Combine treatment and control from a ChIP-seq experiment into a single track.
    // At each genomic location, the number of binned reads from the treatment
    // experiment is divided by the number of control reads. To avoid division by
    // zero, a pseudocount is added to both treatment and control. The parameter
    // d determines the extension of reads.
    pub fn normalize(
        &mut self,
        control  : &dyn Track,
        c1       : f64,
        c2       : f64,
        log_scale: bool,
    ) -> Result<(), Box<dyn Error>> {

        if c1 <= 0.0 || c2 <= 0.0 {
            return Err("pseudocounts must be strictly positive".into());
        }

        for name in self.track.get_seq_names() {
            let mut seq1 = self.track.get_sequence_mut(&name)?;
            let seq2 = match control.get_sequence(&name) {
                Ok(seq) => seq,
                Err(_)  => continue,
            };

            for i in 0..seq1.n_bins() {
                let value = if log_scale {
                    ((seq1.at_bin(i) + c1) / (seq2.at_bin(i) + c2) * c2 / c1).ln()
                } else {
                    (seq1.at_bin(i) + c1) / (seq2.at_bin(i) + c2) * c2 / c1
                };

                seq1.set_bin(i, value);
            }
        }

        Ok(())
    }

    pub fn quantile_normalize_to_counts(&mut self, x: Vec<f64>, y: Vec<usize>) -> Result<(), Box<dyn Error>> {
        let mut map_in: HashMap<OrderedFloat, usize> = HashMap::new();
        let mut map_tr: HashMap<OrderedFloat, f64  > = HashMap::new();

        // Mapping values to count occurrences
        self.map(&mut |_seqname : &str, _position, value : f64| {
            if !value.is_nan() {
                *map_in.entry(OrderedFloat(value)).or_insert(0) += 1;
            }
            value
        })?;

        let dist_ref = CumDist::from_counts(x, y);
        let dist_in  = CumDist::new(map_in);

        if dist_ref.x.is_empty() {
            return Ok(());
        }

        // Set the first value to keep data on the same range
        map_tr.insert(OrderedFloat(dist_in.x[0]), dist_ref.x[0]);

        let mut i = 1;
        let mut j = 1;
        while i < dist_ref.x.len() {
            let p_ref = dist_ref.y[i] as f64 / dist_ref.num() as f64;
            while j < dist_in.x.len() {
                let p_in = dist_in.y[j] as f64 / dist_in.num() as f64;
                if p_in > p_ref {
                    break;
                }
                // Map input x_j to reference x_i
                map_tr.insert(OrderedFloat(dist_in.x[j]), dist_ref.x[i]);
                j += 1;
            }
            i += 1;
        }

        // Applying the transformation
        self.map(&mut |_seqname : &str, _position, value : f64| {
            if value.is_nan() {
                value
            } else {
                *map_tr.get(&OrderedFloat(value)).unwrap_or(&value)
            }
        })?;

        Ok(())
    }

    pub fn quantile_normalize(&mut self, track_ref: &dyn Track) -> Result<(), Box<dyn Error>> {
        let mut map_ref = HashMap::new();
        let mut map_in  = HashMap::new();
        let mut map_tr  = HashMap::new();

        // Map `track_ref` to `map_ref`
        GenericTrack::wrap(track_ref).map(&mut |_ : &str, _, value : f64| {
            if !value.is_nan() {
                *map_ref.entry(OrderedFloat(value)).or_insert(0) += 1;
            }
        })?;

        // Map `self` to `map_in`
        self.map(&mut |_ : &str, _, value : f64| {
            if !value.is_nan() {
                *map_in.entry(OrderedFloat(value)).or_insert(0) += 1;
            }
            value
        })?;

        let dist_ref = CumDist::new(map_ref);
        let dist_in  = CumDist::new(map_in);

        if dist_ref.x.is_empty() {
            return Ok(());
        }

        // Set the first value to keep data on the same range
        map_tr.insert(OrderedFloat(dist_in.x[0]), dist_ref.x[0]);

        let mut j = 1;
        for i in 1..dist_ref.x.len() {
            let p_ref = dist_ref.y[i] as f64 / dist_ref.num() as f64;
            while j < dist_in.x.len() {
                let p_in = dist_in.y[j] as f64 / dist_in.num() as f64;
                if p_in > p_ref {
                    break;
                }
                map_tr.insert(OrderedFloat(dist_in.x[j]), dist_ref.x[i]);
                j += 1;
            }
        }

        // Apply the mapping to normalize `track`
        self.map(&mut |_ : &str, _, value : f64| {
            if value.is_nan() {
                value
            } else {
                *map_tr.get(&OrderedFloat(value)).unwrap_or(&value)
            }
        })?;

        Ok(())
    }

    // Smoothen track data with an adaptive window method. For each region the smallest window
    // size among windowSizes is selected which contains at least minCounts counts. If the
    // minimum number of counts is not reached, the larges window size is selected.
    pub fn smoothen(&mut self, min_counts: f64, window_sizes: Vec<usize>) -> Result<(), Box<dyn Error>> {
        if window_sizes.is_empty() {
            return Ok(());
        }
        
        let mut window_sizes = window_sizes;
        window_sizes.sort(); // Sort window sizes in ascending order
        
        let offset1 = div_int_up  (window_sizes[0] - 1, 2);
        let offset2 = div_int_down(window_sizes[0] - 1, 2);
        let nw      = window_sizes.len(); // Number of window sizes
        
        // Loop over sequences
        for name in self.track.get_seq_names() {

            let mut seq = self.track.get_sequence_mut(&name)?;
            let nbins = seq.n_bins();
            let mut rst = vec![f64::NEG_INFINITY; nbins];
            
            // Loop over sequence bins
            for i in offset1..(nbins - offset2) {
                let mut counts : f64 = f64::NEG_INFINITY;
                let mut wsize  : i64 = -1;
                
                for k in 0..nw {
                    let mut from = i as isize - div_int_up  (window_sizes[k] - 1, 2) as isize;
                    let mut to   = i as isize + div_int_down(window_sizes[k] - 1, 2) as isize;
                    
                    if from < 0 {
                        to += -from;
                    }
                    
                    if to >= nbins as isize {
                        from -= to - (nbins as isize - 1);
                        to = nbins as isize - 1;
                    }
                    
                    let from = std::cmp::max(0, from) as usize;
                    let to   = std::cmp::min(nbins - 1, to as usize);
                    
                    counts = 0.0;
                    for j in from..=to {
                        counts += seq.at_bin(j);
                    }
                    wsize = (to - from + 1) as i64;
                    
                    if counts >= min_counts {
                        break;
                    }
                }
                
                if wsize != -1 {
                    rst[i] = counts / wsize as f64;
                }
            }
            
            for i in 0..nbins {
                seq.set_bin(i, rst[i]);
            }
        }
        
        Ok(())
    }

    pub fn map<F>(&mut self, mut f: F) -> Result<(), Box<dyn Error>>
    where
        F: FnMut(&str, usize, f64) -> f64,
    {
        let bin_size = self.track.get_bin_size();

        for name in self.track.get_seq_names() {

            let mut seq = self.track.get_sequence_mut(&name)?;
            
            for i in 0..seq.n_bins() {
                // Call the function `f` with name, index and value
                let v = f(&name, i * bin_size, seq.at_bin(i));

                seq.set_bin(i, v);
            }
        }

        Ok(())
    }

    pub fn window_map<F>(
        &mut self,
        track      : &dyn Track,
        window_size: usize,
        mut f      : F,
    ) -> Result<(), Box<dyn Error>>
    where
        F: FnMut(&str, usize, &[f64]) -> f64,
    {
        if window_size == 0 {
            return Err(Box::new(InvalidWindowSizeError));
        }

        let mut v    = vec![f64::NAN; window_size];
        let bin_size = self.track.get_bin_size();

        if self.track.get_bin_size() != track.get_bin_size() {
            return Err(Box::new(BinSizeMismatchError));
        }

        for name in self.track.get_seq_names() {
            let mut seq1 = self.track.get_sequence_mut(&name)?;
            let     seq2 = track.get_sequence(&name)?;

            if seq1.n_bins() != seq2.n_bins() {
                return Err(Box::new(SequenceLengthMismatchError(name)));
            }

            for i in 0..seq2.n_bins() {
                for j in 0..window_size {
                    let k = i as isize - (window_size / 2) as isize + j as isize;
                    v[j] = if k < 0 || k >= seq2.n_bins() as isize {
                        f64::NAN
                    } else {
                        seq2.at_bin(k as usize)
                    };
                }
                seq1.set_bin(i, f(&name, i * bin_size, &v));
            }
        }

        Ok(())
    }

    pub fn map_list<F>(
        &mut self,
        tracks: &[&dyn Track],
        mut f : F,
    ) -> Result<(), Box<dyn Error>>
    where
        F: FnMut(&str, usize, &[f64]) -> f64,
    {
        if tracks.is_empty() {
            return Ok(());
        }
    
        let n        = tracks.len();
        let bin_size = self.track.get_bin_size();
        let mut v    = vec![f64::NAN; n];

        // Check bin sizes
        for t in tracks.iter() {
            if bin_size != t.get_bin_size() {
                return Err(Box::new(BinSizeMismatchError));
            }
        }

        for name in self.track.get_seq_names() {

            let mut dst       = self.track.get_sequence_mut(&name)?;
            let mut sequences = Vec::new();

            // Collect source sequences
            for (k, t) in tracks.iter().enumerate() {
                if let Ok(seq) = t.get_sequence(&name) {
                    if seq.n_bins() != dst.n_bins() {
                        return Err(Box::new(SequenceLengthMismatchError(format!(
                            "sequence `{}` in track `{}` has invalid length (`{}` instead of `{}`)",
                            name,
                            k,
                            seq.n_bins(),
                            dst.n_bins()
                        ))));
                    }
                    sequences.push(seq);
                }
            }

            // Reduce length of v if some tracks are missing a sequence
            let v_len = sequences.len();
            v.truncate(v_len);

            // Loop over sequence
            for i in 0..dst.n_bins() {
                // Copy values to local vector
                for (j, seq) in sequences.iter().enumerate() {
                    v[j] = seq.at_bin(i);
                }
                // Apply function
                dst.set_bin(i, f(&name, i * bin_size, &v));
            }
        }
    
        Ok(())
    }

    pub fn window_map_list<F>(
        &mut self,
        tracks     : &[&dyn Track],
        window_size: usize,
        mut f      : F,
    ) -> Result<(), Box<dyn Error>>
    where
        F: FnMut(&str, usize, &[Vec<f64>]) -> f64,
    {
        if tracks.is_empty() {
            return Ok(());
        }
        if window_size == 0 {
            return Err(Box::new(InvalidWindowSizeError));
        }
    
        let n        = tracks.len();
        let bin_size = self.track.get_bin_size();
        let mut v: Vec<Vec<f64>> = vec![vec![f64::NAN; window_size]; n];
    
        // Check bin sizes
        for t in tracks.iter() {
            if bin_size != t.get_bin_size() {
                return Err(Box::new(BinSizeMismatchError));
            }
        }

        for name in self.track.get_seq_names() {

            let mut dst       = self.track.get_sequence_mut(&name)?;
            let mut sequences = Vec::new();

            // Collect source sequences
            for (k, t) in tracks.iter().enumerate() {
                if let Ok(seq) = t.get_sequence(&name) {
                    if seq.n_bins() != dst.n_bins() {
                        return Err(Box::new(SequenceLengthMismatchError(format!(
                            "sequence `{}` in track `{}` has invalid length (`{}` instead of `{}`)",
                            name,
                            k,
                            seq.n_bins(),
                            dst.n_bins()
                        ))));
                    }
                    sequences.push(seq);
                }
            }

            // Reduce length of v if some tracks are missing a sequence
            let v_len = sequences.len();
            v.truncate(v_len);

            // Loop over sequence
            for i in 0..dst.n_bins() {
                // Copy values to local vector
                for (j, seq) in sequences.iter().enumerate() {
                    for k in 0..window_size {
                        let t = i as isize - (window_size / 2) as isize + k as isize;
                        if t < 0 || t >= seq.n_bins() as isize {
                            v[j][k] = f64::NAN;
                        } else {
                            v[j][k] = seq.at_bin(t as usize);
                        }
                    }
                }
                // Apply function
                dst.set_bin(i, f(&name, i * bin_size, &v));
            }
        }
    
        Ok(())
    }

}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct ReadOutOfRangeError(Read);

impl fmt::Display for ReadOutOfRangeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "read {:?} is out of range", self.0)
    }
}

impl Error for ReadOutOfRangeError {}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct InvalidWindowSizeError;

impl fmt::Display for InvalidWindowSizeError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid window size")
    }
}

impl Error for InvalidWindowSizeError {}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct BinSizeMismatchError;

impl fmt::Display for BinSizeMismatchError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "bin sizes do not match")
    }
}
impl Error for BinSizeMismatchError {}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct SequenceLengthMismatchError(String);

impl fmt::Display for SequenceLengthMismatchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Error for SequenceLengthMismatchError {}
