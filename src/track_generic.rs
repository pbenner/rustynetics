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

    /// Adds a single read to the coverage track by incrementing the value of each bin that
    /// overlaps with the read. If the read is single-end, it is extended in the 3' direction
    /// to have a length of `d`. This behavior mimics the `extsize` parameter in macs2.
    /// Reads are not extended if `d` is zero.
    ///
    /// # Arguments
    ///
    /// * `read` - A reference to the read to be added. Contains sequence name and position data.
    /// * `d` - The extension size. If greater than zero, the read will be extended in the 3' direction;
    ///         otherwise, the read is used as is.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the read was successfully added to the track.
    /// If the read's position is out of range, it returns an error.
    ///
    /// # Errors
    ///
    /// This function will return an error if the read's position falls outside of the track's bin range.
    /// Specifically, a `ReadOutOfRangeError` is returned if the read cannot be mapped to any valid bins.
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

    /// Adds a single read to the coverage track by calculating and adding the fraction of overlap
    /// between the read and each bin. If the read is single-end, it is extended in the 3' direction
    /// to have a length of `d`, similar to the `extsize` parameter in macs2. If `d` is zero, the read
    /// is not extended.
    ///
    /// The amount added to each bin is proportional to the fraction of the bin that overlaps with the
    /// read. For example, if a bin partially overlaps with the read, only the fraction of overlap will
    /// be added to the bin's value.
    ///
    /// # Arguments
    ///
    /// * `read` - A reference to the read to be added. Contains sequence name and position data.
    /// * `d` - The extension size. If greater than zero, the read will be extended in the 3' direction;
    ///         otherwise, the read is used as is.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the read was successfully added to the track.
    /// If the read's position is out of range, an error is returned.
    ///
    /// # Errors
    ///
    /// This function returns an error if the read's position is outside of the valid bin range.
    /// Specifically, a `ReadOutOfRangeError` is returned if the read cannot be mapped to any valid bins.
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

    /// Adds a single read to the coverage track by calculating and adding the number of overlapping
    /// nucleotides between the read and each bin. If the read is single-end, it is extended in the 3' direction
    /// to have a length of `d`, similar to the `extsize` parameter in macs2. If `d` is zero, the read is not extended.
    ///
    /// The value added to each bin corresponds to the total number of nucleotides that overlap between the read
    /// and the bin. For instance, if a read spans several bins, the function will increment each bin by the number
    /// of nucleotides that fall within that bin.
    ///
    /// # Arguments
    ///
    /// * `read` - A reference to the read to be added. Contains sequence name and position data.
    /// * `d` - The extension size. If greater than zero, the read will be extended in the 3' direction;
    ///         otherwise, the read is used as is.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the read was successfully added to the track.
    /// If the read's position is out of range, an error is returned.
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

    /// Adds multiple reads to the coverage track, applying an extension in the 3' direction for single-end reads
    /// according to the given extension size `d`. This behavior is similar to the `extsize` parameter in macs2.
    /// If `d` is zero, the reads are not extended.
    ///
    /// The method of how bins are incremented depends on the `method` parameter:
    ///
    /// - `"default"` or `"simple"`: Increments the value of each bin that overlaps the read by 1.
    /// - `"overlap"`: Increments the value of each overlapping bin by the number of nucleotides that overlap the bin.
    /// - `"mean overlap"`: Increments the value of each overlapping bin by the fraction of overlapping nucleotides within the bin.
    ///
    /// # Arguments
    ///
    /// * `reads` - An iterator over reads to be added to the track.
    /// * `d` - The extension size. If greater than zero, each read will be extended in the 3' direction.
    ///         If zero, reads are used as is.
    /// * `method` - A string specifying how the coverage should be computed. Possible values are:
    ///     - `"default"` or `"simple"`: Increments the value of overlapping bins by 1.
    ///     - `"overlap"`: Increments overlapping bins by the number of overlapping nucleotides.
    ///     - `"mean overlap"`: Increments overlapping bins by the fraction of overlapping nucleotides.
    ///
    /// # Returns
    ///
    /// Returns the number of reads successfully added to the track.
    ///
    /// # Panics
    ///
    /// This function will panic if an invalid `method` is provided (i.e., any value other than `"default"`, `"overlap"`,
    /// or `"mean overlap"`).
    ///
    /// # Errors
    ///
    /// The function will return an error internally if a read's position is out of the track's valid bin range.
    ///
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

    /// Combines treatment and control tracks from a ChIP-seq experiment into a single normalized track.
    /// At each genomic location, the number of binned reads from the treatment track is divided by the number
    /// of control reads, and a pseudocount is added to both treatment and control values to avoid division by zero.
    ///
    /// The normalization is performed across all sequences in the track, and the result is stored in the treatment track.
    /// The method supports an optional logarithmic transformation of the normalized values.
    ///
    /// # Arguments
    ///
    /// * `control` - A reference to the control track against which the treatment track will be normalized.
    /// * `c1` - The pseudocount added to the treatment track to avoid division by zero. Must be strictly positive.
    /// * `c2` - The pseudocount added to the control track to avoid division by zero. Must be strictly positive.
    /// * `log_scale` - If `true`, the result is transformed to the natural logarithm of the ratio; otherwise, the raw ratio is used.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the normalization is successful.
    ///
    /// # Errors
    ///
    /// Returns an error if either `c1` or `c2` are non-positive pseudocount values, as pseudocounts must be strictly positive.
    /// Additionally, an error may occur if sequences cannot be retrieved from either the treatment or control track.
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

    /// Smoothens the track data using an adaptive window method. For each bin, the function selects
    /// the smallest window size from the provided `window_sizes` that contains at least `min_counts` counts
    /// within the window. If no window size can satisfy the minimum count, the largest window size is used.
    ///
    /// This method loops over the track's sequences and applies the smoothing operation to each bin,
    /// adjusting the values based on the selected window size.
    ///
    /// # Arguments
    ///
    /// * `min_counts` - The minimum number of counts required in a window for it to be considered valid.
    /// * `window_sizes` - A list of window sizes to select from. The function tries to use the smallest window
    ///                    size that contains at least `min_counts` counts. If none of the windows satisfy
    ///                    this condition, the largest window size is applied.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the smoothing operation is completed successfully.
    ///
    /// # Errors
    ///
    /// Returns an error if there are issues retrieving sequences or updating bins in the track.
    ///
    /// # Panics
    ///
    /// This function does not panic. However, if `window_sizes` is empty, the function will immediately return
    /// without modifying the track.
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

    /// Applies a user-defined function to sliding windows of data from two tracks, modifying the bins of the current track (`self`) based on the values from the provided `track`.
    /// The user-defined function `f` is called on each window of data and should return a new value for the corresponding bin in the current track.
    ///
    /// The function operates over each sequence in the tracks, applying the sliding window approach to each bin. It ensures that the bin sizes of the two tracks match and that
    /// the lengths of their sequences are consistent. If the sequences differ in length or if the bin sizes do not match, an error is returned.
    ///
    /// # Arguments
    ///
    /// * `track` - A reference to the second track from which windowed values are retrieved for comparison or computation.
    /// * `window_size` - The size of the window (number of bins) for the sliding window operation. Each window contains values from `track` centered around the current bin.
    /// * `f` - A user-defined function that takes the sequence name (`&str`), the genomic position of the current bin (`usize`), and the array of values in the window (`&[f64]`).
    ///         The function should return a new value to be set in the current track for that bin.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the windowed mapping is applied successfully to all sequences and bins.
    ///
    /// # Errors
    ///
    /// * Returns an error if `window_size` is zero.
    /// * Returns an error if the bin sizes of the two tracks do not match.
    /// * Returns an error if the sequences in the two tracks have different numbers of bins.
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

    /// Applies a user-defined function to corresponding bins across multiple tracks, modifying the bins of the current track (`self`) based on the values from the provided list of tracks.
    /// The function loops over all bins in all tracks and calls the user-defined function `f` to compute a new value for each bin in the current track.
    ///
    /// This function checks that the bin sizes of all tracks match and ensures that each sequence has the same number of bins across all tracks before applying the function.
    ///
    /// # Arguments
    ///
    /// * `tracks` - A slice of references to the tracks that will be used for the computation. Each track should have the same bin size as the current track.
    /// * `f` - A user-defined function that takes the sequence name (`&str`), the genomic position (`usize`), and an array of bin values (`&[f64]`) from the tracks. The function should return a new value to be set in the current track for that bin.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the function is applied successfully to all sequences and bins.
    ///
    /// # Errors
    ///
    /// * Returns an error if the list of tracks is empty.
    /// * Returns an error if the bin sizes of any track do not match the bin size of the current track.
    /// * Returns an error if any sequence in the provided tracks has a different number of bins than the corresponding sequence in the current track.
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

    /// Applies a user-defined function to a sliding window of bins across multiple tracks, modifying the bins of the current track (`self`).
    /// The function processes each bin by considering the values from a window of surrounding bins in each track. The size of the sliding window is specified by `window_size`.
    /// The user-defined function `f` computes a new value for each bin in the current track, based on the corresponding windows from the other tracks.
    ///
    /// # Arguments
    ///
    /// * `tracks` - A slice of references to the tracks that will be used for the computation. Each track should have the same bin size as the current track.
    /// * `window_size` - The size of the sliding window used to extract bins from the other tracks for each bin in the current track. Must be greater than 0.
    /// * `f` - A user-defined function that takes the sequence name (`&str`), the genomic position (`usize`), and a slice of windows (`&[Vec<f64>]`) from the tracks. Each window is a `Vec<f64>` representing the values of the bins within the window for a specific track. The function should return a new value to be set in the current track for that bin.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the function is applied successfully to all sequences and bins.
    ///
    /// # Errors
    ///
    /// * Returns an error if `window_size` is 0.
    /// * Returns an error if the bin sizes of any track do not match the bin size of the current track.
    /// * Returns an error if any sequence in the provided tracks has a different number of bins than the corresponding sequence in the current track.
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
