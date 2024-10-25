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

use std::error::Error;

use crate::granges::GRanges;
use crate::meta::MetaData;
use crate::track_generic::GenericTrack;

/* -------------------------------------------------------------------------- */

impl<'a> GenericTrack<'a> {

    /// Converts the track data into genomic ranges (GRanges) format.
    ///
    /// This function iterates over the sequences in the track and identifies contiguous intervals of equal values. For each interval, it records the start and end positions along with the corresponding value. The result is a `GRanges` object that contains the genomic ranges for each sequence along with their associated values.
    ///
    /// # Arguments
    ///
    /// * `name` - A string slice that represents the name of the meta column that contains the corresponding track values.
    ///
    /// # Returns
    ///
    /// Returns a `Result<GRanges, Box<dyn Error>>`. On success, it returns a `GRanges` object containing the genomic intervals and associated values. On failure, it returns an error box containing information about the failure.
    ///
    /// # Errors
    ///
    /// This function will return an error if:
    /// - There is a failure while retrieving the sequence from the track.
    /// - Any other errors related to the metadata or GRanges construction occur.
    pub fn granges(
        &self,
        name : &str
    ) -> Result<GRanges, Box<dyn Error>> {

        let bin_size     = self.track.get_bin_size();
        let mut seqnames = Vec::new();
        let mut from     = Vec::new();
        let mut to       = Vec::new();
        let mut values   = Vec::new();

        for seq_name in self.track.get_seq_names() {

            let sequence = match self.track.get_sequence(&seq_name) {
                Ok(seq) => seq,
                Err(e) => return Err(e),
            };

            if sequence.n_bins() == 0 {
                continue;
            }

            // Current values
            let mut c_from = 0;
            let mut c_to   = bin_size;
            let mut c_val  = sequence.at_bin(0);

            for i in 1..sequence.n_bins() {

                let v = sequence.at_bin(i);

                if v != c_val {

                    seqnames.push(seq_name.clone());
                    from    .push(c_from);
                    to      .push(c_to);
                    values  .push(c_val);

                    c_from = c_to;
                    c_to   = c_from + bin_size;
                    c_val  = v;
                } else {
                    c_to  += bin_size;
                }
            }

            // Append last result
            seqnames.push(seq_name.clone());
            from    .push(c_from);
            to      .push(c_to);
            values  .push(c_val);
        }

        let mut r = GRanges::new(seqnames, from, to, vec![]);

        r.meta.add(name, MetaData::FloatArray(values))?;

        Ok(r)
    }

}
