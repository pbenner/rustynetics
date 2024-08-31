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

use crate::granges::GRanges;
use crate::meta::MetaData;
use crate::track_generic::GenericTrack;

/* -------------------------------------------------------------------------- */

impl<'a> GenericTrack<'a> {

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
                    from.push(c_from);
                    to.push(c_to);
                    values.push(c_val);

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

        r.meta.add_meta(name, MetaData::FloatArray(values))?;

        Ok(r)
    }

}
