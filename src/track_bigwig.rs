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

use std::fs::File;
use std::io::{Read, Seek};
use std::error::Error;

use crate::bigwig::BigWigReader;
use crate::genome::Genome;
use crate::granges_row::GRangesRow;
use crate::track_statistics::BinSummaryStatistics;
use crate::track::{Track, TrackSequence};

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct BigWigTrack<R: Read + Seek> {
    name        : String,
    bin_size    : usize,
    bin_overlap : i32,
    bwr         : BigWigReader<R>,
    bin_sum_stat: BinSummaryStatistics,
    init        : f64,
}

/* -------------------------------------------------------------------------- */

impl<R: Read + Seek> BigWigTrack<R> {

    pub fn new(reader: R, name: String, f: BinSummaryStatistics, bin_size: i32, bin_overlap: i32, init: f64) -> Result<Self, Box<dyn Error>> {

        let bwr      = BigWigReader::new(reader)?;
        let bin_size = if bin_size == 0 {
            bwr.get_bin_size()?
        } else {
            bin_size
        };

        Ok(Self {
            name,
            bin_size,
            bin_overlap,
            bwr,
            bin_sum_stat: f,
            init,
        })
    }

    pub fn open(&mut self, filename: &str, name: &str, f: BinSummaryStatistics, bin_size: i32, bin_overlap: i32, init: f64) -> Result<(), Box<dyn Error>> {
        let file = File::open(filename)?;
        let lazy_track = LazyTrack::new(file, name.to_string(), f, bin_size, bin_overlap, init)?;
        self.lazy_track = lazy_track;
        Ok(())
    }

}

/* -------------------------------------------------------------------------- */

impl<R: Read + Seek> Track for BigWigTrack<R> {

    // Access methods
    fn get_bin_size(&self) -> usize {
        self.bin_size
    }

    fn get_name(&self) -> String {
        self.name.clone()
    }

    fn get_seq_names(&self) -> Vec<String> {
        self.bwr.genome().seqnames.clone()
    }

    fn get_genome(&self) -> &Genome {
        self.bwr.genome()
    }

    fn get_sequence(&self, query: &str) -> Result<TrackSequence, Box<dyn Error>> {
        let (seq, bin_size) = self.bwr.query_sequence(query, self.bin_sum_stat, self.bin_size, self.bin_overlap, self.init)?;
        Ok(TrackSequence { seq, bin_size })
    }

    fn get_slice(&self, r: &GRangesRow) -> Result<Vec<f64>, Box<dyn Error>> {
        let bin_from = r.range().from / self.bin_size;
        let bin_to   = r.range().to   / self.bin_size;

        let mut seq = vec![0.0; (bin_to - bin_from) as usize];
        
        for item in self.bwr.query(&r.seqname, r.range().from, r.range().to, self.bin_size) {
            if let Err(err) = item {
                return Err(Box::new(err));
            }
            let s = item.unwrap();

            for i in (s.data.from..s.data.to).step_by(self.bin_size as usize) {
                let j = (i - r.range().from as i32) / (self.bin_size as i32);
                if (j as usize) < seq.len() {
                    seq[j as usize] = s.data.statistics.sum / s.data.statistics.valid;
                }
            }
        }
        Ok(seq)
    }

}
