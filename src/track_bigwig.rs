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

use std::io::{Read, Seek, Write};
use std::{rc::Rc, cell::RefCell};
use std::error::Error;
use std::fs::File;

use crate::bbi::{BBI_MAX_ZOOM_LEVELS, BBI_RES_INCREMENT};
use crate::bigwig::{BigWigReader, BigWigWriter, BigWigParameters, OptionBigWig};
use crate::bigwig::BigWigFile;
use crate::genome::Genome;
use crate::granges_row::GRangesRow;
use crate::netfile::NetFile;
use crate::track_statistics::BinSummaryStatistics;
use crate::track::{Track, TrackSequence};
use crate::track_generic::GenericTrack;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct BigWigTrack<R: Read + Seek> {
    name        : String,
    bin_size    : usize,
    bin_overlap : usize,
    // BigWigReader must be mutable, while we want to implement Track, which
    // only allows &self references for most methods. Wrap into RefCell
    bwr         : RefCell<BigWigReader<R>>,
    bin_sum_stat: BinSummaryStatistics,
    genome      : Genome,
    init        : f64,
}

/* -------------------------------------------------------------------------- */

impl BigWigTrack<NetFile> {

    pub fn new(filename : &str, name: String, f: BinSummaryStatistics, bin_size: usize, bin_overlap: usize, init: f64) -> Result<Self, Box<dyn Error>> {

        let mut bwr  = BigWigFile::new_reader(filename)?;
        let genome   = bwr.genome().clone();
        let bin_size = if bin_size == 0 {
            bwr.get_bin_size()?
        } else {
            bin_size
        };

        Ok(Self {
            name,
            bin_size,
            bin_overlap,
            bwr :RefCell::new(bwr),
            bin_sum_stat: f,
            genome,
            init,
        })
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
        self.bwr.borrow().genome().seqnames.clone()
    }

    fn get_genome(&self) -> &Genome {
        &self.genome
    }

    fn get_sequence(&self, query: &str) -> Result<TrackSequence, Box<dyn Error>> {
        let (seq, bin_size) = self.bwr.borrow_mut().query_sequence(query, self.bin_sum_stat, self.bin_size, self.bin_overlap, self.init)?;
        let rc_seq = Rc::new(RefCell::new(seq));
        Ok(TrackSequence::new(rc_seq, bin_size))
    }

    fn get_slice(&self, r: &GRangesRow) -> Result<Vec<f64>, Box<dyn Error>> {
        let bin_from = r.range().from / self.bin_size;
        let bin_to   = r.range().to   / self.bin_size;

        let mut seq = vec![0.0; (bin_to - bin_from) as usize];
        
        for item in self.bwr.borrow_mut().query(&r.seqname(), r.range().from, r.range().to, self.bin_size) {
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

/* -------------------------------------------------------------------------- */

impl<'a> GenericTrack<'a> {

    fn bigwig_automatic_reduction_levels(&self, items_per_slot : usize) -> Vec<i32> {

        let c = (BBI_RES_INCREMENT as usize) * self.track.get_bin_size();
        let mut n = Vec::new();
        let mut l = 0;

        // Get length of longest track
        for &length in &self.track.get_genome().lengths {
            if length / self.track.get_bin_size() > l {
                l = length / self.track.get_bin_size();
            }
        }

        // Initial zoom level
        let mut r = std::cmp::max(100, c);

        // Compute number of zoom levels
        while n.len() <= BBI_MAX_ZOOM_LEVELS {
            if l / r > items_per_slot {
                n.push(r as i32);
                r *= c;
            } else {
                break;
            }
        }

        n
    }

    pub fn write_bigwig<W: Write + Seek>(&self, writer: &mut W, mut parameters_arg: Vec<OptionBigWig>) -> Result<(), Box<dyn Error>> {

        // We need idems_per_slot to compute default zoom levels
        let mut items_per_slot = BigWigParameters::default().items_per_slot;
        // Variable for checking if reduction levels are given
        let mut has_reduction_levels = false;

        for p in &parameters_arg {
            if let OptionBigWig::ItemsPerSlot(x) = *p {
                items_per_slot = x;
            }
            if let OptionBigWig::ReductionLevels(_) = *p {
                has_reduction_levels = true;
            }
        }
        // If no reduction levels are given, compute a default set of zoom levels
        if has_reduction_levels == false {
            parameters_arg.push(OptionBigWig::ReductionLevels(self.bigwig_automatic_reduction_levels(items_per_slot)));
        }

        // Create new BigWig writer
        let mut bww = BigWigWriter::new(writer, self.track.get_genome().clone(), parameters_arg)?;

        // Write data
        for name in self.track.get_seq_names() {
            let sequence = self.track.get_sequence(&name)?;
            bww.write(&name, &sequence.clone_as_vec(), self.track.get_bin_size())?;
        }

        bww.write_index()?;

        // Write zoomed data
        for (i, &reduction_level) in bww.parameters().reduction_levels.clone().iter().enumerate() {
            bww.start_zoom_data(i)?;
            for name in self.track.get_seq_names() {
                let sequence = self.track.get_sequence(&name)?;
                bww.write_zoom(&name, &sequence.clone_as_vec(), self.track.get_bin_size(), reduction_level as usize, i)?;
            }
            bww.write_index_zoom(i)?;
        }

        bww.close()?;

        Ok(())
    }

    pub fn export_bigwig(&self, filename: &str, params: Vec<OptionBigWig>) -> Result<(), Box<dyn Error>> {
        let mut file = File::create(filename)?;
        self.write_bigwig(&mut file, params)?;
        Ok(())
    }

}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::fs;

    use crate::track_bigwig::{BigWigFile, OptionBigWig};
    use crate::track_generic::GenericTrack;
    use crate::track_simple::SimpleTrack;
    use crate::genome::Genome;

    #[test]
    fn test_track_bigwig_1() {

        let filename = "tests/test_bigwig_tmp.bw";
        let nan = f64::NAN;

        let seq_1 = vec![0.0,0.0,0.0,nan,4.5,5.6,0.0,7.8,8.9,0.0];
        let seq_2 = vec![0.1,1.2,2.3,3.4,4.5,5.6,0.0,0.0,8.9,9.0,0.1,1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.0];
        let seq_3 = vec![nan,nan,nan,nan,4.5,5.6,nan,nan,nan,nan];

        let sequences = vec![seq_1.clone(), seq_2.clone(), seq_3.clone()];
        let seqnames  = vec!["test1", "test2", "test3"].into_iter().map(|x| { x.to_string() }).collect();
        let lengths   = vec![100, 200, 100];
        let genome    = Genome::new(seqnames, lengths);

        let track = SimpleTrack::new("track_name".to_string(), sequences, genome, 10).unwrap();

        let params = vec![
            OptionBigWig::ReductionLevels(vec![20])
        ];

        if let Err(e) = GenericTrack::wrap(&track).export_bigwig(filename, params) {
            println!("{}", e);
        }

        let result = BigWigFile::new_reader(filename);

        assert!(result.is_ok());

        if let Ok(mut bw) = result {

            assert_eq!(bw.query("test1", 0, 100, 10).count(),  9);
            assert_eq!(bw.query("test2", 0, 200, 10).count(), 20);
            assert_eq!(bw.query("test3", 0, 100, 10).count(),  2);

            for item in bw.query("test1", 0, 100, 10) {

                let result = item.unwrap();
                let i      = result.data.from as usize / 10;

                assert_eq!(result.data_type, 3);
                assert_eq!(result.data.from, (i as i32)*10);
                assert_eq!(result.data.to  , (i as i32)*10+10);

                assert!((result.data.statistics.sum - seq_1[i]).abs() < 1e-4);

            }

            for (i, item) in bw.query("test2", 0, 200, 10).enumerate() {

                let result = item.unwrap();

                assert_eq!(result.data_type, 3);
                assert_eq!(result.data.from, (i as i32)*10);
                assert_eq!(result.data.to  , (i as i32)*10+10);

                assert!((result.data.statistics.sum - seq_2[i]).abs() < 1e-4);

            }

            for item in bw.query("test3", 0, 100, 10) {

                let result = item.unwrap();
                let i      = result.data.from as usize / 10;

                assert_eq!(result.data_type, 2);

                assert!((result.data.statistics.sum - seq_3[i]).abs() < 1e-4);

            }

            for item in bw.query("test1", 0, 100, 20) {

                let result  = item.unwrap();
                let i       = result.data.from as usize / 10;
                let mut val = 0.0;
                let mut sum = 0.0;

                if !seq_1[i].is_nan() {
                    sum += seq_1[i+0]; val += 1.0;
                }
                if !seq_1[i+1].is_nan() {
                    sum += seq_1[i+1]; val += 1.0;
                }

                assert_eq!(result.data_type, 1);

                assert!((result.data.statistics.valid - val).abs() < 1e-4);
                assert!((result.data.statistics.sum   - sum).abs() < 1e-4);

            }

            for item in bw.query("test2", 0, 200, 20) {

                let result = item.unwrap();
                let i      = result.data.from as usize / 10;

                assert_eq!(result.data_type, 1);

                assert!((result.data.statistics.valid - 2.0).abs() < 1e-4);
                assert!((result.data.statistics.sum - (seq_2[i] + seq_2[i+1])).abs() < 1e-4);

            }

            for item in bw.query("test3", 0, 100, 20) {

                let result = item.unwrap();
                let i      = result.data.from as usize / 10;

                assert_eq!(result.data_type, 1);

                assert!((result.data.statistics.valid - 2.0).abs() < 1e-4);
                assert!((result.data.statistics.sum - (seq_3[i] + seq_3[i+1])).abs() < 1e-4);

            }

        }

        assert!(fs::remove_file(filename).is_ok());

    }
}
