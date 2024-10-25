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

use std::fs::File;
use std::error::Error;
use std::io::{self, Read, Seek, Write};

use crate::bigwig::{BigWigReader, OptionBigWig};
use crate::track_generic::GenericTrack;
use crate::track_simple::SimpleTrack;
use crate::track_statistics::BinSummaryStatistics;

/* -------------------------------------------------------------------------- */

impl SimpleTrack {

    pub fn read_bigwig<R: Read + Seek>(
        &mut self,
        reader: &mut R,
        name  : &str,
        f     : BinSummaryStatistics,
        mut bin_size: usize,
        bin_overlap: usize,
        init: f64,
    ) -> Result<(), Box<dyn Error>> {

        let mut bwr = BigWigReader::new(reader)?;
        let mut sequences: Vec<Vec<f64>> = Vec::new();
        let genome = bwr.genome().clone();

        for seqname in genome.seqnames {
            let (s, b) = bwr.query_sequence(&seqname, f, bin_size, bin_overlap, init)?;
            if bin_size == 0 {
                bin_size = b;
            }
            sequences.push(s);
        }

        let tmp = SimpleTrack::new(name.to_string(), sequences, bwr.genome().clone(), bin_size)?;
        *self = tmp;

        Ok(())
    }

    pub fn import_bigwig(
        &mut self,
        filename   : &str,
        name       : &str,
        s          : BinSummaryStatistics,
        bin_size   : usize,
        bin_overlap: usize,
        init       : f64,
    ) -> Result<(), Box<dyn Error>> {

        let mut f = File::open(filename)?;

        Ok(self.read_bigwig(&mut f, name, s, bin_size, bin_overlap, init)
            .map_err(|e| {
                io::Error::new(io::ErrorKind::Other, format!("importing bigWig file from `{}` failed: {}", filename, e))
            }
        )?)

    }

    pub fn write_bigwig<W: Write + Seek>(
        &self,
        writer: &mut W,
        params: Vec<OptionBigWig>
    ) -> Result<(), Box<dyn Error>> {

        GenericTrack { track: self }.write_bigwig(writer, params)

    }

    pub fn export_bigwig(
        &self,
        filename: &str,
        params  : Vec<OptionBigWig>
    ) -> Result<(), Box<dyn Error>> {

        let mut f = File::create(filename)?;

        Ok(self.write_bigwig(&mut f, params).map_err(|e| {
            io::Error::new(io::ErrorKind::Other, format!("exporting bigWig file to `{}` failed: {}", filename, e))
        })?)

    }

}
