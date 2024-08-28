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

use std::io::{self, Read, Seek, SeekFrom};
use std::result::Result;
use std::error::Error;

use async_stream::stream;
use futures::executor::block_on_stream;
use futures::executor::BlockingStream;
use futures_core::stream::Stream;
use futures::StreamExt;

use byteorder::{ByteOrder, ReadBytesExt, LittleEndian};

use crate::genome::Genome;
use crate::bbi::BbiQueryType;
use crate::bbi::BbiFile;
use crate::netfile::NetFile;

/* -------------------------------------------------------------------------- */

const BIGWIG_MAGIC : u32 = 0x888FFC26;

/* -------------------------------------------------------------------------- */

pub fn is_bigwig_file(filename: &str) -> Result<bool, Box<dyn Error>> {

    let mut file = NetFile::open(filename)?;

    let magic = file.read_u32::<LittleEndian>()?;

    Ok(BIGWIG_MAGIC == magic)

}

/* -------------------------------------------------------------------------- */

// Struct for BigWigParameters
#[derive(Clone, Debug)]
pub struct BigWigParameters {
    pub block_size      : usize,
    pub items_per_slot  : usize,
    pub reduction_levels: Option<Vec<i32>>,
}

/* -------------------------------------------------------------------------- */

impl Default for BigWigParameters {
    // Default constructor for BigWigParameters
    fn default() -> Self {
        BigWigParameters {
            block_size      : 256,
            items_per_slot  : 1024,
            reduction_levels: None,
        }
    }
}

/* -------------------------------------------------------------------------- */

pub struct BigWigReader<R: Read + Seek> {
    reader: R,
    bwf   : BbiFile,
    genome: Genome,
}

/* -------------------------------------------------------------------------- */

pub type BigWigFile = NetFile;

/* -------------------------------------------------------------------------- */

impl BigWigFile {

    pub fn new_reader(filename : &str) -> Result<BigWigReader<Self>, Box<dyn Error>> {

        let file = NetFile::open(filename)?;

        BigWigReader::new(file)

    }

}

/* -------------------------------------------------------------------------- */

impl<R: Read + Seek> BigWigReader<R> {

    pub fn new(mut reader: R) -> Result<Self, Box<dyn Error>> {
        let bwf = BigWigReader::<R>::open_bwf(&mut reader)?;

        let r = BigWigReader {
            reader: reader,
            bwf   : bwf,
            genome: Genome::default(),
        };

        Ok(r.initialize::<LittleEndian>()?)
    }

    fn open_bwf(reader: &mut R) -> io::Result<BbiFile> {

        let mut bwf = BbiFile::default();

        reader.seek(SeekFrom::Start(0))?;

        bwf.open::<LittleEndian, R>(reader, BIGWIG_MAGIC)?;

        Ok(bwf)
    }

    fn initialize<E: ByteOrder>(mut self) -> io::Result<Self> {

        self.genome.seqnames = vec![String::new(); self.bwf.chrom_data.keys.len()];
        self.genome.lengths  = vec![0; self.bwf.chrom_data.keys.len()];

        for i in 0..self.bwf.chrom_data.keys.len() {
            if self.bwf.chrom_data.values[i].len() != 8 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "invalid chromosome list",
                ));
            }

            let idx = (&self.bwf.chrom_data.values[i][0..4]).read_u32::<E>()? as usize;
             
            if idx >= self.bwf.chrom_data.keys.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "invalid chromosome index",
                ));
            }
            self.genome.seqnames[idx] = String::from_utf8_lossy(&self.bwf.chrom_data.keys[i]).trim_end_matches('\x00').to_string();
            self.genome.lengths [idx] = (&self.bwf.chrom_data.values[i][4..8]).read_u32::<E>()? as usize;
        }

        Ok(self)
    }

    pub fn query_stream<'a>(
        &'a mut self,
        seq_regex: &'a str,
        from     : usize,
        to       : usize,
        bin_size : usize
    ) -> impl Stream<Item = io::Result<BbiQueryType>> + 'a {

        stream! {

            let re = regex::Regex::new(&format!("^{}$", seq_regex)).unwrap();

            for seqname in &self.genome.seqnames {
                if !re.is_match(seqname) {
                    continue;
                }
                if let Some(idx) = self.genome.get_idx(seqname) {

                    let mut iterator = self.bwf.query_stream::<LittleEndian, R>(&mut self.reader, idx as u32, from as u32, to as u32, bin_size as u32);

                    while let Some(item) = iterator.next().await {

                        yield item;

                    }
                }
            }
        }
    }

    pub fn query<'a>(
        &'a mut self,
        seq_regex: &'a str,
        from     : usize,
        to       : usize,
        bin_size : usize,
    ) -> BlockingStream<impl Stream<Item = io::Result<BbiQueryType>> + 'a> {

        let s = Box::pin(self.query_stream(seq_regex, from, to, bin_size));

        block_on_stream(s)
    }

    pub fn genome(&self) -> &Genome {
        &self.genome
    }
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use byteorder::LittleEndian;
    use approx::assert_relative_eq;

    use crate::bigwig::BigWigFile;
    use crate::netfile::NetFile;

    #[test]
    fn test_bigwig_1() {

        let result =  BigWigFile::new_reader("tests/test_bigwig_1.bw");

        assert!(result.is_ok());

        if let Ok(mut bw) = result {

            assert_eq!(bw.genome().len(), 2);

            assert_eq!(bw.genome().seqnames[0], "test1");
            assert_eq!(bw.genome().seqnames[1], "test2");

            let mut sum_id   = 0;
            let mut sum_from = 0;
            let mut sum_to   = 0;
            let mut sum_min  = 0.0;
            let mut sum_max  = 0.0;

            for result in bw.bwf.query::<LittleEndian, NetFile>(&mut bw.reader, 0, 0, 100, 1) {
                if let Ok(item) = result {
                    sum_id   += item.data.chrom_id;
                    sum_from += item.data.from;
                    sum_to   += item.data.to;
                    sum_min  += item.data.statistics.min;
                    sum_max  += item.data.statistics.max;
                }
            }

            assert_eq!(sum_id  ,   0);
            assert_eq!(sum_from, 450);
            assert_eq!(sum_to  , 550);

            assert_relative_eq!(sum_min, 7.0, epsilon = f32::EPSILON as f64);
            assert_relative_eq!(sum_max, 7.0, epsilon = f32::EPSILON as f64);

        }
    }
}
