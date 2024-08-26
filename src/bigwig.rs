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

use std::io::{Read, Seek};
use std::io;

use async_stream::stream;
use futures::executor::block_on_stream;
use futures::executor::BlockingStream;
use futures_core::stream::Stream;

use byteorder::{ByteOrder, ReadBytesExt, BigEndian, LittleEndian};

use crate::genome::Genome;
use crate::bbi::RVertex;
use crate::bbi::BbiQueryType;
use crate::bbi::BbiFile;

/* -------------------------------------------------------------------------- */

const BIGWIG_MAGIC : u32 = 0x888FFC26;

/* -------------------------------------------------------------------------- */

// Struct for BigWigParameters
#[derive(Debug)]
struct BigWigParameters {
    block_size      : usize,
    items_per_slot  : usize,
    reduction_levels: Option<Vec<i32>>,
}

/* -------------------------------------------------------------------------- */

impl BigWigParameters {
    // Default constructor for BigWigParameters
    pub fn default() -> Self {
        BigWigParameters {
            block_size      : 256,
            items_per_slot  : 1024,
            reduction_levels: None,
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(PartialEq)]
enum BigWigOrder {
    LE,
    BE,
}

/* -------------------------------------------------------------------------- */

pub struct BigWigReader<R: Read + Seek> {
    reader: R,
    bwf   : BbiFile,
    genome: Genome,
    order : BigWigOrder
}

/* -------------------------------------------------------------------------- */

impl<R: Read + Seek> BigWigReader<R> {

    fn open(reader: &mut R) -> io::Result<(BbiFile, BigWigOrder)> {

        let mut bwf = BbiFile::new();

        // Try to open with little endian
        let r1 = bwf.open::<LittleEndian, R>(reader, BIGWIG_MAGIC);

        if let Err(err) = r1 {
            if err.kind() != io::ErrorKind::InvalidData || err.to_string() != "Invalid magic number" {
                return Err(err)
            }
        } else {
            return Ok((bwf, BigWigOrder::LE))
        }

        let r2 = bwf.open::<BigEndian, R>(reader, BIGWIG_MAGIC);

        if let Err(err) = r2 {
            return Err(err)
        } else {
            return Ok((bwf, BigWigOrder::BE))
        }
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

    pub fn new(mut reader: R) -> io::Result<Self> {
        let (bwf, order) = BigWigReader::<R>::open(&mut reader)?;

        let r = BigWigReader {
            reader: reader,
            bwf   : bwf,
            genome: Genome::default(),
            order : order,
        };

        match r.order {
            BigWigOrder::LE => r.initialize::<LittleEndian>(),
            BigWigOrder::BE => r.initialize::<BigEndian   >(),
        }
    }
    /*
    fn read_blocks(&mut self) -> std::sync::mpsc::Receiver<BigWigReaderType> {
        let (tx, rx) = channel();
        let bwf = self.bwf.clone();
        let mut reader = self.reader.clone();

        thread::spawn(move || {
            Self::fill_channel(tx, &mut reader, &bwf.index.root).unwrap();
        });

        rx
    }

    fn fill_channel(tx: Sender<BigWigReaderType>, reader: &mut R, vertex: &RVertex) -> Result<(), Box<dyn Error>> {
        if vertex.is_leaf != 0 {
            for i in 0..vertex.n_children as usize {
                match vertex.read_block(reader, i) {
                    Ok(block) => tx.send(BigWigReaderType { block: Some(block), error: None }).unwrap(),
                    Err(err)  => tx.send(BigWigReaderType { block: None, error: Some(Box::new(err)) }).unwrap(),
                }
            }
        } else {
            for i in 0..vertex.n_children as usize {
                Self::fill_channel(tx.clone(), reader, &vertex.children[i])?;
            }
        }
        Ok(())
    }
    */
    pub fn query<'a>(
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

                    let mut iterator = match self.order {
                        BigWigOrder::LE => self.bwf.query_iterator::<LittleEndian, R>(&mut self.reader, idx as u32, from as u32, to as u32, bin_size as u32),
                        BigWigOrder::BE => self.bwf.query_iterator::<BigEndian   , R>(&mut self.reader, idx as u32, from as u32, to as u32, bin_size as u32),
                    };
                    while let Some(item) = iterator.next() {

                        yield item;

                    }
                }
            }
        }
    }

    pub fn query_iterator<'a>(
        &'a mut self,
        seq_regex: &'a str,
        from     : usize,
        to       : usize,
        bin_size : usize,
    ) -> BlockingStream<impl Stream<Item = io::Result<BbiQueryType>> + 'a> {

        let s = Box::pin(self.query(seq_regex, from, to, bin_size));
        //pin_mut!(s);

        block_on_stream(s)

    }
}
