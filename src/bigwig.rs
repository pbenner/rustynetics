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
use std::fs::File;
use std::io::{self, Read, Write, Seek, SeekFrom};
use std::result::Result;
use std::error::Error;
use std::collections::BTreeMap;

use async_stream::stream;
use futures::executor::block_on_stream;
use futures::executor::BlockingStream;
use futures_core::stream::Stream;
use futures::StreamExt;

use byteorder::{ByteOrder, ReadBytesExt, LittleEndian};

use crate::genome::Genome;
use crate::bbi::{BbiHeader, BbiFile, BbiQueryType, BbiHeaderZoom, BbiSummaryRecord, BbiSummaryStatistics, RTree, RVertex, RVertexGenerator};
use crate::bbi::{BBI_TYPE_FIXED, BBI_TYPE_VARIABLE, BBI_TYPE_BED_GRAPH};
use crate::netfile::NetFile;
use crate::track_statistics::BinSummaryStatistics;
use crate::utility::div_int_down;

/* -------------------------------------------------------------------------- */

const BIGWIG_MAGIC : u32 = 0x888FFC26;

/* -------------------------------------------------------------------------- */

pub fn is_bigwig_file(filename: &str) -> Result<bool, Box<dyn Error>> {

    let mut file = NetFile::open(filename)?;

    let magic = file.read_u32::<LittleEndian>()?;

    Ok(BIGWIG_MAGIC == magic)

}

/* -------------------------------------------------------------------------- */

pub enum OptionBigWig {
    BlockSize(usize),
    ItemsPerSlot(usize),
    ReductionLevels(Vec<i32>),
}

/* -------------------------------------------------------------------------- */

// Struct for BigWigParameters
#[derive(Clone, Debug)]
pub struct BigWigParameters {
    pub block_size      : usize,
    pub items_per_slot  : usize,
    pub reduction_levels: Vec<i32>,
}

/* -------------------------------------------------------------------------- */

impl BigWigParameters {
    pub fn insert_option(&mut self, option: OptionBigWig) {
        match option {
            OptionBigWig::BlockSize(x)       => self.block_size       = x,
            OptionBigWig::ItemsPerSlot(x)    => self.items_per_slot   = x,
            OptionBigWig::ReductionLevels(x) => self.reduction_levels = x,
        }
    }
}

/* -------------------------------------------------------------------------- */

impl Default for BigWigParameters {
    // Default constructor for BigWigParameters
    fn default() -> Self {
        BigWigParameters {
            block_size      : 256,
            items_per_slot  : 1024,
            reduction_levels: vec![],
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct BigWigSummaryRecord {
    pub chrom     : String,
    pub from      : i32,
    pub to        : i32,
    pub statistics: BbiSummaryStatistics,
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for BigWigSummaryRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "(chrom={}, from={}, to={}, statistics={})",
            self.chrom,
            self.from,
            self.to,
            self.statistics)
    }
}

/* -------------------------------------------------------------------------- */

impl From<(BbiSummaryRecord, &Genome)> for BigWigSummaryRecord {
    // Required method
    fn from(value: (BbiSummaryRecord, &Genome)) -> Self {
        let record = value.0;
        BigWigSummaryRecord{
            chrom     : value.1.seqnames[record.chrom_id as usize].clone(),
            from      : record.from,
            to        : record.to,
            statistics: record.statistics,
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct BigWigQueryType {
    pub data     : BigWigSummaryRecord,
    pub data_type: u8,
}

/* -------------------------------------------------------------------------- */

impl From<(BbiQueryType, &Genome)> for BigWigQueryType {
    // Required method
    fn from(value: (BbiQueryType, &Genome)) -> Self {
        BigWigQueryType{
            data     : (value.0.data, value.1).into(),
            data_type: value.0.data_type,
        }
    }
}

/* -------------------------------------------------------------------------- */

impl<'a> fmt::Display for BigWigQueryType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let type_str = match self.data_type {
            BBI_TYPE_FIXED     => "fixed",
            BBI_TYPE_VARIABLE  => "variable",
            BBI_TYPE_BED_GRAPH => "bedgraph",
            _                  => panic!("internal error"),
        };
        write!(f, "(data={}, type={})",
            self.data, type_str)
    }
}

/* -------------------------------------------------------------------------- */

pub enum BigWigFile {}

/* -------------------------------------------------------------------------- */

impl BigWigFile {

    pub fn new_reader(filename : &str) -> Result<BigWigReader<NetFile>, Box<dyn Error>> {

        let file = NetFile::open(filename)?;

        BigWigReader::new(file)

    }

    pub fn new_writer(filename : &str, genome: Genome, parameters: Vec<OptionBigWig>) -> Result<BigWigWriter<File>, Box<dyn Error>> {

        let file = File::create(filename)?;

        BigWigWriter::new(file, genome, parameters)

    }

}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct BigWigReader<R: Read + Seek> {
    reader: R,
    bwf   : BbiFile,
    genome: Genome,
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
    ) -> impl Stream<Item = io::Result<BigWigQueryType>> + 'a {

        stream! {

            let re = regex::Regex::new(&format!("^{}$", seq_regex)).unwrap();

            for seqname in &self.genome.seqnames {
                if !re.is_match(seqname) {
                    continue;
                }
                if let Some(idx) = self.genome.get_idx(seqname) {

                    let mut iterator = self.bwf.query_stream::<LittleEndian, R>(&mut self.reader, idx as u32, from as u32, to as u32, bin_size as u32);

                    while let Some(item) = iterator.next().await {

                        match item {
                            Ok (r) => yield Ok(BigWigQueryType::from((r, &self.genome))),
                            Err(e) => yield Err(e)
                        }

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
    ) -> BlockingStream<impl Stream<Item = io::Result<BigWigQueryType>> + 'a> {

        let s = Box::pin(self.query_stream(seq_regex, from, to, bin_size));

        block_on_stream(s)
    }

    pub fn genome(&self) -> &Genome {
        &self.genome
    }

    pub fn header(&self) -> &BbiHeader {
        &self.bwf.header
    }

    pub fn query_slice(
        &mut self,
        seqname     : &str,
        from        : usize,
        to          : usize,
        f           : BinSummaryStatistics,
        mut bin_size: usize,
        bin_overlap : usize,
        init        : f64
    ) -> Result<(Vec<f64>, usize), Box<dyn Error>> {

        // We don't want to use regular expressions here, otherwise our sequence may come
        // from multiple chromosomes
        let id = self.genome.get_idx(seqname).ok_or(
            std::io::Error::new(std::io::ErrorKind::InvalidData, format!("Sequence '{}' not found", seqname))
        )?;

        let mut r: Vec<BbiSummaryRecord> = vec![];

        // A bin_size of 0 means that the raw data is returned as is
        if bin_size == 0 {
            for item in self.bwf.query::<LittleEndian, R>(&mut self.reader, id as u32, from as u32, to as u32, bin_size as u32) {
                if let Err(err) = item {
                    return Err(Box::new(err));
                }
                if let Ok(record) = item {
                    // Try to determine bin_size from the first record (this most likely fails for bedGraph files)
                    if bin_size == 0 {
                        if record.data_type == BBI_TYPE_BED_GRAPH {
                            return Err("failed to determine bin-size for bigWig file: data has type bedGraph".into());
                        }
                        bin_size = (record.data.to - record.data.from) as usize;
                        r = vec![BbiSummaryRecord::default(); div_int_down(to - from, bin_size)];
                    }
                    for idx in (record.data.from as usize / bin_size)..(record.data.to as usize / bin_size) {
                        if (idx) < r.len() {
                            r[idx] = record.data;
                        }
                    }
                }
            }
        } else {
            r = vec![BbiSummaryRecord::default(); div_int_down(to - from, bin_size)];
            for item in self.bwf.query::<LittleEndian, R>(&mut self.reader, id as u32, from as u32, to as u32, bin_size as u32) {
                if let Err(err) = item {
                    return Err(Box::new(err));
                }
                if let Ok(record) = item {
                    let i_from = if record.data.from as usize >= from {
                        record.data.from as usize - from
                    } else {
                        0
                    };
                    let i_to = if record.data.to as usize >= from {
                        record.data.to as usize - from
                    } else {
                        0
                    };
                    for idx in (i_from / bin_size)..(i_to / bin_size) {
                        if idx < r.len() {
                            r[idx] = record.data;
                        }
                    }
                }
            }
        }

        // Convert summary records to sequence
        let mut s = vec![init; r.len()];
        if bin_overlap != 0 {
            let mut t = BbiSummaryRecord::default();
            for i in 0..s.len() {
                t.reset();

                let i_from = if i >= bin_overlap {
                    i - bin_overlap
                } else {
                    0
                };
                let i_to = i + bin_overlap;

                for j in i_from..=i_to {
                    if j >= s.len() {
                        break;
                    }
                    let j = j as usize;
                    if r[j].statistics.valid > 0.0 {
                        t.add_record(&r[j]);
                    }
                }
                if t.statistics.valid > 0.0 {
                    s[i] = f(t.statistics.sum, t.statistics.sum_squares, t.statistics.min, t.statistics.max, t.statistics.valid);
                }
            }
        } else {
            for (i, t) in r.iter().enumerate() {
                if t.statistics.valid > 0.0 {
                    s[i] = f(t.statistics.sum, t.statistics.sum_squares, t.statistics.min, t.statistics.max, t.statistics.valid);
                }
            }
        }

        Ok((s, bin_size))
    }

    pub fn query_sequence(&mut self, seqregex: &str, f: BinSummaryStatistics, bin_size: usize, bin_overlap: usize, init: f64) -> Result<(Vec<f64>, usize), Box<dyn Error>> {
        let seqlength = self.genome.seq_length(seqregex)?;
        self.query_slice(seqregex, 0, seqlength, f, bin_size, bin_overlap, init)
    }

    pub fn get_bin_size(&mut self) -> Result<usize, Box<dyn Error>> {
        let mut bin_size = 0;
        for r in self.query(".*", 0, usize::MAX, bin_size) {
            if let Err(e) = r {
                return Err(Box::new(e));
            }
            if let Ok(record) = r {
                if record.data_type == BBI_TYPE_BED_GRAPH {
                    return Err("failed to determine bin-size for bigWig file: data has type bedGraph".into());
                }
                bin_size = (record.data.to - record.data.from) as usize;
                break;
            }
        }
        Ok(bin_size)
    }

}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct BigWigWriter<W: Write + Seek> {
    writer    : W,
    bwf       : BbiFile,
    genome    : Genome,
    parameters: BigWigParameters,
    generator : RVertexGenerator,
    leaves    : BTreeMap<usize, Vec<RVertex>>,
}

/* Some private utility functions
 * -------------------------------------------------------------------------- */

impl<W: Write + Seek> BigWigWriter<W> {

    fn use_fixed_step(&self, sequence: &Vec<f64>) -> bool {
        let n = sequence.iter().filter(|&&x| x.is_nan()).count();
        n < sequence.len() / 2
    }

    fn get_leaves_sorted(&self) -> Vec<RVertex> {
        let mut indices: Vec<_> = self.leaves.keys().cloned().collect();
        indices.sort_unstable();
        
        indices.iter().flat_map(|idx| self.leaves.get(idx).unwrap().clone()).collect()
    }

    fn reset_leaf_map(&mut self) {
        self.leaves.clear();
    }

    pub fn parameters(&self) -> &BigWigParameters {
        return &self.parameters;
    }

}

/* -------------------------------------------------------------------------- */

impl<W: Write + Seek> BigWigWriter<W> {
    pub fn new(writer: W, genome: Genome, parameters_arg: Vec<OptionBigWig>) -> Result<Self, Box<dyn Error>> {

        let mut parameters = BigWigParameters::default();

        for parameter in parameters_arg {
            parameters.insert_option(parameter);
        }

        let mut bwf = BbiFile::default();
        bwf.header.magic = BIGWIG_MAGIC;

        let mut bww = BigWigWriter {
            writer,
            bwf,
            genome,
            parameters: parameters.clone(),
            generator: RVertexGenerator::new(parameters.block_size, parameters.items_per_slot)?,
            leaves: BTreeMap::new(),
        };

        for reduction_level in &parameters.reduction_levels {

            let mut header = BbiHeaderZoom::default();
            header.reduction_level = *reduction_level as u32;

            bww.bwf.header.zoom_headers.push(header);
        }

        bww.bwf.header.zoom_levels         = parameters.reduction_levels.len() as u16;
        bww.bwf.index_zoom                 = vec![RTree::default(); parameters.reduction_levels.len()];
        bww.bwf.header.uncompress_buf_size = 1;
        bww.bwf.chrom_data.value_size      = 8;

        bww.bwf.create::<LittleEndian, W>(&mut bww.writer)?;

        Ok(bww)
    }

    pub fn write(&mut self, seqname: &str, sequence: &Vec<f64>, bin_size: usize) -> Result<(), Box<dyn Error>> {
        let idx = self.genome.get_idx(seqname).ok_or(
            std::io::Error::new(std::io::ErrorKind::InvalidData, format!("Sequence '{}' not found", seqname))
        )?;
        let mut n      = 0;
        let fixed_step = self.use_fixed_step(&sequence);

        for mut tmp in self.generator.generate::<LittleEndian>(idx, sequence, bin_size, 0, fixed_step) {
            for i in 0..tmp.vertex.n_children as usize {
                tmp.vertex.write_block::<LittleEndian, W>(&mut self.writer, &mut self.bwf, i, &tmp.blocks[i])?;
                n += 1;
            }
            self.leaves.entry(idx).or_default().push(tmp.vertex);
        }

        for v in sequence {
            self.bwf.header.summary_add_value(*v, bin_size as u64);
        }
        self.bwf.header.n_blocks += n as u64;

        Ok(())
    }

    pub fn write_zoom(&mut self, seqname: &str, sequence: &Vec<f64>, bin_size: usize, reduction_level: usize, i: usize) -> Result<(), Box<dyn Error>> {
        let idx = self.genome.get_idx(seqname).ok_or(
            std::io::Error::new(std::io::ErrorKind::InvalidData, format!("Sequence '{}' not found", seqname))
        )?;
        let mut n = 0;

        for mut tmp in self.generator.generate::<LittleEndian>(idx, sequence, bin_size, reduction_level, true) {
            for i in 0..tmp.vertex.n_children as usize {
                tmp.vertex.write_block::<LittleEndian, W>(&mut self.writer, &mut self.bwf, i, &tmp.blocks[i])?;
                n += 1;
            }
            self.leaves.entry(idx).or_default().push(tmp.vertex);
        }

        self.bwf.header.zoom_headers[i].n_blocks += n as u32;

        Ok(())
    }

    pub fn write_index(&mut self) -> Result<(), Box<dyn Error>> {
        let mut tree = RTree::default();
        tree.block_size = self.parameters.block_size as u32;
        tree.n_items_per_slot = self.parameters.items_per_slot as u32;
        
        let leaves = self.get_leaves_sorted();
        tree.build_tree(leaves)?;

        self.reset_leaf_map();
        self.bwf.index = tree;
        Ok(self.bwf.write_index::<LittleEndian, W>(&mut self.writer)?)
    }

    pub fn write_index_zoom(&mut self, i: usize) -> Result<(), Box<dyn Error>> {
        let mut tree = RTree::default();
        tree.block_size = self.parameters.block_size as u32;
        tree.n_items_per_slot = self.parameters.items_per_slot as u32;

        let leaves = self.get_leaves_sorted();
        tree.build_tree(leaves)?;

        self.reset_leaf_map();
        self.bwf.index_zoom[i] = tree;
        Ok(self.bwf.write_index_zoom::<LittleEndian, W>(&mut self.writer, i)?)
    }

    pub fn start_zoom_data(&mut self, i: usize) -> Result<(), Box<dyn Error>> {
        let offset = self.writer.seek(SeekFrom::Current(0))?;
        self.bwf.header.zoom_headers[i].data_offset = offset as u64;

        self.writer.write_all(&self.bwf.header.zoom_headers[i].n_blocks.to_le_bytes())?;
        self.bwf.header.zoom_headers[i].write_offsets::<LittleEndian, W>(&mut self.writer)?;

        Ok(())
    }

    pub fn close(&mut self) -> Result<(), Box<dyn Error>> {
        for name in &self.genome.seqnames {
            if self.bwf.chrom_data.key_size < (name.len() + 1) as u32 {
                self.bwf.chrom_data.key_size = (name.len() + 1) as u32;
            }
        }

        for name in &self.genome.seqnames {
            let mut key = vec![0; self.bwf.chrom_data.key_size as usize];
            let mut value = vec![0; self.bwf.chrom_data.value_size as usize];
            key[..name.len()].copy_from_slice(name.as_bytes());

            let idx = self.genome.get_idx(name).ok_or(
                std::io::Error::new(std::io::ErrorKind::InvalidData, format!("Sequence '{}' not found", name))
            )?;
            value[ ..4].copy_from_slice(&(idx as u32).to_le_bytes());
            value[4..8].copy_from_slice(&(self.genome.lengths[idx as usize] as u32).to_le_bytes());

            self.bwf.chrom_data.add(key, value)?;
        }

        self.bwf.write_chrom_list::<LittleEndian, W>(&mut self.writer)?;
        self.bwf.header.write_n_blocks::<LittleEndian, W>(&mut self.writer)?;

        for i in 0..self.bwf.header.zoom_headers.len() {
            self.bwf.header.zoom_headers[i].write_n_blocks::<LittleEndian, W>(&mut self.writer)?;
        }

        self.writer.seek(SeekFrom::End(0))?;
        self.bwf.header.write_summary::<LittleEndian, W>(&mut self.writer)?;
        self.writer.write_all(&self.bwf.header.magic.to_le_bytes())?;

        Ok(())
    }
}

/* Utility functions
 * -------------------------------------------------------------------------- */

pub fn bigwig_read_genome<R: Read + Seek>(file: R) -> Result<Genome, Box<dyn Error>> {
    let reader = BigWigReader::new(file)?;
    Ok(reader.genome().clone())
}

pub fn bigwig_import_genome(filename: &str) -> Result<Genome, Box<dyn Error>> {
    let file = NetFile::open(filename)?;
    let genome = bigwig_read_genome(file)
        .map_err(|err| format!("Importing genome from `{}` failed: {}", filename, err))?;
    Ok(genome)
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

            assert_relative_eq!(sum_min, 49.5, epsilon = 1e-6);
            assert_relative_eq!(sum_max, 49.5, epsilon = 1e-6);

        }
    }
}
