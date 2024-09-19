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
use std::io::{self, BufRead, BufReader, Read};
use std::error::Error;

use async_stream::stream;
use byteorder::{LittleEndian, ReadBytesExt};
use futures::executor::block_on_stream;
use futures_core::stream::Stream;
use futures::StreamExt;

use crate::bgzf::BgzfReader;
use crate::genome::Genome;
use crate::netfile::NetFile;
use crate::range::Range;
use crate::reads;

/* -------------------------------------------------------------------------- */

// Represents a BAM sequence
#[derive(Clone, Debug, Default)]
pub struct BamSeq(Vec<u8>);

/* -------------------------------------------------------------------------- */

impl fmt::Display for BamSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let t = b"=ACMGRSVTWYHKDBN";
        for &byte in &self.0 {
            let b1 = byte >> 4;
            let b2 = byte & 0xf;
            write!(f, "{}", t[b1 as usize] as char)?;
            if b2 != 0 {
                write!(f, "{}", t[b2 as usize] as char)?;
            }
        }
        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

// Represents BAM quality scores
#[derive(Clone, Debug, Default)]
pub struct BamQual(Vec<u8>);

/* -------------------------------------------------------------------------- */

impl fmt::Display for BamQual {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for &byte in &self.0 {
            write!(f, "{}", (byte + 33) as char)?;
        }
        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

// Represents a BAM auxiliary data field
#[derive(Clone, Debug, Default)]
pub struct BamAuxiliary {
    pub tag  : [u8; 2],
    pub value: BamAuxValue,
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub enum BamAuxValue {
    A         (u8),
    C         (i8),
    CUnsigned (u8),
    S         (i16),
    SUnsigned (u16),
    I         (i32),
    IUnsigned (u32),
    F         (f32),
    D         (f64),
    Z         (String),
    H         (String),
    BInt8     (Vec<i8>),
    BUint8    (Vec<u8>),
    BInt16    (Vec<i16>),
    BUint16   (Vec<u16>),
    BInt32    (Vec<i32>),
    BUint32   (Vec<u32>),
    BFloat32  (Vec<f32>),
    None      (),
}

/* -------------------------------------------------------------------------- */

impl Default for BamAuxValue {
    fn default() -> Self {
        Self::None()
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for BamAuxiliary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}:", self.tag[0] as char, self.tag[1] as char)?;
        match &self.value {
            BamAuxValue::A(v)         => write!(f, "{}", *v as char),
            BamAuxValue::C(v)         => write!(f, "{}", v),
            BamAuxValue::CUnsigned(v) => write!(f, "{}", v),
            BamAuxValue::S(v)         => write!(f, "{}", v),
            BamAuxValue::SUnsigned(v) => write!(f, "{}", v),
            BamAuxValue::I(v)         => write!(f, "{}", v),
            BamAuxValue::IUnsigned(v) => write!(f, "{}", v),
            BamAuxValue::F(v)         => write!(f, "{}", v),
            BamAuxValue::D(v)         => write!(f, "{}", v),
            BamAuxValue::Z(v)         => write!(f, "{}", v),
            BamAuxValue::H(v)         => write!(f, "{}", v),
            BamAuxValue::BInt8(v)     => write!(f, "{:?}", v),
            BamAuxValue::BUint8(v)    => write!(f, "{:?}", v),
            BamAuxValue::BInt16(v)    => write!(f, "{:?}", v),
            BamAuxValue::BUint16(v)   => write!(f, "{:?}", v),
            BamAuxValue::BInt32(v)    => write!(f, "{:?}", v),
            BamAuxValue::BUint32(v)   => write!(f, "{:?}", v),
            BamAuxValue::BFloat32(v)  => write!(f, "{:?}", v),
            BamAuxValue::None()       => panic!("internal error"),
        }
    }
}

/* -------------------------------------------------------------------------- */

impl BamAuxiliary {
    fn read<R: BufRead>(reader: &mut R) -> io::Result<(u64, Self)> {

        let mut tag = [0; 2];
        let mut n   = 0 as u64;

        reader.read_exact(&mut tag)?; n += 2;

        let value_type = reader.read_u8()?; n += 1;

        let value = match value_type {
            b'A' => {n += 1; BamAuxValue::A(reader.read_u8()?)},
            b'c' => {n += 1; BamAuxValue::C(reader.read_i8()?)},
            b'C' => {n += 1; BamAuxValue::CUnsigned(reader.read_u8()?)},
            b's' => {n += 2; BamAuxValue::S(reader.read_i16::<LittleEndian>()?)},
            b'S' => {n += 2; BamAuxValue::SUnsigned(reader.read_u16::<LittleEndian>()?)},
            b'i' => {n += 4; BamAuxValue::I(reader.read_i32::<LittleEndian>()?)},
            b'I' => {n += 4; BamAuxValue::IUnsigned(reader.read_u32::<LittleEndian>()?)},
            b'f' => {n += 4; BamAuxValue::F(reader.read_f32::<LittleEndian>()?)},
            b'd' => {n += 8; BamAuxValue::D(reader.read_f64::<LittleEndian>()?)},
            b'Z' => {
                let mut buffer : Vec<u8> = Vec::new();
                reader.read_until(0, &mut buffer)?; n += buffer.len() as u64;
                buffer.pop(); // Remove the trailing null byte
                BamAuxValue::Z(String::from_utf8_lossy(&buffer).to_string())
            }
            b'H' => {
                let mut buffer : Vec<u8>  = Vec::new();
                reader.read_until(0, &mut buffer)?; n += buffer.len() as u64;
                buffer.pop(); // Remove the trailing null byte
                BamAuxValue::H(buffer.iter().map(|b| format!("{:X}", b)).collect::<String>())
            }
            b'B' => {
                let array_type = reader.read_u8()?; n += 1;
                let array_len  = reader.read_i32::<LittleEndian>()?; n += 4;
                match array_type {
                    b'c' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_i8()?); n += 1;
                        }
                        BamAuxValue::BInt8(vec)
                    }
                    b'C' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_u8()?); n += 1;
                        }
                        BamAuxValue::BUint8(vec)
                    }
                    b's' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_i16::<LittleEndian>()?); n += 2;
                        }
                        BamAuxValue::BInt16(vec)
                    }
                    b'S' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_u16::<LittleEndian>()?); n += 2;
                        }
                        BamAuxValue::BUint16(vec)
                    }
                    b'i' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_i32::<LittleEndian>()?); n += 4;
                        }
                        BamAuxValue::BInt32(vec)
                    }
                    b'I' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_u32::<LittleEndian>()?); n += 4;
                        }
                        BamAuxValue::BUint32(vec)
                    }
                    b'f' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_f32::<LittleEndian>()?); n += 4;
                        }
                        BamAuxValue::BFloat32(vec)
                    }
                    _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid array type")),
                }
            }
            _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid auxiliary type")),
        };

        Ok((n, BamAuxiliary { tag, value }))
    }
}

/* -------------------------------------------------------------------------- */

// Represents BAM flags
#[derive(Clone, Debug, Default)]
pub struct BamFlag(u16);

/* -------------------------------------------------------------------------- */

impl BamFlag {
    pub fn bit(&self, i: u8) -> bool {
        (self.0 >> i) & 1 == 1
    }

    pub fn read_paired(&self) -> bool {
        self.bit(0)
    }

    pub fn read_mapped_proper_paired(&self) -> bool {
        self.bit(1)
    }

    pub fn unmapped(&self) -> bool {
        self.bit(2)
    }

    pub fn mate_unmapped(&self) -> bool {
        self.bit(3)
    }

    pub fn reverse_strand(&self) -> bool {
        self.bit(4)
    }

    pub fn mate_reverse_strand(&self) -> bool {
        self.bit(5)
    }

    pub fn first_in_pair(&self) -> bool {
        self.bit(6)
    }

    pub fn second_in_pair(&self) -> bool {
        self.bit(7)
    }

    pub fn secondary_alignment(&self) -> bool {
        self.bit(8)
    }

    pub fn not_passing_filters(&self) -> bool {
        self.bit(9)
    }

    pub fn duplicate(&self) -> bool {
        self.bit(10)
    }
}

/* -------------------------------------------------------------------------- */

// Represents a BAM CIGAR string
#[derive(Clone, Debug, Default)]
pub struct BamCigar(Vec<u32>);

/* -------------------------------------------------------------------------- */

impl fmt::Display for BamCigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for cigar_block in self.parse_cigar() {
            write!(f, "{}{}", cigar_block.n, cigar_block.type_)?;
        }
        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

impl BamCigar {
    fn alignment_length(&self) -> i32 {
        let mut length = 0;
        for cigar_block in self.parse_cigar() {
            match cigar_block.type_ {
                'M' | 'D' | 'N' | '=' | 'X' => length += cigar_block.n,
                _ => {}
            }
        }
        length
    }

    fn parse_cigar(&self) -> impl Iterator<Item = CigarBlock> + '_ {
        let types = b"MIDNSHP=X";
        self.0.iter().map(move |&c| {
            let n = c >> 4;
            let t = types[(c & 0xf) as usize] as char;
            CigarBlock { n: n as i32, type_: t }
        })
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Debug, Default)]
pub struct CigarBlock {
    pub n    : i32,
    pub type_: char,
}

/* -------------------------------------------------------------------------- */

// Represents a BAM header
#[derive(Clone, Debug, Default)]
pub struct BamHeader {
    pub text_length: i32,
    pub text       : String,
    pub n_ref      : i32,
}

/* -------------------------------------------------------------------------- */

// Represents a BAM block
#[derive(Clone, Debug, Default)]
pub struct BamBlock {
    pub ref_id       : i32,
    pub position     : i32,
    pub bin          : u16,
    pub mapq         : u8,
    pub rname_len    : u8,
    pub flag         : BamFlag,
    pub n_cigar_op   : u16,
    pub l_seq        : i32,
    pub next_ref_id  : i32,
    pub next_position: i32,
    pub tlen         : i32,
    pub read_name    : String,
    pub cigar        : BamCigar,
    pub seq          : BamSeq,
    pub qual         : BamQual,
    pub auxiliary    : Vec<BamAuxiliary>,
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, Default)]
pub struct BamReaderType1 {
    pub block: BamBlock,
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, Default)]
pub struct BamReaderType2 {
    pub block1: BamBlock,
    pub block2: BamBlock,
}

/* -------------------------------------------------------------------------- */

#[derive(Debug, Default)]
pub struct BamReaderOptions {
    read_name     : bool,
    read_cigar    : bool,
    read_sequence : bool,
    read_auxiliary: bool,
    read_qual     : bool,
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct BamReader<R: Read> {
    options : BamReaderOptions,
    header  : BamHeader,
    genome  : Genome,
    reader  : BufReader<BgzfReader<R>>,
}

/* -------------------------------------------------------------------------- */

impl<R: Read> BamReader<R> {
    pub fn new(reader: R, options: Option<BamReaderOptions>) -> io::Result<Self> {
        let mut bam_reader = BamReader {
            options : options.unwrap_or_default(),
            genome  : Genome::default(),
            header  : BamHeader::default(),
            reader  : BufReader::new(BgzfReader::new(reader)?),
        };

        // Default options
        bam_reader.options.read_name      = true;
        bam_reader.options.read_cigar     = true;
        bam_reader.options.read_sequence  = true;
        bam_reader.options.read_auxiliary = true;
        bam_reader.options.read_qual      = true;

        let mut magic = [0; 4];
        bam_reader.reader.read_exact(&mut magic)?;

        if &magic != b"BAM\x01" {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "not a BAM file"));
        }

        bam_reader.header.text_length = bam_reader.reader.read_i32::<LittleEndian>()?;
        let mut text_bytes = vec![0; bam_reader.header.text_length as usize];
        bam_reader.reader.read_exact(&mut text_bytes)?;
        bam_reader.header.text = String::from_utf8(text_bytes).unwrap();

        bam_reader.header.n_ref = bam_reader.reader.read_i32::<LittleEndian>()?;
        for _ in 0..bam_reader.header.n_ref {
            let length_name = bam_reader.reader.read_i32::<LittleEndian>()?;
            let mut name_bytes = vec![0; length_name as usize];
            bam_reader.reader.read_exact(&mut name_bytes)?;
            let length_seq = bam_reader.reader.read_i32::<LittleEndian>()?;
            bam_reader.genome.add_sequence(
                String::from_utf8(name_bytes).unwrap().trim_matches('\0').to_string(),
                length_seq as usize,
            ).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        }

        Ok(bam_reader)
    }
}

/* -------------------------------------------------------------------------- */

impl<R: Read> BamReader<R> {

    fn read_single_end_stream<'a>(
        &'a mut self
    ) -> impl Stream<Item = io::Result<BamReaderType1>> + 'a {

        stream! {

            let mut block_size: i32;
            let mut flag_nc   : u32;
            let mut bin_mq_nl : u32;

            let mut block = BamBlock::default();

            loop {
                let mut buf = Vec::new();

                block_size = match self.reader.read_i32::<LittleEndian>() {
                    Ok (v) => v,
                    Err(e) => { yield Err(e); return; }
                };
                block.ref_id = match self.reader.read_i32::<LittleEndian>() {
                    Ok (v) => v,
                    Err(e) => { yield Err(e); return; }
                };
                block.position = match self.reader.read_i32::<LittleEndian>() {
                    Ok (v) => v,
                    Err(e) => { yield Err(e); return; }
                };

                bin_mq_nl = match self.reader.read_u32::<LittleEndian>() {
                    Ok (v) => v,
                    Err(e) => { yield Err(e); return; }
                };    
                block.bin       = ((bin_mq_nl >>   16) & 0xffff) as u16;
                block.mapq      = ((bin_mq_nl >>    8) & 0xff  ) as u8;
                block.rname_len =  (bin_mq_nl  & 0xff) as u8;

                flag_nc = match self.reader.read_u32::<LittleEndian>() {
                    Ok (v) => v,
                    Err(e) => { yield Err(e); return; }
                };
                block.flag       = BamFlag((flag_nc >> 16) as u16);
                block.n_cigar_op = (flag_nc & 0xffff) as u16;
    
                block.l_seq = match self.reader.read_i32::<LittleEndian>() {
                    Ok (v) => v,
                    Err(e) => { yield Err(e); return; }
                };
                block.next_ref_id = match self.reader.read_i32::<LittleEndian>() {
                    Ok (v) => v,
                    Err(e) => { yield Err(e); return; }
                };
                block.next_position = match self.reader.read_i32::<LittleEndian>() {
                    Ok (v) => v,
                    Err(e) => { yield Err(e); return; }
                };
                block.tlen = match self.reader.read_i32::<LittleEndian>() {
                    Ok (v) => v,
                    Err(e) => { yield Err(e); return; }
                };
    
                // Parse the read name
                loop {
                    match self.reader.read_u8() {
                        Ok (b) if b == 0 => {
                            block.read_name = String::from_utf8(buf.clone()).unwrap();
                            break;
                        }
                        Ok (b) => buf.push(b),
                        Err(e) => { yield Err(e); return; }
                    }
                }
    
                // Parse CIGAR block
                if self.options.read_cigar {
                    block.cigar = BamCigar(Vec::with_capacity(block.n_cigar_op as usize));
                    for _ in 0..block.n_cigar_op {
                        if let Err(e) = self.reader.read_u32::<LittleEndian>().map(|v| block.cigar.0.push(v)) {
                            yield Err(e); return;
                        }
                    }
                } else {
                    for _ in 0..block.n_cigar_op {
                        if let Err(e) = io::copy(&mut self.reader, &mut io::sink()).map(|_| ()) {
                            yield Err(e); return;
                        }
                    }
                }
    
                // Parse sequence
                if self.options.read_sequence {
                    let seq_len = (block.l_seq + 1) / 2;
                    block.seq = BamSeq(vec![0; seq_len as usize]);
                    if let Err(e) = self.reader.read_exact(&mut block.seq.0) {
                        yield Err(e); return;
                    }
                } else {
                    if let Err(e) = io::copy(&mut self.reader, &mut io::sink()).map(|_| ()) {
                        yield Err(e); return;
                    }
                }
    
                // Parse qual block
                if self.options.read_qual {
                    block.qual = BamQual(vec![0; block.l_seq as usize]);
                    if let Err(e) = self.reader.read_exact(&mut block.qual.0) {
                        yield Err(e); return;
                    }
                } else {
                    if let Err(e) = io::copy(&mut self.reader, &mut io::sink()).map(|_| ()) {
                        yield Err(e); return;
                    }
                }
    
                // Read auxiliary data
                let mut position = (8 * 4 + block.rname_len as usize + 4 * block.n_cigar_op as usize
                    + (block.l_seq as usize + 1) / 2 + block.l_seq as usize) as i32;
    
                if self.options.read_auxiliary {
                    while position < block_size {
                        match BamAuxiliary::read(&mut self.reader) {
                            Ok ((bytes_read, aux)) => {
                                block.auxiliary.push(aux);
                                position += bytes_read as i32;
                            }
                            Err(e) => { yield Err(e); return; }
                        }
                    }
                } else {
                    if let Err(e) = io::copy(&mut self.reader, &mut io::sink()).map(|_| ()) {
                        yield Err(e); return;
                    }
                }

                yield Ok(BamReaderType1{
                    block: block.clone()
                });

            }
        }
    }

    fn read_paired_end_stream<'a>(
        &'a mut self
    ) -> impl Stream<Item = io::Result<BamReaderType2>> + 'a {

        self.options.read_name = true;

        stream! {

            let mut cache : std::collections::HashMap<String, BamBlock> = std::collections::HashMap::new();

            let mut iterator = Box::pin(self.read_single_end_stream());

            while let Some(item) = iterator.next().await {

                match item {

                    Err(e) => yield Err(e),
                    Ok (r) => {

                        let block1 = r.block;

                        if block1.flag.read_paired() {
                            if let Some(block2) = cache.remove(&block1.read_name) {

                                let paired_block = if block1.position < block2.position {
                                    BamReaderType2 {
                                        block1: block1,
                                        block2: block2,
                                    }
                                } else {
                                    BamReaderType2 {
                                        block1: block2,
                                        block2: block1,
                                    }
                                };

                                yield Ok(paired_block);

                            } else {

                                cache.insert(block1.read_name.clone(), block1);

                            }
                        }
                    }
                }
            }
        }
    }

    pub fn read_simple_stream<'a>(
        &'a mut self,
        join_pairs: bool,
        paired_end_strand_specific: bool
    ) -> impl Stream<Item = io::Result<reads::Read>> + 'a {

        let genome = self.genome.clone();

        stream!{

            self.options.read_cigar = true;

            let mut iterator = Box::pin(self.read_paired_end_stream());

            while let Some(item) = iterator.next().await {

                match item {

                    Err(e) => {yield Err(e); break},
                    Ok (r) => {

                        if r.block1.flag.read_paired() && join_pairs {
                            if r.block1.flag.unmapped() || !r.block1.flag.read_mapped_proper_paired() {
                                continue;
                            }
                            if r.block2.flag.unmapped() || !r.block2.flag.read_mapped_proper_paired() {
                                continue;
                            }

                            let seqname    = genome.seqnames[r.block1.ref_id as usize].clone();
                            let from       = r.block1.position;
                            let to         = r.block2.position + r.block2.cigar.alignment_length() as i32;
                            let mut strand = b'*';
                            let duplicate  = r.block1.flag.duplicate() || r.block2.flag.duplicate();
                            let mapq       = std::cmp::min(r.block1.mapq as i32, r.block2.mapq as i32);

                            if from < 0 {
                                yield Err(io::Error::new(io::ErrorKind::InvalidData, format!("invalid position detected: from={}", from)));
                            }

                            if paired_end_strand_specific {
                                if r.block1.flag.second_in_pair() {
                                    strand = if r.block1.flag.reverse_strand() { b'-' } else { b'+' };
                                } else {
                                    strand = if r.block2.flag.reverse_strand() { b'-' } else { b'+' };
                                }
                            }

                            yield Ok(reads::Read {
                                seqname   : seqname,
                                range     : Range::new(from as usize, to as usize),
                                strand    : strand as char,
                                mapq      : mapq as i64,
                                duplicate : duplicate,
                                paired_end: true,
                            });

                        } else if !r.block1.flag.unmapped() {

                            let seqname   = genome.seqnames[r.block1.ref_id as usize].clone();
                            let from      = r.block1.position;
                            let to        = r.block1.position + r.block1.cigar.alignment_length() as i32;
                            let strand    = if r.block1.flag.reverse_strand() { b'-' } else { b'+' };
                            let mapq      = r.block1.mapq as i32;
                            let duplicate = r.block1.flag.duplicate();
                            let paired    = r.block1.flag.read_paired();

                            yield Ok(reads::Read {
                                seqname   : seqname,
                                range     : Range::new(from as usize, to as usize),
                                strand    : strand as char,
                                mapq      : mapq as i64,
                                duplicate : duplicate,
                                paired_end: paired,
                            });
                        }
                    }
                }
            }
        }
    }
}

/* -------------------------------------------------------------------------- */

impl<R: Read> BamReader<R> {

    pub fn read_single_end<'a>(&'a mut self) -> impl Iterator<Item = io::Result<BamReaderType1>> + 'a {

        let s = Box::pin(self.read_single_end_stream());

        block_on_stream(s)

    }

    pub fn read_paired_end<'a>(&'a mut self) -> impl Iterator<Item = io::Result<BamReaderType2>> + 'a {

        let s = Box::pin(self.read_paired_end_stream());

        block_on_stream(s)

    }

    pub fn read_simple<'a>(
        &'a mut self,
        join_pairs: bool,
        paired_end_strand_specific: bool
    ) -> impl Iterator<Item = io::Result<reads::Read>> + 'a {

        let s = Box::pin(self.read_simple_stream(join_pairs, paired_end_strand_specific));

        block_on_stream(s)

    }
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct BamFile {
    pub reader: BamReader<NetFile>,
}

/* -------------------------------------------------------------------------- */

impl BamFile {
    pub fn open(filename: &str, options: Option<BamReaderOptions>) -> Result<Self, Box<dyn Error>> {
        let file   = NetFile::open(filename)?;
        let reader = BamReader::new(file, options)?;

        Ok(BamFile {
            reader : reader,
        })
    }
}

/* -------------------------------------------------------------------------- */

pub fn bam_read_genome<R: Read>(reader: R) -> Result<Genome, Box<dyn Error>> {
    let bam_reader = BamReader::<R>::new(reader, Some(BamReaderOptions::default()))?;
    Ok(bam_reader.genome)
}

/* -------------------------------------------------------------------------- */

pub fn bam_import_genome(filename: &str) -> Result<Genome, Box<dyn Error>> {
    let file = NetFile::open(filename)?;
    Ok(bam_read_genome(file)?)
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::bam::BamFile;

    #[test]
    fn test_bam_1() {

        let result =  BamFile::open("src/bam_test.1.bam", None);

        assert!(result.is_ok());

        if let Ok(bam) = result {

            let genome = bam.reader.genome;

            assert_eq!(genome.len(), 2);

            assert_eq!(genome.seqnames[0], "ref");
            assert_eq!(genome.seqnames[1], "ref2");

            assert_eq!(genome.lengths[0], 45);
            assert_eq!(genome.lengths[1], 40);

        }
    }

    #[test]
    fn test_bam_2() {

        let result =  BamFile::open("src/bam_test.1.bam", None);

        assert!(result.is_ok());

        let mut bam = result.unwrap();

        for item in bam.reader.read_simple(true, true) {
            if item.is_ok() {
                println!("{}", item.unwrap());
            }
        }

    }


}
