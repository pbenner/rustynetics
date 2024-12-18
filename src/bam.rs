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
use std::io::{self, Read};
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
use crate::read;
use crate::utility_io::{read_until_null, skip_n_bytes};

/* -------------------------------------------------------------------------- */

// Represents a BAM sequence
#[derive(Clone, Debug, Default)]
pub struct BamSeq(pub Vec<u8>);

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
pub struct BamQual(pub Vec<u8>);

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
    fn read<R: Read>(reader: &mut R) -> io::Result<(u64, Self)> {

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
                let buffer : Vec<u8> = read_until_null(reader)?; n += (buffer.len() + 1) as u64;
                BamAuxValue::Z(String::from_utf8_lossy(&buffer).to_string())
            }
            b'H' => {
                let buffer : Vec<u8> = read_until_null(reader)?; n += (buffer.len() + 1) as u64;
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
pub struct BamFlag(pub u16);

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
pub struct BamCigar(pub Vec<u32>);

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
    pub fn alignment_length(&self) -> usize {
        let mut length = 0;
        for cigar_block in self.parse_cigar() {
            match cigar_block.type_ {
                'M' | 'D' | 'N' | '=' | 'X' => length += cigar_block.n as usize,
                _ => {}
            }
        }
        length
    }

    pub fn parse_cigar(&self) -> impl Iterator<Item = CigarBlock> + '_ {
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

/// Options to configure BAM file reading behavior.
///
/// This struct allows users to specify which fields to read from the BAM file,
/// enabling optimized data access based on specific requirements.
///
/// # Fields
/// - `read_name`: If `true`, reads the name of each read.
/// - `read_cigar`: If `true`, reads the CIGAR string for each read.
/// - `read_sequence`: If `true`, reads the sequence of each read.
/// - `read_auxiliary`: If `true`, reads auxiliary fields in the BAM file.
/// - `read_qual`: If `true`, reads the quality scores of each read.
#[derive(Copy, Clone, Debug, Default)]
pub struct BamReaderOptions {
    pub read_name     : bool,
    pub read_cigar    : bool,
    pub read_sequence : bool,
    pub read_auxiliary: bool,
    pub read_qual     : bool,
}

/* -------------------------------------------------------------------------- */

/// A reader for BAM files, supporting various options for reading specific data fields.
///
/// The `BamReader` struct allows for efficient reading and handling of BAM files,
/// providing flexibility in the amount of data loaded based on `BamReaderOptions`.
///
/// # Type Parameters
/// - `R`: A `Read` type that provides access to the BAM file data.
///
/// # Fields
/// - `options`: The `BamReaderOptions` specifying which fields to read.
/// - `header`: The header information from the BAM file, stored as a `BamHeader`.
/// - `genome`: The genomic reference information associated with the BAM file.
/// - `reader`: The underlying `BgzfReader` for reading compressed BAM data.
#[derive(Debug)]
pub struct BamReader<R: Read> {
    options : BamReaderOptions,
    header  : BamHeader,
    genome  : Genome,
    reader  : BgzfReader<R>,
}

/* -------------------------------------------------------------------------- */

impl<R: Read> BamReader<R> {
    /// Creates a new `BamReader` with optional `BamReaderOptions`.
    ///
    /// # Arguments
    /// * `reader` - A reader that implements `Read` for reading BAM files.
    /// * `options_arg` - Optional `BamReaderOptions` to specify read preferences.
    ///
    /// If no options are provided, defaults are used to read all sections
    /// (name, CIGAR, sequence, auxiliary, and quality).
    ///
    /// # Errors
    /// Returns an `io::Error` if the file does not start with the "BAM\x01" magic
    /// bytes or if header data cannot be read.
    pub fn new(reader: R, options_arg: Option<BamReaderOptions>) -> io::Result<Self> {

        let mut bam_reader = BamReader {
            options : options_arg.unwrap_or_default(),
            genome  : Genome::default(),
            header  : BamHeader::default(),
            reader  : BgzfReader::new(reader)?,
        };

        if options_arg.is_none() {
            // Default options
            bam_reader.options.read_name      = true;
            bam_reader.options.read_cigar     = true;
            bam_reader.options.read_sequence  = true;
            bam_reader.options.read_auxiliary = true;
            bam_reader.options.read_qual      = true;
        }

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

    /// Returns a reference to the parsed genome data.
    ///
    /// # Returns
    /// A reference to the `Genome` structure populated from BAM header information.
    pub fn get_genome(&self) -> &Genome {
        &self.genome
    }
}

/* -------------------------------------------------------------------------- */

impl<R: Read> BamReader<R> {

    /// Reads single-end reads from the BAM file as a stream.
    ///
    /// # Returns
    /// An asynchronous stream of `io::Result<BamReaderType1>` where each item
    /// represents a single read.
    ///
    /// Yields `io::Error` if any data reading/parsing fails.
    pub fn read_single_end_stream<'a>(
        &'a mut self
    ) -> impl Stream<Item = io::Result<BamReaderType1>> + 'a {

        stream! {

            let mut block_size: i32;
            let mut flag_nc   : u32;
            let mut bin_mq_nl : u32;

            loop {
                let mut block = BamBlock::default();
                let mut buf   = Vec::new();

                block_size = match self.reader.read_i32::<LittleEndian>() {
                    Ok (v) => v,
                    Err(e) => {
                        // No more reads available, exiting
                        if e.kind() == std::io::ErrorKind::UnexpectedEof {
                            return;
                        }
                        yield Err(e); return;
                    }
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
                    if let Err(e) = skip_n_bytes(&mut self.reader, (block.n_cigar_op * 4) as usize) {
                        yield Err(e); return;
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
                    let seq_len = (block.l_seq + 1) / 2;

                    if let Err(e) = skip_n_bytes(&mut self.reader, seq_len as usize) {
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
                    if let Err(e) = skip_n_bytes(&mut self.reader, block.l_seq as usize) {
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
                    if let Err(e) = skip_n_bytes(&mut self.reader, (block_size - position) as usize) {
                        yield Err(e); return;
                    }
                }

                yield Ok(BamReaderType1{
                    block
                });

            }
        }
    }

    /// Reads paired-end reads from the BAM file as a stream.
    ///
    /// # Returns
    /// An asynchronous stream of `io::Result<BamReaderType2>` where each item
    /// represents a pair of reads.
    ///
    /// Yields `io::Error` if any data reading/parsing fails.
    pub fn read_paired_end_stream<'a>(
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

    /// Reads simplified read data as a stream, with options to join pairs and
    /// use strand-specific behavior.
    ///
    /// # Arguments
    /// * `join_pairs` - If true, pairs reads from paired-end sequencing.
    /// * `paired_end_strand_specific` - If true, applies strand-specific behavior
    ///   for paired-end reads.
    ///
    /// # Returns
    /// An asynchronous stream of `io::Result<read::Read>` for each read.
    pub fn read_simple_stream<'a>(
        &'a mut self,
        join_pairs: bool,
        paired_end_strand_specific: bool
    ) -> impl Stream<Item = io::Result<read::Read>> + 'a {

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

                            yield Ok(read::Read {
                                seqname   : seqname,
                                range     : Range::new(from as usize, to as usize),
                                strand    : strand as char,
                                mapq      : mapq   as i64,
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

                            yield Ok(read::Read {
                                seqname   : seqname,
                                range     : Range::new(from as usize, to as usize),
                                strand    : strand as char,
                                mapq      : mapq   as i64,
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

    /// Reads single-end reads synchronously.
    ///
    /// # Returns
    /// An iterator over `io::Result<BamReaderType1>`, where each item
    /// represents a single read.
    pub fn read_single_end<'a>(&'a mut self) -> impl Iterator<Item = io::Result<BamReaderType1>> + 'a {

        let s = Box::pin(self.read_single_end_stream());

        block_on_stream(s)

    }

    /// Reads paired-end reads synchronously.
    ///
    /// # Returns
    /// An iterator over `io::Result<BamReaderType2>`, where each item
    /// represents a pair of reads.
    pub fn read_paired_end<'a>(&'a mut self) -> impl Iterator<Item = io::Result<BamReaderType2>> + 'a {

        let s = Box::pin(self.read_paired_end_stream());

        block_on_stream(s)

    }

    /// Reads simple read data synchronously, with options to join pairs and
    /// use strand-specific behavior.
    ///
    /// # Arguments
    /// * `join_pairs` - If true, pairs reads from paired-end sequencing.
    /// * `paired_end_strand_specific` - If true, applies strand-specific behavior
    ///   for paired-end reads.
    ///
    /// # Returns
    /// An iterator over `io::Result<read::Read>`.
    pub fn read_simple<'a>(
        &'a mut self,
        join_pairs: bool,
        paired_end_strand_specific: bool
    ) -> impl Iterator<Item = io::Result<read::Read>> + 'a {

        let s = Box::pin(self.read_simple_stream(join_pairs, paired_end_strand_specific));

        block_on_stream(s)

    }
}

/* -------------------------------------------------------------------------- */

/// `BamFile` is a wrapper around `BamReader`, designed to read BAM files from
/// either local files or HTTP sources.
#[derive(Debug)]
pub struct BamFile {
    pub reader: BamReader<NetFile>,
}

/* -------------------------------------------------------------------------- */

impl BamFile {
    /// Opens a BAM file by filename with optional `BamReaderOptions`.
    ///
    /// # Arguments
    /// * `filename` - The file path to the BAM file.
    /// * `options` - Optional `BamReaderOptions` to specify read preferences.
    ///
    /// # Errors
    /// Returns an error if the file cannot be opened or if the reader setup fails.
    pub fn open(filename: &str, options: Option<BamReaderOptions>) -> Result<Self, Box<dyn Error>> {
        let file   = NetFile::open(filename)?;
        let reader = BamReader::new(file, options)?;

        Ok(BamFile {
            reader : reader,
        })
    }
}

/* -------------------------------------------------------------------------- */

/// Reads the genome information from a BAM file.
///
/// # Arguments
///
/// * `reader` - A data source implementing `Read`, potentially from a local or remote BAM file.
///
/// # Returns
///
/// Returns a `Genome` with chromosome names and lengths parsed from the BAM header.
///
/// # Errors
///
/// Will return an error if the BAM file header cannot be parsed correctly.
pub fn bam_read_genome<R: Read>(reader: R) -> Result<Genome, Box<dyn Error>> {
    let bam_reader = BamReader::<R>::new(reader, Some(BamReaderOptions::default()))?;
    Ok(bam_reader.genome)
}

/* -------------------------------------------------------------------------- */

/// Imports genome information from a BAM file located at a given path or URL.
///
/// # Arguments
///
/// * `filename` - Path to a local BAM file or a URL of a remote BAM file.
///
/// # Returns
///
/// A `Result` containing the `Genome` if successful, or an error if the file could not be opened or parsed.
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
    fn test_bam_genome() {
        let result = BamFile::open("tests/test_bam_1.bam", None);
        assert!(result.is_ok());
        if let Ok(bam) = result {
            let genome = bam.reader.genome;
            assert_eq!(genome.len(), 2);
            assert_eq!(genome.seqnames[0], "ref");
            assert_eq!(genome.seqnames[1], "ref2");
            assert_eq!(genome.lengths [0], 45);
            assert_eq!(genome.lengths [1], 40);
        }
    }

    #[test]
    fn test_bam_read_simple() {

        let result = BamFile::open("tests/test_bam_2.bam", None);

        assert!(result.is_ok());

        let mut bam = result.unwrap();
        let mut cnt = 0;

        for item in bam.reader.read_simple(true, true) {
            if !item.is_ok() {
                break;
            }
            cnt += 1;
        }

        assert_eq!(cnt, 2335);
    }


}
