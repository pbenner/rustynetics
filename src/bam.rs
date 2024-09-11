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
use std::io::{self, BufRead};
use byteorder::{LittleEndian, ReadBytesExt};

/* -------------------------------------------------------------------------- */

// Represents a BAM sequence
#[derive(Debug)]
struct BamSeq(Vec<u8>);

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
#[derive(Debug)]
struct BamQual(Vec<u8>);

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
#[derive(Debug)]
struct BamAuxiliary {
    tag  : [u8; 2],
    value: BamAuxValue,
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
enum BamAuxValue {
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
        }
    }
}

/* -------------------------------------------------------------------------- */

impl BamAuxiliary {
    fn read<R: BufRead>(reader: &mut R) -> io::Result<Self> {
        let mut tag = [0; 2];
        reader.read_exact(&mut tag)?;

        let value_type = reader.read_u8()?;
        let value = match value_type {
            b'A' => BamAuxValue::A(reader.read_u8()?),
            b'c' => BamAuxValue::C(reader.read_i8()?),
            b'C' => BamAuxValue::CUnsigned(reader.read_u8()?),
            b's' => BamAuxValue::S(reader.read_i16::<LittleEndian>()?),
            b'S' => BamAuxValue::SUnsigned(reader.read_u16::<LittleEndian>()?),
            b'i' => BamAuxValue::I(reader.read_i32::<LittleEndian>()?),
            b'I' => BamAuxValue::IUnsigned(reader.read_u32::<LittleEndian>()?),
            b'f' => BamAuxValue::F(reader.read_f32::<LittleEndian>()?),
            b'd' => BamAuxValue::D(reader.read_f64::<LittleEndian>()?),
            b'Z' => {
                let mut buffer = Vec::new();
                reader.read_until(0, &mut buffer)?;
                buffer.pop(); // Remove the trailing null byte
                BamAuxValue::Z(String::from_utf8(buffer)?)
            }
            b'H' => {
                let mut buffer = Vec::new();
                reader.read_until(0, &mut buffer)?;
                buffer.pop(); // Remove the trailing null byte
                BamAuxValue::H(buffer.iter().map(|b| format!("{:X}", b)).collect::<String>())
            }
            b'B' => {
                let array_type = reader.read_u8()?;
                let array_len = reader.read_i32::<LittleEndian>()?;
                match array_type {
                    b'c' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_i8()?);
                        }
                        BamAuxValue::BInt8(vec)
                    }
                    b'C' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_u8()?);
                        }
                        BamAuxValue::BUint8(vec)
                    }
                    b's' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_i16::<LittleEndian>()?);
                        }
                        BamAuxValue::BInt16(vec)
                    }
                    b'S' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_u16::<LittleEndian>()?);
                        }
                        BamAuxValue::BUint16(vec)
                    }
                    b'i' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_i32::<LittleEndian>()?);
                        }
                        BamAuxValue::BInt32(vec)
                    }
                    b'I' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_u32::<LittleEndian>()?);
                        }
                        BamAuxValue::BUint32(vec)
                    }
                    b'f' => {
                        let mut vec = Vec::with_capacity(array_len as usize);
                        for _ in 0..array_len {
                            vec.push(reader.read_f32::<LittleEndian>()?);
                        }
                        BamAuxValue::BFloat32(vec)
                    }
                    _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid array type")),
                }
            }
            _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid auxiliary type")),
        };
        Ok(BamAuxiliary { tag, value })
    }
}

/* -------------------------------------------------------------------------- */

// Represents BAM flags
#[derive(Debug)]
struct BamFlag(u16);

/* -------------------------------------------------------------------------- */

impl BamFlag {
    fn bit(&self, i: u8) -> bool {
        (self.0 >> i) & 1 == 1
    }

    fn read_paired(&self) -> bool {
        self.bit(0)
    }

    fn read_mapped_proper_paired(&self) -> bool {
        self.bit(1)
    }

    fn unmapped(&self) -> bool {
        self.bit(2)
    }

    fn mate_unmapped(&self) -> bool {
        self.bit(3)
    }

    fn reverse_strand(&self) -> bool {
        self.bit(4)
    }

    fn mate_reverse_strand(&self) -> bool {
        self.bit(5)
    }

    fn first_in_pair(&self) -> bool {
        self.bit(6)
    }

    fn second_in_pair(&self) -> bool {
        self.bit(7)
    }

    fn secondary_alignment(&self) -> bool {
        self.bit(8)
    }

    fn not_passing_filters(&self) -> bool {
        self.bit(9)
    }

    fn duplicate(&self) -> bool {
        self.bit(10)
    }
}

/* -------------------------------------------------------------------------- */

// Represents a BAM CIGAR string
#[derive(Debug)]
struct BamCigar(Vec<u32>);

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

#[derive(Debug)]
struct CigarBlock {
    n    : i32,
    type_: char,
}

/* -------------------------------------------------------------------------- */

// Represents a BAM header
#[derive(Debug)]
struct BamHeader {
    text_length: i32,
    text       : String,
    n_ref      : i32,
}

/* -------------------------------------------------------------------------- */

// Represents a BAM block
#[derive(Debug)]
struct BamBlock {
    ref_id       : i32,
    position     : i32,
    bin          : u16,
    mapq         : u8,
    rname_len    : u8,
    flag         : BamFlag,
    n_cigar_op   : u16,
    l_seq        : i32,
    next_ref_id  : i32,
    next_position: i32,
    tlen         : i32,
    read_name    : String,
    cigar        : BamCigar,
    seq          : BamSeq,
    qual         : BamQual,
    auxiliary    : Vec<BamAuxiliary>,
}

/* -------------------------------------------------------------------------- */

#[derive(Default)]
struct BamReaderOptions {
    read_name     : bool,
    read_cigar    : bool,
    read_sequence : bool,
    read_auxiliary: bool,
    read_qual     : bool,
}

/* -------------------------------------------------------------------------- */

#[derive(Default)]
struct BamReader {
    options    : BamReaderOptions,
    header     : BamHeader,
    genome     : Genome,
    bgzf_reader: BgzfReader,
}

/* -------------------------------------------------------------------------- */

#[derive(Default)]
struct BamReaderType1 {
    bam_block: BamBlock,
    error    : Option<io::Error>,
}

/* -------------------------------------------------------------------------- */

#[derive(Default)]
struct BamReaderType2 {
    block1: BamBlock,
    block2: BamBlock,
    error : Option<io::Error>,
}

/* -------------------------------------------------------------------------- */

impl BamReader {
    pub fn new<R: Read>(reader: R, options: Option<BamReaderOptions>) -> io::Result<Self> {
        let mut bam_reader = BamReader {
            options: options.unwrap_or_default(),
            ..Default::default()
        };

        // Default options
        bam_reader.options.read_name = true;
        bam_reader.options.read_cigar = true;
        bam_reader.options.read_sequence = true;
        bam_reader.options.read_auxiliary = true;
        bam_reader.options.read_qual = true;

        let mut magic = [0; 4];
        bam_reader.bgzf_reader = BgzfReader::new(reader)?;
        bam_reader.bgzf_reader.read_exact(&mut magic)?;

        if &magic != b"BAM\1" {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "not a BAM file"));
        }

        bam_reader.header.text_length = bam_reader.bgzf_reader.read_i32::<LittleEndian>()?;
        let mut text_bytes = vec![0; bam_reader.header.text_length as usize];
        bam_reader.bgzf_reader.read_exact(&mut text_bytes)?;
        bam_reader.header.text = String::from_utf8(text_bytes).unwrap();

        bam_reader.header.n_ref = bam_reader.bgzf_reader.read_i32::<LittleEndian>()?;
        for _ in 0..bam_reader.header.n_ref {
            let length_name = bam_reader.bgzf_reader.read_i32::<LittleEndian>()?;
            let mut name_bytes = vec![0; length_name as usize];
            bam_reader.bgzf_reader.read_exact(&mut name_bytes)?;
            let length_seq = bam_reader.bgzf_reader.read_i32::<LittleEndian>()?;
            bam_reader.genome.add_sequence(
                String::from_utf8(name_bytes).unwrap().trim_matches('\0').to_string(),
                length_seq as usize,
            );
        }

        Ok(bam_reader)
    }

    pub fn read_single_end(&mut self) -> impl Iterator<Item = BamReaderType1> + '_ {
        let mut channel = Vec::new();
        self.read_single_end_into_channel(&mut channel);
        channel.into_iter()
    }

    fn read_single_end_into_channel(&mut self, channel: &mut Vec<BamReaderType1>) {
        let mut block_size: i32;
        let mut flag_nc: u32;
        let mut bin_mq_nl: u32;

        let mut block = BamReaderType1::default();
        let mut block_reserve = BamReaderType1::default();

        loop {
            if let Err(e) = self.bgzf_reader.read_i32_into::<LittleEndian>(&mut block_size) {
                if e.kind() == io::ErrorKind::UnexpectedEof {
                    return;
                }
                channel.push(BamReaderType1 {
                    error: Some(e),
                    ..Default::default()
                });
                return;
            }

            if let Err(e) = self.bgzf_reader.read_i32_into::<LittleEndian>(&mut block.bam_block.ref_id) {
                channel.push(BamReaderType1 {
                    error: Some(e),
                    ..Default::default()
                });
                return;
            }

            // (continue similarly for other fields...)
            // Populate the BamBlock structure

            // Add the block to the channel
            channel.push(block.clone());
            std::mem::swap(&mut block, &mut block_reserve);
        }
    }

    pub fn read_paired_end(&mut self) -> impl Iterator<Item = BamReaderType2> + '_ {
        let mut channel = Vec::new();
        self.read_paired_end_into_channel(&mut channel);
        channel.into_iter()
    }

    fn read_paired_end_into_channel(&mut self, channel: &mut Vec<BamReaderType2>) {
        let mut cache = std::collections::HashMap::new();
        self.options.read_name = true;

        for r in self.read_single_end() {
            if let Some(e) = r.error {
                channel.push(BamReaderType2 {
                    error: Some(e),
                    ..Default::default()
                });
                continue;
            }
            let block1 = r.bam_block;

            if block1.flag.read_paired() {
                if let Some(block2) = cache.remove(&block1.read_name) {
                    let paired_block = if block1.position < block2.position {
                        BamReaderType2 {
                            block1,
                            block2,
                            ..Default::default()
                        }
                    } else {
                        BamReaderType2 {
                            block1: block2,
                            block2: block1,
                            ..Default::default()
                        }
                    };
                    channel.push(paired_block);
                } else {
                    cache.insert(block1.read_name.clone(), block1);
                }
            }
        }
    }

    pub fn read_simple(&mut self, join_pairs: bool, paired_end_strand_specific: bool) -> impl Iterator<Item = Read> + '_ {
        self.options.read_cigar = true;
        let mut channel = Vec::new();

        for r in self.read_paired_end() {
            if let Some(e) = r.error {
                break;
            }

            if r.block1.flag.read_paired() && join_pairs {
                if r.block1.flag.unmapped() || !r.block1.flag.read_mapped_proper_paired() {
                    continue;
                }
                if r.block2.flag.unmapped() || !r.block2.flag.read_mapped_proper_paired() {
                    continue;
                }

                let seqname    = self.genome.seqnames[r.block1.ref_id as usize].clone();
                let from       = r.block1.position;
                let to         = r.block2.position + r.block2.cigar.alignment_length() as i32;
                let mut strand = b'*';
                let duplicate  = r.block1.flag.duplicate() || r.block2.flag.duplicate();
                let mapq       = std::cmp::min(r.block1.mapq as i32, r.block2.mapq as i32);

                if paired_end_strand_specific {
                    if r.block1.flag.second_in_pair() {
                        strand = if r.block1.flag.reverse_strand() { b'-' } else { b'+' };
                    } else {
                        strand = if r.block2.flag.reverse_strand() { b'-' } else { b'+' };
                    }
                }

                channel.push(Read {
                    range: GRange {
                        seqname,
                        range: Range { from, to },
                        strand,
                    },
                    mapq,
                    duplicate,
                    paired: true,
                });
            } else if !r.block1.flag.unmapped() {
                let seqname   = self.genome.seqnames[r.block1.ref_id as usize].clone();
                let from      = r.block1.position;
                let to        = r.block1.position + r.block1.cigar.alignment_length() as i32;
                let strand    = if r.block1.flag.reverse_strand() { b'-' } else { b'+' };
                let mapq      = r.block1.mapq as i32;
                let duplicate = r.block1.flag.duplicate();
                let paired    = r.block1.flag.read_paired();

                channel.push(Read {
                    range: GRange {
                        seqname,
                        range: Range { from, to },
                        strand,
                    },
                    mapq,
                    duplicate,
                    paired,
                });
            }
        }

        channel.into_iter()
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Default)]
struct BamFile {
    bam_reader: BamReader,
    file      : Option<File>,
}

/* -------------------------------------------------------------------------- */

impl BamFile {
    pub fn open<P: AsRef<Path>>(filename: P, options: Option<BamReaderOptions>) -> io::Result<Self> {
        let file   = File::open(filename)?;
        let reader = BamReader::new(BufReader::new(&file), options)?;

        Ok(BamFile {
            bam_reader: reader,
            file: Some(file),
        })
    }

    pub fn close(&mut self) -> io::Result<()> {
        if let Some(file) = self.file.take() {
            file.sync_all()?;
        }
        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

pub fn bam_read_genome<R: Read>(reader: R) -> io::Result<Genome> {
    let bam_reader = BamReader::new(reader, Some(BamReaderOptions::default()))?;
    Ok(bam_reader.genome)
}

/* -------------------------------------------------------------------------- */

pub fn bam_import_genome<P: AsRef<Path>>(filename: P) -> io::Result<Genome> {
    let file = File::open(filename)?;
    bam_read_genome(BufReader::new(file))
}
