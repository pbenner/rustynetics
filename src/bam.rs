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
    n: i32,
    type_: char,
}

/* -------------------------------------------------------------------------- */

// Represents a BAM header
#[derive(Debug)]
struct BamHeader {
    text_length: i32,
    text: String,
    n_ref: i32,
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
