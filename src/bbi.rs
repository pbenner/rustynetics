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

use std::io::{self, Cursor, Read, Seek, SeekFrom, Write};

use std::f32;
use std::f64;

use byteorder::{ByteOrder, ReadBytesExt, WriteBytesExt, BigEndian, LittleEndian};

use flate2::Compression;
use flate2::write::ZlibEncoder;
use flate2::read::ZlibDecoder;

use std::sync::mpsc::{channel, Receiver};

/* -------------------------------------------------------------------------- */

const CIRTREE_MAGIC      : u32   = 0x78ca8c91;
const IDX_MAGIC          : u32   = 0x2468ace0;
const BBI_MAX_ZOOM_LEVELS: usize = 10;
const BBI_RES_INCREMENT  : u32   = 4;
const BBI_TYPE_FIXED     : u8    = 3;
const BBI_TYPE_VARIABLE  : u8    = 2;
const BBI_TYPE_BED_GRAPH : u8    = 1;

/* -------------------------------------------------------------------------- */

fn file_read_at<T: Read + Seek>(file: &mut T, offset: u64, data: &mut [u8]) -> io::Result<Vec<u8>> {
    let current_position = file.seek(SeekFrom::Current(0))?;
    file.seek(SeekFrom::Start(offset))?;
    let buffer = Vec::new();
    file.read_exact(data)?;
    file.seek(SeekFrom::Start(current_position))?;
    Ok(buffer)
}

fn file_write_at<T: Write + Seek>(file: &mut T, offset: u64, data: &[u8]) -> io::Result<()> {
    let current_position = file.seek(SeekFrom::Current(0))?;
    file.seek(SeekFrom::Start(offset))?;
    file.write_all(data)?;
    file.seek(SeekFrom::Start(current_position))?;
    Ok(())
}

fn uncompress_slice(data: &[u8]) -> io::Result<Vec<u8>> {
    let mut decoder = ZlibDecoder::new(data);
    let mut buffer = Vec::new();
    decoder.read_to_end(&mut buffer)?;
    Ok(buffer)
}

fn compress_slice(data: &[u8]) -> io::Result<Vec<u8>> {
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::best());
    encoder.write_all(data)?;
    let compressed_data = encoder.finish()?;
    Ok(compressed_data)
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Copy, Default, Debug)]
struct BbiZoomRecord {
    chrom_id   : u32,
    start      : u32,
    end        : u32,
    valid      : u32,
    min        : f32,
    max        : f32,
    sum        : f32,
    sum_squares: f32,
}

/* -------------------------------------------------------------------------- */

impl BbiZoomRecord {
    fn add_value(&mut self, x: f64) {
        if x.is_nan() {
            return;
        }
        if self.min.is_nan() || self.min > x as f32 {
            self.min = x as f32;
        }
        if self.max.is_nan() || self.max < x as f32 {
            self.max = x as f32;
        }
        self.valid       += 1;
        self.sum         += x as f32;
        self.sum_squares += (x * x) as f32;
    }

    fn read<E: ByteOrder, T: Read + Seek>(&mut self, reader: &mut T) -> io::Result<()> {
        self.chrom_id    = reader.read_u32::<E>()?;
        self.start       = reader.read_u32::<E>()?;
        self.end         = reader.read_u32::<E>()?;
        self.valid       = reader.read_u32::<E>()?;
        self.min         = reader.read_f32::<E>()?;
        self.max         = reader.read_f32::<E>()?;
        self.sum         = reader.read_f32::<E>()?;
        self.sum_squares = reader.read_f32::<E>()?;
        Ok(())
    }

    fn write<E: ByteOrder, T: Write + Seek>(&self, writer: &mut T) -> io::Result<()> {
        writer.write_u32::<E>(self.chrom_id)?;
        writer.write_u32::<E>(self.start)?;
        writer.write_u32::<E>(self.end)?;
        writer.write_u32::<E>(self.valid)?;
        writer.write_f32::<E>(self.min)?;
        writer.write_f32::<E>(self.max)?;
        writer.write_f32::<E>(self.sum)?;
        writer.write_f32::<E>(self.sum_squares)?;
        Ok(())
    }

    fn write_buffer<E: ByteOrder>(&self, buffer : &mut Vec<u8>) -> io::Result<()> {
        let mut cursor = Cursor::new(buffer);

        self.write::<E, _>(&mut cursor)
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
struct BbiSummaryStatistics {
    valid      : f64,
    min        : f64,
    max        : f64,
    sum        : f64,
    sum_squares: f64,
}

/* -------------------------------------------------------------------------- */

impl BbiSummaryStatistics {

    fn new() -> Self {
        BbiSummaryStatistics {
            valid      : 0.0,
            min        : f64::INFINITY,
            max        : f64::NEG_INFINITY,
            sum        : 0.0,
            sum_squares: 0.0,
        }
    }

    fn reset(&mut self) {
        self.valid       = 0.0;
        self.min         = f64::INFINITY;
        self.max         = f64::NEG_INFINITY;
        self.sum         = 0.0;
        self.sum_squares = 0.0;
    }

    fn add_value(&mut self, x: f64) {
        if x.is_nan() {
            return;
        }
        self.valid       += 1.0;
        self.min          = self.min.min(x);
        self.max          = self.max.max(x);
        self.sum         += x;
        self.sum_squares += x * x;
    }

    fn add(&mut self, other: &BbiSummaryStatistics) {
        self.valid += other.valid;
        self.min = self.min.min(other.min);
        self.max = self.max.max(other.max);
        self.sum += other.sum;
        self.sum_squares += other.sum_squares;
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
struct BbiSummaryRecord {
    chrom_id  : i32,
    from      : i32,
    to        : i32,
    statistics: BbiSummaryStatistics,
}

/* -------------------------------------------------------------------------- */

impl BbiSummaryRecord {
    fn new() -> BbiSummaryRecord {
        BbiSummaryRecord {
            chrom_id  : -1,
            from      :  0,
            to        :  0,
            statistics: BbiSummaryStatistics::new(),
        }
    }

    fn reset(&mut self) {
        self.chrom_id = -1;
        self.from = 0;
        self.to = 0;
        self.statistics.reset();
    }

    fn add_record(&mut self, other: &BbiSummaryRecord) {
        if self.chrom_id == -1 {
            self.chrom_id = other.chrom_id;
            self.from = other.from;
            self.to = other.to;
        }
        if self.to < other.from {
            self.statistics.valid += (other.from - self.to) as f64;
            if self.statistics.min > 0.0 {
                self.statistics.min = 0.0;
            }
            if self.statistics.max < 0.0 {
                self.statistics.max = 0.0;
            }
        }
        self.to = other.to;
        self.statistics.add(&other.statistics);
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
struct BbiDataHeader {
    chrom_id  : u32,
    start     : u32,
    end       : u32,
    step      : u32,
    span      : u32,
    kind      : u8,
    reserved  : u8,
    item_count: u16,
}

/* -------------------------------------------------------------------------- */

impl BbiDataHeader {

    fn new() -> Self {
        BbiDataHeader {
            chrom_id  : 0,
            start     : 0,
            end       : 0,
            step      : 0,
            span      : 0,
            kind      : 0,
            reserved  : 0,
            item_count: 0,
        }
    }

    fn read_buffer<E: ByteOrder>(&mut self, buffer: &[u8]) {
        let mut cursor = Cursor::new(buffer);

        self.chrom_id   = cursor.read_u32::<E>().unwrap();
        self.start      = cursor.read_u32::<E>().unwrap();
        self.end        = cursor.read_u32::<E>().unwrap();
        self.step       = cursor.read_u32::<E>().unwrap();
        self.span       = cursor.read_u32::<E>().unwrap();
        self.kind       = cursor.read_u8      ().unwrap();
        self.reserved   = cursor.read_u8      ().unwrap();
        self.item_count = cursor.read_u16::<E>().unwrap();
    }

    fn write_buffer<E: ByteOrder>(&self, buffer: &mut [u8]) {
        let mut cursor = Cursor::new(buffer);

        cursor.write_u32::<E>(self.chrom_id  ).unwrap();
        cursor.write_u32::<E>(self.start     ).unwrap();
        cursor.write_u32::<E>(self.end       ).unwrap();
        cursor.write_u32::<E>(self.step      ).unwrap();
        cursor.write_u32::<E>(self.span      ).unwrap();
        cursor.write_u8      (self.kind      ).unwrap();
        cursor.write_u8      (self.reserved  ).unwrap();
        cursor.write_u16::<E>(self.item_count).unwrap();
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone)]
struct BbiRawBlockDecoderIterator<'a> {
    decoder : &'a BbiRawBlockDecoder<'a>,
    position: usize
}

/* -------------------------------------------------------------------------- */

impl<'a> BbiRawBlockDecoderIterator<'a> {
 
    fn read<E: ByteOrder>(&self) -> BbiRawBlockDecoderType { 

        let mut r = BbiRawBlockDecoderType::new();

        match self.decoder.header.kind {
            BBI_TYPE_BED_GRAPH => {
                r.read_bed_graph::<E>(&self.decoder.header, &self.decoder.buffer[self.position-12..self.position])
            }
            BBI_TYPE_VARIABLE => {
                r.read_variable::<E>(&self.decoder.header, &self.decoder.buffer[self.position-8..self.position])
            }
            BBI_TYPE_FIXED => {
                r.read_fixed::<E>(&self.decoder.header, &self.decoder.buffer[self.position-4..self.position], self.position-4)
            }
            _ => panic!("Unsupported block type"),
        }
        r
    }

}

impl<'a> Iterator for BbiRawBlockDecoderIterator<'a> {

    type Item = Self;

    fn next(&mut self) -> Option<Self::Item> {

        if self.position >= self.decoder.buffer.len() {
            return None;
        }

        match self.decoder.header.kind {
            BBI_TYPE_BED_GRAPH => {
                self.position += 12;
            }
            BBI_TYPE_VARIABLE => {
                self.position += 8;
            }
            BBI_TYPE_FIXED => {
                self.position += 4;
            }
            _ => panic!("Unsupported block type"),
        }
        Some(self.clone())
    }
}

/* -------------------------------------------------------------------------- */

struct BbiRawBlockDecoderType(BbiSummaryRecord);

/* -------------------------------------------------------------------------- */

impl BbiRawBlockDecoderType {

    fn new() -> BbiRawBlockDecoderType {
        BbiRawBlockDecoderType(
            BbiSummaryRecord::new()
        )
    }

    fn read_fixed<E: ByteOrder>(&mut self, header: &BbiDataHeader, buffer: &[u8], i: usize) -> () {
        self.chrom_id               = header.chrom_id as i32;
        self.from                   = (header.start + (i as u32 / 4) * header.step) as i32;
        self.to                     = self.from + header.span as i32;
        self.statistics.valid       = 1.0;
        self.statistics.sum         = f32::from_bits(E::read_u32(&buffer[0..4])) as f64;
        self.statistics.sum_squares = self.statistics.sum * self.statistics.sum;
        self.statistics.min         = self.statistics.sum;
        self.statistics.max         = self.statistics.sum;
    }

    fn read_variable<E: ByteOrder>(&mut self, header: &BbiDataHeader, buffer: &[u8]) -> () {
        self.chrom_id               = header.chrom_id as i32;
        self.from                   = E::read_u32(&buffer[0..4]) as i32;
        self.to                     = self.from + header.span as i32;
        self.statistics.valid       = 1.0;
        self.statistics.sum         = f32::from_bits(E::read_u32(&buffer[4..8])) as f64;
        self.statistics.sum_squares = self.statistics.sum * self.statistics.sum;
        self.statistics.min         = self.statistics.sum;
        self.statistics.max         = self.statistics.sum;
    }

    fn read_bed_graph<E: ByteOrder>(&mut self, header: &BbiDataHeader, buffer: &[u8]) -> () {
        self.chrom_id               = header.chrom_id as i32;
        self.from                   = E::read_u32(&buffer[0..4]) as i32;
        self.to                     = E::read_u32(&buffer[4..8]) as i32;
        self.statistics.valid       = 1.0;
        self.statistics.sum         = f32::from_bits(E::read_u32(&buffer[8..12])) as f64;
        self.statistics.sum_squares = self.statistics.sum * self.statistics.sum;
        self.statistics.min         = self.statistics.sum;
        self.statistics.max         = self.statistics.sum;
    }
}

impl std::ops::Deref for BbiRawBlockDecoderType {
    type Target = BbiSummaryRecord;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for BbiRawBlockDecoderType {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

/* -------------------------------------------------------------------------- */

struct BbiRawBlockDecoder<'a> {
    header: BbiDataHeader,
    buffer: &'a [u8],
}

/* -------------------------------------------------------------------------- */

impl<'a> BbiRawBlockDecoder<'a> {
    fn new<E: ByteOrder>(buffer: &'a [u8]) -> io::Result<BbiRawBlockDecoder<'a>> {
        if buffer.len() < 24 {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "block length is shorter than 24 bytes"));
        }
        let mut decoder = BbiRawBlockDecoder {
            header: BbiDataHeader::new(),
            buffer,
        }; 
        decoder.header.read_buffer::<E>(buffer);
        match decoder.header.kind {
            BBI_TYPE_BED_GRAPH => {
                if buffer.len() % 12 != 0 {
                    return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "bedGraph data block has invalid length"));
                }
            }
            BBI_TYPE_VARIABLE => {
                if buffer.len() % 8 != 0 {
                    return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "variable step data block has invalid length"));
                }
            }
            BBI_TYPE_FIXED => {
                if buffer.len() % 4 != 0 {
                    return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "fixed step data block has invalid length"));
                }
            }
            _ => return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "unsupported block type")),
        }
        Ok(decoder)
    }
    fn decode(&self) -> BbiRawBlockDecoderIterator {
        let iterator = BbiRawBlockDecoderIterator {
            decoder : self,
            position: 0
        };
        iterator
    }
}

/* -------------------------------------------------------------------------- */

struct BbiZoomBlockDecoderType(BbiSummaryRecord);

/* -------------------------------------------------------------------------- */

impl BbiZoomBlockDecoderType {

    const LENGTH : usize = 52;

    fn new() -> BbiZoomBlockDecoderType {
        BbiZoomBlockDecoderType(
            BbiSummaryRecord::new()
        )
    }

    fn read_buffer<E: ByteOrder>(&mut self, buffer : &[u8]) -> io::Result<()> {
        let mut cursor = Cursor::new(buffer);

        self.chrom_id               = cursor.read_i32::<E>()?;
        self.from                   = cursor.read_i32::<E>()?;
        self.to                     = cursor.read_i32::<E>()?;
        self.statistics.valid       = cursor.read_f64::<E>()?;
        self.statistics.min         = cursor.read_f64::<E>()?;
        self.statistics.max         = cursor.read_f64::<E>()?;
        self.statistics.sum         = cursor.read_f64::<E>()?;
        self.statistics.sum_squares = cursor.read_f64::<E>()?;
        Ok(())
    }
}

impl std::ops::Deref for BbiZoomBlockDecoderType {
    type Target = BbiSummaryRecord;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for BbiZoomBlockDecoderType {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

/* -------------------------------------------------------------------------- */

struct BbiZoomBlockDecoder<'a> {
    buffer: &'a [u8],
}

/* -------------------------------------------------------------------------- */

impl<'a> BbiZoomBlockDecoder<'a> {
    fn new(buffer: &'a [u8]) -> BbiZoomBlockDecoder<'a> {
        BbiZoomBlockDecoder {
            buffer,
        }
    }

    fn decode(&self) -> BbiZoomBlockDecoderIterator {
        let iterator = BbiZoomBlockDecoderIterator {
            decoder : self,
            position: 0
        };
        iterator
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone)]
struct BbiZoomBlockDecoderIterator<'a> {
    decoder : &'a BbiZoomBlockDecoder<'a>,
    position: usize
}

/* -------------------------------------------------------------------------- */

impl<'a> BbiZoomBlockDecoderIterator<'a> {

    fn read<E: ByteOrder>(&mut self) -> io::Result<BbiZoomBlockDecoderType> {

        let mut r = BbiZoomBlockDecoderType::new();

        let start = self.position - BbiZoomBlockDecoderType::LENGTH;
        let end   = self.position;

        if let Err(v) = r.read_buffer::<E>(&self.decoder.buffer[start..end]) {
            Err(v)
        } else {
            Ok(r)
        }
    }

}

/* -------------------------------------------------------------------------- */

impl<'a> Iterator for BbiZoomBlockDecoderIterator<'a> {

    type Item = Self;

    fn next(&mut self) -> Option<Self::Item> {

        if self.position >= self.decoder.buffer.len() {
            return None
        }

        self.position += BbiZoomBlockDecoderType::LENGTH;

        Some(self.clone())
    }
}

/* Result type of both the raw and zoom block encoder
 * -------------------------------------------------------------------------- */

#[derive(Default)]
struct BbiBlockEncoderType {
    from : usize,
    to   : usize,
    block: Vec<u8>,
}

/* -------------------------------------------------------------------------- */

#[derive(Clone)]
struct BbiRawBlockEncoder {
    items_per_slot: usize,
    fixed_step    : bool,
}

/* -------------------------------------------------------------------------- */

impl BbiRawBlockEncoder {
    fn new(items_per_slot: usize, fixed_step: bool) -> Self {
        BbiRawBlockEncoder {
            items_per_slot,
            fixed_step,
        }
    }

    fn encode_variable<E: ByteOrder>(&self, buffer: &mut [u8], position: u32, value: f64) {
        E::write_u32_into(&[position, value.to_bits() as u32], &mut buffer[0..8]);
    }

    fn encode_fixed<E: ByteOrder>(&self, buffer: &mut [u8], value: f64) {
        E::write_u32(buffer, value.to_bits() as u32);
    }

    fn encode(&self, chrom_id: usize, sequence: Vec<f64>, bin_size: usize) -> BbiRawBlockEncoderIterator {
        BbiRawBlockEncoderIterator {
            encoder: Box::new(self.clone()),
            chrom_id,
            sequence,
            bin_size,
            position: 0,
            seqbuf  : Vec::new(),
            header  : BbiDataHeader::new(),
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone)]
struct BbiRawBlockEncoderIterator {
    encoder : Box<BbiRawBlockEncoder>,
    chrom_id: usize,
    sequence: Vec<f64>,
    bin_size: usize,
    position: usize,
    seqbuf  : Vec<f64>,
    header  : BbiDataHeader,
}

/* -------------------------------------------------------------------------- */

impl BbiRawBlockEncoderIterator {

    fn write<E: ByteOrder>(&self) -> io::Result<BbiBlockEncoderType> {

        let mut buffer = Cursor::new(Vec::new());
        let mut tmp    = vec![0u8; 24];

        self.header.write_buffer::<E>(&mut tmp);
        buffer.write_all(&tmp).unwrap();

        if self.encoder.fixed_step {
            for entry in &self.seqbuf {
                self.encoder.encode_fixed::<E>(&mut tmp, *entry);
                buffer.write_all(&tmp[..4])?;
            }
        } else {
            for entry in &self.seqbuf {
                self.encoder.encode_variable::<E>(&mut tmp, self.header.end, *entry);
                buffer.write_all(&tmp[..8])?;
            }
        }

        Ok(BbiBlockEncoderType{
            from : self.header.start as usize,
            to   : self.header.end   as usize,
            block: buffer.into_inner(),
        })
    }

}

/* -------------------------------------------------------------------------- */

impl Iterator for BbiRawBlockEncoderIterator {

    type Item = Self;

    fn next(&mut self) -> Option<Self::Item> {

        while self.position < self.sequence.len() && self.sequence[self.position].is_nan() {
            self.position += 1;
        }

        self.seqbuf.clear();

        self.header = BbiDataHeader {
            chrom_id  : self.chrom_id as u32,
            start     : (self.bin_size * self.position) as u32,
            end       : (self.bin_size * self.position) as u32,
            step      : self.bin_size as u32,
            span      : self.bin_size as u32,
            item_count: 0,
            kind      : if self.encoder.fixed_step { 3 } else { 2 },
            reserved  : 0,
        };

        if self.encoder.fixed_step {
            while self.position < self.sequence.len() {
                if self.sequence[self.position].is_nan() {
                    while self.position < self.sequence.len() && self.sequence[self.position].is_nan() {
                        self.position += 1;
                    }
                    break;
                }
                self.seqbuf.push(self.sequence[self.position]);
                self.header.item_count += 1;
                self.header.end        += self.header.step;
                if self.header.item_count as usize == self.encoder.items_per_slot {
                    self.position += 1;
                    break;
                }
                self.position += 1;
            }
        } else {
            while self.position < self.sequence.len() {
                if !self.sequence[self.position].is_nan() {
                    self.seqbuf.push(self.sequence[self.position]);
                    self.header.item_count += 1;
                    self.header.end         = (self.bin_size * self.position) as u32 + self.header.step;
                }
                if self.header.item_count as usize == self.encoder.items_per_slot {
                    self.position += 1;
                    break;
                }
                self.position += 1;
            }
        }

        if self.seqbuf.len() > 0 {
            Some(self.clone())
        }
        else {
            None
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone)]
struct BbiZoomBlockEncoder {
    items_per_slot : usize,
    tmp            : Vec<u8>,
    reduction_level: usize,
}

/* -------------------------------------------------------------------------- */

#[derive(Clone)]
struct BbiZoomBlockEncoderIterator {
    encoder : Box<BbiZoomBlockEncoder>,
    chrom_id: usize,
    sequence: Vec<f64>,
    bin_size: usize,
    position: usize,
    records : Vec<BbiZoomRecord>,
    from    : i64,
    to      : i64,
    count   : usize,
}

/* -------------------------------------------------------------------------- */

impl BbiZoomBlockEncoder {

    fn new(items_per_slot: usize, reduction_level: usize) -> Self {
        BbiZoomBlockEncoder {
            items_per_slot,
            tmp: Vec::new(),
            reduction_level,
        }
    }

    fn encode(&self, chrom_id: usize, sequence: Vec<f64>, bin_size: usize) -> BbiZoomBlockEncoderIterator {
        BbiZoomBlockEncoderIterator {
            encoder  : Box::new(self.clone()),
            chrom_id : chrom_id,
            sequence : sequence,
            bin_size : bin_size,
            position : 0,
            records  : vec![],
            from     : 0,
            to       : 0,
            count    : 0,
        }
    }
}

/* -------------------------------------------------------------------------- */

impl BbiZoomBlockEncoderIterator {

    fn write<E: ByteOrder>(&mut self) -> io::Result<BbiBlockEncoderType> {

        let mut buffer = Vec::new();

        for record in &self.records {
            record.write_buffer::<E>(&mut buffer)?;
        }

        Ok(BbiBlockEncoderType{
            from : self.from as usize,
            to   : self.to   as usize,
            block: buffer,
        })
    }
}

/* -------------------------------------------------------------------------- */

impl Iterator for BbiZoomBlockEncoderIterator {

    type Item = Self;

    fn next(&mut self) -> Option<Self::Item> {

        let n = (self.encoder.reduction_level + self.bin_size - 1) / self.bin_size;

        self.from  = -1;
        self.to    = -1;
        self.count =  0;

        self.records.clear();

        for p in (self.position..self.bin_size * self.sequence.len()).step_by(self.encoder.reduction_level) {

            let i = p / self.bin_size;
            let mut record = BbiZoomRecord::default();

            record.chrom_id = self.chrom_id as u32;
            record.start    = p as u32;
            record.end      = (p + self.encoder.reduction_level) as u32;
            record.min      = f32::NAN;
            record.max      = f32::NAN;

            if record.end > (self.bin_size * self.sequence.len()) as u32 {
                record.end = (self.bin_size * self.sequence.len()) as u32;
            }

            for j in 0..n {
                if i + j < self.sequence.len() {
                    record.add_value(self.sequence[i + j]);
                }
            }

            if record.valid > 0 {

                self.records.push(record);

                if self.from == -1 {
                    self.from = record.start as i64;
                }
                self.to     = record.end as i64;
                self.count += 1;
            }

            if self.count == self.encoder.items_per_slot || p + self.encoder.reduction_level >= self.bin_size * self.sequence.len() {
  
                self.position = p + self.encoder.reduction_level;
                    
                return Some(self.clone())
            }
        }
        None
    }
}

/* -------------------------------------------------------------------------- */

struct BVertex {
    is_leaf : u8,
    keys    : Vec<Vec<u8>>,
    values  : Vec<Vec<u8>>,
    children: Vec<BVertex>,
}

/* -------------------------------------------------------------------------- */

impl BVertex {
    fn build_tree(&mut self, data: &BData, from: usize, to: usize, level: i32) -> Result<usize, String> {

        let mut i = 0;

        if level == 0 {
            self.is_leaf = 1;
            while i < data.items_per_block as usize && from + i < to {
                if data.keys[from + i].len() != data.key_size as usize {
                    return Err(format!("key number `{}` has invalid size", i));
                }
                if data.values[from + i].len() != data.value_size as usize {
                    return Err(format!("value number `{}` has invalid size", i));
                }
                self.keys  .push(data.keys  [from + i].clone());
                self.values.push(data.values[from + i].clone());
                i += 1;
            }
        } else {
            self.is_leaf = 0;
            while i < data.items_per_block as usize && from + i < to {
                self.keys.push(data.keys[from + i].clone());
                let mut child = BVertex {
                    is_leaf: 0,
                    keys    : Vec::new(),
                    values  : Vec::new(),
                    children: Vec::new(),
                };
                let j = child.build_tree(data, from + i, to, level - 1)?;
                self.children.push(child);
                i += j as usize;
            }
        }
        Ok(i)
    }

    fn write_leaf<E: ByteOrder, W: Write>(&self, writer: &mut W) -> io::Result<()> {

        let padding = 0u8;
        let n_vals  = self.keys.len() as u16;

        writer.write_u8(self.is_leaf)?;
        writer.write_u8(padding)?;
        writer.write_u16::<E>(n_vals)?;
        for i in 0..self.keys.len() {
            writer.write_all(&self.keys[i])?;
            writer.write_all(&self.values[i])?;
        }
        Ok(())
    }

    fn write_index<E: ByteOrder, W: Write + Seek>(&self, writer: &mut W) -> io::Result<()> {

        let is_leaf = 0u8;
        let padding = 0u8;
        let n_vals  = self.keys.len() as u16;

        let mut offsets = Vec::new();

        writer.write_u8(is_leaf)?;
        writer.write_u8(padding)?;
        writer.write_u16::<byteorder::LittleEndian>(n_vals)?;

        for i in 0..self.keys.len() {
            writer .write_all(&self.keys[i])?;
            offsets.push(writer.seek(io::SeekFrom::Current(0))?);
            writer .write_u64::<byteorder::LittleEndian>(0)?;
        }
        for i in 0..self.keys.len() {
            let offset = writer.seek(io::SeekFrom::Current(0))? as u64;
            writer.seek(io::SeekFrom::Start(offsets[i]))?;
            writer.write_u64::<E>(offset)?;
            writer.seek(io::SeekFrom::Start(offset))?;
            self.children[i].write::<E, W>(writer)?;
        }
        Ok(())
    }

    fn write<E: ByteOrder, W: Write + Seek>(&self, writer: &mut W) -> io::Result<()> {
        if self.is_leaf != 0 {
            self.write_leaf::<E, W>(writer)
        } else {
            self.write_index::<E, W>(writer)
        }
    }
}

/* -------------------------------------------------------------------------- */

struct BTree {
    key_size       : u32,
    value_size     : u32,
    items_per_block: u32,
    item_count     : u64,
    root           : BVertex,
}

/* -------------------------------------------------------------------------- */

impl BTree {

    fn new(data: &BData) -> Self {
        let mut tree = BTree {
            key_size       : data.key_size,
            value_size     : data.value_size,
            items_per_block: data.items_per_block,
            item_count: data.item_count,
            root: BVertex {
                is_leaf : 0,
                keys    : Vec::new(),
                values  : Vec::new(),
                children: Vec::new(),
            },
        };
        if data.item_count == 1 {
            tree.root.build_tree(data, 0, data.item_count as usize, 0).unwrap();
        } else {
            let d = ((data.item_count as f64).log(data.items_per_block as f64).ceil()) as i32;
            tree.root.build_tree(data, 0, data.item_count as usize, d - 1).unwrap();
        }
        tree
    }

    fn write<E: ByteOrder, W: Write+Seek>(&self, writer: &mut W) -> io::Result<()> {
        let magic = CIRTREE_MAGIC;

        writer.write_u32::<E>(magic)?;
        writer.write_u32::<E>(self.items_per_block)?;
        writer.write_u32::<E>(self.key_size)?;
        writer.write_u32::<E>(self.value_size)?;
        writer.write_u64::<E>(self.item_count)?;
        writer.write_u64::<E>(0)?;
        self.root.write::<E, W>(writer)?;
        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

struct BData {
    key_size       : u32,
    value_size     : u32,
    items_per_block: u32,
    item_count     : u64,
    keys           : Vec<Vec<u8>>,
    values         : Vec<Vec<u8>>,
    ptr_keys       : Vec<i64>,
    ptr_values     : Vec<i64>,
}

/* -------------------------------------------------------------------------- */

impl BData {
    fn new() -> Self {
        BData {
            key_size        : 0,
            value_size      : 0,
            items_per_block : 0,
            item_count      : 0,
            keys            : Vec::new(),
            values          : Vec::new(),
            ptr_keys        : Vec::new(),
            ptr_values      : Vec::new(),
        }
    }

    fn add(&mut self, key: Vec<u8>, value: Vec<u8>) -> Result<(), String> {
        if key.len() as u32 != self.key_size {
            return Err("BData.Add(): key has invalid length".to_string());
        }
        if value.len() as u32 != self.value_size {
            return Err("BData.Add(): value has invalid length".to_string());
        }
        self.keys.push(key);
        self.values.push(value);
        self.items_per_block += 1;
        self.item_count += 1;
        Ok(())
    }

    fn read_vertex_leaf<E: ByteOrder, R: Read + Seek>(&mut self, file: &mut R) -> io::Result<()> {
        let n_vals = file.read_u16::<E>()?;
        for _ in 0..n_vals {
            let mut key = vec![0; self.key_size as usize];
            let mut value = vec![0; self.value_size as usize];

            let ptr_key = file.seek(SeekFrom::Current(0))?;
            file.read_exact(&mut key)?;
            let ptr_value = file.seek(SeekFrom::Current(0))?;
            file.read_exact(&mut value)?;

            self.keys.push(key);
            self.values.push(value);
            self.ptr_keys.push(ptr_key as i64);
            self.ptr_values.push(ptr_value as i64);
        }
        Ok(())
    }

    fn read_vertex_index<E: ByteOrder, R: Read + Seek>(&mut self, file: &mut R) -> io::Result<()> {
        let n_vals = file.read_u16::<E>()?;
        for _ in 0..n_vals {
            let mut key = vec![0; self.key_size as usize];
            let position = file.read_u64::<E>()?;

            file.read_exact(&mut key)?;

            // save current position and jump to child vertex
            let current_position = file.seek(SeekFrom::Current(0))?;
            file.seek(SeekFrom::Start(position as u64))?;
            self.read_vertex::<E, R>(file)?;
            file.seek(SeekFrom::Start(current_position))?;
        }
        Ok(())
    }

    fn read_vertex<E: ByteOrder, R: Read + Seek>(&mut self, file: &mut R) -> io::Result<()> {
        let is_leaf = file.read_u8()?;
        file.read_u8()?; // padding
        if is_leaf != 0 {
            self.read_vertex_leaf::<E, R>(file)
        } else {
            self.read_vertex_index::<E, R>(file)
        }
    }

    fn read<E: ByteOrder, R: Read + Seek>(&mut self, file: &mut R) -> io::Result<()> {
        let magic = file.read_u32::<E>()?;
        if magic != CIRTREE_MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid tree"));
        }

        self.items_per_block = file.read_u32::<E>()?;
        self.key_size = file.read_u32::<E>()?;
        self.value_size = file.read_u32::<E>()?;
        self.item_count = file.read_u64::<E>()?;

        file.read_u32::<E>()?; // padding
        file.read_u32::<E>()?; // padding

        self.read_vertex::<E, R>(file)
    }

    fn write<E: ByteOrder, W: Write+Seek>(&self, file: &mut W) -> io::Result<()> {
        let tree = BTree::new(self);
        tree.write::<E, W>(file)
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Default)]
struct BbiHeaderZoom {
    reduction_level : u32,
    reserved        : u32,
    data_offset     : u64,
    index_offset    : u64,
    n_blocks        : u32,
    ptr_data_offset : u64,
    ptr_index_offset: u64,
}

/* -------------------------------------------------------------------------- */

impl BbiHeaderZoom {
    fn read<E: ByteOrder, R: Read + Seek>(&mut self, file: &mut R) -> io::Result<()> {
        let mut buf = [0u8; 4];

        self.reduction_level = file.read_u32::<E>()?;
        self.reserved = file.read_u32::<E>()?;

        self.ptr_data_offset = file.seek(SeekFrom::Current(0))?;
        self.data_offset = file.read_u64::<E>()?;

        self.ptr_index_offset = file.seek(SeekFrom::Current(0))?;
        self.index_offset = file.read_u64::<E>()?;

        file_read_at(file, self.data_offset, &mut buf)?;
        self.n_blocks = E::read_u32(&buf);

        Ok(())
    }

    fn write<E: ByteOrder, W: Write + Seek>(&mut self, file: &mut W) -> io::Result<()> {
        file.write_u32::<E>(self.reduction_level)?;
        file.write_u32::<E>(self.reserved)?;

        self.ptr_data_offset = file.seek(SeekFrom::Current(0))?;
        file.write_u64::<E>(self.data_offset)?;

        self.ptr_index_offset = file.seek(SeekFrom::Current(0))?;
        file.write_u64::<E>(self.index_offset)?;

        Ok(())
    }

    fn write_offsets<E: ByteOrder, W: Write + Seek>(&self, file: &mut W) -> io::Result<()> {
        let mut buf = [0u8; 8];

        if self.ptr_data_offset != 0 {
            E::write_u64(&mut buf, self.data_offset);

            file_write_at(file, self.ptr_data_offset, &buf)?;
        }
        if self.ptr_index_offset != 0 {
            E::write_u64(&mut buf, self.index_offset);

            file_write_at(file, self.ptr_index_offset, &buf)?;
        }
        Ok(())
    }

    fn write_n_blocks<E: ByteOrder, W: Write + Seek>(&self, file: &mut W) -> io::Result<()> {
        let mut buf = [0u8; 4];

        E::write_u32(&mut buf, self.n_blocks);

        file_write_at(file, self.data_offset, &buf)
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Default)]
struct BbiHeader {
    magic                  : u32,
    version                : u16,
    zoom_levels            : u16,
    ct_offset              : u64,
    data_offset            : u64,
    index_offset           : u64,
    field_count            : u16,
    defined_field_count    : u16,
    sql_offset             : u64,
    summary_offset         : u64,
    uncompress_buf_size    : u32,
    extension_offset       : u64,
    n_bases_covered        : u64,
    min_val                : f64,
    max_val                : f64,
    sum_data               : f64,
    sum_squares            : f64,
    zoom_headers           : Vec<BbiHeaderZoom>,
    n_blocks               : u64,
    ptr_ct_offset          : u64,
    ptr_data_offset        : u64,
    ptr_index_offset       : u64,
    ptr_sql_offset         : u64,
    ptr_summary_offset     : u64,
    ptr_uncompress_buf_size: u64,
    ptr_extension_offset   : u64,
}

/* -------------------------------------------------------------------------- */

impl BbiHeader {
    fn new() -> Self {
        BbiHeader {
            version: 4,
            min_val: f64::NAN,
            max_val: f64::NAN,
            ..Default::default()
        }
    }

    fn summary_add_value(&mut self, x: f64, n: i32) {
        if x.is_nan() {
            return;
        }

        if self.min_val.is_nan() || self.min_val > x {
            self.min_val = x;
        }

        if self.max_val.is_nan() || self.max_val < x {
            self.max_val = x;
        }

        self.n_bases_covered += n as u64;
        self.sum_data        += x;
        self.sum_squares     += x * x;
    }

    fn read_order<E: ByteOrder, R: Read + Seek>(&mut self, file: &mut R) -> io::Result<()> {

        let mut buf = [0u8; 4];

        self.version                 = file.read_u16::<E>()?;
        self.zoom_levels             = file.read_u16::<E>()?;
        self.ptr_ct_offset           = file.seek(SeekFrom::Current(0))?;
        self.ct_offset               = file.read_u64::<E>()?;
        self.ptr_data_offset         = file.seek(SeekFrom::Current(0))?;
        self.data_offset             = file.read_u64::<E>()?;
        self.ptr_index_offset        = file.seek(SeekFrom::Current(0))?;
        self.index_offset            = file.read_u64::<E>()?;
        self.field_count             = file.read_u16::<E>()?;
        self.defined_field_count     = file.read_u16::<E>()?;
        self.ptr_sql_offset          = file.seek(SeekFrom::Current(0))?;
        self.sql_offset              = file.read_u64::<E>()?;
        self.ptr_summary_offset      = file.seek(SeekFrom::Current(0))?;
        self.summary_offset          = file.read_u64::<E>()?;
        self.ptr_uncompress_buf_size = file.seek(SeekFrom::Current(0))?;
        self.uncompress_buf_size     = file.read_u32::<E>()?;
        self.ptr_extension_offset    = file.seek(SeekFrom::Current(0))?;
        self.extension_offset        = file.read_u64::<E>()?;

        self.zoom_headers = Vec::with_capacity(self.zoom_levels as usize);
        for _ in 0..self.zoom_levels {
            let mut zoom_header = BbiHeaderZoom::default();
            zoom_header.read::<E, R>(file)?;
            self.zoom_headers.push(zoom_header);
        }

        if self.summary_offset > 0 {
            file.seek(SeekFrom::Start(self.summary_offset))?;
            self.n_bases_covered = file.read_u64::<E>()?;
            self.min_val         = file.read_f64::<E>()?;
            self.max_val         = file.read_f64::<E>()?;
            self.sum_data        = file.read_f64::<E>()?;
            self.sum_squares     = file.read_f64::<E>()?;
        }

        file_read_at(file, self.data_offset, &mut buf)?;
        self.n_blocks = E::read_u64(&buf);

        Ok(())
    }

    fn read<R: Read + Seek>(&mut self, file: &mut R, magic: u32) -> io::Result<()> {

        // Read magin number
        self.magic = file.read_u32::<LittleEndian>()?;

        if self.magic == magic {

            self.read_order::<LittleEndian, R>(file)?;

        } else {

            file.seek(SeekFrom::Current(-4))?;
            self.magic = file.read_u32::<BigEndian>()?;

            if self.magic != magic {
                return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid magic number"));
            }

            self.read_order::<BigEndian, R>(file)?;

        }

        Ok(())
    }

    fn write_offsets<E: ByteOrder, W: Write + Seek>(&self, file: &mut W) -> io::Result<()> {

        let mut buf = [0u8; 8];

        if self.ptr_ct_offset != 0 {
            E::write_u64(&mut buf, self.ct_offset);
            file_write_at(file, self.ptr_ct_offset, &buf)?;
        }
        if self.ptr_data_offset != 0 {
            E::write_u64(&mut buf, self.data_offset);
            file_write_at(file, self.ptr_data_offset, &buf)?;
        }
        if self.ptr_index_offset != 0 {
            E::write_u64(&mut buf, self.index_offset);
            file_write_at(file, self.ptr_index_offset, &buf)?;
        }
        if self.ptr_sql_offset != 0 {
            E::write_u64(&mut buf, self.sql_offset);
            file_write_at(file, self.ptr_sql_offset, &buf)?;
        }
        if self.ptr_extension_offset != 0 {
            E::write_u64(&mut buf, self.extension_offset);
            file_write_at(file, self.ptr_extension_offset, &buf)?;
        }
        Ok(())
    }

    fn write_uncompress_buf_size<E: ByteOrder, W: Write + Seek>(&self, file: &mut W) -> std::io::Result<()> {

        let mut buf = [0u8; 4];

        if self.ptr_uncompress_buf_size != 0 {
            E::write_u32(&mut buf, self.uncompress_buf_size);
            file_write_at(file, self.ptr_uncompress_buf_size, &buf)?;
        }
        Ok(())
    }

    fn write<E: ByteOrder, W: Write + Seek>(&mut self, file: &mut W) -> std::io::Result<()> {
        file.write_u32::<E>(self.magic)?;
        file.write_u16::<E>(self.version)?;
        file.write_u16::<E>(self.zoom_levels)?;

        self.ptr_ct_offset = file.seek(SeekFrom::Current(0))?;
        file.write_u64::<E>(self.ct_offset)?;

        self.ptr_data_offset = file.seek(SeekFrom::Current(0))?;
        file.write_u64::<E>(self.data_offset)?;

        self.ptr_index_offset = file.seek(SeekFrom::Current(0))?;
        file.write_u64::<E>(self.index_offset)?;

        file.write_u16::<E>(self.field_count)?;
        file.write_u16::<E>(self.defined_field_count)?;

        self.ptr_sql_offset = file.seek(SeekFrom::Current(0))?;
        file.write_u64::<E>(self.sql_offset)?;

        self.ptr_summary_offset = file.seek(SeekFrom::Current(0))?;
        file.write_u64::<E>(self.summary_offset)?;

        self.ptr_uncompress_buf_size = file.seek(SeekFrom::Current(0))?;
        file.write_u32::<E>(self.uncompress_buf_size)?;

        self.ptr_extension_offset = file.seek(SeekFrom::Current(0))?;
        file.write_u64::<E>(self.extension_offset)?;

        for zoom_header in &mut self.zoom_headers {
            zoom_header.write::<E, W>(file)?;
        }

        if self.summary_offset > 0 {
            file.seek(SeekFrom::Start(self.summary_offset))?;
            file.write_u64::<E>(self.n_bases_covered)?;
            file.write_f64::<E>(self.min_val)?;
            file.write_f64::<E>(self.max_val)?;
            file.write_f64::<E>(self.sum_data)?;
            file.write_f64::<E>(self.sum_squares)?;
        }

        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Default)]
struct RTree {
    block_size      : u32,
    n_items         : u64,
    chr_idx_start   : u32,
    base_start      : u32,
    chr_idx_end     : u32,
    base_end        : u32,
    idx_size        : u64,
    n_items_per_slot: u32,
    root            : Option<Box<RVertex>>,
    ptr_idx_size    : i64,
}

/* -------------------------------------------------------------------------- */

impl RTree {
    fn new() -> Self {
        RTree {
            block_size: 256,
            n_items_per_slot: 1024,
            ..Default::default()
        }
    }

    fn is_nil(&self) -> bool {
        self.block_size == 0
    }

    fn read<E: ByteOrder, W: Read + Seek>(&mut self, file: &mut W) -> io::Result<()> {
        let magic = file.read_u32::<E>()?;
        if magic != IDX_MAGIC {
            return Err(io::Error::other("invalid bbi tree"));
        }

        self.block_size       = file.read_u32::<E>()?;
        self.n_items          = file.read_u64::<E>()?;
        self.chr_idx_start    = file.read_u32::<E>()?;
        self.base_start       = file.read_u32::<E>()?;
        self.chr_idx_end      = file.read_u32::<E>()?;
        self.base_end         = file.read_u32::<E>()?;
        self.ptr_idx_size     = file.seek(SeekFrom::Current(0))? as i64;
        self.idx_size         = file.read_u64::<E>()?;
        self.n_items_per_slot = file.read_u32::<E>()?;

        file.read_u32::<E>()?; // Padding

        let mut root = Box::new(RVertex::default());
        root.read::<E, W>(file)?;
        self.root = Some(root);

        Ok(())
    }

    fn write_size<E: ByteOrder, W: Write + Seek>(&self, file: &mut W) -> io::Result<()> {
        file.seek(SeekFrom::Start(self.ptr_idx_size as u64))?;
        file.write_u64::<E>(self.idx_size)?;
        Ok(())
    }

    fn write<E: ByteOrder, W: Write + Seek>(&mut self, file: &mut W) -> io::Result<()> {
        let offset_start = file.seek(SeekFrom::Current(0))?;

        file.write_u32::<E>(IDX_MAGIC)?;
        file.write_u32::<E>(self.block_size)?;
        file.write_u64::<E>(self.n_items)?;
        file.write_u32::<E>(self.chr_idx_start)?;
        file.write_u32::<E>(self.base_start)?;
        file.write_u32::<E>(self.chr_idx_end)?;
        file.write_u32::<E>(self.base_end)?;

        self.ptr_idx_size = file.seek(SeekFrom::Current(0))? as i64;
        file.write_u64::<E>(self.idx_size)?;
        file.write_u32::<E>(self.n_items_per_slot)?;

        file.write_u32::<E>(0)?; // Padding

        if let Some(ref mut root) = self.root {
            root.write::<E, W>(file)?;
        }

        let offset_end = file.seek(SeekFrom::Current(0))?;
        self.idx_size = offset_end - offset_start;

        self.write_size::<E, W>(file)?;

        Ok(())
    }

    fn build_tree_rec(&self, mut leaves: Vec<Box<RVertex>>, level: usize) -> (Option<Box<RVertex>>, Vec<Box<RVertex>>) {
        let mut v = Box::new(RVertex::default());
        let n = leaves.len();

        if n == 0 {
            return (None, leaves);
        }

        if level == 0 {
            let n = n.min(self.block_size as usize);
            v.n_children = n as u16;
            v.children = leaves.drain(0..n).collect();
        } else {
            for _ in 0..self.block_size as usize {
                if leaves.is_empty() {
                    break;
                }
                let (vertex, remaining_leaves) = self.build_tree_rec(leaves, level - 1);
                if let Some(vertex) = vertex {
                    v.n_children += 1;
                    v.children.push(vertex);
                }
                leaves = remaining_leaves;
            }
        }

        for child in &v.children {
            v.chr_idx_start.push(child.chr_idx_start[0]);
            v.chr_idx_end.push(child.chr_idx_end[child.n_children as usize - 1]);
            v.base_start.push(child.base_start[0]);
            v.base_end.push(child.base_end[child.n_children as usize - 1]);
        }

        (Some(v), leaves)
    }

    fn build_tree(&mut self, leaves: Vec<Box<RVertex>>) -> io::Result<()> {
        if leaves.is_empty() {
            return Ok(());
        }

        if leaves.len() == 1 {
            self.root = Some(leaves.into_iter().next().unwrap());
        } else {
            let depth = ((leaves.len() as f64).ln() / (self.block_size as f64).ln()).ceil() as usize;
            let (root, remaining_leaves) = self.build_tree_rec(leaves, depth - 1);

            if !remaining_leaves.is_empty() {
                return Err(io::Error::other("internal error"));
            }

            self.root = root;
        }

        if let Some(ref root) = self.root {
            self.chr_idx_start = root.chr_idx_start[0];
            self.chr_idx_end = root.chr_idx_end[root.n_children as usize - 1];
            self.base_start = root.base_start[0];
            self.base_end = root.base_end[root.n_children as usize - 1];
        }

        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Default, Debug)]
struct RVertex {
    is_leaf        : u8,
    n_children     : u16,
    chr_idx_start  : Vec<u32>,
    base_start     : Vec<u32>,
    chr_idx_end    : Vec<u32>,
    base_end       : Vec<u32>,
    data_offset    : Vec<u64>,
    sizes          : Vec<u64>,
    children       : Vec<Box<RVertex>>,
    ptr_data_offset: Vec<i64>,
    ptr_sizes      : Vec<i64>,
}

/* -------------------------------------------------------------------------- */

impl RVertex {
    fn read_block<R: Read + Seek>(&self, reader: &mut R, bwf: &BbiFile, i: usize) -> io::Result<Vec<u8>> {
        let mut block = vec![0u8; self.sizes[i] as usize];

        reader.seek(SeekFrom::Start(self.data_offset[i]))?;
        reader.read_exact(&mut block)?;

        if bwf.header.uncompress_buf_size != 0 {
            block = uncompress_slice(&block)?;
        }

        Ok(block)
    }

    fn write_block<E: ByteOrder, W: Write + Seek>(&mut self, writer: &mut W, bwf: &mut BbiFile, i: usize, mut block: Vec<u8>) -> io::Result<()> {
        if bwf.header.uncompress_buf_size != 0 {
            if block.len() as u32 > bwf.header.uncompress_buf_size {
                bwf.header.uncompress_buf_size = block.len() as u32;
                bwf.header.write_uncompress_buf_size::<E, W>(writer)?;
            }
            block = compress_slice(&block)?;
        }

        let offset = writer.seek(SeekFrom::Current(0))?;
        self.data_offset[i] = offset;
        if self.ptr_data_offset[i] != 0 {
            writer.seek(SeekFrom::Start(self.ptr_data_offset[i] as u64))?;
            writer.write_u64::<E>(self.data_offset[i])?;
        }

        writer.seek(SeekFrom::Start(offset))?;
        writer.write_all(&block)?;

        self.sizes[i] = block.len() as u64;
        if self.ptr_sizes[i] != 0 {
            writer.seek(SeekFrom::Start(self.ptr_sizes[i] as u64))?;
            writer.write_u64::<E>(self.sizes[i])?;
        }

        Ok(())
    }

    fn read<E: ByteOrder, R: Read + Seek>(&mut self, file: &mut R) -> io::Result<()> {
        let mut padding = [0u8; 1];

        self.is_leaf    = file.read_u8()?;
        file.read_exact(&mut padding)?;
        self.n_children = file.read_u16::<E>()?;

        self.chr_idx_start  .resize(self.n_children as usize, 0);
        self.base_start     .resize(self.n_children as usize, 0);
        self.chr_idx_end    .resize(self.n_children as usize, 0);
        self.base_end       .resize(self.n_children as usize, 0);
        self.data_offset    .resize(self.n_children as usize, 0);
        self.ptr_data_offset.resize(self.n_children as usize, 0);

        if self.is_leaf != 0 {
            self.sizes.resize(self.n_children as usize, 0);
            self.ptr_sizes.resize(self.n_children as usize, 0);
        }

        for i in 0..self.n_children as usize {
            self.chr_idx_start[i] = file.read_u32::<E>()?;
            self.base_start   [i] = file.read_u32::<E>()?;
            self.chr_idx_end  [i] = file.read_u32::<E>()?;
            self.base_end     [i] = file.read_u32::<E>()?;

            let offset = file.seek(SeekFrom::Current(0))?;
            self.ptr_data_offset[i] = offset as i64;
            self.data_offset    [i] = file.read_u64::<E>()?;

            if self.is_leaf != 0 {
                let offset = file.seek(SeekFrom::Current(0))?;
                self.ptr_sizes[i] = offset as i64;
                self.sizes    [i] = file.read_u64::<E>()?;
            }
        }

        if self.is_leaf == 0 {
            for i in 0..self.n_children as usize {
                file.seek(SeekFrom::Start(self.data_offset[i]))?;
                let mut child = Box::new(RVertex::default());
                child.read::<E, R>(file)?;
                self.children.push(child);
            }
        }
        assert_eq!(self.children.len(), self.n_children as usize);

        Ok(())
    }

    fn write<E: ByteOrder, W: Write + Seek>(&mut self, file: &mut W) -> std::io::Result<()> {
        if self.data_offset.len() != self.n_children as usize {
            self.data_offset.resize(self.n_children as usize, 0);
        }
        if self.sizes.len() != self.n_children as usize {
            self.sizes.resize(self.n_children as usize, 0);
        }
        if self.ptr_data_offset.len() != self.n_children as usize {
            self.ptr_data_offset.resize(self.n_children as usize, 0);
        }
        if self.ptr_sizes.len() != self.n_children as usize {
            self.ptr_sizes.resize(self.n_children as usize, 0);
        }

        file.write_u8(self.is_leaf)?;
        file.write_u8(0)?; // padding
        file.write_u16::<E>(self.n_children)?;

        for i in 0..self.n_children as usize {
            file.write_u32::<E>(self.chr_idx_start[i])?;
            file.write_u32::<E>(self.base_start[i])?;
            file.write_u32::<E>(self.chr_idx_end[i])?;
            file.write_u32::<E>(self.base_end[i])?;

            let offset = file.seek(SeekFrom::Current(0))?;
            self.ptr_data_offset[i] = offset as i64;
            file.write_u64::<E>(self.data_offset[i])?;

            let offset = file.seek(SeekFrom::Current(0))?;
            self.ptr_sizes[i] = offset as i64;

            if self.is_leaf != 0 {
                file.write_u64::<E>(self.sizes[i])?;
            }
        }

        if self.is_leaf == 0 {
            for i in 0..self.n_children as usize {
                let offset = file.seek(SeekFrom::Current(0))?;
                self.data_offset[i] = offset;
                file.write_u64::<E>(self.data_offset[i])?;
                self.children[i].write::<E, W>(file)?;
            }
        }

        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

struct RVertexGenerator {
    block_size: usize,
    items_per_slot: usize,
}

struct RVertexGeneratorType {
    vertex: RVertex,
    blocks: Vec<Vec<u8>>,
}

impl RVertexGenerator {
    fn new(block_size: usize, items_per_slot: usize) -> Result<Self, String> {
        if block_size == 0 {
            return Err(format!("invalid block size `{}`", block_size));
        }
        if items_per_slot == 0 {
            return Err(format!("invalid items per slot `{}`", items_per_slot));
        }
        Ok(RVertexGenerator {
            block_size,
            items_per_slot,
        })
    }

    fn generate<E: ByteOrder>(self, chrom_id: usize, sequence: Vec<f64>, bin_size: usize, reduction_level: usize, fixed_step: bool) -> Receiver<RVertexGeneratorType> {
        let (tx, rx) = channel();

        std::thread::spawn(move || {
            self.generate_impl::<E>(tx, chrom_id, sequence, bin_size, reduction_level, fixed_step).unwrap();
        });

        rx
    }

    fn generate_zoom<E: ByteOrder>(&self, tx: std::sync::mpsc::Sender<RVertexGeneratorType>, chrom_id: usize, sequence: Vec<f64>, bin_size: usize, reduction_level: usize) -> Result<(), String> {

        let encoder = BbiZoomBlockEncoder::new(
            self.items_per_slot, reduction_level
        );

        let mut vertex = RVertex::default();
        vertex.is_leaf = 1;
        let mut blocks = Vec::new();

        for mut item in encoder.encode(chrom_id, sequence, bin_size) {
            let chunk = item.write::<E>().unwrap();

            if vertex.n_children as usize == self.block_size {
                tx.send(RVertexGeneratorType { vertex, blocks }).unwrap();
                vertex = RVertex::default();
                vertex.is_leaf = 1;
                blocks = Vec::new();
            }
            vertex.chr_idx_start  .push(chrom_id   as u32);
            vertex.chr_idx_end    .push(chrom_id   as u32);
            vertex.base_start     .push(chunk.from as u32);
            vertex.base_end       .push(chunk.to   as u32);
            vertex.data_offset    .push(0);
            vertex.sizes          .push(0);
            vertex.ptr_data_offset.push(0);
            vertex.ptr_sizes      .push(0);
            vertex.n_children += 1;

            blocks.push(chunk.block);
        }

        if vertex.n_children != 0 {
            tx.send(RVertexGeneratorType { vertex, blocks }).unwrap();
        }

        Ok(())
    }

    fn generate_raw<E: ByteOrder>(&self, tx: std::sync::mpsc::Sender<RVertexGeneratorType>, chrom_id: usize, sequence: Vec<f64>, bin_size: usize, fixed_step: bool) -> Result<(), String> {

        let encoder = BbiRawBlockEncoder::new(
            self.items_per_slot, fixed_step
        );

        let mut vertex = RVertex::default();
        vertex.is_leaf = 1;
        let mut blocks = Vec::new();

        for item in encoder.encode(chrom_id, sequence, bin_size) {
            let chunk = item.write::<E>().unwrap();

            if vertex.n_children as usize == self.block_size {
                tx.send(RVertexGeneratorType { vertex, blocks }).unwrap();
                vertex = RVertex::default();
                vertex.is_leaf = 1;
                blocks = Vec::new();
            }
            vertex.chr_idx_start  .push(chrom_id   as u32);
            vertex.chr_idx_end    .push(chrom_id   as u32);
            vertex.base_start     .push(chunk.from as u32);
            vertex.base_end       .push(chunk.to   as u32);
            vertex.data_offset    .push(0);
            vertex.sizes          .push(0);
            vertex.ptr_data_offset.push(0);
            vertex.ptr_sizes      .push(0);
            vertex.n_children += 1;

            blocks.push(chunk.block);
        }

        if vertex.n_children != 0 {
            tx.send(RVertexGeneratorType { vertex, blocks }).unwrap();
        }

        Ok(())
    }

    fn generate_impl<E: ByteOrder>(&self, tx: std::sync::mpsc::Sender<RVertexGeneratorType>, chrom_id: usize, sequence: Vec<f64>, bin_size: usize, reduction_level: usize, fixed_step: bool) -> Result<(), String> {
        if reduction_level > bin_size {
            self.generate_zoom::<E>(tx, chrom_id, sequence, bin_size, reduction_level)
        } else {
            self.generate_raw::<E>(tx, chrom_id, sequence, bin_size, fixed_step)
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
struct RTreeTraverser<'a> {
    chrom_id: i32,                            // Chromosome ID for the query
    from    : i32,                            // Start of the region query
    to      : i32,                            // End of the region query
    stack   : Vec<RTreeTraverserType<'a>>,    // Stack for keeping track of the current tree position
}

#[derive(Debug)]
struct RTreeTraverserType<'a> {
    vertex: &'a RVertex,  // A reference to the current vertex
    idx: usize,           // Index within the vertex
}

impl<'a> RTreeTraverser<'a> {
    fn new(tree: &'a RTree, chrom_id: i32, from: i32, to: i32) -> Self {
        let mut traverser = RTreeTraverser {
            chrom_id,
            from,
            to,
            stack: Vec::new(),
        };
        // Push the root of the tree onto the stack and initiate traversal
        traverser.stack.push(RTreeTraverserType { vertex: &tree.root, idx: 0 });
        traverser
    }
}

impl<'a> Iterator for RTreeTraverser<'a> {

    type Item = &'a RTreeTraverserType<'a>;

    fn next(&mut self) -> Option<Self::Item> {

        if self.stack.is_empty() {
            return None;
        }

        // Loop over the stack until we find a new position or the stack is empty
        while let Some(mut top) = self.stack.pop() {
            let vertex = top.vertex;
            for i in top.idx..vertex.n_children as usize {
                // If chromosome index start is greater than the query, stop searching this node
                if vertex.chr_idx_start[i] > self.chrom_id {
                    continue;
                }

                // Check if this is the correct chromosome
                if self.chrom_id >= vertex.chr_idx_start[i] && self.chrom_id <= vertex.chr_idx_end[i] {
                    // Check region on the chromosome
                    if vertex.chr_idx_start[i] == vertex.chr_idx_end[i] {
                        // Check if query region is ahead or past the current region
                        if vertex.base_end[i] <= self.from {
                            continue;
                        }
                        if vertex.base_start[i] >= self.to {
                            break;
                        }
                    }

                    // Push current position incremented by one leaf
                    self.stack.push(RTreeTraverserType {
                        vertex,
                        idx: i + 1,
                    });

                    // If this is a non-leaf vertex, traverse its children
                    if vertex.is_leaf == 0 {
                        self.stack.push(RTreeTraverserType {
                            vertex: &vertex.children[i],
                            idx: 0,
                        });
                        break;
                    } else {
                        // Save result and exit
                        return Some(RTreeTraverserType {
                            vertex,
                            idx: i,
                        });
                    }
                }
            }
        }
    }
}

/* -------------------------------------------------------------------------- */

struct BbiQueryType {
    bbi_summary_record: BbiSummaryRecord,
    data_type         : u8,
    quit              : Box<dyn FnOnce()>,
    error             : Option<std::io::Error>,
}

/* -------------------------------------------------------------------------- */

impl BbiQueryType {
    pub fn new(quit: Box<dyn FnOnce()>) -> Self {
        BbiQueryType {
            bbi_summary_record: BbiSummaryRecord::new(),
            data_type: 0,
            quit,
            error: None,
        }
    }
}

/* -------------------------------------------------------------------------- */

struct BbiFile {
    header    : BbiHeader,
    chrom_data: BData,
    index     : RTree,
    index_zoom: Vec<RTree>,
}

/* -------------------------------------------------------------------------- */

impl BbiFile {
    pub fn new() -> Self {
        BbiFile {
            header    : BbiHeader::new(),
            chrom_data: BData::new(),
            index     : RTree::new(),
            index_zoom: vec![],
        }
    }

    pub fn estimate_size(&self, offset: i64, init: usize) -> usize {
        let mut n = i64::MAX;
        let mut k : i64;        

        k = self.header.ct_offset as i64;
        if offset < k as i64 && k < n {
            n = k;
        }

        k = self.header.data_offset as i64;
        if offset < k && k < n {
            n = k;
        }

        k = self.header.index_offset as i64;
        if offset < k && k < n {
            n = k;
        }

        for zoom_header in &self.header.zoom_headers {
            k = zoom_header.index_offset as i64;
            if offset < k && k < n {
                n = k;
            }

            k = zoom_header.data_offset as i64;
            if offset < k && k < n {
                n = k;
            }
        }

        if n == i64::MAX {
            init
        } else {
            (n - offset) as usize
        }
    }

    pub fn read_index<E: ByteOrder, R: Read + Seek>(&mut self, reader: &mut R) -> io::Result<()> {
        let n = self.estimate_size(self.header.index_offset as i64, 1024);
        reader.seek(SeekFrom::Start(self.header.index_offset))?;
        self.index.read::<E, R>(reader)
    }

    pub fn read_zoom_index<E: ByteOrder, R: Read + Seek>(&mut self, reader: &mut R, i: usize) -> io::Result<()> {
        let n = self.estimate_size(self.header.zoom_headers[i].index_offset as i64, 1024);
        reader.seek(SeekFrom::Start(self.header.zoom_headers[i].index_offset))?;
        self.index_zoom[i].read::<E, R>(reader)
    }

    pub fn query_zoom<E: ByteOrder, R: Read + Seek>(
        &mut self,
        reader  : &mut R,
        channel : &mut Vec<BbiQueryType>,
        done    : &mut bool,
        zoom_idx: usize,
        chrom_id: i32,
        from    : i32,
        to      : i32,
        bin_size: i32,
    ) -> bool {
        if self.index_zoom[zoom_idx].is_nil() {
            if let Err(err) = self.read_zoom_index(reader, zoom_idx) {
                channel.push(BbiQueryType {
                    error: Some(err),
                    ..BbiQueryType::new(Box::new(|| {}))
                });
                return false;
            }
        }

        let mut traverser = RTreeTraverser::new(&self.index_zoom[zoom_idx], chrom_id, from, to);
        let mut result = BbiQueryType::new(Box::new(|| *done = true));

        for r in traverser {

            match r.vertex.read_block::<E>(reader, self, r.idx) {
                Err(err) => {
                    channel.push(BbiQueryType {
                        error: Some(err),
                        ..BbiQueryType::new(Box::new(|| {}))
                    });
                },
                Ok(block) => {
                    let mut decoder = BbiZoomBlockDecoder::new(block);

                    for item in decoder.decode() {

                        match item.read::<E>() {
                            Err(err) => {
                                channel.push(BbiQueryType {
                                    error: Some(err),
                                    ..BbiQueryType::new(Box::new(|| {}))
                                });
                            },
                            Ok(record) => {
                                if record.chrom_id != chrom_id || record.from < from || record.to > to {
                                    continue;
                                }
        
                                if result.bbi_summary_record.chrom_id == -1 {
                                    result.bbi_summary_record.chrom_id = record.chrom_id;
                                    result.bbi_summary_record.from     = record.from;
                                    result.bbi_summary_record.to       = record.from;
                                    result.data_type                   = BBI_TYPE_BED_GRAPH;
                                }
        
                                if result.bbi_summary_record.to - result.bbi_summary_record.from >= bin_size
                                    || result.bbi_summary_record.from + bin_size < record.from
                                {
                                    if result.bbi_summary_record.from != result.bbi_summary_record.to {
                                        if *done {
                                            return false;
                                        } else {
                                            channel.push(result);
                                            result = BbiQueryType::new(Box::new(|| *done = true));
                                        }
                                    }
                                }
        
                                result.bbi_summary_record.add_record(&record);        
                            }
                        }
                    }
                }
            }
            traverser.next();
        }

        if result.bbi_summary_record.chrom_id != -1 {
            channel.push(result);
        }
        true
    }

    pub fn query_raw<E: ByteOrder, R: Read + Seek>(
        &mut self,
        reader  : &mut R,
        channel : &mut Vec<BbiQueryType>,
        done    : &mut bool,
        chrom_id: i32,
        from    : i32,
        to      : i32,
        bin_size: i32,
    ) -> bool {
        if self.index.is_nil() {
            if let Err(err) = self.read_index(reader) {
                channel.push(BbiQueryType {
                    error: Some(err),
                    ..BbiQueryType::new(Box::new(|| {}))
                });
                return false;
            }
        }

        let mut traverser = RTreeTraverser::new(&self.index, chrom_id, from, to);
        let mut result = BbiQueryType::new(Box::new(|| *done = true));

        for r in traverser {

            match r.vertex.read_block::<E>(reader, self, r.idx) {
                Err(err) => {
                    channel.push(BbiQueryType {
                        error: Some(err),
                        ..BbiQueryType::new(Box::new(|| {}))
                    });
                },
                Ok(block) => {
                    let mut decoder = BbiRawBlockDecoder::new(block).unwrap();

                    for item in decoder.decode() {

                        let record = item.read::<E>();

                        if record.chrom_id != chrom_id || record.from < from || record.to > to {
                            continue;
                        }

                        if result.bbi_summary_record.chrom_id == -1 {
                            result.bbi_summary_record.chrom_id = record.chrom_id;
                            result.bbi_summary_record.from     = record.from;
                            result.bbi_summary_record.to       = record.from;
                            result.data_type                   = decoder.header.kind;
                        }

                        if result.bbi_summary_record.to - result.bbi_summary_record.from >= bin_size
                            || result.bbi_summary_record.from + bin_size < record.from
                        {
                            if result.bbi_summary_record.from != result.bbi_summary_record.to {
                                if *done {
                                    return false;
                                } else {
                                    channel.push(result);
                                    result = BbiQueryType::new(Box::new(|| *done = true));
                                }
                            }
                        }

                        result.bbi_summary_record.add_record(&record);
                    }
                }
            }
            traverser.next();
        }

        if result.bbi_summary_record.chrom_id != -1 {
            channel.push(result);
        }
        true
    }

    pub fn query<E: ByteOrder, R: Read + Seek>(
        &mut self,
        reader  : &mut R,
        channel : &mut Vec<BbiQueryType>,
        chrom_id: i32,
        from    : i32,
        to      : i32,
        bin_size: i32,
    ) -> bool {
        let mut done = false;

        if bin_size != 0 {
            let from = (from / bin_size) * bin_size;
            let to = ((to + bin_size - 1) / bin_size) * bin_size;

            let mut zoom_idx = -1;
            for (i, zoom_header) in self.header.zoom_headers.iter().enumerate() {
                if bin_size >= zoom_header.reduction_level as i32
                    && bin_size % zoom_header.reduction_level as i32 == 0
                {
                    zoom_idx = i as i32;
                    break;
                }
            }

            if zoom_idx != -1 {
                return self.query_zoom::<E, R>(reader, channel, &mut done, zoom_idx as usize, chrom_id, from, to, bin_size);
            }
        }

        self.query_raw::<E, R>(reader, channel, &mut done, chrom_id, from, to, bin_size)
    }
}
