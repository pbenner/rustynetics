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

use std::io::{Read, Seek, Write};
use std::io::Cursor;
use std::io::SeekFrom;

use std::f32;
use std::f64;

use byteorder::{ByteOrder, ReadBytesExt, WriteBytesExt, BigEndian, LittleEndian};

use flate2::Compression;
use flate2::write::ZlibEncoder;
use flate2::read::ZlibDecoder;

/* -------------------------------------------------------------------------- */

const CIRTREE_MAGIC   : u32   = 0x78ca8c91;
const IDX_MAGIC       : u32   = 0x2468ace0;
const BbiMaxZoomLevels: usize = 10;
const BbiResIncrement : u32   = 4;
const BbiTypeFixed    : u8    = 3;
const BbiTypeVariable : u8    = 2;
const BbiTypeBedGraph : u8    = 1;

/* -------------------------------------------------------------------------- */

fn file_read_at<T: Read + Seek>(file: &mut T, offset: u64) -> Result<Vec<u8>, std::io::Error> {
    let current_position = file.seek(SeekFrom::Current(0))?;
    file.seek(SeekFrom::Start(offset))?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    file.seek(SeekFrom::Start(current_position))?;
    Ok(buffer)
}

fn file_write_at<T: Write + Seek>(file: &mut T, offset: u64, data: &[u8]) -> Result<(), std::io::Error> {
    let current_position = file.seek(SeekFrom::Current(0))?;
    file.seek(SeekFrom::Start(offset))?;
    file.write_all(data)?;
    file.seek(SeekFrom::Start(current_position))?;
    Ok(())
}

fn uncompress_slice(data: &[u8]) -> Result<Vec<u8>, std::io::Error> {
    let mut decoder = ZlibDecoder::new(data);
    let mut buffer = Vec::new();
    decoder.read_to_end(&mut buffer)?;
    Ok(buffer)
}

fn compress_slice(data: &[u8]) -> Result<Vec<u8>, std::io::Error> {
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::best());
    encoder.write_all(data)?;
    let compressed_data = encoder.finish()?;
    Ok(compressed_data)
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
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

    fn read<E: ByteOrder, T: Read + Seek>(&mut self, reader: &mut T) -> Result<(), std::io::Error> {
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

    fn write<E: ByteOrder, T: Write + Seek>(&self, writer: &mut T) -> Result<(), std::io::Error> {
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
    fn reset(&mut self) {
        self.valid = 0.0;
        self.min = f64::INFINITY;
        self.max = f64::NEG_INFINITY;
        self.sum = 0.0;
        self.sum_squares = 0.0;
    }

    fn add_value(&mut self, x: f64) {
        if x.is_nan() {
            return;
        }
        self.valid += 1.0;
        self.min = self.min.min(x);
        self.max = self.max.max(x);
        self.sum += x;
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
            chrom_id: -1,
            from: 0,
            to: 0,
            statistics: BbiSummaryStatistics {
                valid: 0.0,
                min: f64::INFINITY,
                max: f64::NEG_INFINITY,
                sum: 0.0,
                sum_squares: 0.0,
            },
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

#[derive(Debug)]
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
    fn read_buffer<E: ByteOrder>(&mut self, buffer: &[u8]) {
        let mut cursor = Cursor::new(buffer);

        self.chrom_id   = cursor.read_u32::<E>().unwrap();
        self.start      = cursor.read_u32::<E>().unwrap();
        self.end        = cursor.read_u32::<E>().unwrap();
        self.step       = cursor.read_u32::<E>().unwrap();
        self.span       = cursor.read_u32::<E>().unwrap();
        self.kind       = cursor.read_u8().unwrap();
        self.reserved   = cursor.read_u8().unwrap();
        self.item_count = cursor.read_u16::<E>().unwrap();
    }

    fn write_buffer<E: ByteOrder>(&self, buffer: &mut [u8]) {
        let mut cursor = Cursor::new(buffer);

        cursor.write_u32::<E>(self.chrom_id).unwrap();
        cursor.write_u32::<E>(self.start).unwrap();
        cursor.write_u32::<E>(self.end).unwrap();
        cursor.write_u32::<E>(self.step).unwrap();
        cursor.write_u32::<E>(self.span).unwrap();
        cursor.write_u8(self.kind).unwrap();
        cursor.write_u8(self.reserved).unwrap();
        cursor.write_u16::<E>(self.item_count).unwrap();
    }
}

/* -------------------------------------------------------------------------- */

trait BbiBlockDecoderIterator {
    fn get (&self) -> &BbiBlockDecoderType;
    fn ok  (&self) -> bool;
    fn next(&mut self);
}

/* -------------------------------------------------------------------------- */

struct BbiRawBlockDecoderIterator<'a> {
    decoder: &'a BbiRawBlockDecoder<'a>,
    i: usize,
    record: BbiBlockDecoderType,
}

/* -------------------------------------------------------------------------- */

trait BbiBlockDecoder {
    fn decode(&mut self) -> Option<Box<BbiBlockDecoderIterator>>;
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
struct BbiBlockDecoderType {
    record: BbiSummaryRecord,
}

/* -------------------------------------------------------------------------- */

struct BbiRawBlockDecoder<'a> {
    header: BbiDataHeader,
    buffer: &'a [u8],
    order: byteorder::LittleEndian,
}

/* -------------------------------------------------------------------------- */

impl<'a> BbiRawBlockDecoder<'a> {
    fn new(buffer: &'a [u8], order: byteorder::LittleEndian) -> Result<BbiRawBlockDecoder<'a>, std::io::Error> {
        if buffer.len() < 24 {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "block length is shorter than 24 bytes"));
        }
        let mut decoder = BbiRawBlockDecoder {
            header: BbiDataHeader {
                chrom_id: 0,
                start: 0,
                end: 0,
                step: 0,
                span: 0,
                kind: 0,
                reserved: 0,
                item_count: 0,
            },
            buffer,
            order,
        };
        decoder.header.read_buffer(buffer, order);
        match decoder.header.kind {
            BbiTypeBedGraph => {
                if buffer.len() % 12 != 0 {
                    return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "bedGraph data block has invalid length"));
                }
            }
            BbiTypeVariable => {
                if buffer.len() % 8 != 0 {
                    return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "variable step data block has invalid length"));
                }
            }
            BbiTypeFixed => {
                if buffer.len() % 4 != 0 {
                    return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "fixed step data block has invalid length"));
                }
            }
            _ => return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "unsupported block type")),
        }
        Ok(decoder)
    }

    fn read_fixed(&self, i: usize) -> BbiBlockDecoderType {
        let mut record = BbiBlockDecoderType {
            record: BbiSummaryRecord::new(),
        };
        record.record.chrom_id               = self.header.chrom_id as i32;
        record.record.from                   = (self.header.start + (i as u32 / 4) * self.header.step) as i32;
        record.record.to                     = record.record.from + self.header.span as i32;
        record.record.statistics.valid       = 1.0;
        record.record.statistics.sum         = f32::from_bits(self.order.read_u32(&self.buffer[i..i + 4]).unwrap()) as f64;
        record.record.statistics.sum_squares = record.record.statistics.sum * record.record.statistics.sum;
        record.record.statistics.min         = record.record.statistics.sum;
        record.record.statistics.max         = record.record.statistics.sum;
        record
    }

    fn read_variable(&self, i: usize) -> BbiBlockDecoderType {
        let mut record = BbiBlockDecoderType {
            record: BbiSummaryRecord::new(),
        };
        record.record.chrom_id               = self.header.chrom_id as i32;
        record.record.from                   = self.order.read_u32(&self.buffer[i..i + 4]).unwrap() as i32;
        record.record.to                     = record.record.from + self.header.span as i32;
        record.record.statistics.valid       = 1.0;
        record.record.statistics.sum         = f32::from_bits(self.order.read_u32(&self.buffer[i + 4..i + 8]).unwrap()) as f64;
        record.record.statistics.sum_squares = record.record.statistics.sum * record.record.statistics.sum;
        record.record.statistics.min         = record.record.statistics.sum;
        record.record.statistics.max         = record.record.statistics.sum;
        record
    }

    fn read_bed_graph(&self, i: usize) -> BbiBlockDecoderType {
        let mut record = BbiBlockDecoderType {
            record: BbiSummaryRecord::new(),
        };
        record.record.chrom_id = self.header.chrom_id as i32;
        record.record.from = self.order.read_u32(&self.buffer[i..i + 4]).unwrap() as i32;
        record.record.to = self.order.read_u32(&self.buffer[i + 4..i + 8]).unwrap() as i32;
        record.record.statistics.valid = 1.0;
        record.record.statistics.sum = f32::from_bits(self.order.read_u32(&self.buffer[i + 8..i + 12]).unwrap()) as f64;
        record.record.statistics.sum_squares = record.record.statistics.sum * record.record.statistics.sum;
        record.record.statistics.min = record.record.statistics.sum;
        record.record.statistics.max = record.record.statistics.sum;
        record
    }
}

impl<'a> BbiBlockDecoder for BbiRawBlockDecoder<'a> {
    fn decode(&mut self) -> Option<Box<dyn BbiBlockDecoderIterator + 'a>> {
        if self.header.item_count == 0 {
            return None;
        }
        Some(Box::new(BbiRawBlockDecoderIterator {
            decoder: self,
            i: 0,
            record: BbiBlockDecoderType {
                record: BbiSummaryRecord::new(),
            },
        }))
    }
}

impl<'a> BbiBlockDecoderIterator for BbiRawBlockDecoderIterator<'a> {
    fn get(&self) -> &BbiBlockDecoderType {
        &self.record
    }

    fn ok(&self) -> bool {
        self.i < self.decoder.buffer.len()
    }

    fn next(&mut self) {
        if self.i >= self.decoder.buffer.len() {
            return;
        }
        match self.decoder.header.kind {
            BbiTypeBedGraph => {
                self.record = self.decoder.read_bed_graph(self.i);
                self.i += 12;
            }
            BbiTypeVariable => {
                self.record = self.decoder.read_variable(self.i);
                self.i += 8;
            }
            BbiTypeFixed => {
                self.record = self.decoder.read_fixed(self.i);
                self.i += 4;
            }
            _ => panic!("Unsupported block type"),
        }
    }
}
