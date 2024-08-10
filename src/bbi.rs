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

/* -------------------------------------------------------------------------- */

const CIRTREE_MAGIC      : u32   = 0x78ca8c91;
const IDX_MAGIC          : u32   = 0x2468ace0;
const BBI_MAX_ZOOM_LEVELS: usize = 10;
const BBI_RES_INCREMENT  : u32   = 4;
const BBI_TYPE_FIXED     : u8    = 3;
const BBI_TYPE_VARIABLE  : u8    = 2;
const BBI_TYPE_BED_GRAPH : u8    = 1;

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
    fn new<E: ByteOrder>(buffer: &'a [u8]) -> Result<BbiRawBlockDecoder<'a>, std::io::Error> {
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

    fn read_buffer<E: ByteOrder>(&mut self, buffer : &[u8]) -> Result<(), io::Error> {
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
 
    fn read<E: ByteOrder>(&mut self) -> Result<BbiZoomBlockDecoderType, io::Error> { 

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
        let tree = BTree::new(self); // Assuming a `BTree` struct similar to `BData`
        tree.write::<E, W>(file)
    }
}
