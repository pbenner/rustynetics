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

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::error::Error;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use list_comprehension_macro::comp;

use crate::range::Range;
use crate::granges::GRanges;
use crate::granges_error::MissingColumn;
use crate::meta::MetaData;
use crate::error::ArgumentError;

/* Write GRanges to BED files
 * -------------------------------------------------------------------------- */

impl GRanges {

    pub fn write_bed3<W: Write>(&self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        for i in 0..self.num_rows() {
            write!(writer, "{}\t{}\t{}\n", self.seqnames[i], self.ranges[i].from, self.ranges[i].to)?;
        }
        Ok(())
    }

    pub fn write_bed6<W: Write>(&self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        let name  = self.meta.get_column_str("name").ok_or(
            Box::new(MissingColumn("name".to_string()))
        )?;
        let score = self.meta.get_column_int("score").ok_or(
            Box::new(MissingColumn("score".to_string()))
        )?;
        for i in 0..self.num_rows() {
            write!(writer, "{}\t{}\t{}\t{}\t{}\t{}\n", self.seqnames[i], self.ranges[i].from, self.ranges[i].to, name[i], score[i], self.strand.get(i).unwrap_or(&'.'))?;
        }
        Ok(())
    }

    pub fn write_bed9<W: Write>(&self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        let name  = self.meta.get_column_str("name").ok_or(
            Box::new(MissingColumn("name".to_string()))
        )?;
        let score = self.meta.get_column_int("score").ok_or(
            Box::new(MissingColumn("score".to_string()))
        )?;
        let item_rgb = match self.meta.get_column_str("itemRgb") {
            Some(v) => v.clone(),
            None    => vec![String::from("0,0,0"); self.num_rows()]
        };
        let thick_start = match self.meta.get_column_int("thickStart") {
            Some(v) => v.clone(),
            None    => comp![ self.ranges[i].from as i64 for i in 0..self.num_rows() ]
        };
        let thick_end   = match self.meta.get_column_int("thickEnd") {
            Some(v) => v.clone(),
            None    => comp![ self.ranges[i].to   as i64 for i in 0..self.num_rows() ]
        };

        for i in 0..self.num_rows() {
            write!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", self.seqnames[i], self.ranges[i].from, self.ranges[i].to, name[i], score[i], self.strand.get(i).unwrap_or(&'.'), thick_start[i], thick_end[i], item_rgb[i])?;
        }
        Ok(())
    }

    pub fn write_bed<W: Write>(&self, writer: &mut W, columns: usize) -> Result<(), Box<dyn Error>> {
        match columns {
            3 => self.write_bed3(writer),
            6 => self.write_bed6(writer),
            9 => self.write_bed9(writer),
            _ => Err(Box::new(ArgumentError("Invalid number of columns".to_string()))),
        }
    }

}

/* Export GRanges to BED files
 * -------------------------------------------------------------------------- */

impl GRanges {

    pub fn export_bed3(&self, filename: &str, compress: bool) -> Result<(), Box<dyn Error>> {
        let file = File::create(filename)?;
        let mut writer: Box<dyn Write> = if compress {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };
        self.write_bed3(&mut writer)?;
        writer.flush()?;
        Ok(())
    }

    pub fn export_bed6(&self, filename: &str, compress: bool) -> Result<(), Box<dyn Error>> {
        let file = File::create(filename)?;
        let mut writer: Box<dyn Write> = if compress {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };
        self.write_bed6(&mut writer)?;
        writer.flush()?;
        Ok(())
    }

    pub fn export_bed9(&self, filename: &str, compress: bool) -> Result<(), Box<dyn Error>> {
        let file = File::create(filename)?;
        let mut writer: Box<dyn Write> = if compress {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };
        self.write_bed9(&mut writer)?;
        writer.flush()?;
        Ok(())
    }

    pub fn export_bed(&mut self, filename: &str, columns: usize, compress: bool) -> Result<(), Box<dyn Error>> {
        let file = File::create(filename)?;
        let mut writer: Box<dyn Write> = if compress {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };
        self.write_bed(&mut writer, columns)?;
        Ok(())
    }

}

/* Bufead GRanges from BED files
 * -------------------------------------------------------------------------- */

impl GRanges {

    pub fn bufread_bed3<R: Read + BufRead>(&mut self, reader: &mut R) -> Result<(), Box<dyn Error>> {
        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 3 {
                return Err(Box::new(io::Error::new(io::ErrorKind::InvalidInput, "Bed file must have at least 3 columns".to_string())));
            }
            let from = fields[1].parse::<usize>().unwrap();
            let to   = fields[2].parse::<usize>().unwrap();
            self.seqnames.push(fields[0].to_string());
            self.ranges  .push(Range::new(from, to));
            self.strand  .push('*');
            line.clear();
        }
        Ok(())
    }


    pub fn bufread_bed6<R: Read + BufRead>(&mut self, reader: &mut R) -> Result<(), Box<dyn Error>> {
        let mut line  = String::new();
        let mut name  = Vec::new();
        let mut score = Vec::new();
        while reader.read_line(&mut line)? > 0 {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 6 {
                return Err(Box::new(io::Error::new(io::ErrorKind::InvalidInput, "Bed file must have at least 6 columns".to_string())));
            }
            let from   = fields[1].parse::<usize>().unwrap();
            let to     = fields[2].parse::<usize>().unwrap();
            let strand = fields[5].chars().next().unwrap_or('.');
            self.seqnames.push(fields[0].to_string());
            self.ranges  .push(Range::new(from, to));
            self.strand  .push(strand);
            name .push(fields[3].to_string());
            score.push(fields[4].parse::<i64>().unwrap());
            line.clear();
        }
        self.meta.add("name" , MetaData::StringArray(name))?;
        self.meta.add("score", MetaData::IntArray   (score))?;
        Ok(())
    }

    pub fn bufread_bed9<R: Read + BufRead>(&mut self, reader_: &mut R) -> Result<(), Box<dyn Error>> {
        let mut reader = BufReader::new(reader_);
        let mut line   = String::new();
        let mut name   = Vec::new();
        let mut score  = Vec::new();
        let mut thick_start = Vec::new();
        let mut thick_end   = Vec::new();
        let mut item_rgb    = Vec::new();
        while reader.read_line(&mut line)? > 0 {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() > 0 && (fields[0] == "track" || fields[0] == "browser") {
                line.clear();
                continue;
            }
            if fields.len() < 9 {
                return Err(Box::new(io::Error::new(io::ErrorKind::InvalidInput, "Bed file must have at least 9 columns".to_string())));
            }
            let from        = fields[1].parse::<usize>().unwrap();
            let to          = fields[2].parse::<usize>().unwrap();
            let strand      = fields[5].chars().next().unwrap_or('.');
            self.seqnames.push(fields[0].to_string());
            self.ranges  .push(Range::new(from, to));
            self.strand  .push(strand);
            name .push(fields[3].to_string());
            score.push(fields[4].parse::<i64>().unwrap());
            thick_start.push(fields[6].parse::<i64>().unwrap());
            thick_end  .push(fields[7].parse::<i64>().unwrap());
            item_rgb   .push(fields[8].to_string());
            line.clear();
        }
        self.meta.add("name"      , MetaData::StringArray(name ))?;
        self.meta.add("score"     , MetaData::IntArray   (score))?;
        self.meta.add("thickStart", MetaData::IntArray   (thick_start))?;
        self.meta.add("thickEnd"  , MetaData::IntArray   (thick_end  ))?;
        self.meta.add("itemRgb"   , MetaData::StringArray(item_rgb   ))?;
        Ok(())
    }

}

/* Read GRanges from BED files
 * -------------------------------------------------------------------------- */

 impl GRanges {
    pub fn read_bed3<R: Read>(&mut self, reader: &mut R) -> Result<(), Box<dyn Error>> {
        self.bufread_bed3(&mut BufReader::new(reader))
    }

    pub fn read_bed6<R: Read>(&mut self, reader: &mut R) -> Result<(), Box<dyn Error>> {
        self.bufread_bed6(&mut BufReader::new(reader))
    }

    pub fn read_bed9<R: Read>(&mut self, reader: &mut R) -> Result<(), Box<dyn Error>> {
        self.bufread_bed9(&mut BufReader::new(reader))
    }

    pub fn read_bed<R: Read>(&mut self, reader: &mut R, columns: usize) -> Result<(), Box<dyn Error>> {
        match columns {
            3 => self.read_bed3(reader),
            6 => self.read_bed6(reader),
            9 => self.read_bed9(reader),
            _ => Err(Box::new(ArgumentError("Invalid number of columns".to_string()))),
        }
    }

}

/* Import GRanges from BED files
 * -------------------------------------------------------------------------- */

impl GRanges {

    pub fn import_bed3(&mut self, filename: &str, compress: bool) -> Result<(), Box<dyn Error>> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        self.bufread_bed3(&mut reader)?;
        Ok(())
    }

    pub fn import_bed6(&mut self, filename: &str, compress: bool) -> Result<(), Box<dyn Error>> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        self.read_bed6(&mut reader)?;
        Ok(())
    }

    pub fn import_bed9(&mut self, filename: &str, compress: bool) -> Result<(), Box<dyn Error>> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        self.read_bed9(&mut reader)?;
        Ok(())
    }

    pub fn import_bed(&mut self, filename: &str, columns: usize, compress: bool) -> Result<(), Box<dyn Error>> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        self.read_bed(&mut reader, columns)?;
        Ok(())
    }

}
