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

use std::io::{BufRead, BufReader, Write};
use std::fs::File;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use list_comprehension_macro::comp;

use crate::range::Range;
use crate::granges::GRanges;
use crate::error::Error;
use crate::meta::MetaData;

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn write_bed3(&self, writer: &mut dyn Write) -> Result<(), Error> {
        for i in 0..self.num_rows() {
            write!(writer, "{}\t{}\t{}\n", self.seqnames[i], self.ranges[i].from, self.ranges[i].to)?;
        }
        Ok(())
    }

    pub fn export_bed3(&self, filename: &str, compress: bool) -> Result<(), Error> {
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

    pub fn write_bed6(&self, writer: &mut dyn Write) -> Result<(), Error> {
        let name  = match self.meta.get_column_str(&String::from("name" )) {
            Some(v) => v,
            None    => return Err(Error::Generic("Meta data does not contain a column `name'".to_string()))
        };
        let score = match self.meta.get_column_int(&String::from("score")) {
            Some(v) => v,
            None    => return Err(Error::Generic("Meta data does not contain a column `score'".to_string()))
        };
        for i in 0..self.num_rows() {
            write!(writer, "{}\t{}\t{}\t{}\t{}\t{}\n", self.seqnames[i], self.ranges[i].from, self.ranges[i].to, name[i], score[i], self.strand.get(i).unwrap_or(&'.'))?;
        }
        Ok(())
    }

    pub fn export_bed6(&self, filename: &str, compress: bool) -> Result<(), Error> {
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

    pub fn write_bed9(&self, writer: &mut dyn Write) -> Result<(), Error> {
        let name  = self.meta.get_column_str("name").ok_or(
            Error::from("Meta data does not contain a column `name'")
            )?;
        let score = self.meta.get_column_int("score").ok_or(
            Error::from("Meta data does not contain a column `score'")
            )?;
        let item_rgb = match self.meta.get_column_str("itemRgb") {
            Some(v) => v.clone(),
            None    => vec![String::from("0,0,0"); self.num_rows()]
        };
        let thick_start = match self.meta.get_column_int("thickStart") {
            Some(v) => v.clone(),
            None    => comp![ self.ranges[i].from as i64 for i in 0..self.num_rows() ]
        };
        let thick_end   = match self.meta.get_column_int("thickStart") {
            Some(v) => v.clone(),
            None    => comp![ self.ranges[i].to   as i64 for i in 0..self.num_rows() ]
        };

        for i in 0..self.num_rows() {
            write!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", self.seqnames[i], self.ranges[i].from, self.ranges[i].to, name[i], score[i], self.strand.get(i).unwrap_or(&'.'), thick_start[i], thick_end[i], item_rgb[i])?;
        }
        Ok(())
    }

    pub fn export_bed9(&self, filename: &str, compress: bool) -> Result<(), Error> {
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

    pub fn read_bed3(&mut self, reader: &mut dyn BufRead) -> Result<(), Error> {
        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            let fields: Vec<&str> = line.trim().split('\t').collect();
            if fields.len() < 3 {
                return Err(Error::Generic("Bed file must have at least 3 columns".to_string()));
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

    pub fn import_bed3(&mut self, filename: &str, compress: bool) -> Result<(), Error> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        self.read_bed3(&mut reader)?;
        Ok(())
    }

    pub fn read_bed6(&mut self, reader: &mut dyn BufRead) -> Result<(), Error> {
        let mut line  = String::new();
        let mut name  = Vec::new();
        let mut score = Vec::new();
        while reader.read_line(&mut line)? > 0 {
            let fields: Vec<&str> = line.trim().split('\t').collect();
            if fields.len() < 6 {
                return Err(Error::Generic("Bed file must have at least 6 columns".to_string()));
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
        self.meta.add_meta("name" , MetaData::StringArray(name))?;
        self.meta.add_meta("score", MetaData::IntArray   (score))?;
        Ok(())
    }

    pub fn import_bed6(&mut self, filename: &str, compress: bool) -> Result<(), Error> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        self.read_bed6(&mut reader)?;
        Ok(())
    }

    pub fn read_bed9(&mut self, reader: &mut dyn BufRead) -> Result<(), Error> {
        let mut line  = String::new();
        let mut name  = Vec::new();
        let mut score = Vec::new();
        let mut thick_start = Vec::new();
        let mut thick_end   = Vec::new();
        let mut item_rgb    = Vec::new();
        while reader.read_line(&mut line)? > 0 {
            let fields: Vec<&str> = line.trim().split('\t').collect();
            if fields.len() < 9 {
                return Err(Error::Generic("Bed file must have at least 9 columns".to_string()));
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
        self.meta.add_meta("name"      , MetaData::StringArray(name ))?;
        self.meta.add_meta("score"     , MetaData::IntArray   (score))?;
        self.meta.add_meta("thickStart", MetaData::IntArray   (thick_start))?;
        self.meta.add_meta("thickEnd"  , MetaData::IntArray   (thick_end  ))?;
        self.meta.add_meta("item_rgb"  , MetaData::StringArray(item_rgb   ))?;
        Ok(())
    }

    pub fn import_bed9(&mut self, filename: &str, compress: bool) -> Result<(), Error> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        self.read_bed9(&mut reader)?;
        Ok(())
    }

    pub fn read_bed(&mut self, reader: &mut dyn BufRead, columns: usize) -> Result<(), Error> {
        match columns {
            3 => self.read_bed3(reader),
            6 => self.read_bed6(reader),
            9 => self.read_bed9(reader),
            _ => Err(Error::Generic("Invalid number of columns".to_string())),
        }
    }

    pub fn import_bed(&mut self, filename: &str, columns: usize, compress: bool) -> Result<(), Error> {
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
