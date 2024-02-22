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

use std::io::{self, BufRead, BufReader, Write};
use std::fs::File;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use list_comprehension_macro::comp;
use std::convert::TryInto;

use crate::range::Range;
use crate::granges::{GRanges};

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn write_bed3(&self, writer: &mut dyn Write) -> io::Result<()> {
        for i in 0..self.length() {
            write!(writer, "{}\t{}\t{}\n", self.seqnames[i], self.ranges[i].from, self.ranges[i].to)?;
        }
        Ok(())
    }

    pub fn export_bed3(&self, filename: &str, compress: bool) -> io::Result<()> {
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

    pub fn write_bed6(&self, writer: &mut dyn Write) -> io::Result<()> {
        let name  = self.meta.get_column_str(&String::from("name" )).unwrap();
        let score = self.meta.get_column_int(&String::from("score")).unwrap();

        for i in 0..self.length() {
            let r = write!(writer, "{}\t{}\t{}\t{}\t{}\t{}\n", self.seqnames[i], self.ranges[i].from, self.ranges[i].to, name[i], score[i], self.strand.get(i).unwrap_or(&'.'));

            if r.is_err() {
                return r
            }
        }
        Ok(())
    }

    pub fn export_bed6(&self, filename: &str, compress: bool) -> io::Result<()> {
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

    pub fn write_bed9(&self, writer: &mut dyn Write) -> io::Result<()> {
        let name  = self.meta.get_column_str(&String::from("name" )).unwrap();
        let score = self.meta.get_column_int(&String::from("score")).unwrap();
        let item_rgb = match self.meta.get_column_str(&String::from("itemRgb")) {
            Some(v) => v,
            None    => &vec![String::from("0,0,0"); self.length()]
        };
        let thick_start = self.meta.get_column_int(&String::from("thickStart")).unwrap_or(&comp![ self.ranges[i].from.try_into().unwrap() for i in 0..self.length() ]);
        let thick_end   = self.meta.get_column_int(&String::from("thickEnd"  )).unwrap_or(&comp![ self.ranges[i].to  .try_into().unwrap() for i in 0..self.length() ]);

        for i in 0..self.length() {

            let r = write!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", self.seqnames[i], self.ranges[i].from, self.ranges[i].to, name[i], score[i], self.strand.get(i).unwrap_or(&'.'), thick_start[i], thick_end[i], item_rgb[i]);
            
            if r.is_err() {
                return r
            }
        }
        Ok(())
    }

    pub fn export_bed9(&self, filename: &str, compress: bool) -> io::Result<()> {
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

    pub fn read_bed3(&mut self, reader: &mut dyn BufRead) -> io::Result<()> {
        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            let fields: Vec<&str> = line.trim().split('\t').collect();
            if fields.len() < 3 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "Bed file must have at least 3 columns"));
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

    pub fn import_bed3(&mut self, filename: &str, compress: bool) -> io::Result<()> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        self.read_bed3(&mut reader)?;
        Ok(())
    }

    pub fn read_bed6(&mut self, reader: &mut dyn BufRead) -> io::Result<()> {
        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            let fields: Vec<&str> = line.trim().split('\t').collect();
            if fields.len() < 6 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "Bed file must have at least 6 columns"));
            }
            let from   = fields[1].parse::<usize>().unwrap();
            let to     = fields[2].parse::<usize>().unwrap();
            let name   = fields[3].to_string();
            let score  = fields[4].parse::<usize>().unwrap();
            let strand = fields[5].chars().next().unwrap_or('.');
            self.seqnames.push(fields[0].to_string());
            self.ranges.push(Range::new(from, to));
            self.strand.push(strand);
            //self.add_meta("name" .to_string(), vec![name]);
            //self.add_meta("score".to_string(), vec![score.to_string()]);
            line.clear();
        }
        Ok(())
    }

    pub fn import_bed6(&mut self, filename: &str, compress: bool) -> io::Result<()> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        self.read_bed6(&mut reader)?;
        Ok(())
    }

    pub fn read_bed9(&mut self, reader: &mut dyn BufRead) -> io::Result<()> {
        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            let fields: Vec<&str> = line.trim().split('\t').collect();
            if fields.len() < 9 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "Bed file must have at least 9 columns"));
            }
            let from        = fields[1].parse::<usize>().unwrap();
            let to          = fields[2].parse::<usize>().unwrap();
            let name        = fields[3].to_string();
            let score       = fields[4].parse::<usize>().unwrap();
            let strand      = fields[5].chars().next().unwrap_or('.');
            let thick_start = fields[6].parse::<usize>().unwrap();
            let thick_end   = fields[7].parse::<usize>().unwrap();
            let item_rgb    = fields[8].to_string();
            self.seqnames.push(fields[0].to_string());
            self.ranges  .push(Range::new(from, to));
            self.strand  .push(strand);
            //self.add_meta("name".to_string(), vec![name]);
            //self.add_meta("score".to_string(), vec![score.to_string()]);
            //self.add_meta("thickStart".to_string(), vec![thick_start.to_string()]);
            //self.add_meta("thickEnd".to_string(), vec![thick_end.to_string()]);
            //self.add_meta("itemRgb".to_string(), vec![item_rgb]);
            line.clear();
        }
        Ok(())
    }

    pub fn import_bed9(&mut self, filename: &str, compress: bool) -> io::Result<()> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        self.read_bed9(&mut reader)?;
        Ok(())
    }

    pub fn read_bed(&mut self, reader: &mut dyn BufRead, columns: usize) -> io::Result<()> {
        match columns {
            3 => self.read_bed3(reader),
            6 => self.read_bed6(reader),
            9 => self.read_bed9(reader),
            _ => Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid number of columns")),
        }
    }

    pub fn import_bed(&mut self, filename: &str, columns: usize, compress: bool) -> io::Result<()> {
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
