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

    /// Writes the content of the `GRanges` object in BED3 format (first 3 columns: chromosome, start, end) to the provided writer.
    ///
    /// # Arguments
    /// - `writer`: Mutable reference to an object implementing the `Write` trait, where the BED3 data will be written.
    ///
    /// # Errors
    /// Returns an error if writing to the `writer` fails.
    pub fn write_bed3<W: Write>(&self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        for i in 0..self.num_rows() {
            write!(writer, "{}\t{}\t{}\n", self.seqnames[i], self.ranges[i].from, self.ranges[i].to)?;
        }
        Ok(())
    }

    /// Writes the content of the `GRanges` object in BED6 format (first 6 columns: chromosome, start, end, name, score, strand) to the provided writer.
    /// Requires "name" and "score" metadata columns in `GRanges`.
    ///
    /// # Arguments
    /// - `writer`: Mutable reference to an object implementing the `Write` trait, where the BED6 data will be written.
    ///
    /// # Errors
    /// Returns an error if writing fails or if required metadata columns are missing.
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

    /// Writes the content of the `GRanges` object in BED9 format (first 9 columns) to the provided writer.
    /// Includes optional metadata columns for "name", "score", "thickStart", "thickEnd", and "itemRgb".
    ///
    /// # Arguments
    /// - `writer`: Mutable reference to an object implementing the `Write` trait, where the BED9 data will be written.
    ///
    /// # Errors
    /// Returns an error if writing fails or if required metadata columns are missing.
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

    /// Generalized BED writer that writes the `GRanges` object in either BED3, BED6, or BED9 format, depending on the number of columns specified.
    ///
    /// # Arguments
    /// - `writer`: Mutable reference to an object implementing the `Write` trait.
    /// - `columns`: Number of BED columns to output, must be 3, 6, or 9.
    ///
    /// # Errors
    /// Returns an error if writing fails or if the specified number of columns is invalid.
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

    /// Writes BED3 format to a specified file, with optional gzip compression.
    ///
    /// # Arguments
    /// - `filename`: Path to the output file.
    /// - `compress`: If true, the output is gzip-compressed.
    ///
    /// # Errors
    /// Returns an error if file creation, writing, or compression fails.
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

    /// Writes BED6 format to a specified file, with optional gzip compression.
    /// Requires "name" and "score" metadata columns in `GRanges`.
    ///
    /// # Arguments
    /// - `filename`: Path to the output file.
    /// - `compress`: If true, the output is gzip-compressed.
    ///
    /// # Errors
    /// Returns an error if file creation, writing, or compression fails.
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

    /// Writes BED9 format to a specified file, with optional gzip compression.
    /// Requires additional metadata columns ("name", "score", "thickStart", "thickEnd", "itemRgb").
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

    /// Exports the `GRanges` object to a BED file with the specified number of columns (3, 6, or 9),
    /// with optional gzip compression. This function is a flexible exporter, supporting different
    /// BED formats by delegating to `export_bed3`, `export_bed6`, or `export_bed9` as needed.
    ///
    /// # Arguments
    /// - `filename`: Path to the output file where the BED data will be saved.
    /// - `columns`: Number of BED columns to output, must be 3, 6, or 9.
    /// - `compress`: If true, the output file is gzip-compressed.
    ///
    /// # Errors
    /// Returns an error if file creation, writing, or compression fails,
    /// or if an invalid number of columns is specified.
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

    /// Reads a BED3 file into the `GRanges` object. Only processes first 3 columns; ignores metadata lines and optional fields.
    ///
    /// # Arguments
    /// - `reader`: Reference to an object implementing `Read` and `BufRead`.
    ///
    /// # Errors
    /// Returns an error if reading fails or if rows contain fewer than 3 columns.
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

    /// Reads a BED6 file into the `GRanges` object, including "name" and "score" columns as metadata.
    /// Additional columns and metadata lines (e.g., "track", "browser") are ignored.
    ///
    /// # Errors
    /// Returns an error if reading fails, rows contain fewer than 6 columns, or if metadata columns are missing.
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

    /// Reads a BED9 file into the `GRanges` object, supporting all 9 columns and relevant metadata fields.
    /// Skips "track" and "browser" metadata lines.
    ///
    /// # Errors
    /// Returns an error if reading fails, rows contain fewer than 9 columns, or if required metadata columns are missing.
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

    /// Imports a BED3 file from disk, optionally gzip-compressed.
    ///
    /// # Arguments
    /// - `filename`: Path to the input BED file.
    /// - `compress`: If true, the input file is treated as gzip-compressed.
    ///
    /// # Errors
    /// Returns an error if file reading or decompression fails.
    pub fn read_bed3<R: Read>(&mut self, reader: &mut R) -> Result<(), Box<dyn Error>> {
        self.bufread_bed3(&mut BufReader::new(reader))
    }

    /// Imports a BED6 file from disk, optionally gzip-compressed. Loads "name" and "score" metadata columns.
    ///
    /// # Errors
    /// Returns an error if file reading, decompression fails, or if metadata columns are missing.
    pub fn read_bed6<R: Read>(&mut self, reader: &mut R) -> Result<(), Box<dyn Error>> {
        self.bufread_bed6(&mut BufReader::new(reader))
    }

    /// Imports a BED9 file from disk, optionally gzip-compressed, and loads all supported metadata columns.
    pub fn read_bed9<R: Read>(&mut self, reader: &mut R) -> Result<(), Box<dyn Error>> {
        self.bufread_bed9(&mut BufReader::new(reader))
    }

    /// Reads a BED file with specified columns (3, 6, or 9) into the `GRanges` object.
    ///
    /// # Arguments
    /// - `reader`: Reference to an object implementing `Read`.
    /// - `columns`: Number of columns expected in the file, must be 3, 6, or 9.
    ///
    /// # Errors
    /// Returns an error if reading fails or if `columns` is an unsupported value.
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

    /// Imports a BED3 file (with 3 mandatory columns) into the `GRanges` object,
    /// with optional gzip compression.
    ///
    /// # Arguments
    /// - `filename`: Path to the input BED file to read.
    /// - `compress`: If true, reads the file as gzip-compressed.
    ///
    /// # Errors
    /// Returns an error if file reading or decompression fails, or if the file format is invalid.
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

    /// Imports a BED6 file (with 6 mandatory columns) into the `GRanges` object,
    /// with optional gzip compression. The function extracts additional metadata for `name`, `score`,
    /// and `strand` columns from the file.
    ///
    /// # Arguments
    /// - `filename`: Path to the input BED file to read.
    /// - `compress`: If true, reads the file as gzip-compressed.
    ///
    /// # Errors
    /// Returns an error if file reading or decompression fails, or if the file format is invalid.
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

    /// Imports a BED9 file (with 9 mandatory columns) into the `GRanges` object,
    /// with optional gzip compression. The function extracts additional metadata for `name`,
    /// `score`, `strand`, `thickStart`, `thickEnd`, and `itemRgb` columns from the file.
    ///
    /// # Arguments
    /// - `filename`: Path to the input BED file to read.
    /// - `compress`: If true, reads the file as gzip-compressed.
    ///
    /// # Errors
    /// Returns an error if file reading or decompression fails, or if the file format is invalid.
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

    /// Imports a BED file with a specified number of columns (3, 6, or 9) into the `GRanges` object,
    /// with optional gzip compression. This method dispatches to `import_bed3`, `import_bed6`, or
    /// `import_bed9` based on the `columns` argument.
    ///
    /// # Arguments
    /// - `filename`: Path to the input BED file to read.
    /// - `columns`: Number of BED columns to read; must be 3, 6, or 9.
    /// - `compress`: If true, reads the file as gzip-compressed.
    ///
    /// # Errors
    /// Returns an error if file reading, decompression, or column specification fails, or if the file format is invalid.
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
