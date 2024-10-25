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

use std::any::Any;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::fs::File;

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;

use crate::granges::GRanges;
use crate::granges_table_reader::GRangesTableReader;
use crate::granges_table_writer::GRangesTableWriter;
use crate::meta_table_reader::MetaTableReader;

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct OptionPrintScientific(pub bool);

#[derive(Debug)]
pub struct OptionPrintStrand(pub bool);

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn bufread_table<R: BufRead>(&mut self, mut buf_reader: R, names: &[&str], types: &[&str]) -> io::Result<()> {
        let mut mreader = MetaTableReader   ::new(names, types);
        let mut greader = GRangesTableReader::new();

        let mut line = String::new();

        // Read first line as header
        buf_reader.read_line(&mut line)?;
        greader.read_header(&line)?;
        mreader.read_header(&line)?;

        line.clear();

        let mut line_counter = 0;

        while buf_reader.read_line(&mut line)? > 0 {

            if line.is_empty() {
                continue;
            }
            greader.read_line(&line, line_counter)?;
            mreader.read_line(&line, line_counter)?;

            line_counter += 1;

            line.clear();
        }
        greader.push(self);
        mreader.push(&mut self.meta);

        Ok(())
    }

    pub fn write_table<W: Write>(&self, writer: &mut W, args: &[&dyn Any]) -> io::Result<()> {

        let mut use_scientific = false;
        let mut use_strand     = false;

        for arg in args {
            if let Some(option) = arg.downcast_ref::<OptionPrintScientific>() {
                use_scientific = option.0;
            }
            if let Some(option) = arg.downcast_ref::<OptionPrintStrand>() {
                use_strand = option.0;
            }
        }

        let mut gwriter = GRangesTableWriter::new(self, use_scientific, use_strand);

        let meta_str = self.meta.print_table(args);
        let mut meta_reader = BufReader::new(meta_str.as_bytes());

        gwriter.determine_widths()?;
        gwriter.write_header(writer, &mut meta_reader)?;

        for i in 0..self.num_rows() {
            gwriter.write_row(writer, &mut meta_reader, i)?;
        }
        Ok(())
    }

}

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn read_table<R: Read>(&mut self, reader: R, names: &[&str], types: &[&str]) -> io::Result<()> {

        let buf_reader = BufReader::new(reader);

        self.bufread_table(buf_reader, names, types)

    }

}

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn import_table(&mut self, filename: &str, names: &[&str], types: &[&str], compress: bool) -> io::Result<()> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        self.bufread_table(&mut reader, names, types)?;

        Ok(())
    }

    pub fn export_table(&self, filename: &str, compress: bool, args: &[&dyn Any]) -> io::Result<()> {
        let file = File::create(filename)?;
        let mut writer: Box<dyn Write> = if compress {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };

        self.write_table(&mut writer, args)?;

        Ok(())
    }

}
