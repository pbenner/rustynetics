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

use crate::meta::Meta;
use crate::meta_table_reader::MetaTableReader;
use crate::meta_table_writer::MetaTableWriter;
use crate::granges_table::OptionPrintScientific;

/* -------------------------------------------------------------------------- */

impl Meta {

    pub fn write_table<W: Write>(&self, writer: &mut W, args: &[&dyn Any]) -> io::Result<()> {
        let mut use_scientific = false;
        for arg in args {
            if let Some(option) = arg.downcast_ref::<OptionPrintScientific>() {
                use_scientific = option.0;
            }
        }
        let mut mwriter = MetaTableWriter::new(self, use_scientific);

        mwriter.determine_widths()?;
        mwriter.write_header(writer)?;

        for i in 0..self.num_rows() {
            mwriter.write_row(writer, i)?;
        }
        Ok(())
    }

    pub fn read_table<R: Read>(&mut self, reader: R, names: &[&str], types: &[&str]) -> io::Result<()> {

        let mut mreader = MetaTableReader::new(names, types);

        let mut buf_reader = BufReader::new(reader);

        // Parse header
        let mut line = String::new();
        buf_reader.read_line(&mut line)?;
        mreader.read_header(&line)?;

        line.clear();

        let mut line_counter = 0;

        while buf_reader.read_line(&mut line)? > 0 {
 
            if line.is_empty() {
                continue;
            }
            mreader.read_line(&line, line_counter)?;

            line.clear();

            line_counter += 1;
        }
        mreader.push(self);

        Ok(())
    }

    pub fn print_table(&self, args: &[&dyn Any]) -> String {
        let mut buffer = Vec::new();
        {
            let mut writer = io::Cursor::new(&mut buffer);
            self.write_table(&mut writer, args).unwrap();
        }
        String::from_utf8(buffer).unwrap()
    }
}

