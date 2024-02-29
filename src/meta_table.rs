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

use std::any::Any;
use std::io::{self, BufRead, BufReader, Read, Write};

use crate::meta::Meta;
use crate::meta_table_reader::MetaTableReader;
use crate::meta_table_writer::MetaTableWriter;

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct OptionPrintScientific {
    pub value: bool,
}

/* -------------------------------------------------------------------------- */

impl Meta {

    pub fn write_table<W: Write>(&self, writer: &mut W, args: &[&dyn Any]) -> io::Result<()> {
        let mut use_scientific = false;
        for arg in args {
            if let Some(option) = arg.downcast_ref::<OptionPrintScientific>() {
                use_scientific = option.value;
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

