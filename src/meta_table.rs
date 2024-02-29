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
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};

use crate::meta::Meta;
use crate::meta::MetaData;
use crate::meta_table_reader::MetaTableReader;

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct OptionPrintScientific {
    pub value: bool,
}

/* -------------------------------------------------------------------------- */

impl Meta {

    fn print_meta_cell_slice<W: Write>(&self, writer: &mut W, widths: &[usize], i: usize, j: usize, data: &MetaData, use_scientific: bool) -> io::Result<()> {
        let mut tmp_buffer = Vec::new();
        {
            let mut tmp_writer = io::Cursor::new(&mut tmp_buffer);
            match data {
                MetaData::StringMatrix(v) => {
                    for (k, s) in v[i].iter().enumerate() {
                        if k != 0 {
                            write!(tmp_writer, ",")?;
                        }
                        write!(tmp_writer, "{}", s)?;
                    }
                }
                MetaData::FloatMatrix(v) => {
                    for (k, f) in v[i].iter().enumerate() {
                        if k != 0 {
                            write!(tmp_writer, ",")?;
                        }
                        if use_scientific {
                            write!(tmp_writer, "{:e}", f)?;
                        } else {
                            write!(tmp_writer, "{}", f)?;
                        }
                    }
                }
                MetaData::IntMatrix(v) => {
                    for (k, i) in v[i].iter().enumerate() {
                        if k != 0 {
                            write!(tmp_writer, ",")?;
                        }
                        write!(tmp_writer, "{}", i)?;
                    }
                }
                _ => unreachable!(),
            }
        }
        write!(writer, " {:width$}", String::from_utf8(tmp_buffer).unwrap(), width = widths[j] - 1)
    }

    fn print_meta_cell<W: Write>(&self, writer: &mut W, widths: &[usize], i: usize, j: usize, use_scientific: bool) -> io::Result<()> {
        match &self.meta_data[j] {
            MetaData::StringArray(v) => {
                write!(writer, " {:width$}", v[i], width = widths[j] - 1)
            }
            MetaData::FloatArray(v) => {
                if use_scientific {
                    write!(writer, " {:width$e}", v[i], width = widths[j] - 1)
                } else {
                    write!(writer, " {:width$}", v[i], width = widths[j] - 1)
                }
            }
            MetaData::IntArray(v) => {
                write!(writer, " {:width$}", v[i], width = widths[j] - 1)
            }
            MetaData::RangeArray(v) => {
                let s = format!("[{},{})", v[i].from, v[i].to);
                write!(writer, " {:width$}", s, width = widths[j] - 1)
            }
            _ => self.print_meta_cell_slice(writer, widths, i, j, &self.meta_data[j], use_scientific),
        }
    }

    fn print_meta_row<W: Write>(&self, writer: &mut W, widths: &[usize], i: usize, use_scientific: bool) -> io::Result<()> {
        if i != 0 {
            writeln!(writer)?;
        }
        for j in 0..self.num_cols() {
            self.print_meta_cell(writer, widths, i, j, use_scientific)?;
        }
        Ok(())
    }

    fn meta_update_max_widths(&self, i: usize, widths: &mut [usize], use_scientific: bool) -> io::Result<()> {
        for j in 0..self.num_cols() {
            let mut tmp_buffer = Vec::new();
            {
                let mut tmp_writer = BufWriter::new(&mut tmp_buffer);
                self.print_meta_cell(&mut tmp_writer, widths, i, j, use_scientific)?;
                tmp_writer.flush()?;
            }
            let width = tmp_buffer.len();
            if width > widths[j] {
                widths[j] = width;
            }
        }
        Ok(())
    }

    fn print_meta_header<W: Write>(&self, writer: &mut W, widths: &[usize]) -> io::Result<()> {
        for j in 0..self.num_cols() {
            write!(writer, " {:width$}", self.meta_name[j], width = widths[j] - 1)?;
        }
        writeln!(writer)
    }

    pub fn write_table<W: Write>(&self, writer: &mut W, args: &[&dyn Any]) -> io::Result<()> {
        let mut use_scientific = false;
        for arg in args {
            if let Some(option) = arg.downcast_ref::<OptionPrintScientific>() {
                use_scientific = option.value;
            }
        }

        let mut widths = vec![0; self.num_cols()];
        for j in 0..self.num_cols() {
            widths[j] = self.meta_name[j].len()+1;
        }

        for i in 0..self.num_rows() {
            self.meta_update_max_widths(i, &mut widths, use_scientific)?;
        }
        self.print_meta_header(writer, &widths)?;

        for i in 0..self.num_rows() {
            self.print_meta_row(writer, &widths, i, use_scientific)?;
        }
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
}

