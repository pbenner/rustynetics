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

/* -------------------------------------------------------------------------- */

use std::io;
use std::io::BufWriter;
use std::io::Write;

use crate::error::Error;
use crate::meta::Meta;
use crate::meta::MetaData;

/* -------------------------------------------------------------------------- */

impl Meta {

    fn update_max_widths(&self, i: usize, widths: &mut [usize], use_scientific: bool) -> io::Result<()> {
        for j in 0..self.num_cols() {
            let mut tmp_buffer = Vec::new();
            {
                let mut tmp_writer = BufWriter::new(&mut tmp_buffer);
                self.print_cell(&mut tmp_writer, widths, i, j, use_scientific)?;
                tmp_writer.flush()?;
            }
            let width = tmp_buffer.len();
            if width > widths[j] {
                widths[j] = width;
            }
        }
        Ok(())
    }

    fn print_cell_slice(&self,
                        writer: &mut dyn Write,
                        widths: &[usize],
                        i: usize,
                        j: usize,
                        data: &MetaData,
                        use_scientific: bool)
        -> io::Result<()> {
        let mut buffer = Vec::new();

        match data {
            MetaData::StringMatrix(v) => {
                for k in 0..v[i].len() {
                    write!(buffer, " {}", v[i][k])?;
                }
            }
            MetaData::FloatMatrix(v) => {
                if use_scientific {
                    for k in 0..v[i].len() {
                        write!(buffer, " {:e}", v[i][k])?;
                    }
                } else {
                    for k in 0..v[i].len() {
                        write!(buffer, " {:.2}", v[i][k])?;
                    }
                }
            }
            MetaData::IntMatrix(v) => {
                for k in 0..v[i].len() {
                    write!(buffer, " {}", v[i][k])?;
                }
            }
            _ => unreachable!(),
        }
        write!(writer, "{:width$}" , String::from_utf8(buffer).unwrap(), width = widths[j] - 1)
    }

    fn print_cell(&self,
                  writer: &mut dyn Write,
                  widths: &[usize],
                  i: usize,
                  j: usize,
                  use_scientific: bool)
        -> io::Result<()> {
        match &self.meta_data[j] {
            MetaData::StringArray(v) => write!(writer, " {:width$}", v[i], width = widths[j] - 1),
            MetaData::FloatArray(v)  => {
                if use_scientific {
                    write!(writer, " {:width$e}" , v[i], width = widths[j] - 1)
                } else {
                    write!(writer, " {:width$.2}", v[i], width = widths[j] - 1)
                }
            }
            MetaData::IntArray(v)   => write!(writer, " {:width$}", v[i], width = widths[j] - 1),
            MetaData::RangeArray(v) => {
                write!(writer, " {:width$}", v[i], width = widths[j] - 1)
            },
            _ => self.print_cell_slice(writer, widths, i, j, &self.meta_data[j], use_scientific),
        }
    }

    fn print_header(&self, writer: &mut dyn Write, widths: &[usize]) -> io::Result<()> {
        for j in 0..self.num_cols() {
            write!(writer, " {:width$}", self.meta_name[j], width = widths[j] - 1)?;
        }
        writeln!(writer)
    }

    fn print_row(&self, writer: &mut dyn Write, widths: &[usize], i: usize, use_scientific: bool) -> io::Result<()> {
        if i != 0 {
            writeln!(writer)?;
        }
        for j in 0..self.num_cols() {
            self.print_cell(writer, widths, i, j, use_scientific)?;
        }
        Ok(())
    }

    fn print_all(&self, writer: &mut dyn Write, widths : &[usize], n: usize, use_scientific: bool) -> io::Result<()> {
        if self.num_rows() <= n + 1 {
            for i in 0..self.num_rows() {
                self.print_row(writer, &widths, i, use_scientific)?;
            }
        } else {
            // Print first n/2 rows
            for i in 0..n / 2 {
                self.print_row(writer, &widths, i, use_scientific)?;
            }
            // Print gap
            writeln!(writer)?;
            for j in 0..self.num_cols() {
                write!(writer, " {:>width$}", "...", width = widths[j] - 1)?;
            }
            // Print last n/2 rows
            for i in self.num_rows() - n / 2..self.num_rows() {
                self.print_row(writer, &widths, i, use_scientific)?;
            }
        }
        Ok(())
    }

    fn print_pretty(&self, n: usize, use_scientific: bool) -> io::Result<Vec<u8>> {

        let mut buffer = Vec::new();

        let mut widths = vec![0; self.num_cols()];
        for j in 0..self.num_cols() {
            let width = format!(" {}", self.meta_name[j]).len();
            widths[j] = width;
        }

        for i in 0..self.num_rows() {
            self.update_max_widths(i, &mut widths, use_scientific)?;
        }

        self.print_header(&mut buffer, &widths)?;
        self.print_all   (&mut buffer, &widths, n, use_scientific)?;

        Ok(buffer)
    }

    pub fn pretty_string(&self, n: usize, use_scientific: bool) -> Result<String, Error> {
        let r = self.print_pretty(n, use_scientific)?;
        let s = match String::from_utf8(r) {
            Ok (v) => v,
            Err(_) => panic!("internal error")
        };
        Ok(s)
    }

}
