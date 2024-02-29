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

    fn print_pretty<W: Write>(&self, writer: &mut W, n: usize, use_scientific: bool) -> io::Result<()> {

        let mut widths = vec![0; self.num_cols()];
        for j in 0..self.num_cols() {
            widths[j] = self.meta_name[j].len();
        }

        if self.num_rows() <= n + 1 {
            for i in 0..self.num_rows() {
                update_max_widths(self, i, &mut widths, use_scientific)?;
            }
        } else {
            // Print first n/2 rows
            for i in 0..n / 2 {
                update_max_widths(self, i, &mut widths, use_scientific)?;
            }
            // Print last n/2 rows
            for i in self.num_rows() - n / 2..self.num_rows() {
                update_max_widths(self, i, &mut widths, use_scientific)?;
            }
        }

        print_header(self, writer, &widths)?;
        print_all   (self, writer, &widths, n, use_scientific)?;

        Ok(())
    }

    pub fn format_pretty(&self, n: usize, use_scientific: bool) -> Result<String, Error> {
        let mut buffer = Vec::new();
        {
            let mut writer = BufWriter::new(&mut buffer);

            self.print_pretty(&mut writer, n, use_scientific)?;

            writer.flush().unwrap();
        }
        let s = match String::from_utf8(buffer) {
            Ok (v) => v,
            Err(_) => panic!("internal error")
        };
        Ok(s)
    }

}

/* -------------------------------------------------------------------------- */

fn update_max_widths(meta: &Meta, i: usize, widths: &mut [usize], use_scientific: bool) -> io::Result<()> {
    for j in 0..meta.num_cols() {
        let mut tmp_buffer = Vec::new();
        {
            let mut tmp_writer = BufWriter::new(&mut tmp_buffer);
            print_cell(meta, &mut tmp_writer, widths, i, j, use_scientific)?;
            tmp_writer.flush()?;
        }
        let width = tmp_buffer.len()-1;
        if width > widths[j] {
            widths[j] = width;
        }
    }
    Ok(())
}

fn print_cell_slice<W: Write>(writer: &mut W,
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
    write!(writer, "{:>width$}" , String::from_utf8(buffer).unwrap(), width = widths[j]+1)
}

fn print_cell<W: Write>(meta: &Meta,
              writer: &mut W,
              widths: &[usize],
              i: usize,
              j: usize,
              use_scientific: bool)
    -> io::Result<()> {
    match &meta.meta_data[j] {
        MetaData::StringArray(v) => write!(writer, " {:>width$}", v[i], width = widths[j]),
        MetaData::FloatArray(v)  => {
            if use_scientific {
                write!(writer, " {:>width$e}" , v[i], width = widths[j])
            } else {
                write!(writer, " {:>width$.2}", v[i], width = widths[j])
            }
        }
        MetaData::IntArray(v)   => write!(writer, " {:>width$}", v[i], width = widths[j]),
        MetaData::RangeArray(v) => {
            write!(writer, " {:>width$}", v[i], width = widths[j])
        },
        _ => print_cell_slice(writer, widths, i, j, &meta.meta_data[j], use_scientific),
    }
}

fn print_header<W: Write>(meta: &Meta, writer: &mut W, widths: &[usize]) -> io::Result<()> {
    for j in 0..meta.num_cols() {
        write!(writer, " {:>width$}", meta.meta_name[j], width = widths[j])?;
    }
    writeln!(writer)
}

fn print_row<W: Write>(meta: &Meta, writer: &mut W, widths: &[usize], i: usize, use_scientific: bool) -> io::Result<()> {
    if i != 0 {
        writeln!(writer)?;
    }
    for j in 0..meta.num_cols() {
        print_cell(meta, writer, widths, i, j, use_scientific)?;
    }
    Ok(())
}

fn print_all<W: Write>(meta: &Meta, writer: &mut W, widths : &[usize], n: usize, use_scientific: bool) -> io::Result<()> {
    if meta.num_rows() <= n + 1 {
        for i in 0..meta.num_rows() {
            print_row(meta, writer, &widths, i, use_scientific)?;
        }
    } else {
        // Print first n/2 rows
        for i in 0..n / 2 {
            print_row(meta, writer, &widths, i, use_scientific)?;
        }
        // Print gap
        writeln!(writer)?;
        for j in 0..meta.num_cols() {
            write!(writer, " {:>width$}", "...", width = widths[j])?;
        }
        // Print last n/2 rows
        for i in meta.num_rows() - n / 2..meta.num_rows() {
            print_row(meta, writer, &widths, i, use_scientific)?;
        }
    }
    Ok(())
}
