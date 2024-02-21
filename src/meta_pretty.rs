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

use std::fmt::Write;

use crate::error::Error;
use crate::meta::Meta;
use crate::meta::MetaData;

/* -------------------------------------------------------------------------- */

impl Meta {

    fn write_pretty(&self, n: usize, args: &[impl Into<OptionPrintScientific>]) -> String {
        let use_scientific = args
            .iter()
            .any(|arg| matches!(arg.into(), OptionPrintScientific { value: true }));

        let mut buffer = String::new();

        let print_cell_slice = |writer: &mut dyn Write,
                                widths: &[usize],
                                i: usize,
                                j: usize,
                                data: &MetaData|
         -> fmt::Result {
            match data {
                MetaData::StringSlice(v) => {
                    for k in 0..v[i].len() {
                        write!(writer, " {}", v[i][k])?;
                    }
                }
                MetaData::FloatSlice(v) => {
                    if use_scientific {
                        for k in 0..v[i].len() {
                            write!(writer, " {:e}", v[i][k])?;
                        }
                    } else {
                        for k in 0..v[i].len() {
                            write!(writer, " {:.2}", v[i][k])?;
                        }
                    }
                }
                MetaData::IntSlice(v) => {
                    for k in 0..v[i].len() {
                        write!(writer, " {}", v[i][k])?;
                    }
                }
                _ => panic!("invalid meta data"),
            }
            Ok(())
        };

        let print_cell = |writer: &mut dyn Write,
                          widths: &[usize],
                          i: usize,
                          j: usize|
         -> fmt::Result {
            match &self.meta_data[j] {
                MetaData::String(v) => write!(writer, " {:width$}", v[i], width = widths[j] - 1),
                MetaData::Float(v) => {
                    if use_scientific {
                        write!(writer, " {:width$e}", v[i], width = widths[j] - 1)
                    } else {
                        write!(writer, " {:width$.2}", v[i], width = widths[j] - 1)
                    }
                }
                MetaData::Int(v) => write!(writer, " {:width$d}", v[i], width = widths[j] - 1),
                MetaData::Range(v) => write!(writer, " [{}, {}]", v[i].from, v[i].to),
                _ => print_cell_slice(writer, widths, i, j, &self.meta_data[j]),
            }
        };

        let print_row = |writer: &mut dyn Write, widths: &[usize], i: usize| -> fmt::Result {
            if i != 0 {
                writeln!(writer)?;
            }
            for j in 0..self.meta_length() {
                print_cell(writer, widths, i, j)?;
            }
            Ok(())
        };

        let update_max_widths = |i: usize, widths: &mut [usize]| -> fmt::Result {
            for j in 0..self.meta_length() {
                let mut tmp_buffer = String::new();
                let mut tmp_writer = std::io::BufWriter::new(&mut tmp_buffer);
                print_cell(&mut tmp_writer, widths, i, j)?;
                tmp_writer.flush()?;
                let width = tmp_buffer.len();
                if width > widths[j] {
                    widths[j] = width;
                }
            }
            Ok(())
        };

        let print_header = |writer: &mut dyn Write, widths: &[usize]| -> fmt::Result {
            for j in 0..self.meta_length() {
                write!(writer, " {:width$}", self.meta_name[j], width = widths[j] - 1)?;
            }
            writeln!(writer)?;
            Ok(())
        };

        let apply_rows = |f1: &mut dyn FnMut(usize) -> fmt::Result,
                          f2: &mut dyn FnMut() -> fmt::Result|
         -> fmt::Result {
            if self.length() <= n + 1 {
                for i in 0..self.length() {
                    f1(i)?;
                }
            } else {
                for i in 0..n / 2 {
                    f1(i)?;
                }
                f2()?;
                for i in self.length() - n / 2..self.length() {
                    f1(i)?;
                }
            }
            Ok(())
        };

        let mut widths = vec![0; self.meta_length()];
        for j in 0..self.meta_length() {
            let width = format!(" {}", self.meta_name[j]).len();
            widths[j] = width;
        }

        apply_rows(
            &mut |i| update_max_widths(i, &mut widths),
            &mut || Ok(()),
        )?;

        print_header(&mut buffer, &widths)?;
        apply_rows(
            &mut |i| print_row(&mut buffer, &widths, i),
            &mut || {
                writeln!(buffer)?;
                for j in 0..self.meta_length() {
                    write!(buffer, " {:width$}", "...", width = widths[j] - 1)?;
                }
                Ok(())
            },
        )?;

        buffer
    }

    fn print_pretty(&self, n: usize, args: &[impl Into<OptionPrintScientific>]) -> String {
        self.write_pretty(n, args)
    }

}

struct OptionPrintScientific {
    value: bool,
}
