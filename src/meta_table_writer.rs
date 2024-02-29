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

 use std::io;
 use std::io::{BufWriter, Write};

 use crate::meta::{Meta, MetaData};
 
 /* -------------------------------------------------------------------------- */
 
 pub struct MetaTableWriter<'a> {
    meta          : &'a Meta,
    widths        : Vec<usize>,
    use_scientific: bool
 }
 
/* -------------------------------------------------------------------------- */
 
impl<'a> MetaTableWriter<'a> {

    pub fn new(meta: &'a Meta, use_scientific: bool) -> Self {
        MetaTableWriter{meta: meta, widths: vec![0; meta.num_cols()], use_scientific: use_scientific}
    }

    pub fn determine_widths(&mut self) -> io::Result<()> {
        for j in 0..self.meta.num_cols() {
            self.widths[j] = self.meta.meta_name[j].len()+1;
        }
        for i in 0..self.meta.num_rows() {
            meta_update_max_widths(self.meta, i, &mut self.widths, self.use_scientific)?;
        }
        Ok(())
    }

    pub fn write_header<W: Write>(&self, writer: &mut W) -> io::Result<()> {

        for j in 0..self.meta.num_cols() {
            write!(writer, " {:width$}", self.meta.meta_name[j], width = self.widths[j] - 1)?;
        }
        writeln!(writer)
    
    }

    pub fn write_row<W: Write>(&self, writer: &mut W, i: usize) -> io::Result<()> {
        if i != 0 {
            writeln!(writer)?;
        }
        for j in 0..self.meta.num_cols() {
            print_meta_cell(self.meta, writer, i, j, &self.widths, self.use_scientific)?;
        }
        Ok(())
    }    

}

/* -------------------------------------------------------------------------- */

fn meta_update_max_widths(meta: &Meta, i: usize, widths: &mut [usize], use_scientific: bool) -> io::Result<()> {
    for j in 0..meta.num_cols() {
        let mut tmp_buffer = Vec::new();
        {
            let mut tmp_writer = BufWriter::new(&mut tmp_buffer);
            print_meta_cell(meta, &mut tmp_writer, i, j, widths, use_scientific)?;
            tmp_writer.flush()?;
        }
        let width = tmp_buffer.len();
        if width > widths[j] {
            widths[j] = width;
        }
    }
    Ok(())
}

fn print_meta_cell_slice<W: Write>(data: &MetaData, writer: &mut W, i: usize, j: usize, widths: &[usize], use_scientific: bool) -> io::Result<()> {
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

fn print_meta_cell<W: Write>(meta: &Meta, writer: &mut W, i: usize, j: usize, widths: &[usize], use_scientific: bool) -> io::Result<()> {
    match &meta.meta_data[j] {
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
        _ => print_meta_cell_slice(&meta.meta_data[j], writer, i, j, widths, use_scientific),
    }
}
