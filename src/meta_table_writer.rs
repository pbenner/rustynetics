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
            update_max_widths(self.meta, i, &mut self.widths, self.use_scientific)?;
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
            write_cell(self.meta, writer, i, j, &self.widths, self.use_scientific)?;
        }
        Ok(())
    }    

}

/* -------------------------------------------------------------------------- */

fn update_max_widths(meta: &Meta, i: usize, widths: &mut [usize], use_scientific: bool) -> io::Result<()> {
    for j in 0..meta.num_cols() {
        let mut tmp_buffer = Vec::new();
        {
            let mut tmp_writer = BufWriter::new(&mut tmp_buffer);
            write_cell(meta, &mut tmp_writer, i, j, widths, use_scientific)?;
            tmp_writer.flush()?;
        }
        let width = tmp_buffer.len()-1;
        if width > widths[j] {
            widths[j] = width;
        }
    }
    Ok(())
}

fn write_cell_slice<W: Write>(data: &MetaData, writer: &mut W, i: usize, j: usize, widths: &[usize], use_scientific: bool) -> io::Result<()> {
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

fn write_cell<W: Write>(meta: &Meta, writer: &mut W, i: usize, j: usize, widths: &[usize], use_scientific: bool) -> io::Result<()> {
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
        _ => write_cell_slice(&meta.meta_data[j], writer, i, j, widths, use_scientific),
    }
}
