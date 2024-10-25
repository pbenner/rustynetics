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
 use std::io::{BufRead, Write};

 use crate::granges::GRanges;
 
 /* -------------------------------------------------------------------------- */
 
 pub struct GRangesTableWriter<'a> {
    granges       : &'a GRanges,
    widths        : Vec<usize>,
    use_scientific: bool,
    use_strand    : bool,
 }
 
/* -------------------------------------------------------------------------- */
 
impl<'a> GRangesTableWriter<'a> {

    pub fn new(granges: &'a GRanges, use_scientific: bool, use_strand: bool) -> Self {
        GRangesTableWriter{
            granges       : granges,
            widths        : vec![8, 4, 2, 6],
            use_scientific: use_scientific,
            use_strand    : use_strand,
        }
    }

    pub fn determine_widths(&mut self) -> io::Result<()> {
        for i in 0..self.granges.num_rows() {
            update_max_widths(self.granges, i, &mut self.widths, self.use_scientific)?;
        }
        Ok(())
    }

    pub fn write_header<R: BufRead, W: Write>(&self, writer: &mut W, meta_reader: &mut R) -> io::Result<()> {
        write_header(writer, &self.widths, self.use_strand)?;
        write_row_meta(self.granges, writer, meta_reader)?;
        writeln!(writer)    
    }

    pub fn write_row<R: BufRead, W: Write>(&self, writer: &mut W, meta_reader: &mut R, i: usize) -> io::Result<()> {
        write_row     (self.granges, writer, i, &self.widths, self.use_strand)?;
        write_row_meta(self.granges, writer, meta_reader)?;
        writeln!(writer)    
    }

}

/* -------------------------------------------------------------------------- */

fn update_max_widths(granges: &GRanges, i: usize, widths: &mut [usize], strand: bool) -> io::Result<()> {
    let seqname_width = granges.seqnames[i].len();
    if seqname_width > widths[0] {
        widths[0] = seqname_width;
    }
    let from_width = granges.ranges[i].from.to_string().len();
    if from_width > widths[1] {
        widths[1] = from_width;
    }
    let to_width = granges.ranges[i].to.to_string().len();
    if to_width > widths[2] {
        widths[2] = to_width;
    }
    if strand {
        let strand_width = granges.strand[i].to_string().len();
        if strand_width > widths[3] {
            widths[3] = strand_width;
        }
    }
    Ok(())
}

fn write_header<W: Write>(writer: &mut W, widths: &[usize], strand: bool) -> io::Result<()> {
    if strand {
        write!(writer,
            "{:width0$} {:width1$} {:width2$} {:width3$}",
            "seqnames", "from", "to", "strand",
            width0=widths[0], width1=widths[1], width2=widths[2], width3=widths[3])?;
    } else {
        write!(writer,
            "{:width0$} {:width1$} {:width2$}",
            "seqnames", "from", "to",
            width0=widths[0], width1=widths[1], width2=widths[2])?;
    }
    Ok(())
}

fn write_row_meta(granges: &GRanges, writer: &mut dyn Write, meta_reader: &mut dyn BufRead) -> io::Result<()> {
    if granges.meta.num_cols() > 0 {
        let mut line = String::new();
        meta_reader.read_line(&mut line)?;
        write!(writer, " | {}", line.trim_end_matches('\n'))?;
    }
    Ok(())
}

fn write_row<W: Write>(granges: &GRanges, writer: &mut W, i: usize, widths: &[usize], strand: bool) -> io::Result<()> {
    if strand {
        write!(writer,
            "{:width0$} {:width1$} {:width2$} {:width3$}",
            granges.seqnames[i], granges.ranges[i].from, granges.ranges[i].to, granges.strand[i],
            width0=widths[0], width1=widths[1], width2=widths[2], width3=widths[3])?;
    } else {
        write!(writer,
            "{:width0$} {:width1$} {:width2$}",
            granges.seqnames[i], granges.ranges[i].from, granges.ranges[i].to,
            width0=widths[0], width1=widths[1], width2=widths[2])?;
    }
    Ok(())
}
