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
        for j in 0..self.granges.meta.num_cols() {
            write!(writer, " {:width$}", self.granges.meta.meta_name[j], width = self.widths[j] - 1)?;
        }
        write_row_meta(self.granges, writer, meta_reader)?;
        writeln!(writer)    
    }

    pub fn write_row<R: BufRead, W: Write + ?Sized>(&self, writer: &mut W, meta_reader: &mut R, i: usize) -> io::Result<()> {
        if i != 0 {
            writeln!(writer)?;
        }
        write_row     (self.granges, writer, i, &self.widths, self.use_strand)?;
        write_row_meta(self.granges, writer, meta_reader)?;
        Ok(())
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

fn write_header<W: Write + ?Sized>(writer: &mut W, widths: &[usize], strand: bool) -> io::Result<()> {
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
        write!(writer, " | ")?;
        let mut line = String::new();
        meta_reader.read_line(&mut line)?;
        write!(writer, "{}", line.trim_end_matches('\n'))?;
    }
    Ok(())
}

fn write_row<W: Write + ?Sized>(granges: &GRanges, writer: &mut W, i: usize, widths: &[usize], strand: bool) -> io::Result<()> {
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
