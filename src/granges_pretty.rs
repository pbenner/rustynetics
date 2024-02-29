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
use std::io::BufReader;
use std::io::BufWriter;
use std::io::BufRead;
use std::io::Write;

use crate::granges::GRanges;
use crate::error::Error;

/* -------------------------------------------------------------------------- */

impl GRanges {

    fn print_pretty(&self, writer: &mut dyn Write, n: usize) -> io::Result<()> {
        let meta_str = format!("{}", self.meta);
        let mut meta_reader = BufReader::new(meta_str.as_bytes());

        let mut widths : [usize; 5] = [1, 8, 1, 1, 6];

        for i in 0..self.num_rows() {
            update_max_widths(self, i, &mut widths);
        }

        let widths_row    = [widths[0], widths[1], widths[2], widths[3], widths[4]];
        let widths_header = [widths[0], widths[1], widths[2] + widths[3] + 4, widths[4]];

        write_header(self, writer, &mut meta_reader, &widths_header)?;
        write_all   (self, writer, &mut meta_reader, &widths_header, &widths_row, n)?;

        Ok(())
    }

    pub fn pretty_string(&self, n: usize) -> Result<String, Error> {
        let mut buffer = Vec::new();
        {
            let mut writer = BufWriter::new(&mut buffer);

            self.print_pretty(&mut writer, n)?;

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

fn update_max_width(widths: &mut [usize], j: usize, args: String) {
    let width = args.len();
    if width > widths[j] {
        widths[j] = width;
    }
}

fn update_max_widths(granges: &GRanges, i: usize, widths: &mut [usize; 5]) {
    update_max_width(widths, 0, (i + 1)                 .to_string());
    update_max_width(widths, 1, granges.seqnames[i]     .to_string());
    update_max_width(widths, 2, granges.ranges  [i].from.to_string());
    update_max_width(widths, 3, granges.ranges  [i].to  .to_string());
    update_max_width(widths, 4, granges.strand  [i]     .to_string());
}

fn write_meta_row(granges: &GRanges, writer: &mut dyn Write, reader: &mut dyn BufRead) -> io::Result<()> {
    if granges.meta.num_cols() > 0 {
        write!(writer, " |")?;
        let mut line = String::new();
        reader.read_line(&mut line)?;
        write!(writer, "{}", line.trim_end_matches('\n'))?;
    }
    Ok(())
}

fn write_header(granges: &GRanges, writer: &mut dyn Write, meta_reader: &mut dyn BufRead, widths: &[usize]) -> io::Result<()> {
    write!(writer,
        "{:width0$} {:width1$} {:width2$} {:width3$}",
        "", "seqnames", "ranges", "strand",
        width0=widths[0], width1=widths[1], width2=widths[2], width3=widths[3])?;
    write_meta_row(granges, writer, meta_reader)?;
    Ok(())
}

fn write_row(granges: &GRanges, writer: &mut dyn Write, meta_reader: &mut dyn BufRead, widths: &[usize], i: usize) -> io::Result<()> {
    writeln!(writer)?;
    write!(writer,
        "{:width0$} {:width1$} [{:width2$}, {:width3$}) {:width4$}",
        i + 1, granges.seqnames[i], granges.ranges[i].from, granges.ranges[i].to, granges.strand[i],
        width0=widths[0], width1=widths[1], width2=widths[2], width3=widths[3], width4=widths[4])?;
    write_meta_row(granges, writer, meta_reader)?;
    Ok(())
}

fn write_all(granges: &GRanges, writer: &mut dyn Write, meta_reader: &mut dyn BufRead, widths_header : &[usize], widths_row: &[usize], n: usize) -> io::Result<()> {
    if granges.num_rows() <= n + 1 {
        for i in 0..granges.num_rows() {
            write_row(granges, writer, meta_reader, &widths_row, i)?;
        }
    } else {
        // Print first n/2 rows
        for i in 0..n / 2 {
            write_row(granges, writer, meta_reader, &widths_row, i)?;
        }
        // Print gap
        writeln!(writer)?;
        write!(writer, 
            "{:width0$} {:width1$} {:width2$} {:width3$}",
            "", "...", "...", "",
            width0=widths_header[0], width1=widths_header[1], width2=widths_header[2], width3=widths_header[3])?;
        // Print last n/2 rows
        for i in granges.num_rows() - n / 2..granges.num_rows() {
            write_row(granges, writer, meta_reader, &widths_row, i)?;
        }
    }
    Ok(())
}
