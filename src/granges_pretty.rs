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

//use std::io::BufReader;
use std::io::BufWriter;
use std::io::Result;
use std::io::Write;

use crate::granges::{GRanges};

/* -------------------------------------------------------------------------- */

impl GRanges {

    fn update_max_width(&self, widths: &mut [usize], j: usize, args: String) {
        let width = args.len();
        if width > widths[j] {
            widths[j] = width;
        }
    }

    fn update_max_widths(&self, i: usize, widths: &mut [usize; 5]) {
        self.update_max_width(widths, 0, (i + 1)              .to_string());
        self.update_max_width(widths, 1, self.seqnames[i]     .to_string());
        self.update_max_width(widths, 2, self.ranges  [i].from.to_string());
        self.update_max_width(widths, 3, self.ranges  [i].to  .to_string());
        self.update_max_width(widths, 4, self.strand  [i]     .to_string());
    }

    fn print_header(&self, writer: &mut dyn Write, widths: &[usize]) -> Result<()> {
        write!(writer,
            "{:width0$} {:width1$} {:width2$} {:width3$}",
            "", "seqnames", "ranges", "strand",
            width0=widths[0], width1=widths[1], width2=widths[2], width3=widths[3])?;
        //print_meta_row(writer)?;
        Ok(())
    }

    fn print_row(&self, writer: &mut dyn Write, widths: &[usize], i: usize) -> Result<()> {
        writeln!(writer)?;
        write!(writer,
            "{:width0$} {:width1$} [{:width2$}, {:width3$}) {:width4$}",
            i + 1, self.seqnames[i], self.ranges[i].from, self.ranges[i].to, self.strand[i],
            width0=widths[0], width1=widths[1], width2=widths[2], width3=widths[3], width4=widths[4])?;
        //print_meta_row(writer)?;
        Ok(())
    }

    fn print_all(&self, writer: &mut dyn Write, widths_header : &[usize], widths_row: &[usize], n: usize) -> Result<()> {
        if self.length() <= n + 1 {
            for i in 0..self.length() {
                self.print_row(writer, &widths_row, i)?;
            }
        } else {
            // Print first n/2 rows
            for i in 0..n / 2 {
                self.print_row(writer, &widths_row, i)?;
            }
            // Print gap
            writeln!(writer)?;
            write!(writer, 
                "{:width0$} {:width1$} {:width2$} {:width3$}",
                "", "...", "...", "",
                width0=widths_header[0], width1=widths_header[1], width2=widths_header[2], width3=widths_header[3])?;
            // Print last n/2 rows
            for i in self.length() - n / 2..self.length() {
                self.print_row(writer, &widths_row, i)?;
            }
        }
        Ok(())
    }

    fn write_pretty(&self, writer: &mut dyn Write, n: usize) -> Result<()> {
        //let meta_str = format!("{}", self.meta);
        //let meta_reader = BufReader::new(meta_str.as_bytes());

        //let print_meta_row = |writer: &mut dyn Write| -> Result<()> {
        //    if self.meta_length() != 0 {
        //        write!(writer, " | ")?;
        //        let mut line = String::new();
        //        meta_reader.read_line(&mut line)?;
        //        write!(writer, "{}", line)?;
        //    }
        //    Ok(())
        //};

        let mut widths : [usize; 5] = [1, 8, 1, 1, 6];

        for i in 0..self.length() {
            self.update_max_widths(i, &mut widths);
        }

        let widths_row    = [widths[0], widths[1], widths[2], widths[3], widths[4]];
        let widths_header = [widths[0], widths[1], widths[2] + widths[3] + 4, widths[4]];

        self.print_header(writer, &widths_header)?;
        self.print_all   (writer, &widths_header, &widths_row, n)?;

        Ok(())
    }

    pub fn print_pretty(&self, n: usize) -> String {
        let mut buffer = Vec::new();
        {
            let mut writer = BufWriter::new(&mut buffer);

            if let Err(_) = self.write_pretty(&mut writer, n) {
                return String::new();
            }

            writer.flush().unwrap();
        }

        String::from_utf8(buffer).unwrap()
    }
}
