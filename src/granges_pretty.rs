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

use std::io::{self, BufRead, BufReader, BufWriter, Write};

use crate::granges::GRanges;

/* -------------------------------------------------------------------------- */

impl GRanges {

    /// Prints the contents of the `GRanges` instance in a pretty-printed format to the provided writer.
    ///
    /// This method formats the `GRanges` data into a human-readable string and writes it to the specified
    /// output stream. It calculates the maximum widths of each column to ensure proper alignment in the output.
    ///
    /// # Arguments
    /// - `writer`: A mutable reference to a writer that implements the `Write` trait, where the formatted output will be written.
    /// - `n`: The maximum number of rows to print. This allows for controlling the output size and is useful for
    ///         displaying a limited number of rows.
    ///
    /// # Returns
    /// An `io::Result<()>`, which will be `Ok(())` if the operation succeeds or an error if the writing fails.
    fn print_pretty<W: Write>(&self, writer: &mut W, n: usize) -> io::Result<()> {
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

    /// Formats the contents of the `GRanges` instance into a pretty-printed string.
    ///
    /// This method collects the formatted representation of the `GRanges` data into a string, which can be
    /// useful for displaying or logging purposes. It internally calls `print_pretty` to handle the formatting.
    ///
    /// # Arguments
    /// - `n`: The maximum number of rows to include in the formatted output.
    ///
    /// # Returns
    /// An `io::Result<String>` containing the formatted string if successful, or an error if the formatting fails.
    pub fn format_pretty(&self, n: usize) -> io::Result<String> {
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

fn write_meta_row<W: Write>(granges: &GRanges, writer: &mut W, reader: &mut dyn BufRead) -> io::Result<()> {
    if granges.meta.num_cols() > 0 {
        write!(writer, " |")?;
        let mut line = String::new();
        reader.read_line(&mut line)?;
        write!(writer, "{}", line.trim_end_matches('\n'))?;
    }
    Ok(())
}

fn write_header<W: Write>(granges: &GRanges, writer: &mut W, meta_reader: &mut dyn BufRead, widths: &[usize]) -> io::Result<()> {
    write!(writer,
        "{:width0$} {:width1$} {:width2$} {:width3$}",
        "", "seqnames", "ranges", "strand",
        width0=widths[0], width1=widths[1], width2=widths[2], width3=widths[3])?;
    write_meta_row(granges, writer, meta_reader)?;
    Ok(())
}

fn write_row<W: Write>(granges: &GRanges, writer: &mut W, meta_reader: &mut dyn BufRead, widths: &[usize], i: usize) -> io::Result<()> {
    writeln!(writer)?;
    write!(writer,
        "{:width0$} {:width1$} [{:width2$}, {:width3$}) {:width4$}",
        i + 1, granges.seqnames[i], granges.ranges[i].from, granges.ranges[i].to, granges.strand[i],
        width0=widths[0], width1=widths[1], width2=widths[2], width3=widths[3], width4=widths[4])?;
    write_meta_row(granges, writer, meta_reader)?;
    Ok(())
}

fn write_all<W: Write>(granges: &GRanges, writer: &mut W, meta_reader: &mut dyn BufRead, widths_header : &[usize], widths_row: &[usize], n: usize) -> io::Result<()> {
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
        write_meta_row(granges, writer, meta_reader)?;
        // Print last n/2 rows
        for i in granges.num_rows() - n / 2..granges.num_rows() {
            write_row(granges, writer, meta_reader, &widths_row, i)?;
        }
    }
    Ok(())
}
