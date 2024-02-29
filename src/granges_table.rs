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

use std::any::Any;
use std::io::{self, BufRead, BufReader, Write};
use std::fs::File;

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;

use crate::granges::GRanges;
use crate::granges_table_reader::GRangesTableReader;
use crate::meta_table_reader::MetaTableReader;

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn write_table(&self, writer: &mut dyn Write, strand: bool, args: &[&dyn Any]) -> io::Result<()> {
        let meta_str = self.meta.print_table(args);
        let mut mreader = BufReader::new(meta_str.as_bytes());

        let mut widths = vec![8, 4, 2, 6];
        for i in 0..self.num_rows() {
            update_max_widths(self, i, &mut widths, strand)?;
        }
        write_header(writer, &widths, strand)?;
        self.meta_print_table_row(writer, &mut mreader)?;
        writeln!(writer)?;

        for i in 0..self.num_rows() {
            write_row(self, writer, &widths, i, strand)?;
            self.meta_print_table_row(writer, &mut mreader)?;
            writeln!(writer)?;
        }
        Ok(())
    }

    fn meta_print_table_row(&self, writer: &mut dyn Write, reader: &mut dyn BufRead) -> io::Result<()> {
        if self.meta.num_cols() > 0 {
            write!(writer, " | ")?;
            let mut line = String::new();
            reader.read_line(&mut line)?;
            write!(writer, "{}", line.trim_end_matches('\n'))?;
        }
        Ok(())
    }

    pub fn export_table(&self, filename: &str, strand: bool, compress: bool, args: &[&dyn Any]) -> io::Result<()> {
        let file = File::create(filename)?;
        let mut writer: Box<dyn Write> = if compress {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };

        self.write_table(&mut writer, strand, args)?;

        Ok(())
    }

    fn read_table(&mut self, reader: &mut dyn BufRead, names: &[&str], types: &[&str]) -> io::Result<()> {
        let mut mreader = MetaTableReader   ::new(names, types);
        let mut greader = GRangesTableReader::new();

        let mut buf_reader = BufReader::new(reader);
        let mut line = String::new();

        // Read first line as header
        buf_reader.read_line(&mut line)?;
        greader.read_header(&line)?;
        mreader.read_header(&line)?;

        line.clear();

        let mut line_counter = 0;

        while buf_reader.read_line(&mut line)? > 0 {

            if line.is_empty() {
                continue;
            }
            greader.read_line(&line, line_counter)?;
            mreader.read_line(&line, line_counter)?;

            line_counter += 1;

            line.clear();
        }
        greader.push(self);
        mreader.push(&mut self.meta);

        Ok(())
    }

    pub fn import_table(&mut self, filename: &str, names: &[&str], types: &[&str], compress: bool) -> io::Result<()> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        self.read_table(&mut reader, names, types)?;

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

fn write_header(writer: &mut dyn Write, widths: &[usize], strand: bool) -> io::Result<()> {
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

fn write_row(granges: &GRanges, writer: &mut dyn Write, widths: &[usize], i: usize, strand: bool) -> io::Result<()> {
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
