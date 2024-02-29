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
use std::io::{self, BufRead, BufReader, Seek, Write};
use std::fs::File;
use std::str::FromStr;

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;

use crate::range::Range;
use crate::granges::GRanges;

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn write_table(&self, writer: &mut dyn Write, strand: bool, args: &[&dyn Any]) -> io::Result<()> {
        let meta_str = self.meta.print_table(args);
        let mut meta_reader = BufReader::new(meta_str.as_bytes());

        let mut widths = vec![8, 4, 2, 6];
        for i in 0..self.num_rows() {
            update_max_widths(i, &mut widths, self, strand)?;
        }
        print_header(writer, &widths, strand)?;
        self.meta_print_table_row(writer, &mut meta_reader)?;
        writeln!(writer)?;

        for i in 0..self.num_rows() {
            print_row(writer, &widths, i, self, strand)?;
            self.meta_print_table_row(writer, &mut meta_reader)?;
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
        let mut buf_reader = BufReader::new(reader);
        let mut line = String::new();

        let mut col_seqname = -1;
        let mut col_from    = -1;
        let mut col_to      = -1;
        let mut col_strand  = -1;

        buf_reader.read_line(&mut line)?;
        let fields: Vec<&str> = line.trim().split_whitespace().collect();
        for (i, field) in fields.iter().enumerate() {
            match *field {
                "seqnames" => col_seqname = i as i32,
                "from"     => col_from    = i as i32,
                "to"       => col_to      = i as i32,
                "strand"   => col_strand  = i as i32,
                "start" if col_from == -1 => col_from = i as i32,
                "end"   if col_to   == -1 => col_to   = i as i32,
                _ => (),
            }
        }
        line.clear();

        if col_seqname == -1 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "is missing a seqnames column"));
        }
        if col_from == -1 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "is missing a from column"));
        }
        if col_to == -1 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "is missing a to column"));
        }

        self.seqnames.clear();
        self.ranges  .clear();
        self.strand  .clear();

        let mut line_counter = 0;

        while buf_reader.read_line(&mut line)? > 0 {
            let fields: Vec<&str> = line.trim().split_whitespace().collect();
            if fields.len() < col_seqname as usize || fields.len() < col_from as usize || fields.len() < col_to as usize {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid table"));
            }

            let seqname = fields[col_seqname as usize].to_string();
            let from    = usize::from_str(fields[col_from as usize]).map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("parsing `from' column `{}` failed at line `{}`", col_from + 1, line_counter + 1)))?;
            let to      = usize::from_str(fields[col_to   as usize]).map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("parsing `to' column `{}` failed at line `{}`"  , col_to   + 1, line_counter + 1)))?;

            self.seqnames.push(seqname);
            self.ranges  .push(Range{ from, to });

            if col_strand != -1 {
                if fields.len() < col_strand as usize {
                    return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid table"));
                }
                let strand = fields[col_strand as usize].chars().next().unwrap();
                self.strand.push(strand);
            } else {
                self.strand.push('*');
            }
            line_counter += 1;

            line.clear();
        }
        self.meta.read_table(buf_reader, names, types)?;
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

fn update_max_widths(i: usize, widths: &mut [usize], granges: &GRanges, strand: bool) -> io::Result<()> {
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

fn print_header(writer: &mut dyn Write, widths: &[usize], strand: bool) -> io::Result<()> {
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

fn print_row(writer: &mut dyn Write, widths: &[usize], i: usize, granges: &GRanges, strand: bool) -> io::Result<()> {
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
