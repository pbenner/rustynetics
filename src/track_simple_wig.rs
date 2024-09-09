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

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::{rc::Rc, cell::RefCell};

use crate::track_simple::SimpleTrack;

/* -------------------------------------------------------------------------- */

impl SimpleTrack {

    fn write_wiggle_fixed_step<W: Write>(
        &self,
        writer  : &mut W,
        seqname : &str,
        sequence: &[f64],
    ) -> io::Result<()> {
        let mut gap = true;
        for (i, &value) in sequence.iter().enumerate() {
            if !value.is_nan() {
                if gap {
                    writeln!(
                        writer,
                        "fixedStep chrom={} start={} span={} step={}",
                        seqname,
                        i * self.bin_size + 1,
                        self.bin_size,
                        self.bin_size
                    )?;
                    gap = false;
                }
                writeln!(writer, "{}", value)?;
            } else {
                gap = true;
            }
        }
        Ok(())
    }

    fn write_wiggle_variable_step<W: Write>(
        &self,
        writer  : &mut W,
        seqname : &str,
        sequence: &[f64],
    ) -> io::Result<()> {
        writeln!(
            writer,
            "variableStep chrom={} span={}",
            seqname, self.bin_size
        )?;
        for (i, &value) in sequence.iter().enumerate() {
            if !value.is_nan() {
                writeln!(writer, "{} {}", i * self.bin_size + 1, value)?;
            }
        }
        Ok(())
    }

    pub fn write_wiggle(&self, filename: &str, description: &str) -> io::Result<()> {

        let mut file   = File::create(filename)?;
        let mut writer = io::BufWriter::new(&mut file);

        writeln!(
            writer,
            "track type=wiggle_0 name=\"{}\" description=\"{}\"",
            self.name, description
        )?;

        for (seqname, sequence) in &self.data {
            let n = sequence.borrow().iter().filter(|&&x| x.is_nan() || x == 0.0).count();
            if n >= sequence.borrow().len() / 2 {
                // sparse data track
                self.write_wiggle_variable_step(&mut writer, seqname, &sequence.borrow())?;
            } else {
                // dense data track
                self.write_wiggle_fixed_step(&mut writer, seqname, &sequence.borrow())?;
            }
        }
        Ok(())
    }

}

/* -------------------------------------------------------------------------- */

impl SimpleTrack {

    fn read_wiggle_header(&mut self, scanner: &mut dyn Iterator<Item = String>) -> Result<(), String> {

        let fields: Vec<String> = fields_quoted(&scanner.next().unwrap());

        for field in fields.iter().skip(1) {

            let header_fields: Vec<&str> = field.split('=').collect();
            if header_fields.len() != 2 {
                return Err("invalid declaration line".into());
            }
            match header_fields[0] {
                "name" => self.name = remove_quotes(header_fields[1]),
                "type" => {
                    if remove_quotes(header_fields[1]) != "wiggle_0" {
                        return Err("unsupported wiggle format".into());
                    }
                }
                _ => (),
            }
        }
        Ok(())
    }

    fn read_wiggle_fixed_step(
        &mut self,
        scanner: &mut dyn Iterator<Item = String>,
    ) -> Result<(), String> {

        let fields: Vec<String> = fields_quoted(&scanner.next().unwrap());
        let mut seqname  = String::new();
        let mut position = 0;

        for field in fields.iter().skip(1) {
            let header_fields: Vec<&str> = field.split('=').collect();
            if header_fields.len() != 2 {
                return Err("invalid declaration line".into());
            }
            match header_fields[0] {
                "chrom" => seqname = remove_quotes(header_fields[1]),
                "start" => {
                    let t = header_fields[1].parse::<i64>().map_err(|_| "invalid start value")?;

                    if t <= 0 {
                        return Err("declaration line defines invalid start position".into());
                    }

                    position = self.index((t - 1) as usize);
                }
                "step" => {
                    let t = header_fields[1].parse::<i64>().map_err(|_| "invalid step value")?;
                    if self.bin_size != t as usize {
                        return Err("step sizes do not match the binSize of the track".into());
                    }
                }
                _ => (),
            }
        }

        if seqname.is_empty() {
            return Err("declaration line is missing the chromosome name".into());
        }

        let sequence = self.data.get(&seqname).cloned().unwrap_or_default().borrow().to_vec();

        for line in scanner {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() != 1 {
                break;
            }
            let value = fields[0].parse::<f64>().map_err(|_| "invalid data value")?;
            if let Some(seq) = self.data.get_mut(&seqname) {
                if position < seq.borrow().len() {
                    seq.borrow_mut()[position] = value;
                    position += 1;
                }
            }
        }

        self.data.insert(seqname, Rc::new(RefCell::new(sequence)));
        Ok(())
    }

    fn read_wiggle_variable_step(
        &mut self,
        scanner: &mut dyn Iterator<Item = String>,
    ) -> Result<(), String> {

        let fields: Vec<String> = fields_quoted(&scanner.next().unwrap());
        let mut seqname = String::new();

        for field in fields.iter().skip(1) {
            let header_fields: Vec<&str> = field.split('=').collect();
            if header_fields.len() != 2 {
                return Err("invalid declaration line".into());
            }
            match header_fields[0] {
                "chrom" => seqname = remove_quotes(header_fields[1]),
                "span" => {
                    let t = header_fields[1].parse::<i64>().map_err(|_| "invalid span value")?;
                    if self.bin_size != t as usize {
                        return Err("span does not match the binSize of the track".into());
                    }
                }
                _ => (),
            }
        }
        if seqname.is_empty() {
            return Err("declaration line is missing the chromosome name".into());
        }

        let sequence = self.data.get(&seqname).cloned().unwrap_or_default();

        for line in scanner {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() != 2 {
                break;
            }
            let position = fields[0].parse::<i64>().map_err(|_| "invalid position value")?;
            let value = fields[1].parse::<f64>().map_err(|_| "invalid data value")?;

            if position <= 0 {
                return Err("invalid chromosomal position".into());
            }
            let index = self.index((position - 1) as usize);
            if index < sequence.borrow().len() {
                sequence.borrow_mut()[index] = value;
            }
        }

        self.data.insert(seqname, sequence);
        Ok(())
    }

    pub fn read_wiggle<R: BufRead>(
        &mut self,
        reader: R
    ) -> Result<(), String> {

        let mut scanner = reader.lines().map(|l| l.unwrap());
        let mut header = false;

        while let Some(line) = scanner.next() {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.is_empty() {
                break;
            }

            match fields[0] {
                "track" => {
                    if !header {
                        header = true;
                        self.read_wiggle_header(&mut scanner)?;
                    } else {
                        return Err("file contains more than one track definition line".into());
                    }
                }
                "browser" => continue, // skip any browser options
                "fixedStep" => self.read_wiggle_fixed_step(&mut scanner)?,
                "variableStep" => self.read_wiggle_variable_step(&mut scanner)?,
                _ => return Err("unknown sequence type (i.e., not fixedStep or variableStep)".into()),
            }
        }
        Ok(())
    }

    pub fn import_wiggle(
        &mut self,
        filename: &str
    ) -> Result<(), String> {

        let file = File::open(&filename).map_err(|_| "failed to open file")?;
        let reader = BufReader::new(file);

        if is_gzip(&filename) {
            let gz_reader = flate2::read::GzDecoder::new(reader);
            self.read_wiggle(BufReader::new(gz_reader))
        } else {
            self.read_wiggle(reader)
        }
    }

}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

fn remove_quotes(s: &str) -> String {
    s.trim_matches('"').to_string()
}

fn fields_quoted(s: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut in_quotes = false;
    let mut field = String::new();

    for c in s.chars() {
        if c == '"' {
            in_quotes = !in_quotes;
        }
        if c.is_whitespace() && !in_quotes {
            if !field.is_empty() {
                fields.push(field.clone());
                field.clear();
            }
        } else {
            field.push(c);
        }
    }
    if !field.is_empty() {
        fields.push(field);
    }
    fields
}

fn is_gzip<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref()
        .extension()
        .map_or(false, |ext| ext == "gz" || ext == "gzip")
}
