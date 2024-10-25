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

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::{rc::Rc, cell::RefCell};

use crate::track_simple::SimpleTrack;
use crate::utility::is_gzip;

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
