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
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::str::FromStr;

use crate::meta::Meta;
use crate::meta::MetaData;

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct OptionPrintScientific {
    pub value: bool,
}

/* -------------------------------------------------------------------------- */

impl Meta {

    fn print_meta_cell_slice<W: Write>(&self, writer: &mut W, widths: &[usize], i: usize, j: usize, data: &MetaData, use_scientific: bool) -> io::Result<()> {
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

    fn print_meta_cell<W: Write>(&self, writer: &mut W, widths: &[usize], i: usize, j: usize, use_scientific: bool) -> io::Result<()> {
        match &self.meta_data[j] {
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
                write!(writer, " {:width$}", v[i], width = widths[j] - 1)
            }
            _ => self.print_meta_cell_slice(writer, widths, i, j, &self.meta_data[j], use_scientific),
        }
    }

    fn print_meta_row<W: Write>(&self, writer: &mut W, widths: &[usize], i: usize, use_scientific: bool) -> io::Result<()> {
        if i != 0 {
            writeln!(writer)?;
        }
        for j in 0..self.num_cols() {
            self.print_meta_cell(writer, widths, i, j, use_scientific)?;
        }
        Ok(())
    }

    fn meta_update_max_widths(&self, i: usize, widths: &mut [usize], use_scientific: bool) -> io::Result<()> {
        for j in 0..self.num_cols() {
            let mut tmp_buffer = Vec::new();
            {
                let mut tmp_writer = BufWriter::new(&mut tmp_buffer);
                self.print_meta_cell(&mut tmp_writer, widths, i, j, use_scientific)?;
                tmp_writer.flush()?;
            }
            let width = tmp_buffer.len();
            if width > widths[j] {
                widths[j] = width;
            }
        }
        Ok(())
    }

    fn print_meta_header<W: Write>(&self, writer: &mut W, widths: &[usize]) -> io::Result<()> {
        for j in 0..self.num_cols() {
            write!(writer, " {:width$}", self.meta_name[j], width = widths[j] - 1)?;
        }
        writeln!(writer)
    }

    pub fn write_table<W: Write>(&self, writer: &mut W, args: &[&dyn Any]) -> io::Result<()> {
        let mut use_scientific = false;
        for arg in args {
            if let Some(option) = arg.downcast_ref::<OptionPrintScientific>() {
                use_scientific = option.value;
            }
        }

        let mut widths = vec![0; self.num_cols()];
        for j in 0..self.num_cols() {
            widths[j] = self.meta_name[j].len()+1;
        }
        (0..self.num_cols()).map(|j| widths[j] = self.meta_name[j].len());

        for i in 0..self.num_rows() {
            self.meta_update_max_widths(i, &mut widths, use_scientific)?;
        }
        self.print_meta_header(writer, &widths)?;

        for i in 0..self.num_rows() {
            self.print_meta_row(writer, &widths, i, use_scientific)?;
        }
        Ok(())
    }

    pub fn print_table(&self, args: &[&dyn Any]) -> String {
        let mut buffer = Vec::new();
        {
            let mut writer = io::Cursor::new(&mut buffer);
            self.write_table(&mut writer, args).unwrap();
        }
        String::from_utf8(buffer).unwrap()
    }

    pub fn read_table<R: Read>(&mut self, reader: R, names: &[&str], types: &[&str]) -> io::Result<()> {
        let mut reader = BufReader::new(reader);

        if names.len() != types.len() {
            panic!("invalid arguments");
        }
        let mut meta_reader = MetaTableReader::new(names, types);

        // Parse header
        let mut line = String::new();
        reader.read_line(&mut line)?;

        meta_reader.read_header(&line);

        let mut i = 2;
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            if line.is_empty() {
                continue;
            }
            meta_reader.read_line(&line, i);

            i += 1;
        }

        meta_reader.push(self);

        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

struct MetaTableReader<'a> {
    idx_map  : std::collections::HashMap<&'a str, i32>,
    meta_map : std::collections::HashMap<&'a str, MetaData>
}

/* -------------------------------------------------------------------------- */

impl<'a> MetaTableReader<'a> {

    pub fn new(names: &[&'a str], types: &[&'a str]) -> Self {
        let mut idx_map  = std::collections::HashMap::new();
        let mut meta_map = std::collections::HashMap::new();

        for i in 0..names.len() {
            idx_map.insert(names[i], -1);
            match types[i] {
                "String"      => meta_map.insert(names[i], MetaData::StringArray(Vec::new())),
                "Int"         => meta_map.insert(names[i], MetaData::IntArray(Vec::new())),
                "Float"       => meta_map.insert(names[i], MetaData::FloatArray(Vec::new())),
                "Vec<String>" => meta_map.insert(names[i], MetaData::StringMatrix(Vec::new())),
                "Vec<Int>"    => meta_map.insert(names[i], MetaData::IntMatrix(Vec::new())),
                "Vec<Float>"  => meta_map.insert(names[i], MetaData::FloatMatrix(Vec::new())),
                _ => panic!("invalid types argument"),
            };
        }

        MetaTableReader{
            idx_map : idx_map,
            meta_map: meta_map,
        }
    }

    pub fn read_header(&mut self, line: &String) {

        let fields: Vec<&str> = line.split_whitespace().collect();

        for (i, field) in fields.iter().enumerate() {
            if let Some(idx) = self.idx_map.get_mut(*field) {
                *idx = i as i32;
            }
        }
    }

    pub fn read_line(&mut self, line: &String, i: i32) -> io::Result<()> {

        let fields: Vec<&str> = line.split_whitespace().collect();
        for (name, idx) in &self.idx_map {
            if *idx == -1 {
                continue;
            }
            if *idx >= fields.len() as i32 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid table"));
            }
            match self.meta_map.get_mut(name).unwrap() {
                MetaData::StringArray(v) => {
                    v.push(fields[*idx as usize].to_string());
                }
                MetaData::IntArray(v) => {
                    let value = i64::from_str(fields[*idx as usize]).map_err(|e| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("parsing meta information failed at line `{}`: {}", i, e),
                        )
                    })?;
                    v.push(value);
                }
                MetaData::FloatArray(v) => {
                    if fields[*idx as usize] == "NA" || fields[*idx as usize] == "NaN" {
                        v.push(f64::NAN);
                    } else {
                        let value = f64::from_str(fields[*idx as usize]).map_err(|e| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("parsing meta information failed at line `{}`: {}", i, e),
                            )
                        })?;
                        v.push(value);
                    }
                }
                MetaData::IntMatrix(v) => {
                    let data: Vec<&str> = fields[*idx as usize].split(',').collect();
                    if data.len() == 1 && data[0] == "nil" {
                        v.push(Vec::<i64>::new());
                    } else {
                        let mut entry = Vec::with_capacity(data.len());
                        for d in data {
                            let value = i64::from_str(d).map_err(|e| {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("parsing meta information failed at line `{}`: {}", i, e),
                                )
                            })?;
                            entry.push(value);
                        }
                        v.push(entry);
                    }
                }
                MetaData::FloatMatrix(v) => {
                    let data: Vec<&str> = fields[*idx as usize].split(',').collect();
                    if data.len() == 1 && data[0] == "nil" {
                        v.push(Vec::<f64>::new());
                    } else {
                        let mut entry = Vec::with_capacity(data.len());
                        for d in data {
                            let value = f64::from_str(d).map_err(|e| {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("parsing meta information failed at line `{}`: {}", i, e),
                                )
                            })?;
                            entry.push(value);
                        }
                        v.push(entry);
                    }
                }
                MetaData::StringMatrix(v) => {
                    let data: Vec<&str> = fields[*idx as usize].split(',').collect();
                    if data.len() == 1 && data[0] == "nil" {
                        v.push(Vec::<String>::new());
                    } else {
                        let mut entry = Vec::with_capacity(data.len());
                        for d in data {
                            entry.push(d.to_string());
                        }
                        v.push(entry);
                    }
                }
                _ => unreachable!(),
            }
        }
        Ok(())
    }

    pub fn push(&mut self, meta: &mut Meta) {
        for (name, idx) in &self.idx_map {
            if *idx != -1 {
                if let Err(_) = meta.add_meta(&name, self.meta_map.remove(name).unwrap()) {
                    panic!("internal error")
                }
            }
        }
    }
}
