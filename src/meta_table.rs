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

use std::fmt;
use std::any::Any;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::str::FromStr;

use crate::meta::Meta;

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
struct OptionPrintScientific {
    value: bool,
}

/* -------------------------------------------------------------------------- */

impl Meta {

    fn print_meta_cell_slice<W: Write>(&self, writer: &mut W, widths: &[usize], i: usize, j: usize, data: &dyn fmt::Debug, use_scientific: bool) -> io::Result<()> {
        let mut tmp_buffer = Vec::new();
        {
            let mut tmp_writer = io::Cursor::new(&mut tmp_buffer);
            match data {
                v @ &Vec::<String>::new() => {
                    if v.is_empty() {
                        write!(tmp_writer, "nil")?;
                    }
                    for (k, s) in v.iter().enumerate() {
                        if k != 0 {
                            write!(tmp_writer, ",")?;
                        }
                        write!(tmp_writer, "{}", s)?;
                    }
                }
                v @ &Vec::<f64>::new() => {
                    if v.is_empty() {
                        write!(tmp_writer, "nil")?;
                    }
                    for (k, f) in v.iter().enumerate() {
                        if k != 0 {
                            write!(tmp_writer, ",")?;
                        }
                        if use_scientific {
                            write!(tmp_writer, "{:e}", f)?;
                        } else {
                            write!(tmp_writer, "{:?}", f)?;
                        }
                    }
                }
                v @ &Vec::<i32>::new() => {
                    if v.is_empty() {
                        write!(tmp_writer, "nil")?;
                    }
                    for (k, i) in v.iter().enumerate() {
                        if k != 0 {
                            write!(tmp_writer, ",")?;
                        }
                        write!(tmp_writer, "{}", i)?;
                    }
                }
                _ => panic!("invalid meta data"),
            }
        }
        write!(writer, " {:width$}s", String::from_utf8(tmp_buffer).unwrap(), width = widths[j] - 1)
    }

    fn print_cell<W: Write>(&self, writer: &mut W, widths: &[usize], i: usize, j: usize, use_scientific: bool) -> io::Result<()> {
        match &self.meta_data[j] {
            v @ &Vec::<String>::new() => {
                let format = format!(" %{}s", widths[j] - 1);
                write!(writer, " %{:width$}s", v[i], width = widths[j] - 1)
            }
            v @ &Vec::<f64>::new() => {
                if use_scientific {
                    write!(writer, " %{:width$}e", v[i], width = widths[j] - 1)
                } else {
                    write!(writer, " %{:width$}f", v[i], width = widths[j] - 1)
                }
            }
            v @ &Vec::<i32>::new() => {
                write!(writer, " %{:width$}d", v[i], width = widths[j] - 1)
            }
            _ => print_cell_slice(writer, widths, i, j, &self.meta_data[j]),
        }
    }

    fn print_row<W: Write>(&self, writer: &mut W, widths: &[usize], i: usize) -> io::Result<()> {
        if i != 0 {
            writeln!(writer)?;
        }
        for j in 0..self.meta_length() {
            print_cell(writer, widths, i, j)?;
        }
        Ok(())
    }

    fn update_max_widths<W: Write>(&self, i: usize, widths: &mut [usize]) -> io::Result<()> {
        for j in 0..self.meta_length() {
            let width = print_cell(&mut io::sink(), widths, i, j)?;
            if width > widths[j] {
                widths[j] = width;
            }
        }
        Ok(())
    }

    fn print_header<W: Write>(&self, writer: &mut W, widths: &[usize]) -> io::Result<()> {
        for j in 0..self.meta_length() {
            write!(writer, " %{:width$}s", self.meta_name[j], width = widths[j] - 1)?;
        }
        writeln!(writer)
    }

    fn apply_rows(&self, f1: &mut dyn FnMut(usize) -> io::Result<()>) -> io::Result<()> {
        for i in 0..self.length() {
            f1(i)?;
        }
        Ok(())
    }

    fn write_table<W: Write>(&self, writer: &mut W, header: bool, args: &[&dyn Any]) -> io::Result<()> {
        let mut use_scientific = false;
        for arg in args {
            if let Some(option) = arg.downcast_ref::<OptionPrintScientific>() {
                use_scientific = option.value;
            }
        }

        let mut widths = vec![0; self.meta_length()];
        for j in 0..self.meta_length() {
            let width = write!(io::sink(), " {}", self.meta_name[j])?;
            widths[j] = width;
        }

        apply_rows(&mut |i| update_max_widths(i, &mut widths))?;

        if header {
            print_header(writer, &widths)?;
        }

        apply_rows(&mut |i| print_row(writer, &widths, i))
    }

    fn print_table(&self, header: bool, args: &[&dyn fmt::Debug]) -> String {
        let mut buffer = Vec::new();
        {
            let mut writer = io::Cursor::new(&mut buffer);
            self.write_table(&mut writer, header, args).unwrap();
        }
        String::from_utf8(buffer).unwrap()
    }

    fn read_table<R: Read>(&mut self, reader: R, names: &[String], types: &[String]) -> io::Result<()> {
        let mut reader = BufReader::new(reader);

        if names.len() != types.len() {
            panic!("invalid arguments");
        }

        let mut idx_map = std::collections::HashMap::new();
        let mut meta_map = std::collections::HashMap::new();
        for i in 0..names.len() {
            idx_map.insert(names[i].clone(), -1);
            match types[i].as_str() {
                "Vec<String>" => meta_map.insert(names[i].clone(), Vec::<String>::new()),
                "Vec<i32>" => meta_map.insert(names[i].clone(), Vec::<i32>::new()),
                "Vec<f64>" => meta_map.insert(names[i].clone(), Vec::<f64>::new()),
                "Vec<Vec<String>>" => meta_map.insert(names[i].clone(), Vec::<Vec<String>>::new()),
                "Vec<Vec<i32>>" => meta_map.insert(names[i].clone(), Vec::<Vec<i32>>::new()),
                "Vec<Vec<f64>>" => meta_map.insert(names[i].clone(), Vec::<Vec<f64>>::new()),
                _ => panic!("invalid types argument"),
            };
        }

        let mut line = String::new();
        reader.read_line(&mut line)?;
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 4 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid table"));
        }
        for (i, field) in fields.iter().enumerate() {
            if let Some(idx) = idx_map.get_mut(*field) {
                *idx = i as i32;
            }
        }

        let mut i = 2;
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            if line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split_whitespace().collect();
            for (name, idx) in &idx_map {
                if *idx == -1 {
                    continue;
                }
                if *idx >= fields.len() as i32 {
                    return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid table"));
                }
                match meta_map.get_mut(name).unwrap() {
                    v @ &mut Vec::<String>::new() => {
                        v.push(fields[*idx as usize].to_string());
                    }
                    v @ &mut Vec::<i32>::new() => {
                        let value = i32::from_str(fields[*idx as usize]).map_err(|e| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("parsing meta information failed at line `{}`: {}", i, e),
                            )
                        })?;
                        v.push(value);
                    }
                    v @ &mut Vec::<f64>::new() => {
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
                    v @ &mut Vec::<Vec<i32>>::new() => {
                        let data: Vec<&str> = fields[*idx as usize].split(',').collect();
                        if data.len() == 1 && data[0] == "nil" {
                            v.push(Vec::<i32>::new());
                        } else {
                            let mut entry = Vec::with_capacity(data.len());
                            for d in data {
                                let value = i32::from_str(d).map_err(|e| {
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
                    v @ &mut Vec::<Vec<f64>>::new() => {
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
                    v @ &mut Vec::<Vec<String>>::new() => {
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
            i += 1;
        }

        for (name, idx) in idx_map {
            if idx != -1 {
                self.add_meta(name, meta_map.remove(&name).unwrap());
            }
        }

        Ok(())
    }
}
