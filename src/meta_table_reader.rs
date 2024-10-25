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

use std::io;
use std::str::FromStr;

use crate::meta::{Meta, MetaData};

/* -------------------------------------------------------------------------- */

pub struct MetaTableReader<'a> {
    idx_map  : std::collections::HashMap<&'a str, i32>,
    meta_map : std::collections::HashMap<&'a str, MetaData>
}

/* -------------------------------------------------------------------------- */

impl<'a> MetaTableReader<'a> {

    pub fn new(names: &[&'a str], types: &[&'a str]) -> Self {
        if names.len() != types.len() {
            panic!("invalid arguments");
        }

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

    pub fn read_header(&mut self, line: &String) -> io::Result<()> {

        let fields: Vec<&str> = line.split_whitespace().collect();

        for (i, field) in fields.iter().enumerate() {
            if let Some(idx) = self.idx_map.get_mut(*field) {
                *idx = i as i32;
            }
        }
        Ok(())
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
                if let Err(_) = meta.add(&name, self.meta_map.remove(name).unwrap()) {
                    panic!("internal error")
                }
            }
        }
    }
}
