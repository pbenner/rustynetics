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

/* -------------------------------------------------------------------------- */

use std::fmt;
use std::io::{self, Write};
use std::error::Error;
use std::cmp::Ordering;

use crate::range::Range;
use crate::utility::remove_duplicates_int;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
struct Meta {
    meta_name: Vec<String>,
    meta_data: Vec<Box<dyn MetaData>>,
    rows: usize,
}

impl Meta {
    fn new(names: Vec<String>, data: Vec<Box<dyn MetaData>>) -> Self {
        if names.len() != data.len() {
            panic!("NewMeta(): invalid parameters!");
        }
        let mut meta = Meta {
            meta_name: vec![],
            meta_data: vec![],
            rows: 0,
        };
        for i in 0..names.len() {
            meta.add_meta(names[i].clone(), data[i].clone());
        }
        meta
    }

    fn clone(&self) -> Self {
        let mut result = Meta::new(vec![], vec![]);
        for i in 0..self.meta_length() {
            let cloned_data = self.meta_data[i].clone_data();
            result.add_meta(self.meta_name[i].clone(), cloned_data);
        }
        result
    }

    fn length(&self) -> usize {
        self.rows
    }

    fn meta_length(&self) -> usize {
        self.meta_name.len()
    }

    fn row(&self, i: usize) -> MetaRow {
        MetaRow::new(self.clone(), i)
    }

    fn add_meta(&mut self, name: String, meta: Box<dyn MetaData>) {
        let n = match meta.as_ref().data_type() {
            DataType::String2D(v) => v.len(),
            DataType::String1D(v) => v.len(),
            DataType::Float2D(v) => v.len(),
            DataType::Float1D(v) => v.len(),
            DataType::Int2D(v) => v.len(),
            DataType::Int1D(v) => v.len(),
            DataType::Range1D(v) => v.len(),
        };
        if self.meta_length() > 0 {
            if n != self.rows {
                panic!("AddMeta(): column `{}` has invalid length: expected length of `{}` but column has length `{}`", name, self.rows, n);
            }
        } else {
            self.rows = n;
        }
        self.delete_meta(&name);
        self.meta_data.push(meta);
        self.meta_name.push(name);
    }

    fn delete_meta(&mut self, name: &str) {
        let mut i = 0;
        while i < self.meta_length() {
            if self.meta_name[i] == name {
                self.meta_name.remove(i);
                self.meta_data.remove(i);
            } else {
                i += 1;
            }
        }
    }

    fn rename_meta(&mut self, name_old: &str, name_new: &str) {
        if name_old == name_new {
            return;
        }
        self.delete_meta(name_new);
        for i in 0..self.meta_length() {
            if self.meta_name[i] == name_old {
                self.meta_name[i] = name_new.to_string();
            }
        }
    }

    fn get_meta(&self, name: &str) -> Option<&Box<dyn MetaData>> {
        for i in 0..self.meta_length() {
            if self.meta_name[i] == name {
                return Some(&self.meta_data[i]);
            }
        }
        None
    }

    fn get_meta_str(&self, name: &str) -> Result<Vec<String>, Box<dyn Error>> {
        if let Some(meta) = self.get_meta(name) {
            match meta.as_ref().data_type() {
                DataType::String2D(v) => Ok(v.clone()),
                DataType::String1D(v) => Ok(v.iter().map(|s| s.to_string()).collect()),
                _ => Err(Box::new(MetaError::InvalidMetaType)),
            }
        } else {
            Err(Box::new(MetaError::InvalidMetaName))
        }
    }

    fn get_meta_float(&self, name: &str) -> Result<Vec<f64>, Box<dyn Error>> {
        if let Some(meta) = self.get_meta(name) {
            match meta.as_ref().data_type() {
                DataType::Float2D(v) => Ok(v.clone()),
                DataType::Float1D(v) => Ok(v.clone()),
                _ => Err(Box::new(MetaError::InvalidMetaType)),
            }
        } else {
            Err(Box::new(MetaError::InvalidMetaName))
        }
    }

    fn get_meta_int(&self, name: &str) -> Result<Vec<i32>, Box<dyn Error>> {
        if let Some(meta) = self.get_meta(name) {
            match meta.as_ref().data_type() {
                DataType::Int2D(v) => Ok(v.clone()),
                DataType::Int1D(v) => Ok(v.clone()),
                _ => Err(Box::new(MetaError::InvalidMetaType)),
            }
        } else {
            Err(Box::new(MetaError::InvalidMetaName))
        }
    }

    fn get_meta_range(&self, name: &str) -> Result<Vec<Range>, Box<dyn Error>> {
        if let Some(meta) = self.get_meta(name) {
            match meta.as_ref().data_type() {
                DataType::Range1D(v) => Ok(v.clone()),
                _ => Err(Box::new(MetaError::InvalidMetaType)),
            }
        } else {
            Err(Box::new(MetaError::InvalidMetaName))
        }
    }

    fn append(&self, meta2: &Meta) -> Meta {
        let mut result = Meta::new(vec![], vec![]);
        let m1 = self.clone();
        let m2 = meta2.clone();
        for j in 0..m2.meta_length() {
            let name = &m2.meta_name[j];
            let dat1 = m1.get_meta(name).unwrap();
            let dat2 = &m2.meta_data[j];
            let t = match dat1.as_ref().data_type() {
                DataType::String2D(v) => {
                    let mut t = v.clone();
                    t.extend_from_slice(dat2.as_ref().as_string2d().unwrap());
                    DataType::String2D(t)
                }
                DataType::String1D(v) => {
                    let mut t = v.clone();
                    t.extend_from_slice(dat2.as_ref().as_string1d().unwrap());
                    DataType::String1D(t)
                }
                DataType::Float2D(v) => {
                    let mut t = v.clone();
                    t.extend_from_slice(dat2.as_ref().as_float2d().unwrap());
                    DataType::Float2D(t)
                }
                DataType::Float1D(v) => {
                    let mut t = v.clone();
                    t.extend_from_slice(dat2.as_ref().as_float1d().unwrap());
                    DataType::Float1D(t)
                }
                DataType::Int2D(v) => {
                    let mut t = v.clone();
                    t.extend_from_slice(dat2.as_ref().as_int2d().unwrap());
                    DataType::Int2D(t)
                }
                DataType::Int1D(v) => {
                    let mut t = v.clone();
                    t.extend_from_slice(dat2.as_ref().as_int1d().unwrap());
                    DataType::Int1D(t)
                }
                DataType::Range1D(v) => {
                    let mut t = v.clone();
                    t.extend_from_slice(dat2.as_ref().as_range1d().unwrap());
                    DataType::Range1D(t)
                }
            };
            result.add_meta(name.clone(), Box::new(t));
        }
        result
    }

    fn remove(&self, indices: &[usize]) -> Meta {
        if indices.is_empty() {
            return self.clone();
        }
        let indices = remove_duplicates_int(indices);
        let mut indices = indices.to_vec();
        indices.sort();
        let n = self.length();
        let m = n - indices.len();
        let mut idx = vec![0; m];
        for (i, j, k) in (0..n).zip(0..m).zip(0..indices.len()) {
            while k < indices.len() - 1 && i > indices[k] {
                k += 1;
            }
            if i != indices[k] {
                idx[j] = i;
                j += 1;
            }
        }
        self.subset(&idx)
    }

    fn subset(&self, indices: &[usize]) -> Meta {
        let n = indices.len();
        let m = self.meta_length();
        let mut data = vec![];
        for j in 0..m {
            let meta_data = &self.meta_data[j];
            let t = match meta_data.as_ref().data_type() {
                DataType::String2D(v) => {
                    let mut l = vec![];
                    for i in 0..n {
                        l.push(v[indices[i]].clone());
                    }
                    DataType::String2D(l)
                }
                DataType::String1D(v) => {
                    let mut l = vec![];
                    for i in 0..n {
                        l.push(v[indices[i]].clone());
                    }
                    DataType::String1D(l)
                }
                DataType::Float2D(v) => {
                    let mut l = vec![];
                    for i in 0..n {
                        l.push(v[indices[i]].clone());
                    }
                    DataType::Float2D(l)
                }
                DataType::Float1D(v) => {
                    let mut l = vec![];
                    for i in 0..n {
                        l.push(v[indices[i]].clone());
                    }
                    DataType::Float1D(l)
                }
                DataType::Int2D(v) => {
                    let mut l = vec![];
                    for i in 0..n {
                        l.push(v[indices[i]].clone());
                    }
                    DataType::Int2D(l)
                }
                DataType::Int1D(v) => {
                    let mut l = vec![];
                    for i in 0..n {
                        l.push(v[indices[i]].clone());
                    }
                    DataType::Int1D(l)
                }
                DataType::Range1D(v) => {
                    let mut l = vec![];
                    for i in 0..n {
                        l.push(v[indices[i]].clone());
                    }
                    DataType::Range1D(l)
                }
            };
            data.push(Box::new(t) as Box<dyn MetaData>);
        }
        Meta::new(self.meta_name.clone(), data)
    }

    fn slice(&self, ifrom: usize, ito: usize) -> Meta {
        let n = ito - ifrom;
        let m = self.meta_length();
        let mut data = vec![];
        for j in 0..m {
            let meta_data = &self.meta_data[j];
            let t = match meta_data.as_ref().data_type() {
                DataType::String2D(v) => {
                    let mut l = vec![];
                    for i in ifrom..ito {
                        l.push(v[i].clone());
                    }
                    DataType::String2D(l)
                }
                DataType::String1D(v) => {
                    let mut l = vec![];
                    for i in ifrom..ito {
                        l.push(v[i].clone());
                    }
                    DataType::String1D(l)
                }
                DataType::Float2D(v) => {
                    let mut l = vec![];
                    for i in ifrom..ito {
                        l.push(v[i].clone());
                    }
                    DataType::Float2D(l)
                }
                DataType::Float1D(v) => {
                    let mut l = vec![];
                    for i in ifrom..ito {
                        l.push(v[i].clone());
                    }
                    DataType::Float1D(l)
                }
                DataType::Int2D(v) => {
                    let mut l = vec![];
                    for i in ifrom..ito {
                        l.push(v[i].clone());
                    }
                    DataType::Int2D(l)
                }
                DataType::Int1D(v) => {
                    let mut l = vec![];
                    for i in ifrom..ito {
                        l.push(v[i].clone());
                    }
                    DataType::Int1D(l)
                }
                DataType::Range1D(v) => {
                    let mut l = vec![];
                    for i in ifrom..ito {
                        l.push(v[i].clone());
                    }
                    DataType::Range1D(l)
                }
            };
            data.push(Box::new(t) as Box<dyn MetaData>);
        }
        Meta::new(self.meta_name.clone(), data)
    }

    fn merge(&self, indices: &[usize]) -> Meta {
        let slice_max = |s: &[usize]| -> usize {
            let mut max = 0;
            for &v in s {
                if v > max {
                    max = v;
                }
            }
            max
        };
        let n = slice_max(indices) + 1;
        let m = self.meta_length();
        let mut data = vec![];
        for j in 0..m {
            let meta_data = &self.meta_data[j];
            let t = match meta_data.as_ref().data_type() {
                DataType::String1D(v) => {
                    let mut l = vec![vec![]; n];
                    for i in 0..v.len() {
                        l[indices[i]].push(v[i].clone());
                    }
                    DataType::String2D(l)
                }
                DataType::Float1D(v) => {
                    let mut l = vec![vec![]; n];
                    for i in 0..v.len() {
                        l[indices[i]].push(v[i].clone());
                    }
                    DataType::Float2D(l)
                }
                DataType::Int1D(v) => {
                    let mut l = vec![vec![]; n];
                    for i in 0..v.len() {
                        l[indices[i]].push(v[i].clone());
                    }
                    DataType::Int2D(l)
                }
                _ => panic!("cannot merge {:?}", meta_data.as_ref().data_type()),
            };
            data.push(Box::new(t) as Box<dyn MetaData>);
        }
        Meta::new(self.meta_name.clone(), data)
    }

    fn reduce_string<F>(&mut self, name: &str, name_new: &str, f: F)
    where
        F: Fn(&[String]) -> String,
    {
        if let Some(meta) = self.get_meta(name) {
            match meta.as_ref().data_type() {
                DataType::String2D(v) => {
                    let r = v.iter().map(|s| f(s)).collect();
                    if !name_new.is_empty() && name != name_new {
                        self.add_meta(name_new.to_string(), Box::new(DataType::String1D(r)));
                    } else {
                        self.delete_meta(name);
                        self.add_meta(name.to_string(), Box::new(DataType::String1D(r)));
                    }
                }
                DataType::String1D(v) => {
                    let r = f(v);
                    if !name_new.is_empty() && name != name_new {
                        self.add_meta(name_new.to_string(), Box::new(DataType::String1D(vec![r])));
                    } else {
                        self.delete_meta(name);
                        self.add_meta(name.to_string(), Box::new(DataType::String1D(vec![r])));
                    }
                }
                _ => panic!("invalid meta data"),
            }
        } else {
            panic!("invalid meta name");
        }
    }

    fn reduce_float<F>(&mut self, name: &str, name_new: &str, f: F)
    where
        F: Fn(&[f64]) -> f64,
    {
        if let Some(meta) = self.get_meta(name) {
            match meta.as_ref().data_type() {
                DataType::Float2D(v) => {
                    let r = v.iter().map(|s| f(s)).collect();
                    if !name_new.is_empty() && name != name_new {
                        self.add_meta(name_new.to_string(), Box::new(DataType::Float1D(r)));
                    } else {
                        self.delete_meta(name);
                        self.add_meta(name.to_string(), Box::new(DataType::Float1D(r)));
                    }
                }
                DataType::Float1D(v) => {
                    let r = f(v);
                    if !name_new.is_empty() && name != name_new {
                        self.add_meta(name_new.to_string(), Box::new(DataType::Float1D(vec![r])));
                    } else {
                        self.delete_meta(name);
                        self.add_meta(name.to_string(), Box::new(DataType::Float1D(vec![r])));
                    }
                }
                _ => panic!("invalid meta data"),
            }
        } else {
            panic!("invalid meta name");
        }
    }

    fn reduce_int<F>(&mut self, name: &str, name_new: &str, f: F)
    where
        F: Fn(&[i32]) -> i32,
    {
        if let Some(meta) = self.get_meta(name) {
            match meta.as_ref().data_type() {
                DataType::Int2D(v) => {
                    let r = v.iter().map(|s| f(s)).collect();
                    if !name_new.is_empty() && name != name_new {
                        self.add_meta(name_new.to_string(), Box::new(DataType::Int1D(r)));
                    } else {
                        self.delete_meta(name);
                        self.add_meta(name.to_string(), Box::new(DataType::Int1D(r)));
                    }
                }
                DataType::Int1D(v) => {
                    let r = f(v);
                    if !name_new.is_empty() && name != name_new {
                        self.add_meta(name_new.to_string(), Box::new(DataType::Int1D(vec![r])));
                    } else {
                        self.delete_meta(name);
                        self.add_meta(name.to_string(), Box::new(DataType::Int1D(vec![r])));
                    }
                }
                _ => panic!("invalid meta data"),
            }
        } else {
            panic!("invalid meta name");
        }
    }

    fn sorted_indices(&self, name: &str, reverse: bool) -> Result<Vec<usize>, Box<dyn Error>> {
        let mut l = vec![];
        if let Some(meta) = self.get_meta(name) {
            match meta.as_ref().data_type() {
                DataType::Float1D(v) => {
                    for i in 0..v.len() {
                        l.push((i, v[i]));
                    }
                }
                DataType::Int1D(v) => {
                    for i in 0..v.len() {
                        l.push((i, v[i] as f64));
                    }
                }
                DataType::String1D(v) => {
                    for i in 0..v.len() {
                        l.push((i, i as f64));
                    }
                }
                _ => return Err(Box::new(MetaError::InvalidMetaType)),
            }
        } else {
            return Err(Box::new(MetaError::InvalidMetaName));
        }
        if reverse {
            l.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
        } else {
            l.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));
        }
        let j: Vec<usize> = l.iter().map(|(i, _)| *i).collect();
        Ok(j)
    }

    fn sort(&self, name: &str, reverse: bool) -> Result<Meta, Box<dyn Error>> {
        let j = self.sorted_indices(name, reverse)?;
        Ok(self.subset(&j))
    }

    fn write_pretty<W: Write>(
        &self,
        writer: &mut W,
        n: usize,
        args: &[&dyn std::any::Any],
    ) -> Result<(), Box<dyn Error>> {
        let mut use_scientific = false;
        for arg in args {
            if let Some(a) = arg.downcast_ref::<OptionPrintScientific>() {
                use_scientific = a.value;
            }
        }
        let print_cell_slice = |writer: &mut W, widths: &[usize], i: usize, j: usize, data: &dyn MetaData| -> Result<usize, Box<dyn Error>> {
            let mut tmp_buffer = vec![];
            let mut tmp_writer = io::BufWriter::new(&mut tmp_buffer);
            match data.data_type() {
                DataType::String2D(v) => {
                    for k in 0..v[i].len() {
                        write!(tmp_writer, " {}", v[i][k])?;
                    }
                }
                DataType::Float2D(v) => {
                    if use_scientific {
                        for k in 0..v[i].len() {
                            write!(tmp_writer, " {:e}", v[i][k])?;
                        }
                    } else {
                        for k in 0..v[i].len() {
                            write!(tmp_writer, " {:.6}", v[i][k])?;
                        }
                    }
                }
                DataType::Int2D(v) => {
                    for k in 0..v[i].len() {
                        write!(tmp_writer, " {}", v[i][k])?;
                    }
                }
                _ => panic!("invalid meta data"),
            }
            tmp_writer.flush()?;
            let l = write!(writer, "{:width$}", String::from_utf8_lossy(&tmp_buffer), width = widths[j] - 1)?;
            Ok(l)
        };
        let print_cell = |writer: &mut W, widths: &[usize], i: usize, j: usize| -> Result<usize, Box<dyn Error>> {
            match self.meta_data[j].as_ref().data_type() {
                DataType::String1D(v) => {
                    write!(writer, " {:width$}",  v[i], width = widths[j] - 1)?;
                }
                DataType::Float1D(v) => {
                    if use_scientific {
                        write!(writer, " {:width$e}", v[i], width = widths[j] - 1)?;
                    } else {
                        write!(writer, " {:width$.6}", v[i], width = widths[j] - 1)?;
                    }
                }
                DataType::Int1D(v) => {
                    write!(writer, " {:width$}", v[i], width = widths[j] - 1)?;
                }
                DataType::Range1D(v) => {
                    write!(writer, " [{}, {})", v[i].from, v[i].to)?;
                }
                _ => print_cell_slice(writer, widths, i, j, self.meta_data[j].as_ref())?,
            }
            Ok(0)
        };
        let print_row = |writer: &mut W, widths: &[usize], i: usize| -> Result<(), Box<dyn Error>> {
            if i != 0 {
                writeln!(writer)?;
            }
            for j in 0..self.meta_length() {
                print_cell(writer, widths, i, j)?;
            }
            Ok(())
        };
        let update_max_widths = |i: usize, widths: &mut [usize]| -> Result<(), Box<dyn Error>> {
            for j in 0..self.meta_length() {
                let width = print_cell(&mut io::sink(), widths, i, j)?;
                if width > widths[j] {
                    widths[j] = width;
                }
            }
            Ok(())
        };
        let print_header = |writer: &mut W, widths: &[usize]| -> Result<(), Box<dyn Error>> {
            for j in 0..self.meta_length() {
                write!(writer, " {:width$}", self.meta_name[j], width = widths[j] - 1)?;
            }
            writeln!(writer)?;
            Ok(())
        };
        let apply_rows = |f1: &mut dyn FnMut(usize) -> Result<(), Box<dyn Error>>, f2: &mut dyn FnMut() -> Result<(), Box<dyn Error>>| {
            if self.length() <= n + 1 {
                for i in 0..self.length() {
                    f1(i)?;
                }
            } else {
                for i in 0..n / 2 {
                    f1(i)?;
                }
                f2()?;
                for i in self.length() - n / 2..self.length() {
                    f1(i)?;
                }
            }
            Ok(())
        };
        let mut widths = vec![0; self.meta_length()];
        for j in 0..self.meta_length() {
            let width = write!(io::sink(), " {}", self.meta_name[j])?;
            widths[j] = width;
        }
        apply_rows(
            &mut |i| -> Result<(), Box<dyn Error>> { update_max_widths(i, &mut widths) },
            &mut || -> Result<(), Box<dyn Error>> {
                writeln!(writer)?;
                for j in 0..self.meta_length() {
                    write!(writer, " {:width$}", "...", width = widths[j] - 1)?;
                }
                Ok(())
            },
        )?;
        print_header(writer, &widths)?;
        apply_rows(
            &mut |i| -> Result<(), Box<dyn Error>> { print_row(writer, &widths, i) },
            &mut || -> Result<(), Box<dyn Error>> {
                writeln!(writer)?;
                for j in 0..self.meta_length() {
                    write!(writer, " {:width$}", "...", width = widths[j] - 1)?;
                }
                Ok(())
            },
        )?;
        Ok(())
    }

    fn print_pretty(&self, n: usize, args: &[&dyn std::any::Any]) -> String {
        let mut buffer = vec![];
        let mut writer = io::BufWriter::new(&mut buffer);
        if let Err(_) = self.write_pretty(&mut writer, n, args) {
            return String::new();
        }
        writer.flush().unwrap();
        String::from_utf8_lossy(&buffer).to_string()
    }
}

impl fmt::Display for Meta {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut buffer = vec![];
        let mut writer = io::BufWriter::new(&mut buffer);
        if let Err(_) = self.write_pretty(&mut writer, 10, &OptionPrintScientific { value: false }) {
            return Ok(());
        }
        writer.flush().unwrap();
        write!(f, "{}", String::from_utf8_lossy(&buffer))
    }
}

#[derive(Clone, Debug)]
struct MetaRow {
    meta: Meta,
    idx: usize,
}

impl MetaRow {
    fn new(meta: Meta, idx: usize) -> Self {
        MetaRow { meta, idx }
    }

    fn get_meta(&self, name: &str) -> Option<&dyn MetaData> {
        if let Some(meta) = self.meta.get_meta(name) {
            match meta.as_ref().data_type() {
                DataType::String2D(v) => Some(v[self.idx].as_ref()),
                DataType::String1D(v) => Some(v[self.idx].as_ref()),
                DataType::Float2D(v) => Some(v[self.idx].as_ref()),
                DataType::Float1D(v) => Some(v[self.idx].as_ref()),
                DataType::Int2D(v) => Some(v[self.idx].as_ref()),
                DataType::Int1D(v) => Some(v[self.idx].as_ref()),
                _ => panic!("Row(): invalid type!"),
            }
        } else {
            None
        }
    }

    fn get_meta_str(&self, name: &str) -> Result<String, Box<dyn Error>> {
        if let Some(meta) = self.get_meta(name) {
            match meta.data_type() {
                DataType::String2D(v) => Ok(v[self.idx].clone()),
                DataType::String1D(v) => Ok(v[self.idx].clone()),
                _ => Err(Box::new(MetaError::InvalidMetaType)),
            }
        } else {
            Err(Box::new(MetaError::InvalidMetaName))
        }
    }

    fn get_meta_float(&self, name: &str) -> Result<f64, Box<dyn Error>> {
        if let Some(meta) = self.get_meta(name) {
            match meta.data_type() {
                DataType::Float2D(v) => Ok(v[self.idx]),
                DataType::Float1D(v) => Ok(v[self.idx]),
                _ => Err(Box::new(MetaError::InvalidMetaType)),
            }
        } else {
            Err(Box::new(MetaError::InvalidMetaName))
        }
    }
}

trait MetaData: fmt::Debug {
    fn data_type(&self) -> DataType;
    fn clone_data(&self) -> Box<dyn MetaData>;
    fn as_string2d(&self) -> Option<&[Vec<String>]>;
    fn as_string1d(&self) -> Option<&[String]>;
    fn as_float2d(&self) -> Option<&[Vec<f64>]>;
    fn as_float1d(&self) -> Option<&[f64]>;
    fn as_int2d(&self) -> Option<&[Vec<i32>]>;
    fn as_int1d(&self) -> Option<&[i32]>;
    fn as_range1d(&self) -> Option<&[Range]>;
}

#[derive(Debug)]
enum DataType {
    String2D(Vec<Vec<String>>),
    String1D(Vec<String>),
    Float2D(Vec<Vec<f64>>),
    Float1D(Vec<f64>),
    Int2D(Vec<Vec<i32>>),
    Int1D(Vec<i32>),
    Range1D(Vec<Range>),
}

impl MetaData for DataType {
    fn data_type(&self) -> DataType {
        self.clone()
    }

    fn clone_data(&self) -> Box<dyn MetaData> {
        match self {
            DataType::String2D(v) => Box::new(DataType::String2D(v.clone())),
            DataType::String1D(v) => Box::new(DataType::String1D(v.clone())),
            DataType::Float2D(v) => Box::new(DataType::Float2D(v.clone())),
            DataType::Float1D(v) => Box::new(DataType::Float1D(v.clone())),
            DataType::Int2D(v) => Box::new(DataType::Int2D(v.clone())),
            DataType::Int1D(v) => Box::new(DataType::Int1D(v.clone())),
            DataType::Range1D(v) => Box::new(DataType::Range1D(v.clone())),
        }
    }

    fn as_string2d(&self) -> Option<&[Vec<String>]> {
        match self {
            DataType::String2D(v) => Some(v),
            _ => None,
        }
    }

    fn as_string1d(&self) -> Option<&[String]> {
        match self {
            DataType::String1D(v) => Some(v),
            _ => None,
        }
    }

    fn as_float2d(&self) -> Option<&[Vec<f64>]> {
        match self {
            DataType::Float2D(v) => Some(v),
            _ => None,
        }
    }

    fn as_float1d(&self) -> Option<&[f64]> {
        match self {
            DataType::Float1D(v) => Some(v),
            _ => None,
        }
    }

    fn as_int2d(&self) -> Option<&[Vec<i32>]> {
        match self {
            DataType::Int2D(v) => Some(v),
            _ => None,
        }
    }

    fn as_int1d(&self) -> Option<&[i32]> {
        match self {
            DataType::Int1D(v) => Some(v),
            _ => None,
        }
    }

    fn as_range1d(&self) -> Option<&[Range]> {
        match self {
            DataType::Range1D(v) => Some(v),
            _ => None,
        }
    }
}

#[derive(Debug)]
enum MetaError {
    InvalidMetaType,
    InvalidMetaName,
}

impl fmt::Display for MetaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MetaError::InvalidMetaType => write!(f, "Invalid meta type"),
            MetaError::InvalidMetaName => write!(f, "Invalid meta name"),
        }
    }
}

impl Error for MetaError {}

#[derive(Debug)]
struct OptionPrintScientific {
    value: bool,
}
