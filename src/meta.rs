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
use std::collections::HashMap;
use std::error::Error;
use std::io::{self, Write};

use crate::range::Range;
use crate::utility::remove_duplicates_int;

/* -------------------------------------------------------------------------- */

#[derive(Debug, Clone)]
enum MetaData {
    StringMatrix(Vec<Vec<String>>),
    StringArray(Vec<String>),
    FloatMatrix(Vec<Vec<f64>>),
    FloatArray(Vec<f64>),
    IntMatrix(Vec<Vec<i32>>),
    IntArray(Vec<i32>),
    RangeArray(Vec<Range>),
}

impl MetaData {
    fn clone_data(&self) -> Box<MetaData> {
        match self {
            MetaData::StringMatrix(v) => Box::new(self.clone()),
            MetaData::StringArray(v)  => Box::new(self.clone()),
            MetaData::FloatMatrix(v)  => Box::new(self.clone()),
            MetaData::FloatArray(v)   => Box::new(self.clone()),
            MetaData::IntMatrix(v)    => Box::new(self.clone()),
            MetaData::IntArray(v)     => Box::new(self.clone()),
            MetaData::RangeArray(v)   => Box::new(self.clone()),
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
struct Meta {
    meta_name: Vec<String>,
    meta_data: Vec<Box<MetaData>>,
    rows: usize,
}

/* -------------------------------------------------------------------------- */

impl Meta {
    fn new(names: Vec<String>, data: Vec<Box<MetaData>>) -> Result<Self, Box<dyn Error>> {
        if names.len() != data.len() {
            return Err("NewMeta(): invalid parameters!".into());
        }
        let mut meta = Meta {
            meta_name: Vec::new(),
            meta_data: Vec::new(),
            rows: 0,
        };
        for i in 0..names.len() {
            meta.add_meta(names[i].clone(), data[i].clone())?;
        }
        Ok(meta)
    }

    fn clone(&self) -> Self {
        let mut result = Meta::new(Vec::new(), Vec::new()).unwrap();
        for i in 0..self.meta_length() {
            let cloned_data = self.meta_data[i].clone_data();
            result.add_meta(self.meta_name[i].clone(), cloned_data).unwrap();
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

    fn add_meta(&mut self, name: String, meta: Box<MetaData>) -> Result<(), Box<dyn Error>> {
        let n = match meta.as_ref() {
            MetaData::StringMatrix(v) => v.len(),
            MetaData::StringArray(v) => v.len(),
            MetaData::FloatMatrix(v) => v.len(),
            MetaData::FloatArray(v) => v.len(),
            MetaData::IntMatrix(v) => v.len(),
            MetaData::IntArray(v) => v.len(),
            MetaData::RangeArray(v) => v.len(),
        };
        if self.meta_length() > 0 {
            if n != self.rows {
                return Err(format!(
                    "AddMeta(): column `{}` has invalid length: expected length of `{}` but column has length `{}`",
                    name,
                    self.rows,
                    n
                )
                .into());
            }
        } else {
            self.rows = n;
        }
        self.delete_meta(&name);
        self.meta_data.push(meta);
        self.meta_name.push(name);
        Ok(())
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

    fn get_meta(&self, name: &str) -> Option<&Box<MetaData>> {
        for i in 0..self.meta_length() {
            if self.meta_name[i] == name {
                return Some(&self.meta_data[i]);
            }
        }
        None
    }

    fn get_meta_str(&self, name: &str) -> Vec<String> {
        match self.get_meta(name) {
            Some(meta) => match meta.as_ref() {
                MetaData::StringMatrix(v) => v.clone(),
                MetaData::StringArray(v)  => v.iter().map(|s| v.clone()).collect(),
                _ => Vec::new(),
            },
            None => Vec::new(),
        }
    }

    fn get_meta_float(&self, name: &str) -> Vec<f64> {
        match self.get_meta(name) {
            Some(meta) => match meta.as_ref() {
                MetaData::FloatMatrix(v) => v.clone(),
                MetaData::FloatArray(v) => v.iter().map(|f| *f).collect(),
                _ => Vec::new(),
            },
            None => Vec::new(),
        }
    }

    fn get_meta_int(&self, name: &str) -> Vec<i32> {
        match self.get_meta(name) {
            Some(meta) => match meta.as_ref() {
                MetaData::IntMatrix(v) => v.clone(),
                MetaData::IntArray(v) => v.iter().map(|i| *i).collect(),
                _ => Vec::new(),
            },
            None => Vec::new(),
        }
    }

    fn get_meta_range(&self, name: &str) -> Vec<Range> {
        match self.get_meta(name) {
            Some(meta) => match meta.as_ref() {
                MetaData::RangeArray(v) => v.clone(),
                _ => Vec::new(),
            },
            None => Vec::new(),
        }
    }

    fn append(&self, meta2: &Meta) -> Result<Meta, Box<dyn Error>> {
        let mut result = Meta::new(Vec::new(), Vec::new()).unwrap();
        let m1 = self.clone();
        let m2 = meta2.clone();
        for j in 0..m2.meta_length() {
            let t: Box<MetaData> = match m2.meta_data[j].as_ref() {
                MetaData::StringMatrix(v) => Box::new(v.clone()),
                MetaData::StringArray(v) => Box::new(v.clone()),
                MetaData::FloatMatrix(v) => Box::new(v.clone()),
                MetaData::FloatArray(v) => Box::new(v.clone()),
                MetaData::IntMatrix(v) => Box::new(v.clone()),
                MetaData::IntArray(v) => Box::new(v.clone()),
                MetaData::RangeArray(v) => Box::new(v.clone()),
            };
            let name = m2.meta_name[j].clone();
            result.add_meta(name, t)?;
        }
        Ok(result)
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
        let mut j = 0;
        for i in 0..self.length() {
            while j < indices.len() - 1 && i > indices[j] {
                j += 1;
            }
            if i != indices[j] {
                idx[j] = i;
                j += 1;
            }
        }
        self.subset(&idx)
    }

    fn subset(&self, indices: &[usize]) -> Meta {
        let n = indices.len();
        let m = self.meta_length();
        let mut data: Vec<Box<MetaData>> = Vec::new();
        for j in 0..m {
            let cloned_data: Box<MetaData> = match self.meta_data[j].as_ref() {
                MetaData::StringMatrix(v) => {
                    let mut l = Vec::new();
                    for i in 0..n {
                        l.push(v[indices[i]].clone());
                    }
                    Box::new(l)
                }
                MetaData::StringArray(v) => {
                    let mut l = Vec::new();
                    for i in 0..n {
                        l.push(v[indices[i]].clone());
                    }
                    Box::new(l)
                }
                MetaData::FloatMatrix(v) => {
                    let mut l = Vec::new();
                    for i in 0..n {
                        l.push(v[indices[i]].clone());
                    }
                    Box::new(l)
                }
                MetaData::FloatArray(v) => {
                    let mut l = Vec::new();
                    for i in 0..n {
                        l.push(v[indices[i]]);
                    }
                    Box::new(l)
                }
                MetaData::IntMatrix(v) => {
                    let mut l = Vec::new();
                    for i in 0..n {
                        l.push(v[indices[i]].clone());
                    }
                    Box::new(l)
                }
                MetaData::IntArray(v) => {
                    let mut l = Vec::new();
                    for i in 0..n {
                        l.push(v[indices[i]]);
                    }
                    Box::new(l)
                }
                MetaData::RangeArray(v) => {
                    let mut l = Vec::new();
                    for i in 0..n {
                        l.push(v[indices[i]].clone());
                    }
                    Box::new(l)
                }
            };
            data.push(cloned_data);
        }
        Meta::new(self.meta_name.clone(), data).unwrap()
    }

    fn slice(&self, ifrom: usize, ito: usize) -> Meta {
        let n = ito - ifrom;
        let m = self.meta_length();
        let mut data: Vec<Box<MetaData>> = Vec::new();
        for j in 0..m {
            let cloned_data: Box<MetaData> = match self.meta_data[j].as_ref() {
                MetaData::StringMatrix(v) => {
                    let mut l = Vec::new();
                    for i in ifrom..ito {
                        l.push(v[i].clone());
                    }
                    Box::new(l)
                }
                MetaData::StringArray(v) => {
                    let mut l = Vec::new();
                    for i in ifrom..ito {
                        l.push(v[i].clone());
                    }
                    Box::new(l)
                }
                MetaData::FloatMatrix(v) => {
                    let mut l = Vec::new();
                    for i in ifrom..ito {
                        l.push(v[i].clone());
                    }
                    Box::new(l)
                }
                MetaData::FloatArray(v) => {
                    let mut l = Vec::new();
                    for i in ifrom..ito {
                        l.push(v[i]);
                    }
                    Box::new(l)
                }
                MetaData::IntMatrix(v) => {
                    let mut l = Vec::new();
                    for i in ifrom..ito {
                        l.push(v[i].clone());
                    }
                    Box::new(l)
                }
                MetaData::IntArray(v) => {
                    let mut l = Vec::new();
                    for i in ifrom..ito {
                        l.push(v[i]);
                    }
                    Box::new(l)
                }
                MetaData::RangeArray(v) => {
                    let mut l = Vec::new();
                    for i in ifrom..ito {
                        l.push(v[i].clone());
                    }
                    Box::new(l)
                }
            };
            data.push(cloned_data);
        }
        Meta::new(self.meta_name.clone(), data).unwrap()
    }

    fn merge(&self, indices: &[usize]) -> Result<Meta, Box<dyn Error>> {
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
        let mut data: Vec<Box<MetaData>> = Vec::new();
        for j in 0..m {
            let cloned_data: Box<MetaData> = match self.meta_data[j].as_ref() {
                MetaData::StringArray(v) => {
                    let mut l = vec![Vec::new(); n];
                    for i in 0..v.len() {
                        l[indices[i]].push(v[i].clone());
                    }
                    Box::new(l)
                }
                MetaData::FloatArray(v) => {
                    let mut l = vec![Vec::new(); n];
                    for i in 0..v.len() {
                        l[indices[i]].push(v[i]);
                    }
                    Box::new(l)
                }
                MetaData::IntArray(v) => {
                    let mut l = vec![Vec::new(); n];
                    for i in 0..v.len() {
                        l[indices[i]].push(v[i]);
                    }
                    Box::new(l)
                }
                _ => return Err("cannot merge".into()),
            };
            data.push(cloned_data);
        }
        Meta::new(self.meta_name.clone(), data)
    }

    fn reduce_string<F>(&mut self, name: &str, name_new: &str, f: F)
    where
        F: Fn(&[String]) -> String,
    {
        let t = match self.get_meta(name) {
            Some(meta) => match meta.as_ref() {
                MetaData::StringMatrix(v) => v.clone(),
                MetaData::StringArray(v) => v.iter().map(|s| s.clone()).collect(),
                _ => Vec::new(),
            },
            None => Vec::new(),
        };
        let r = t.iter().map(|v| v.clone()).collect::<Vec<String>>();
        if !name_new.is_empty() && name != name_new {
            self.add_meta(name_new.to_string(), Box::new(r)).unwrap();
        } else {
            self.delete_meta(name);
            self.add_meta(name.to_string(), Box::new(r)).unwrap();
        }
    }

    fn reduce_float<F>(&mut self, name: &str, name_new: &str, f: F)
    where
        F: Fn(&[f64]) -> f64,
    {
        let t = match self.get_meta(name) {
            Some(meta) => match meta.as_ref() {
                MetaData::FloatMatrix(v) => v.clone(),
                MetaData::FloatArray(v) => v.iter().map(|f| *f).collect(),
                _ => Vec::new(),
            },
            None => Vec::new(),
        };
        let r = t.iter().map(|v| *v).collect::<Vec<f64>>();
        if !name_new.is_empty() && name != name_new {
            self.add_meta(name_new.to_string(), Box::new(r)).unwrap();
        } else {
            self.delete_meta(name);
            self.add_meta(name.to_string(), Box::new(r)).unwrap();
        }
    }

    fn reduce_int<F>(&mut self, name: &str, name_new: &str, f: F)
    where
        F: Fn(&[i32]) -> i32,
    {
        let t = match self.get_meta(name) {
            Some(meta) => match meta.as_ref() {
                MetaData::IntMatrix(v) => v.clone(),
                MetaData::IntArray(v) => v.iter().map(|i| *i).collect(),
                _ => Vec::new(),
            },
            None => Vec::new(),
        };
        let r = t.iter().map(|v| *v).collect::<Vec<i32>>();
        if !name_new.is_empty() && name != name_new {
            self.add_meta(name_new.to_string(), Box::new(r)).unwrap();
        } else {
            self.delete_meta(name);
            self.add_meta(name.to_string(), Box::new(r)).unwrap();
        }
    }

    fn sorted_indices(&self, name: &str, reverse: bool) -> Result<Vec<usize>, Box<dyn Error>> {
        let mut l: Vec<(usize, f64, i32, String)> = Vec::new();
        if let Some(meta) = self.get_meta(name) {
            match meta.as_ref() {
                MetaData::FloatArray(v) => {
                    for (i, &f) in v.iter().enumerate() {
                        l.push((i, f, 0, "".to_string()));
                    }
                }
                MetaData::IntArray(v) => {
                    for (i, &j) in v.iter().enumerate() {
                        l.push((i, 0.0, j, "".to_string()));
                    }
                }
                MetaData::StringArray(v) => {
                    for (i, s) in v.iter().enumerate() {
                        l.push((i, 0.0, 0, s.clone()));
                    }
                }
                _ => return Err("Invalid type for sorting!".into()),
            }
        } else {
            return Err("Meta column not found!".into());
        }
        if reverse {
            l.sort_by(|a, b| match (a.1, a.2, a.3) {
                (fa, _, _) => match (b.1, b.2, b.3) {
                    (fb, _, _) => fa.partial_cmp(&fb).unwrap(),
                },
            });
        } else {
            l.sort_by(|a, b| match (a.1, a.2, a.3) {
                (fa, _, _) => match (b.1, b.2, b.3) {
                    (fb, _, _) => fb.partial_cmp(&fa).unwrap(),
                },
            });
        }
        let j: Vec<usize> = l.iter().map(|(i, _, _, _)| *i).collect();
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
            let format = format!(" {:width$}", "", width = widths[j] - 1);
            let l = write!(writer, "{}", format_args!(format, String::from_utf8_lossy(&tmp_buffer)))?;
            Ok(l)
        };
        let print_cell = |writer: &mut W, widths: &[usize], i: usize, j: usize| -> Result<usize, Box<dyn Error>> {
            match self.meta_data[j].as_ref().data_type() {
                DataType::String1D(v) => {
                    let format = format!(" {:width$}", v[i], width = widths[j] - 1);
                    write!(writer, "{}", format_args!(format))?;
                }
                DataType::Float1D(v) => {
                    if use_scientific {
                        let format = format!(" {:width$e}", v[i], width = widths[j] - 1);
                        write!(writer, "{}", format_args!(format))?;
                    } else {
                        let format = format!(" {:width$.6}", v[i], width = widths[j] - 1);
                        write!(writer, "{}", format_args!(format))?;
                    }
                }
                DataType::Int1D(v) => {
                    let format = format!(" {:width$}", v[i], width = widths[j] - 1);
                    write!(writer, "{}", format_args!(format))?;
                }
                DataType::Range1D(v) => {
                    let format = format!(" [{}, {}]", v[i].from, v[i].to);
                    write!(writer, "{}", format)?;
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
                let format = format!(" {:width$}", self.meta_name[j], width = widths[j] - 1);
                write!(writer, "{}", format_args!(format))?;
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
                    let format = format!(" {:width$}", "...", width = widths[j] - 1);
                    write!(writer, "{}", format_args!(format))?;
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
                    let format = format!(" {:width$}", "...", width = widths[j] - 1);
                    write!(writer, "{}", format_args!(format))?;
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

#[derive(Debug, Clone)]
struct MetaRow {
    meta: Meta,
    row: usize,
}

impl MetaRow {
    fn new(meta: Meta, row: usize) -> Self {
        MetaRow { meta, row }
    }

    fn get_meta(&self, name: &str) -> Option<&Box<MetaData>> {
        self.meta.get_meta(name)
    }

    fn get_meta_str(&self, name: &str) -> Vec<String> {
        self.meta.get_meta_str(name)
    }

    fn get_meta_float(&self, name: &str) -> Vec<f64> {
        self.meta.get_meta_float(name)
    }

    fn get_meta_int(&self, name: &str) -> Vec<i32> {
        self.meta.get_meta_int(name)
    }

    fn get_meta_range(&self, name: &str) -> Vec<Range> {
        self.meta.get_meta_range(name)
    }
}
