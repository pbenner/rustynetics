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
use list_comprehension_macro::comp;

use crate::error::Error;
use crate::range::Range;

/* -------------------------------------------------------------------------- */

#[derive(Debug, Clone)]
pub enum MetaData {
    StringMatrix(Vec<Vec<String>>),
    StringArray(Vec<String>),
    FloatMatrix(Vec<Vec<f64>>),
    FloatArray(Vec<f64>),
    IntMatrix(Vec<Vec<i64>>),
    IntArray(Vec<i64>),
    RangeArray(Vec<Range>),
}

/* -------------------------------------------------------------------------- */

impl MetaData {
    pub fn len(&self) -> usize {
        match self {
            MetaData::FloatArray  (v) => v.len(),
            MetaData::IntArray    (v) => v.len(),
            MetaData::StringArray (v) => v.len(),
            MetaData::StringMatrix(v) => v.len(),
            MetaData::FloatMatrix (v) => v.len(),
            MetaData::IntMatrix   (v) => v.len(),
            MetaData::RangeArray  (v) => v.len(),
        }
    }

    pub fn concat(&self, data : &Self) -> Result<Self, Error> {
        match (self, data) {
            (MetaData::FloatArray  (v), MetaData::FloatArray  (w)) => Ok(MetaData::FloatArray  ([v.clone(), w.clone()].concat())),
            (MetaData::IntArray    (v), MetaData::IntArray    (w)) => Ok(MetaData::IntArray    ([v.clone(), w.clone()].concat())),
            (MetaData::StringArray (v), MetaData::StringArray (w)) => Ok(MetaData::StringArray ([v.clone(), w.clone()].concat())),
            (MetaData::StringMatrix(v), MetaData::StringMatrix(w)) => Ok(MetaData::StringMatrix([v.clone(), w.clone()].concat())),
            (MetaData::FloatMatrix (v), MetaData::FloatMatrix (w)) => Ok(MetaData::FloatMatrix ([v.clone(), w.clone()].concat())),
            (MetaData::IntMatrix   (v), MetaData::IntMatrix   (w)) => Ok(MetaData::IntMatrix   ([v.clone(), w.clone()].concat())),
            (MetaData::RangeArray  (v), MetaData::RangeArray  (w)) => Ok(MetaData::RangeArray  ([v.clone(), w.clone()].concat())),
            _ => Err(format!("MetaData types do not match").into())
        }
    }

    pub fn slice(&self, ifrom : usize, ito : usize) -> Self {
        match self {
            MetaData::FloatArray  (v) => MetaData::FloatArray  (comp![ v[i] for i in ifrom..ito ]),
            MetaData::IntArray    (v) => MetaData::IntArray    (comp![ v[i] for i in ifrom..ito ]),
            MetaData::StringArray (v) => MetaData::StringArray (comp![ v[i].clone() for i in ifrom..ito ]),
            MetaData::StringMatrix(v) => MetaData::StringMatrix(comp![ v[i].clone() for i in ifrom..ito ]),
            MetaData::FloatMatrix (v) => MetaData::FloatMatrix (comp![ v[i].clone() for i in ifrom..ito ]),
            MetaData::IntMatrix   (v) => MetaData::IntMatrix   (comp![ v[i].clone() for i in ifrom..ito ]),
            MetaData::RangeArray  (v) => MetaData::RangeArray  (comp![ v[i].clone() for i in ifrom..ito ]),
        }
    }

    pub fn subset(&self, indices : &[usize]) -> Self {
        match self {
            MetaData::FloatArray  (v) => MetaData::FloatArray  (comp![ v[*i] for i in indices ]),
            MetaData::IntArray    (v) => MetaData::IntArray    (comp![ v[*i] for i in indices ]),
            MetaData::StringArray (v) => MetaData::StringArray (comp![ v[*i].clone() for i in indices ]),
            MetaData::StringMatrix(v) => MetaData::StringMatrix(comp![ v[*i].clone() for i in indices ]),
            MetaData::FloatMatrix (v) => MetaData::FloatMatrix (comp![ v[*i].clone() for i in indices ]),
            MetaData::IntMatrix   (v) => MetaData::IntMatrix   (comp![ v[*i].clone() for i in indices ]),
            MetaData::RangeArray  (v) => MetaData::RangeArray  (comp![ v[*i].clone() for i in indices ]),
        }
    }

    pub fn get_float(&self) -> Option<&Vec<f64>> {
        match self {
            MetaData::FloatArray(v) => Some(v),
            _ => None
        }
    }

    pub fn get_int(&self) -> Option<&Vec<i64>> {
        match self {
            MetaData::IntArray(v) => Some(v),
            _ => None
        }
    }

    pub fn get_str(&self) -> Option<&Vec<String>> {
        match self {
            MetaData::StringArray(v) => Some(v),
            _ => None
        }
    }

    pub fn get_float_mut(&mut self) -> Option<&mut Vec<f64>> {
        match self {
            MetaData::FloatArray(v) => Some(v),
            _ => None
        }
    }

    pub fn get_int_mut(&mut self) -> Option<&mut Vec<i64>> {
        match self {
            MetaData::IntArray(v) => Some(v),
            _ => None
        }
    }

    pub fn get_str_mut(&mut self) -> Option<&mut Vec<String>> {
        match self {
            MetaData::StringArray(v) => Some(v),
            _ => None
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct Meta {
    pub meta_name: Vec<String>,
    pub meta_data: Vec<MetaData>,
    rows: usize,
}

/* -------------------------------------------------------------------------- */

impl Meta {
    pub fn new(names: Vec<String>, data: Vec<MetaData>) -> Result<Self, Error> {
        if names.len() != data.len() {
            return Err(format!("Invalid parameters!").into());
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

    pub fn new_empty(rows : usize) -> Meta {

        let meta = Meta {
            meta_name: Vec::new(),
            meta_data: Vec::new(),
            rows: rows,
        };
        meta
    }

    pub fn num_rows(&self) -> usize {
        self.rows
    }

    pub fn num_cols(&self) -> usize {
        self.meta_name.len()
    }

    pub fn append(&self, meta : &Meta) -> Result<Self, Error> {
        let n1 = self.num_rows();
        let n2 = meta.num_rows();

        let mut meta_name = Vec::new();
        let mut meta_data = Vec::new();

        for j in 0..self.num_cols() {
            let meta_col1 = self.meta_data[j].clone();
            let meta_col2 = meta.get_column(&meta.meta_name[j]);

            if meta_col2.is_none() {
                return Err(format!("Column {} not found in meta object", &meta.meta_name[j]).into())
            }
            let meta_col3 = meta_col1.concat(meta_col2.unwrap())?;

            meta_name.push(self.meta_name[j].clone());
            meta_data.push(meta_col3);
        }
        let meta = Meta {
            meta_name: meta_name,
            meta_data: meta_data,
            rows: n1+n2,
        };
        Ok(meta)
    }

    pub fn add_meta(&mut self, name: String, data: MetaData) -> Result<(), Error> {
        let n = data.len();
        if self.meta_name.len() > 0 {
            if n != self.rows {
                return Err(format!("Column '{}' has invalid length: expected length of '{}' but column has length '{}'", name, self.rows, n).into());
            }
        }
        if self.meta_name.len() == 0 {
            self.rows = n;
        }
        self.meta_name.push(name);
        self.meta_data.push(data);
        Ok(())
    }

    pub fn delete_meta(&mut self, name: &String) {
        if let Some(index) = self.meta_name.iter().position(|x| x == name) {
            self.meta_name.remove(index);
            self.meta_data.remove(index);
        }
    }

    pub fn rename_meta(&mut self, name_old: &String, name_new: &String) {
        if name_old == name_new {
            return;
        }
        if let Some(index) = self.meta_name.iter().position(|x| x == name_old) {
            self.meta_name[index] = name_new.to_string();
        }
    }

    pub fn get_column(&self, name: &String) -> Option<&MetaData> {
        self.meta_name.iter().position(|x| x == name).map(|index| &self.meta_data[index])
    }

    pub fn get_column_int(&self, name: &String) -> Option<&Vec<i64>> {
        let r = self.get_column(name)?;
        r.get_int()
    }

    pub fn get_column_float(&self, name: &String) -> Option<&Vec<f64>> {
        let r = self.get_column(name)?;
        r.get_float()
    }

    pub fn get_column_str(&self, name: &String) -> Option<&Vec<String>> {
        let r = self.get_column(name)?;
        r.get_str()
    }

    pub fn get_column_mut(&mut self, name: &String) -> Option<&mut MetaData> {
        self.meta_name.iter().position(|x| x == name).map(move |index| &mut self.meta_data[index])
    }

    pub fn get_column_int_mut(&mut self, name: &String) -> Option<&mut Vec<i64>> {
        let r = self.get_column_mut(name)?;
        r.get_int_mut()
    }

    pub fn get_column_float_mut(&mut self, name: &String) -> Option<&mut Vec<f64>> {
        let r = self.get_column_mut(name)?;
        r.get_float_mut()
    }

    pub fn get_column_str_mut(&mut self, name: &String) -> Option<&mut Vec<String>> {
        let r = self.get_column_mut(name)?;
        r.get_str_mut()
    }

    pub fn slice(&self, ifrom : usize, ito : usize) -> Meta {

        let n = ito-ifrom;
        let m = self.meta_name.len();
        let mut data = Vec::new();

        for j in 0..m {
            data.push(self.meta_data[j].slice(ifrom, ito));
        }

        Meta {
            meta_name: self.meta_name.clone(),
            meta_data: data,
            rows: n,
        }
    }

    pub fn subset(&self, indices: &[usize]) -> Meta {
        let n = indices.len();
        let m = self.meta_name.len();
        let mut data = Vec::new();

        for j in 0..m {
            data.push(self.meta_data[j].subset(indices));
        }

        Meta {
            meta_name: self.meta_name.clone(),
            meta_data: data,
            rows: n,
        }
    }

    pub fn sort(&self, name: &String, reverse: bool) -> Result<Self, Error> {
        let mut indices: Vec<usize> = (0..self.rows).collect();

        if reverse {
            match self.get_column(name).unwrap() {
                MetaData::StringArray(v) => indices.sort_by(|&i, &j| v[j].cmp(&v[i])),
                MetaData::FloatArray (v) => indices.sort_by(|&i, &j| v[j].partial_cmp(&v[i]).unwrap()),
                MetaData::IntArray   (v) => indices.sort_by(|&i, &j| v[j].cmp(&v[i])),
                _ => ()
            }
        }
        else {
            match self.get_column(name).unwrap() {
                MetaData::StringArray(v) => indices.sort_by(|&i, &j| v[i].cmp(&v[j])),
                MetaData::FloatArray (v) => indices.sort_by(|&i, &j| v[i].partial_cmp(&v[j]).unwrap()),
                MetaData::IntArray   (v) => indices.sort_by(|&i, &j| v[i].cmp(&v[j])),
                _ => ()
            }
        }
        Ok(self.subset(&indices))
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for Meta {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.pad(&format!("{}", self.pretty_string(10, false).unwrap()))
    }
}
