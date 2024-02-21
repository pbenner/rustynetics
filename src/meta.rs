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
    IntMatrix(Vec<Vec<i32>>),
    IntArray(Vec<i32>),
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
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct Meta {
    meta_name: Vec<String>,
    meta_data: Vec<MetaData>,
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

        let mut meta = Meta {
            meta_name: Vec::new(),
            meta_data: Vec::new(),
            rows: rows,
        };
        meta
    }

    pub fn add_meta(&mut self, name: String, meta: MetaData) -> Result<Self, Error> {
        let n = meta.len();
        if self.meta_name.len() > 0 {
            if n != self.rows {
                return Err(format!("Column '{}' has invalid length: expected length of '{}' but column has length '{}'", name, self.rows, n).into());
            }
        } else {
            self.rows = n;
        }
        self.delete_meta(&name);
        self.meta_name.push(name);
        self.meta_data.push(meta);
        Ok(())
    }

    pub fn append(&mut self, meta : &Meta) -> Result<(), Error> {
        let m = self.meta_name.len();

        for j in 0..m {
            self.add_meta(meta.meta_name[j], meta.meta_data[j])?;
        }
        Ok(())
    }

    pub fn delete_meta(&mut self, name: &str) {
        if let Some(index) = self.meta_name.iter().position(|x| x == name) {
            self.meta_name.remove(index);
            self.meta_data.remove(index);
        }
    }

    pub fn rename_meta(&mut self, name_old: &str, name_new: &str) {
        if name_old == name_new {
            return;
        }
        if let Some(index) = self.meta_name.iter().position(|x| x == name_old) {
            self.meta_name[index] = name_new.to_string();
        }
    }

    pub fn get_meta(&self, name: &str) -> Option<&MetaData> {
        self.meta_name.iter().position(|x| x == name).map(|index| &self.meta_data[index])
    }

    pub fn slice(&self, ifrom : usize, ito : usize) -> Result<Meta, Error> {

        let n = ito-ifrom;
        let m = self.meta_name.len();
        let mut data = Vec::new();

        for j in 0..m {
            data.push(self.meta_data[j].slice(ifrom, ito));
        }

        Ok(Meta {
            meta_name: self.meta_name.clone(),
            meta_data: data,
            rows: n,
        })
    }

    pub fn subset(&self, indices: &[usize]) -> Result<Meta, Error> {
        let n = indices.len();
        let m = self.meta_name.len();
        let mut data = Vec::new();

        for j in 0..m {
            data.push(self.meta_data[j].subset(indices));
        }

        Ok(Meta {
            meta_name: self.meta_name.clone(),
            meta_data: data,
            rows: n,
        })
    }

    pub fn sort(&self, name: &str, reverse: bool) -> Result<Self, Error> {
        let mut indices: Vec<usize> = (0..self.rows).collect();

        if reverse {
            match self.get_meta(name).unwrap() {
                MetaData::StringArray(v) => indices.sort_by(|&i, &j| v[j].cmp(&v[i])),
                MetaData::FloatArray (v) => indices.sort_by(|&i, &j| v[j].partial_cmp(&v[i]).unwrap()),
                MetaData::IntArray   (v) => indices.sort_by(|&i, &j| v[j].cmp(&v[i])),
                _ => ()
            }
        }
        else {
            match self.get_meta(name).unwrap() {
                MetaData::StringArray(v) => indices.sort_by(|&i, &j| v[i].cmp(&v[j])),
                MetaData::FloatArray (v) => indices.sort_by(|&i, &j| v[i].partial_cmp(&v[j]).unwrap()),
                MetaData::IntArray   (v) => indices.sort_by(|&i, &j| v[i].cmp(&v[j])),
                _ => ()
            }
        }

        self.subset(&indices)
    }
}

/* -------------------------------------------------------------------------- */

//impl fmt::Display for Meta {
//    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//        for i in 0..self.rows {
//            for j in 0..self.meta_name.len() {
//                write!(f, "{}: {}\n", self.meta_name[j], self.meta_data[j][i])?;
//            }
//            write!(f, "\n")?;
//        }
//        Ok(())
//    }
//}
