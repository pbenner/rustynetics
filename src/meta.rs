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
use std::error::Error;
use list_comprehension_macro::comp;

use crate::range::Range;

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

/* -------------------------------------------------------------------------- */

impl MetaData {
    fn len(&self) -> usize {
        self.len()
    }

    fn subset(&self, indices : &[usize]) -> Self {
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

#[derive(Debug)]
struct Meta {
    meta_name: Vec<String>,
    meta_data: Vec<MetaData>,
    rows: usize,
}

/* -------------------------------------------------------------------------- */

impl Meta {
    fn new(names: Vec<String>, data: Vec<MetaData>) -> Result<Meta, Box<dyn Error>> {
        if names.len() != data.len() {
            return Err("Invalid parameters!".into());
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

    fn add_meta(&mut self, name: String, meta: MetaData) -> Result<(), Box<dyn Error>> {
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

    fn delete_meta(&mut self, name: &str) {
        if let Some(index) = self.meta_name.iter().position(|x| x == name) {
            self.meta_name.remove(index);
            self.meta_data.remove(index);
        }
    }

    fn rename_meta(&mut self, name_old: &str, name_new: &str) {
        if name_old == name_new {
            return;
        }
        if let Some(index) = self.meta_name.iter().position(|x| x == name_old) {
            self.meta_name[index] = name_new.to_string();
        }
    }

    fn get_meta(&self, name: &str) -> Option<&MetaData> {
        self.meta_name.iter().position(|x| x == name).map(|index| &self.meta_data[index])
    }

    fn subset(&self, indices: &[usize]) -> Result<Meta, Box<dyn Error>> {
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

    fn sort(&self, name: &str, reverse: bool) -> Result<Meta, Box<dyn Error>> {
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
