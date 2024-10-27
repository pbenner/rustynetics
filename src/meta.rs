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

/* -------------------------------------------------------------------------- */

use std::fmt;
use std::error::Error;

use list_comprehension_macro::comp;

use crate::range::Range;

/* -------------------------------------------------------------------------- */

#[derive(Debug, Clone)]
pub enum MetaData {
    StringMatrix (Vec<Vec<String>>),
    StringArray  (Vec<String>     ),
    FloatMatrix  (Vec<Vec<f64>>   ),
    FloatArray   (Vec<f64>        ),
    IntMatrix    (Vec<Vec<i64>>   ),
    IntArray     (Vec<i64>        ),
    RangeArray   (Vec<Range>      ),
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

    pub fn concat(&self, data : &Self) -> Result<Self, Box<dyn Error>> {
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

    pub fn get_range(&self) -> Option<&Vec<Range>> {
        match self {
            MetaData::RangeArray(v) => Some(v),
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

    pub fn get_range_mut(&mut self) -> Option<&mut Vec<Range>> {
        match self {
            MetaData::RangeArray(v) => Some(v),
            _ => None
        }
    }

}

/* -------------------------------------------------------------------------- */

impl PartialEq for MetaData {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (MetaData::FloatArray  (v), MetaData::FloatArray  (w)) => v == w,
            (MetaData::IntArray    (v), MetaData::IntArray    (w)) => v == w,
            (MetaData::StringArray (v), MetaData::StringArray (w)) => v == w,
            (MetaData::StringMatrix(v), MetaData::StringMatrix(w)) => v == w,
            (MetaData::FloatMatrix (v), MetaData::FloatMatrix (w)) => v == w,
            (MetaData::IntMatrix   (v), MetaData::IntMatrix   (w)) => v == w,
            (MetaData::RangeArray  (v), MetaData::RangeArray  (w)) => v == w,
            _ => false
        }
    }
}

/* -------------------------------------------------------------------------- */

/// A structure for storing metadata associated with genomic ranges in a
/// flexible, column-based format. The `Meta` struct allows storing columns of
/// various data types (integers, floats, strings, and ranges), facilitating
/// metadata storage and manipulation for genomic datasets.
///
/// Each `Meta` instance comprises named metadata columns, where each column can
/// have a different type (handled by `MetaData`). Columns are organized by name
/// in `meta_name` and data in `meta_data`, with rows representing entries.
///
/// # Fields
/// - `meta_name`: A vector of column names, each corresponding to an entry in `meta_data`.
/// - `meta_data`: A vector of `MetaData`, each holding values for a column.
/// - `rows`: The number of rows (entries) in each column.
///
/// # Examples
/// ```
/// use rustynetics::meta::{Meta, MetaData};
///
/// // Creating a new Meta instance with column names and data
/// let meta = Meta::new(
///     vec!["Gene", "Expression"],
///     vec![MetaData::StringArray(vec!["GeneA".to_string(), "GeneB".to_string()]),
///          MetaData::FloatArray(vec![0.85, 1.23])]
/// ).unwrap();
/// ```
///
/// # Usage
/// The `Meta` struct provides methods for adding, appending, retrieving, and
/// managing metadata columns. Metadata can be accessed and modified by column name,
/// and transformations such as subsetting, sorting, or slicing rows are also supported.
///
/// # Note
/// This struct is commonly used in bioinformatics to attach structured annotations
/// to genomic ranges, enabling efficient storage and retrieval of associated metadata.
#[derive(Clone, Debug)]
pub struct Meta {
    pub meta_name: Vec<String>,
    pub meta_data: Vec<MetaData>,
    rows: usize,
}

/* -------------------------------------------------------------------------- */

impl Default for Meta {
    fn default() -> Meta {
        Self {
            meta_name: Vec::new(),
            meta_data: Vec::new(),
            rows: 0,
        }
    }
}

/* -------------------------------------------------------------------------- */

impl Meta {

    /// Creates a new `Meta` instance with the given column names and data.
    ///
    /// # Arguments
    /// - `names`: A vector of names for each column.
    /// - `data`: A vector of `MetaData`, one for each column.
    ///
    /// # Errors
    /// Returns an error if `names` and `data` have mismatched lengths.
    pub fn new(names: Vec<&str>, data: Vec<MetaData>) -> Result<Self, Box<dyn Error>> {
        if names.len() != data.len() {
            return Err(format!("Invalid parameters!").into());
        }
        let mut meta = Meta {
            meta_name: Vec::new(),
            meta_data: Vec::new(),
            rows: 0,
        };
        for i in 0..names.len() {
            meta.add(names[i], data[i].clone())?;
        }
        Ok(meta)
    }

    /// Returns the number of rows in the `Meta` instance.
    pub fn num_rows(&self) -> usize {
        self.rows
    }

    /// Returns the number of columns in the `Meta` instance.
    pub fn num_cols(&self) -> usize {
        self.meta_name.len()
    }

    /// Appends rows from another `Meta` instance to the current instance, matching columns by name.
    ///
    /// This function combines two `Meta` instances by appending rows from `meta` to the existing
    /// rows in `self`. The column names in both `Meta` instances must match exactly.
    ///
    /// # Arguments
    /// - `meta`: A reference to the `Meta` instance to append to `self`.
    ///
    /// # Returns
    /// A new `Meta` instance containing rows from both `self` and `meta`.
    ///
    /// # Errors
    /// Returns an error if:
    /// - The number of rows in `meta` columns does not match the columns in `self`.
    /// - A column name in `meta` is not found in `self`.
    ///
    /// # Example
    /// ```
    /// use rustynetics::meta::{Meta, MetaData};
    ///
    /// let meta1 = Meta::new(
    ///     vec!["Gene", "Expression"],
    ///     vec![MetaData::StringArray(vec!["GeneA".to_string()]),
    ///          MetaData::FloatArray(vec![0.85])]
    /// ).unwrap();
    ///
    /// let meta2 = Meta::new(
    ///     vec!["Gene", "Expression"],
    ///     vec![MetaData::StringArray(vec!["GeneB".to_string()]),
    ///          MetaData::FloatArray(vec![1.23])]
    /// ).unwrap();
    ///
    /// let combined = meta1.append(&meta2).unwrap();
    /// assert_eq!(combined.num_rows(), 2);
    /// ```
    pub fn append(&self, meta : &Meta) -> Result<Self, Box<dyn Error>> {
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

    /// Adds a new column to the `Meta` instance.
    ///
    /// # Arguments
    /// - `name`: The name of the new column.
    /// - `data`: The `MetaData` for the new column.
    ///
    /// # Errors
    /// Returns an error if the length of `data` does not match the existing row count.
    pub fn add(&mut self, name: &str, data: MetaData) -> Result<(), Box<dyn Error>> {
        let n = data.len();
        if self.meta_name.len() > 0 {
            if n != self.rows {
                return Err(format!("Column '{}' has invalid length: expected length of '{}' but column has length '{}'", name, self.rows, n).into());
            }
        }
        if self.meta_name.len() == 0 {
            self.rows = n;
        }
        self.meta_name.push(String::from(name));
        self.meta_data.push(data);
        Ok(())
    }

    /// Deletes a column by name, if it exists.
    ///
    /// # Arguments
    /// - `name`: The name of the column to delete.
    pub fn delete_meta(&mut self, name: &str) {
        if let Some(index) = self.meta_name.iter().position(|x| x == name) {
            self.meta_name.remove(index);
            self.meta_data.remove(index);
        }
    }

    /// Renames a column in `Meta` by name.
    ///
    /// # Arguments
    /// - `name_old`: The existing name of the column.
    /// - `name_new`: The new name for the column.
    pub fn rename_meta(&mut self, name_old: &str, name_new: &str) {
        if name_old == name_new {
            return;
        }
        if let Some(index) = self.meta_name.iter().position(|x| x == name_old) {
            self.meta_name[index] = name_new.to_string();
        }
    }

    /// Retrieves a reference to a column in the `Meta` instance by its name.
    ///
    /// This function searches for a column with the specified name and returns a reference to the
    /// `MetaData` object associated with that column if it exists.
    ///
    /// # Arguments
    /// - `name`: The name of the column to retrieve.
    ///
    /// # Returns
    /// An `Option` containing a reference to the `MetaData` of the specified column. 
    /// Returns `None` if no column with the given name is found.
    ///
    /// # Example
    /// ```
    /// use rustynetics::meta::{Meta, MetaData};
    ///
    /// let meta = Meta::new(
    ///     vec!["Gene", "Expression"],
    ///     vec![MetaData::StringArray(vec!["GeneA".to_string(), "GeneB".to_string()]),
    ///          MetaData::FloatArray(vec![1.0, 2.0])]
    /// ).unwrap();
    ///
    /// if let Some(column) = meta.get_column("Expression") {
    ///     // Process column data here
    /// }
    /// ```
    ///
    /// # Errors
    /// This function does not return an error but rather returns `None` if the specified
    /// column name is not found.
    pub fn get_column(&self, name: &str) -> Option<&MetaData> {
        self.meta_name.iter().position(|x| x == name).map(|index| &self.meta_data[index])
    }

    /// Retrieves a reference to a column containing integer data by its name.
    ///
    /// This function returns a reference to the `Vec<i64>` of integers for the specified column
    /// if the column contains integer data. If the column exists but does not contain integers,
    /// or if the column name is not found, it returns `None`.
    ///
    /// # Arguments
    /// - `name`: The name of the column to retrieve.
    ///
    /// # Returns
    /// An `Option` containing a reference to the integer `Vec<i64>`. Returns `None` if
    /// the column is not found or if it contains non-integer data.
    pub fn get_column_int(&self, name: &str) -> Option<&Vec<i64>> {
        let r = self.get_column(name)?;
        r.get_int()
    }

    /// Retrieves a reference to a column containing floating-point data by its name.
    ///
    /// This function returns a reference to the `Vec<f64>` of floats for the specified column
    /// if the column contains float data. If the column exists but does not contain floats,
    /// or if the column name is not found, it returns `None`.
    ///
    /// # Arguments
    /// - `name`: The name of the column to retrieve.
    ///
    /// # Returns
    /// An `Option` containing a reference to the float `Vec<f64>`. Returns `None` if
    /// the column is not found or if it contains non-float data.
    pub fn get_column_float(&self, name: &str) -> Option<&Vec<f64>> {
        let r = self.get_column(name)?;
        r.get_float()
    }

    /// Retrieves a reference to a column containing string data by its name.
    ///
    /// This function returns a reference to the `Vec<String>` of strings for the specified column
    /// if the column contains string data. If the column exists but does not contain strings,
    /// or if the column name is not found, it returns `None`.
    ///
    /// # Arguments
    /// - `name`: The name of the column to retrieve.
    ///
    /// # Returns
    /// An `Option` containing a reference to the string `Vec<String>`. Returns `None` if
    /// the column is not found or if it contains non-string data.
    pub fn get_column_str(&self, name: &str) -> Option<&Vec<String>> {
        let r = self.get_column(name)?;
        r.get_str()
    }

    /// Retrieves a reference to a column containing range data by its name.
    ///
    /// This function returns a reference to the `Vec<Range>` of ranges for the specified column
    /// if the column contains range data. If the column exists but does not contain ranges,
    /// or if the column name is not found, it returns `None`.
    ///
    /// # Arguments
    /// - `name`: The name of the column to retrieve.
    ///
    /// # Returns
    /// An `Option` containing a reference to the range `Vec<Range>`. Returns `None` if
    /// the column is not found or if it contains non-range data.
    pub fn get_column_range(&self, name: &str) -> Option<&Vec<Range>> {
        let r = self.get_column(name)?;
        r.get_range()
    }

    /// Retrieves a mutable reference to a column in the `Meta` instance by its name.
    ///
    /// This function returns a mutable reference to the `MetaData` object of the specified column,
    /// allowing modifications to the column data if it exists. If the column name is not found,
    /// it returns `None`.
    ///
    /// # Arguments
    /// - `name`: The name of the column to retrieve.
    ///
    /// # Returns
    /// An `Option` containing a mutable reference to the `MetaData` of the specified column.
    /// Returns `None` if no column with the given name is found.
    ///
    /// # Example
    /// ```
    /// use rustynetics::meta::{Meta, MetaData};
    ///
    /// let mut meta = Meta::new(
    ///     vec!["Count"],
    ///     vec![MetaData::IntArray(vec![10, 20, 30])],
    /// ).unwrap();
    ///
    /// if let Some(column) = meta.get_column_mut("Count") {
    ///     if let Some(r) = column.get_int_mut() {
    ///         r[1] = 40;
    ///     }
    /// }
    /// ```
    pub fn get_column_mut(&mut self, name: &str) -> Option<&mut MetaData> {
        self.meta_name.iter().position(|x| x == name).map(move |index| &mut self.meta_data[index])
    }

    /// Retrieves a mutable reference to a column containing integer data by its name.
    ///
    /// This function returns a mutable reference to the `Vec<i64>` of integers for the specified column
    /// if the column exists and contains integer data. If the column exists but does not contain integers,
    /// or if the column name is not found, it returns `None`.
    ///
    /// # Arguments
    /// - `name`: The name of the column to retrieve.
    ///
    /// # Returns
    /// An `Option` containing a mutable reference to the integer `Vec<i64>`. Returns `None` if
    /// the column is not found or if it contains non-integer data.
    pub fn get_column_int_mut(&mut self, name: &str) -> Option<&mut Vec<i64>> {
        let r = self.get_column_mut(name)?;
        r.get_int_mut()
    }

    /// Retrieves a mutable reference to a column containing floating-point data by its name.
    ///
    /// This function returns a mutable reference to the `Vec<f64>` of floats for the specified column
    /// if the column exists and contains float data. If the column exists but does not contain floats,
    /// or if the column name is not found, it returns `None`.
    ///
    /// # Arguments
    /// - `name`: The name of the column to retrieve.
    ///
    /// # Returns
    /// An `Option` containing a mutable reference to the float `Vec<f64>`. Returns `None` if
    /// the column is not found or if it contains non-float data.
    pub fn get_column_float_mut(&mut self, name: &str) -> Option<&mut Vec<f64>> {
        let r = self.get_column_mut(name)?;
        r.get_float_mut()
    }

    /// Retrieves a mutable reference to a column containing string data by its name.
    ///
    /// This function returns a mutable reference to the `Vec<String>` of strings for the specified column
    /// if the column exists and contains string data. If the column exists but does not contain strings,
    /// or if the column name is not found, it returns `None`.
    ///
    /// # Arguments
    /// - `name`: The name of the column to retrieve.
    ///
    /// # Returns
    /// An `Option` containing a mutable reference to the string `Vec<String>`. Returns `None` if
    /// the column is not found or if it contains non-string data.
    pub fn get_column_str_mut(&mut self, name: &str) -> Option<&mut Vec<String>> {
        let r = self.get_column_mut(name)?;
        r.get_str_mut()
    }

    /// Retrieves a mutable reference to a column containing range data by its name.
    ///
    /// This function returns a mutable reference to the `Vec<Range>` of ranges for the specified column
    /// if the column exists and contains range data. If the column exists but does not contain ranges,
    /// or if the column name is not found, it returns `None`.
    ///
    /// # Arguments
    /// - `name`: The name of the column to retrieve.
    ///
    /// # Returns
    /// An `Option` containing a mutable reference to the range `Vec<Range>`. Returns `None` if
    /// the column is not found or if it contains non-range data.
    pub fn get_column_range_mut(&mut self, name: &str) -> Option<&mut Vec<Range>> {
        let r = self.get_column_mut(name)?;
        r.get_range_mut()
    }

    /// Creates a subset of the `Meta` struct by slicing rows between specified indices.
    ///
    /// This function creates a new `Meta` instance containing only the rows from `ifrom` to `ito`.
    /// It clones the specified range of rows for each column without modifying the original `Meta` instance.
    ///
    /// # Arguments
    /// - `ifrom`: The starting row index (inclusive).
    /// - `ito`: The ending row index (exclusive).
    ///
    /// # Returns
    /// A new `Meta` instance containing only the specified slice of rows for each column.
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

    /// Creates a subset of the `Meta` struct based on a list of row indices.
    ///
    /// This function returns a new `Meta` instance containing only the rows specified by `indices`.
    /// Each index in `indices` is cloned into the new `Meta` instance. This method allows for non-sequential
    /// row selections, providing more flexibility for data manipulation.
    ///
    /// # Arguments
    /// - `indices`: A slice of row indices to include in the new `Meta` instance.
    ///
    /// # Returns
    /// A new `Meta` instance containing only the rows specified by `indices`.
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

    /// Sorts the `Meta` struct by a specified column, either in ascending or descending order.
    ///
    /// This function sorts the rows of `Meta` based on the values in the specified column.
    /// The column data type must support ordering (e.g., `MetaData::IntArray`, `MetaData::FloatArray`, `MetaData::StringArray`).
    /// If the `reverse` flag is `true`, the sort is performed in descending order; otherwise, it’s ascending.
    ///
    /// # Arguments
    /// - `name`: The name of the column by which to sort the `Meta` rows.
    /// - `reverse`: A boolean flag indicating whether to sort in descending order (`true`) or ascending (`false`).
    ///
    /// # Returns
    /// A `Result` containing a sorted `Meta` instance if successful, or an error if the specified column name
    /// is not found or has an unsupported data type.
    pub fn sort(&self, name: &str, reverse: bool) -> Result<Self, Box<dyn Error>> {
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

    /// Returns an iterator over the columns in the `Meta` struct.
    ///
    /// This iterator yields pairs of column names and their associated `MetaData` values,
    /// allowing for iteration over each column in the `Meta` instance.
    ///
    /// # Returns
    /// An iterator yielding tuples of (`&String`, `&MetaData`), representing each column's name
    /// and its data, respectively.
    pub fn iter(&self) -> impl Iterator<Item = (&String, &MetaData)> {
        self.meta_name.iter().zip(self.meta_data.iter())
    }
}

/* -------------------------------------------------------------------------- */

impl PartialEq for Meta {
    fn eq(&self, other: &Self) -> bool {
        if self.num_cols() != other.num_cols() {
            return false;
        }
        if self.num_rows() != other.num_rows() {
            return false;
        }
        for j in 0..self.num_cols() {
            let name = &self.meta_name[j];

            let meta_col1 = self.meta_data[j].clone();
            let meta_col2 = other.get_column(name);

            if meta_col2.is_none() {
                return false;
            }
            if meta_col1 != *meta_col2.unwrap() {
                return false;
            }
        }
        true
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for Meta {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.pad(&format!("{}", self.format_pretty(10, false).unwrap()))
    }
}
