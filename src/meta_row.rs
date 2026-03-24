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

use crate::meta::{Meta, MetaData};
use crate::range::Range;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, PartialEq)]
pub enum MetaRowValue<'a> {
    String(&'a str),
    Strings(&'a [String]),
    Float(f64),
    Floats(&'a [f64]),
    Int(i64),
    Ints(&'a [i64]),
    Range(Range),
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Copy, Debug)]
pub struct MetaRow<'a> {
    meta: &'a Meta,
    row: usize,
}

/* -------------------------------------------------------------------------- */

impl<'a> MetaRow<'a> {
    pub fn new(meta: &'a Meta, row: usize) -> Self {
        Self { meta, row }
    }

    pub fn get_meta(&self, name: &str) -> Option<MetaRowValue<'a>> {
        match self.meta.get_column(name)? {
            MetaData::StringArray(v) => Some(MetaRowValue::String(&v[self.row])),
            MetaData::StringMatrix(v) => Some(MetaRowValue::Strings(&v[self.row])),
            MetaData::FloatArray(v) => Some(MetaRowValue::Float(v[self.row])),
            MetaData::FloatMatrix(v) => Some(MetaRowValue::Floats(&v[self.row])),
            MetaData::IntArray(v) => Some(MetaRowValue::Int(v[self.row])),
            MetaData::IntMatrix(v) => Some(MetaRowValue::Ints(&v[self.row])),
            MetaData::RangeArray(v) => Some(MetaRowValue::Range(v[self.row])),
        }
    }

    pub fn get_meta_str(&self, name: &str) -> Option<&'a str> {
        match self.get_meta(name)? {
            MetaRowValue::String(v) => Some(v),
            _ => None,
        }
    }

    pub fn get_meta_strs(&self, name: &str) -> Option<&'a [String]> {
        match self.get_meta(name)? {
            MetaRowValue::Strings(v) => Some(v),
            _ => None,
        }
    }

    pub fn get_meta_float(&self, name: &str) -> Option<f64> {
        match self.get_meta(name)? {
            MetaRowValue::Float(v) => Some(v),
            _ => None,
        }
    }

    pub fn get_meta_floats(&self, name: &str) -> Option<&'a [f64]> {
        match self.get_meta(name)? {
            MetaRowValue::Floats(v) => Some(v),
            _ => None,
        }
    }

    pub fn get_meta_int(&self, name: &str) -> Option<i64> {
        match self.get_meta(name)? {
            MetaRowValue::Int(v) => Some(v),
            _ => None,
        }
    }

    pub fn get_meta_ints(&self, name: &str) -> Option<&'a [i64]> {
        match self.get_meta(name)? {
            MetaRowValue::Ints(v) => Some(v),
            _ => None,
        }
    }

    pub fn get_meta_range(&self, name: &str) -> Option<Range> {
        match self.get_meta(name)? {
            MetaRowValue::Range(v) => Some(v),
            _ => None,
        }
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::meta::{Meta, MetaData};
    use crate::meta_row::MetaRowValue;
    use crate::range::Range;

    #[test]
    fn meta_row_reads_scalar_and_matrix_values() {
        let meta = Meta::new(
            vec!["name", "score", "bins", "window"],
            vec![
                MetaData::StringArray(vec!["gene1".into(), "gene2".into()]),
                MetaData::FloatArray(vec![1.5, 2.5]),
                MetaData::IntMatrix(vec![vec![1, 2], vec![3, 4, 5]]),
                MetaData::RangeArray(vec![Range::new(1, 3), Range::new(5, 8)]),
            ],
        )
        .unwrap();

        let row = meta.row(1);

        assert_eq!(row.get_meta_str("name"), Some("gene2"));
        assert_eq!(row.get_meta_float("score"), Some(2.5));
        assert_eq!(row.get_meta_range("window"), Some(Range::new(5, 8)));
        assert_eq!(row.get_meta_ints("bins"), Some(&[3, 4, 5][..]));

        assert_eq!(
            row.get_meta("bins"),
            Some(MetaRowValue::Ints(&[3, 4, 5][..]))
        );
    }
}
