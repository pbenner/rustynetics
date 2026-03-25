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

use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::ops::{Deref, DerefMut};

use flate2::read::GzDecoder;

use crate::granges::GRanges;
use crate::meta::MetaData;
use crate::utility::is_gzip;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, PartialEq)]
pub struct GPeaks {
    pub granges: GRanges,
}

/* -------------------------------------------------------------------------- */

impl GPeaks {
    pub fn new(
        seqnames: Vec<String>,
        from: Vec<usize>,
        to: Vec<usize>,
        abs_summit: Vec<i64>,
        pileup: Vec<f64>,
        pvalue: Vec<f64>,
        fold_enrichment: Vec<f64>,
        qvalue: Vec<f64>,
    ) -> Self {
        let mut granges = GRanges::new(seqnames, from, to, Vec::new());
        granges
            .meta
            .add("abs_summit", MetaData::IntArray(abs_summit))
            .unwrap();
        granges
            .meta
            .add("pileup", MetaData::FloatArray(pileup))
            .unwrap();
        granges
            .meta
            .add("-log10(pvalue)", MetaData::FloatArray(pvalue))
            .unwrap();
        granges
            .meta
            .add("fold_enrichment", MetaData::FloatArray(fold_enrichment))
            .unwrap();
        granges
            .meta
            .add("-log10(qvalue)", MetaData::FloatArray(qvalue))
            .unwrap();
        Self { granges }
    }

    pub fn abs_summit(&self) -> &Vec<i64> {
        self.granges.meta.get_column_int("abs_summit").unwrap()
    }

    pub fn pileup(&self) -> &Vec<f64> {
        self.granges.meta.get_column_float("pileup").unwrap()
    }

    pub fn pvalue(&self) -> &Vec<f64> {
        self.granges
            .meta
            .get_column_float("-log10(pvalue)")
            .unwrap()
    }

    pub fn fold_enrichment(&self) -> &Vec<f64> {
        self.granges
            .meta
            .get_column_float("fold_enrichment")
            .unwrap()
    }

    pub fn qvalue(&self) -> &Vec<f64> {
        self.granges
            .meta
            .get_column_float("-log10(qvalue)")
            .unwrap()
    }

    pub fn read_xls<R: Read>(reader: R) -> Result<Self, Box<dyn Error>> {
        let reader = BufReader::new(reader);

        let mut saw_header = false;
        let mut seqnames = Vec::new();
        let mut from = Vec::new();
        let mut to = Vec::new();
        let mut abs_summit = Vec::new();
        let mut pileup = Vec::new();
        let mut pvalue = Vec::new();
        let mut fold_enrichment = Vec::new();
        let mut qvalue = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.is_empty() {
                continue;
            }
            if fields[0] == "#" {
                continue;
            }
            if fields.len() != 10 {
                return Err(Box::new(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid peaks xls file",
                )));
            }

            if !saw_header {
                if fields[0] != "chr"
                    || fields[1] != "start"
                    || fields[2] != "end"
                    || fields[4] != "abs_summit"
                    || fields[5] != "pileup"
                    || fields[6] != "-log10(pvalue)"
                    || fields[7] != "fold_enrichment"
                    || fields[8] != "-log10(qvalue)"
                {
                    return Err(Box::new(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "invalid Xls header",
                    )));
                }
                saw_header = true;
                continue;
            }

            seqnames.push(fields[0].to_string());
            from.push(fields[1].parse()?);
            to.push(fields[2].parse::<usize>()? + 1);
            abs_summit.push(fields[4].parse()?);
            pileup.push(fields[5].parse()?);
            pvalue.push(fields[6].parse()?);
            fold_enrichment.push(fields[7].parse()?);
            qvalue.push(fields[8].parse()?);
        }

        Ok(Self::new(
            seqnames,
            from,
            to,
            abs_summit,
            pileup,
            pvalue,
            fold_enrichment,
            qvalue,
        ))
    }

    pub fn import_xls(filename: &str) -> Result<Self, Box<dyn Error>> {
        let file = File::open(filename)?;
        if is_gzip(filename) {
            Self::read_xls(BufReader::new(GzDecoder::new(file)))
        } else {
            Self::read_xls(BufReader::new(file))
        }
    }
}

/* -------------------------------------------------------------------------- */

impl Deref for GPeaks {
    type Target = GRanges;

    fn deref(&self) -> &Self::Target {
        &self.granges
    }
}

impl DerefMut for GPeaks {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.granges
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::io::Cursor;
    use std::path::PathBuf;

    use crate::peaks::GPeaks;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join(name)
    }

    #[test]
    fn import_xls_matches_go_fixture() {
        let peaks = GPeaks::import_xls(fixture("peaks_test.xls").to_str().unwrap()).unwrap();

        assert_eq!(peaks.num_rows(), 16);
        assert_eq!(peaks.abs_summit()[0], 5865);
        assert_eq!(peaks.pileup()[0], 33.0);
    }

    #[test]
    fn invalid_header_is_rejected() {
        let err = GPeaks::read_xls(Cursor::new(
            "chr start end length summit pileup -log10(pvalue) fold_enrichment -log10(qvalue) name\n",
        ))
            .unwrap_err();
        assert!(err.to_string().contains("invalid Xls header"));
    }
}
