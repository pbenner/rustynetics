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

use std::fs::File;
use std::io::{self, Read, Seek};

use crate::bigwig::BigWigReader;
use crate::granges::GRanges;
use crate::meta::MetaData;
use crate::track_statistics::BinSummaryStatistics;

/* -------------------------------------------------------------------------- */

impl GRanges {
    pub fn read_bigwig<R: Read + Seek>(
        &mut self,
        reader: &mut R,
        name: &str,
        f: BinSummaryStatistics,
        mut bin_size: usize,
        bin_overlap: usize,
        init: f64,
        rev_neg_strand: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut bwr = BigWigReader::new(reader)?;
        let mut counts = Vec::with_capacity(self.num_rows());

        for i in 0..self.num_rows() {
            let (mut values, detected_bin_size) = bwr.query_slice(
                &self.seqnames[i],
                self.ranges[i].from,
                self.ranges[i].to,
                f,
                bin_size,
                bin_overlap,
                init,
            )?;

            if bin_size == 0 {
                bin_size = detected_bin_size;
            }
            if rev_neg_strand && self.strand[i] == '-' {
                values.reverse();
            }
            counts.push(values);
        }

        self.meta.delete_meta(name);
        self.meta.add(name, MetaData::FloatMatrix(counts))?;

        Ok(())
    }

    pub fn import_bigwig(
        &mut self,
        filename: &str,
        name: &str,
        f: BinSummaryStatistics,
        bin_size: usize,
        bin_overlap: usize,
        init: f64,
        rev_neg_strand: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::open(filename)?;
        self.read_bigwig(
            &mut file,
            name,
            f,
            bin_size,
            bin_overlap,
            init,
            rev_neg_strand,
        )
        .map_err(|e| {
            Box::new(io::Error::new(
                io::ErrorKind::Other,
                format!("importing bigWig file from `{}` failed: {}", filename, e),
            )) as Box<dyn std::error::Error>
        })
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::granges::GRanges;
    use crate::meta::MetaData;
    use crate::track_statistics::bin_summary_statistics_from_string;

    #[test]
    fn granges_import_bigwig_reads_values_into_float_matrix() {
        let mut granges = GRanges::new(
            vec!["chrY".into(), "chrY".into()],
            vec![1838100, 1838100],
            vec![1838600, 1838300],
            vec!['+', '-'],
        );

        granges
            .import_bigwig(
                "tests/test_bigwig_2.bw",
                "signal",
                bin_summary_statistics_from_string("mean").unwrap(),
                100,
                0,
                0.0,
                true,
            )
            .unwrap();

        match granges.meta.get_column("signal").unwrap() {
            MetaData::FloatMatrix(v) => {
                assert_eq!(v[0], vec![1.0, 1.0, 0.0, 0.0, 0.0]);
                assert_eq!(v[1], vec![1.0, 1.0]);
            }
            _ => panic!("expected FloatMatrix"),
        }
    }
}
