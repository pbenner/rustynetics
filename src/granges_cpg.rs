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

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, Read};

use flate2::read::GzDecoder;
use mysql::from_row;
use mysql::prelude::*;
use mysql::Pool;

use crate::cpg::observed_over_expected_cpg as score_cpg;
use crate::granges::GRanges;
use crate::meta::MetaData;
use crate::utility::is_gzip;

/* -------------------------------------------------------------------------- */

impl GRanges {
    pub fn import_cpg_islands_from_ucsc(genome: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let url = format!("mysql://genome@genome-mysql.cse.ucsc.edu:3306/{}", genome);
        let pool = Pool::new(url.as_str())?;
        let mut conn = pool.get_conn()?;
        let query = "SELECT chrom, chromStart, chromEnd, length, cpgNum, gcNum, perCpg, perGc, obsExp FROM cpgIslandExt";

        let mut seqnames = Vec::new();
        let mut from = Vec::new();
        let mut to = Vec::new();
        let mut length = Vec::new();
        let mut cpg_num = Vec::new();
        let mut gc_num = Vec::new();
        let mut per_cpg = Vec::new();
        let mut per_gc = Vec::new();
        let mut obs_exp = Vec::new();

        let mut result = conn.query_iter(query)?;

        while let Some(result_set) = result.iter() {
            for row in result_set {
                let r: (String, i32, i32, i32, i32, i32, f64, f64, f64) = from_row(row.unwrap());

                seqnames.push(r.0);
                from.push(r.1 as usize);
                to.push(r.2 as usize);
                length.push(r.3 as i64);
                cpg_num.push(r.4 as i64);
                gc_num.push(r.5 as i64);
                per_cpg.push(r.6);
                per_gc.push(r.7);
                obs_exp.push(r.8);
            }
        }

        let mut granges = GRanges::new(seqnames, from, to, Vec::new());
        granges
            .meta
            .add("length", MetaData::IntArray(length))
            .unwrap();
        granges
            .meta
            .add("cpgNum", MetaData::IntArray(cpg_num))
            .unwrap();
        granges
            .meta
            .add("gcNum", MetaData::IntArray(gc_num))
            .unwrap();
        granges
            .meta
            .add("perCpg", MetaData::FloatArray(per_cpg))
            .unwrap();
        granges
            .meta
            .add("perGc", MetaData::FloatArray(per_gc))
            .unwrap();
        granges
            .meta
            .add("obsExp", MetaData::FloatArray(obs_exp))
            .unwrap();

        Ok(granges)
    }

    pub fn read_cpg_islands_table<R: Read>(reader: R) -> Result<Self, Box<dyn std::error::Error>> {
        let mut cpg = GRanges::default();
        cpg.read_table(
            reader,
            &["length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp"],
            &["Int", "Int", "Int", "Float", "Float", "Float"],
        )?;
        Ok(cpg)
    }

    pub fn import_cpg_islands_table(filename: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(filename)?;
        let reader: Box<dyn Read> = if is_gzip(filename) {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        Self::read_cpg_islands_table(reader)
    }

    pub fn observed_over_expected_cpg<B: AsRef<[u8]>>(
        &self,
        genomic_sequences: &HashMap<String, B>,
    ) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
        let mut result = Vec::with_capacity(self.num_rows());

        for i in 0..self.num_rows() {
            let sequence = genomic_sequences.get(&self.seqnames[i]).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("sequence `{}` not found", self.seqnames[i]),
                )
            })?;
            let sequence = sequence.as_ref();

            if self.ranges[i].to > sequence.len() {
                return Err(Box::new(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "range `({}, {})` out of bounds for sequence `{}`",
                        self.ranges[i].from, self.ranges[i].to, self.seqnames[i]
                    ),
                )));
            }

            result.push(score_cpg(&sequence[self.ranges[i].from..self.ranges[i].to]));
        }

        Ok(result)
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::collections::HashMap;
    use std::io::Cursor;

    use crate::granges::GRanges;
    use crate::meta::MetaData;

    #[test]
    fn read_cpg_islands_table_roundtrips_expected_columns() {
        let mut granges = GRanges::new(vec!["chr1".into()], vec![10], vec![20], Vec::new());
        granges
            .meta
            .add("length", MetaData::IntArray(vec![10]))
            .unwrap();
        granges
            .meta
            .add("cpgNum", MetaData::IntArray(vec![2]))
            .unwrap();
        granges
            .meta
            .add("gcNum", MetaData::IntArray(vec![4]))
            .unwrap();
        granges
            .meta
            .add("perCpg", MetaData::FloatArray(vec![20.0]))
            .unwrap();
        granges
            .meta
            .add("perGc", MetaData::FloatArray(vec![40.0]))
            .unwrap();
        granges
            .meta
            .add("obsExp", MetaData::FloatArray(vec![1.5]))
            .unwrap();

        let mut buffer = Vec::new();
        granges.write_table(&mut buffer, &[]).unwrap();

        let parsed = GRanges::read_cpg_islands_table(Cursor::new(buffer)).unwrap();
        assert_eq!(parsed.seqnames, vec!["chr1"]);
        assert_eq!(parsed.meta.get_column_int("length").unwrap(), &vec![10]);
        assert_eq!(parsed.meta.get_column_float("obsExp").unwrap(), &vec![1.5]);
    }

    #[test]
    fn observed_over_expected_cpg_scores_each_range() {
        let granges = GRanges::new(
            vec!["chr1".into(), "chr1".into()],
            vec![0, 2],
            vec![4, 6],
            Vec::new(),
        );
        let mut sequences = HashMap::new();
        sequences.insert("chr1".to_string(), b"ACGCGT".to_vec());

        let scores = granges.observed_over_expected_cpg(&sequences).unwrap();

        assert_eq!(scores, vec![2.0, 2.0]);
    }
}
