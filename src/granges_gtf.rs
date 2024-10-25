
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

use std::fs::{File};
use std::io::{self, BufRead, BufReader, Write};
use std::collections::HashMap;
use std::error::Error;

use flate2::read::GzDecoder;
use regex::Regex;

use crate::granges::GRanges;
use crate::meta::MetaData;
use crate::range::Range;
use crate::utility::is_gzip;

/* -------------------------------------------------------------------------- */

impl GRanges {

    fn parse_gtf_line(line: &str) -> Vec<String> {
        let re = Regex::new(r#""[^"]*"|[^ \t;]+"#).unwrap();
        re.find_iter(line)
            .map(|mat| mat.as_str().to_string())
            .collect()
    }

    fn parse_optional_fields(
        &self,
        fields  : &[String],
        type_map: &HashMap<String, String>,
        gtf_opt : &mut HashMap<String, Vec<String>>,
        gtf_def : &HashMap<String, Option<String>>,
        length  : usize,
    ) -> io::Result<()> {

        for i in (0..fields.len()).step_by(2) {
            let name      = &fields[i];
            let value_str = &fields[i + 1];
            if let Some(_) = type_map.get(name) {
                let entry = gtf_opt.entry(name.clone()).or_insert_with(Vec::new);
                entry.push(value_str.clone());
            }
        }

        for (name, values) in gtf_opt.iter_mut() {
            if values.len() < length {
                if let Some(def_val) = gtf_def.get(name) {
                    if let Some(def) = def_val {
                        values.push(def.clone());
                    } else {
                        return Err(io::Error::new(io::ErrorKind::InvalidData, format!("optional field `{}` is missing at line `{}` with no default", name, length + 1)));
                    }
                }
            }
        }

        Ok(())
    }

}

/* Read/write
 * -------------------------------------------------------------------------- */
 
impl GRanges {

    pub fn read_gtf<R: BufRead>(
        reader   : R,
        opt_names: Vec<&str>,
        opt_types: Vec<&str>,
        defaults : Vec<Option<String>>,
    ) -> Result<Self, Box<dyn Error>> {

        let mut granges  = GRanges::default();

        let mut type_map = HashMap::new();
        let mut gtf_opt  = HashMap::new();
        let mut gtf_def  = HashMap::new();

        let mut source  = vec![];
        let mut feature = vec![];
        let mut score   = vec![];
        let mut frame   = vec![];

        let mut has_score = false;
        let mut has_frame = false;

        for (i, name) in opt_names.iter().enumerate() {
            type_map.insert(String::from(*name), String::from(opt_types[i]));
            if let Some(def) = defaults.get(i) {
                gtf_def.insert(String::from(*name), def.clone());
            }
        }

        for line in reader.lines() {

            let line = line?;
            let fields: Vec<String> = Self::parse_gtf_line(&line);

            if fields.len() < 8 {
                return Err(Box::new(io::Error::new(io::ErrorKind::InvalidData, "file must have at least eight columns")));
            }

            let from = fields[3].parse::<usize>()?;
            let to   = fields[4].parse::<usize>()?;

            granges.seqnames.push(fields[0].clone());
            granges.ranges  .push(Range::new(from, to));
            granges.strand  .push(fields[6].chars().next().unwrap_or('*'));

            source .push(fields[1].clone());
            feature.push(fields[2].clone());

            if fields[5] == "." {
                score.push(0.0);
            } else {
                score.push(fields[5].parse::<f64>()?);
                has_score = true;
            }

            if fields[7] == "." {
                frame.push(-1);
            } else {
                frame.push(fields[7].parse::<i64>()?);
                has_frame = true;
            }

            let optional_fields = &fields[8..];
            granges.parse_optional_fields(optional_fields, &type_map, &mut gtf_opt, &gtf_def, granges.num_rows())?;
        }

        // Add meta data to GRanges
        granges.meta.add("source" , MetaData::StringArray(source))?;
        granges.meta.add("feature", MetaData::StringArray(feature))?;

        if has_score {
            granges.meta.add("score", MetaData::FloatArray(score))?;
        }
        if has_frame {
            granges.meta.add("frame", MetaData::IntArray(frame))?;
        }

        for (name, values) in gtf_opt {
            granges.meta.add(&name, MetaData::StringArray(values))?;
        }

        Ok(granges)
    }

    pub fn write_gtf<W: Write>(&self, writer: W) -> Result<(), Box<dyn Error>> {
        let mut w = io::BufWriter::new(writer);

        let v1 = vec![];
        let v2 = vec![];

        let source = self.meta.get_column_str("source").ok_or(
            Box::new(io::Error::new(io::ErrorKind::InvalidData, "no source column available"))
        )?;
        let feature = self.meta.get_column_str("feature").ok_or(
            Box::new(io::Error::new(io::ErrorKind::InvalidData, "no feature column available"))
        )?;
        let score = self.meta.get_column_float("score").unwrap_or(
            &v1
        );
        let frame = self.meta.get_column_int("frame").unwrap_or(
            &v2
        );

        for i in 0..self.num_rows() {

            let mut printed_tab = false;

            let s = if score.len() == 0 {
                ".".to_string()
            } else {
                score[i].to_string()
            };

            let f = if frame.len() == 0 {
                ".".to_string()
            } else {
                frame[i].to_string()
            };

            write!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.seqnames[i],
                source       [i],
                feature      [i],
                self.ranges  [i].from,
                self.ranges  [i].to,
                s,
                self.strand  [i],
                f,
            )?;

            for (name, item) in self.meta.iter() {
                if name == "source" || name == "feature" || name == "score" || name == "frame" {
                    continue;
                }
                if printed_tab {
                    write!(w, " ")?;
                } else {
                    write!(w, "\t")?;
                    printed_tab = true;
                }

                write!(w, "{} ", name)?;

                match item {
                    MetaData::FloatArray(values) => {
                        write!(w, "\"{}\"", values[i])?;
                    },
                    MetaData::IntArray(values) => {
                        write!(w, "\"{}\"", values[i])?;
                    },
                    MetaData::StringArray(values) => {
                        write!(w, "\"{}\"", values[i])?;
                    },
                    MetaData::FloatMatrix(values) => {
                        write!(w, "\"")?;
                        for (j, v) in values[i].iter().enumerate() {
                            if j != 0 {
                                write!(w, " ")?;
                            }
                            write!(w, "{}", v)?;
                        }
                        write!(w, "\"")?;
                    },
                    MetaData::IntMatrix(values) => {
                        write!(w, "\"")?;
                        for (j, v) in values[i].iter().enumerate() {
                            if j != 0 {
                                write!(w, " ")?;
                            }
                            write!(w, "{}", v)?;
                        }
                        write!(w, "\"")?;
                    },
                    MetaData::StringMatrix(values) => {
                        write!(w, "\"")?;
                        for (j, v) in values[i].iter().enumerate() {
                            if j != 0 {
                                write!(w, " ")?;
                            }
                            write!(w, "{}", v)?;
                        }
                        write!(w, "\"")?;
                    },
                    MetaData::RangeArray(_) => {
                        todo!()
                    },
                }
                write!(w, ";")?;
            }
            writeln!(w)?;
        }

        Ok(())
    }

}


/* Import/export
 * -------------------------------------------------------------------------- */
 
impl GRanges {

    pub fn import_gtf(
        filename : &str,
        opt_names: Vec<&str>,
        opt_types: Vec<&str>,
        opt_def  : Vec<Option<String>>,
    ) -> Result<Self, Box<dyn Error>> {

        let file = File::open(filename)?;
        let reader: Box<dyn BufRead> = if is_gzip(filename) {
            let decoder = GzDecoder::new(file);
            Box::new(BufReader::new(decoder))
        } else {
            Box::new(BufReader::new(file))
        };

        GRanges::read_gtf(reader, opt_names, opt_types, opt_def)

    }

    pub fn export_gtf(&self, filename: &str) -> Result<(), Box<dyn Error>> {
        let file = File::create(filename)?;
        self.write_gtf(file)
    }

}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::fs;

    use crate::granges::GRanges;

    #[test]
    fn test_granges_gtf() {

        let granges = GRanges::import_gtf("src/granges_gtf.gtf",
            vec!["gene_id"], // Names of optional fields
            vec!["str"    ], // Types of optional fields
            vec![]
        ).unwrap();

        assert_eq!(granges.num_rows(), 2);
        assert_eq!(granges.ranges[0].from, 11869);
        assert_eq!(granges.ranges[0].to  , 14409);
        assert_eq!(granges.ranges[1].from, 11870);
        assert_eq!(granges.ranges[1].to  , 14410);

        let source  = granges.meta.get_column_str("source");
        let feature = granges.meta.get_column_str("feature");
        let score   = granges.meta.get_column_float("score");
        let frame   = granges.meta.get_column_int("frame");

        assert!(source.is_some());
        assert!(feature.is_some());
        assert!(score.is_none());
        assert!(frame.is_none());

        assert_eq!(source.unwrap()[0], "transcribed_unprocessed_pseudogene");
        assert_eq!(source.unwrap()[1], "processed_transcript");

        assert_eq!(feature.unwrap()[0], "gene");
        assert_eq!(feature.unwrap()[1], "transcript");

        assert!(granges.export_gtf("src/granges_gtf.gtf.tmp").is_ok());

        assert!(fs::remove_file("src/granges_gtf.gtf.tmp").is_ok());

    }

}
