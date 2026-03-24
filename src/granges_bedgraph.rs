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
use std::io::{self, BufRead, BufReader, Write};
use std::str::FromStr;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::granges::GRanges;
use crate::meta::MetaData;
use crate::utility::is_gzip;

/* -------------------------------------------------------------------------- */

impl GRanges {
    pub fn write_bedgraph<W: Write>(
        &self,
        writer: &mut W,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let values = self
            .meta
            .get_column_float("values")
            .ok_or("values column required for bedGraph export")?;

        for i in 0..self.num_rows() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                self.seqnames[i], self.ranges[i].from, self.ranges[i].to, values[i]
            )?;
        }

        Ok(())
    }

    pub fn export_bedgraph(&self, filename: &str, compress: bool) -> io::Result<()> {
        let file = File::create(filename)?;
        let mut writer: Box<dyn Write> = if compress {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };

        self.write_bedgraph(&mut writer)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))
    }

    pub fn read_bedgraph<R: BufRead>(&mut self, reader: R) -> io::Result<()> {
        let mut seqnames = Vec::new();
        let mut from = Vec::new();
        let mut to = Vec::new();
        let mut values = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split_whitespace().collect();

            if fields.is_empty() {
                continue;
            }
            if fields.len() != 4 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "BedGraph file must have four columns!",
                ));
            }

            seqnames.push(fields[0].to_string());
            from.push(usize::from_str(fields[1]).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidInput, "Invalid integer in column 2")
            })?);
            to.push(usize::from_str(fields[2]).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidInput, "Invalid integer in column 3")
            })?);
            values.push(f64::from_str(fields[3]).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidInput, "Invalid float in column 4")
            })?);
        }

        let mut granges = GRanges::new(seqnames, from, to, Vec::new());
        granges
            .meta
            .add("values", MetaData::FloatArray(values))
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
        *self = granges;

        Ok(())
    }

    pub fn import_bedgraph(&mut self, filename: &str) -> io::Result<()> {
        let file = File::open(filename)?;
        let reader: Box<dyn BufRead> = if is_gzip(filename) {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        self.read_bedgraph(reader)
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::io::{BufReader, Cursor};

    use crate::granges::GRanges;

    #[test]
    fn granges_bedgraph_roundtrip() {
        let data = b"chr1\t10\t20\t1.5\nchr2\t30\t50\t2\n";
        let mut granges = GRanges::default();

        granges
            .read_bedgraph(BufReader::new(Cursor::new(&data[..])))
            .unwrap();

        assert_eq!(granges.seqnames, vec!["chr1", "chr2"]);
        assert_eq!(granges.ranges[0].from, 10);
        assert_eq!(granges.ranges[1].to, 50);
        assert_eq!(
            granges.meta.get_column_float("values").unwrap(),
            &vec![1.5, 2.0]
        );

        let mut buffer = Vec::new();
        granges.write_bedgraph(&mut buffer).unwrap();
        assert_eq!(
            String::from_utf8(buffer).unwrap(),
            "chr1\t10\t20\t1.5\nchr2\t30\t50\t2\n"
        );
    }
}
