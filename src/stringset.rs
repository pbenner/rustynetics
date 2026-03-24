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
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::ops::{Deref, DerefMut, Index};

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::granges::GRanges;
use crate::range::Range;
use crate::utility::is_gzip;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StringSet {
    pub sequences: HashMap<String, Vec<u8>>,
}

/* -------------------------------------------------------------------------- */

impl StringSet {
    pub fn new(seqnames: Vec<String>, sequences: Vec<Vec<u8>>) -> Self {
        if seqnames.len() != sequences.len() {
            panic!("StringSet::new(): invalid parameters");
        }

        let mut data = HashMap::with_capacity(seqnames.len());
        for (name, sequence) in seqnames.into_iter().zip(sequences.into_iter()) {
            data.insert(name, sequence);
        }
        Self { sequences: data }
    }

    pub fn empty() -> Self {
        Self::default()
    }

    pub fn get_slice(&self, name: &str, range: Range) -> Result<&[u8], String> {
        let sequence = self
            .sequences
            .get(name)
            .ok_or_else(|| "StringSet::get_slice(): invalid sequence name".to_string())?;
        let from = range.from.min(sequence.len());
        let to = range.to.min(sequence.len());
        if from != range.from || to != range.to {
            return Err(format!(
                "StringSet::get_slice(): range `({}, {})` out of bounds for sequence `{}`",
                range.from, range.to, name
            ));
        } else {
            Ok(&sequence[from..to])
        }
    }

    pub fn scan(&self, query: &[u8]) -> GRanges {
        if self.sequences.is_empty() || query.is_empty() {
            return GRanges::default();
        }

        let mut seqnames = Vec::new();
        let mut from = Vec::new();
        let mut to = Vec::new();
        let strand = Vec::new();

        for (name, sequence) in &self.sequences {
            extend_scan_hits(&mut seqnames, &mut from, &mut to, name, sequence, query);
        }

        GRanges::new(seqnames, from, to, strand)
    }

    pub fn read_fasta<R: Read>(&mut self, reader: R) -> io::Result<()> {
        let mut reader = BufReader::new(reader);
        let mut line = String::new();
        let mut name = String::new();
        let mut sequence = Vec::new();

        loop {
            line.clear();
            let bytes = reader.read_line(&mut line)?;
            if bytes == 0 {
                break;
            }
            let line = line.trim_end_matches(&['\r', '\n'][..]);
            if line.is_empty() {
                continue;
            }
            if let Some(header) = line.strip_prefix('>') {
                if !name.is_empty() {
                    self.sequences
                        .insert(name.clone(), std::mem::take(&mut sequence));
                }
                name = parse_fasta_name(header).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "StringSet::read_fasta(): invalid FASTA file",
                    )
                })?;
            } else {
                if name.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "StringSet::read_fasta(): invalid FASTA file",
                    ));
                }
                sequence.extend_from_slice(line.as_bytes());
            }
        }

        if !name.is_empty() {
            self.sequences.insert(name, sequence);
        }

        Ok(())
    }

    pub fn import_fasta(&mut self, filename: &str) -> io::Result<()> {
        let file = File::open(filename)?;
        if is_gzip(filename) {
            self.read_fasta(BufReader::new(GzDecoder::new(file)))
        } else {
            self.read_fasta(BufReader::new(file))
        }
    }

    pub fn write_fasta<W: Write>(&self, writer: W) -> io::Result<()> {
        let mut writer = BufWriter::new(writer);
        for (name, sequence) in &self.sequences {
            writeln!(writer, ">{name}")?;
            for chunk in sequence.chunks(80) {
                writer.write_all(chunk)?;
                writer.write_all(b"\n")?;
            }
        }
        writer.flush()
    }

    pub fn export_fasta(&self, filename: &str, compress: bool) -> io::Result<()> {
        let file = File::create(filename)?;
        if compress {
            self.write_fasta(GzEncoder::new(file, Compression::default()))
        } else {
            self.write_fasta(file)
        }
    }
}

/* -------------------------------------------------------------------------- */

impl Deref for StringSet {
    type Target = HashMap<String, Vec<u8>>;

    fn deref(&self) -> &Self::Target {
        &self.sequences
    }
}

impl DerefMut for StringSet {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.sequences
    }
}

impl Index<&str> for StringSet {
    type Output = Vec<u8>;

    fn index(&self, index: &str) -> &Self::Output {
        &self.sequences[index]
    }
}

/* -------------------------------------------------------------------------- */

pub(crate) fn parse_fasta_name(header: &str) -> Option<String> {
    header
        .split(|c: char| c.is_whitespace() || c == '|')
        .find(|field| !field.is_empty())
        .map(|field| field.to_string())
}

pub(crate) fn extend_scan_hits(
    seqnames: &mut Vec<String>,
    from: &mut Vec<usize>,
    to: &mut Vec<usize>,
    name: &str,
    sequence: &[u8],
    query: &[u8],
) {
    if query.len() > sequence.len() {
        return;
    }

    for i in 0..=sequence.len() - query.len() {
        if sequence[i..i + query.len()].eq_ignore_ascii_case(query) {
            seqnames.push(name.to_string());
            from.push(i);
            to.push(i + query.len());
        }
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::path::PathBuf;

    use crate::range::Range;
    use crate::stringset::StringSet;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("gonetics")
            .join(name)
    }

    #[test]
    fn import_fasta_matches_go_fixture_lengths() {
        let mut sequences = StringSet::empty();
        sequences
            .import_fasta(fixture("stringset_test.fa").to_str().unwrap())
            .unwrap();

        assert_eq!(sequences["chr1"].len(), 2011);
        assert_eq!(sequences["chr2"].len(), 1582);
    }

    #[test]
    fn get_slice_reports_out_of_bounds() {
        let sequences = StringSet::new(vec!["chr1".into()], vec![b"acgt".to_vec()]);
        let error = sequences.get_slice("chr1", Range::new(2, 10)).unwrap_err();
        assert!(error.contains("out of bounds"));
    }

    #[test]
    fn scan_finds_case_insensitive_hits() {
        let sequences = StringSet::new(
            vec!["chr1".into(), "chr2".into()],
            vec![b"AcGtacgt".to_vec(), b"tttt".to_vec()],
        );
        let hits = sequences.scan(b"acg");

        assert_eq!(hits.seqnames, vec!["chr1".to_string(), "chr1".to_string()]);
        assert_eq!(hits.ranges[0].from, 0);
        assert_eq!(hits.ranges[1].from, 4);
    }
}
