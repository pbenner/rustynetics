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
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::granges::GRanges;
use crate::range::Range;
use crate::stringset::{extend_scan_hits, parse_fasta_name, StringSet};
use crate::utility::is_gzip;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct OrderedStringSet {
    pub sequences: StringSet,
    pub seqnames: Vec<String>,
}

/* -------------------------------------------------------------------------- */

impl OrderedStringSet {
    pub fn new(seqnames: Vec<String>, sequences: Vec<Vec<u8>>) -> Self {
        if seqnames.len() != sequences.len() {
            panic!("OrderedStringSet::new(): invalid parameters");
        }

        let mut data = StringSet::empty();
        for (name, sequence) in seqnames.iter().cloned().zip(sequences.into_iter()) {
            if data.contains_key(&name) {
                panic!("OrderedStringSet::new(): duplicate sequence name `{name}`");
            }
            data.insert(name, sequence);
        }

        Self {
            sequences: data,
            seqnames,
        }
    }

    pub fn empty() -> Self {
        Self::default()
    }

    pub fn get_slice(&self, name: &str, range: Range) -> Result<&[u8], String> {
        self.sequences.get_slice(name, range)
    }

    pub fn scan(&self, query: &[u8]) -> GRanges {
        if self.sequences.is_empty() || query.is_empty() {
            return GRanges::default();
        }

        let mut seqnames = Vec::new();
        let mut from = Vec::new();
        let mut to = Vec::new();
        let strand = Vec::new();

        for name in &self.seqnames {
            if let Some(sequence) = self.sequences.get(name) {
                extend_scan_hits(&mut seqnames, &mut from, &mut to, name, sequence, query);
            }
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
                    if self.sequences.contains_key(&name) {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("sequence name `{name}` occurred multiple times"),
                        ));
                    }
                    self.sequences
                        .insert(name.clone(), std::mem::take(&mut sequence));
                    self.seqnames.push(name.clone());
                }
                name = parse_fasta_name(header).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "OrderedStringSet::read_fasta(): invalid FASTA file",
                    )
                })?;
            } else {
                if name.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "OrderedStringSet::read_fasta(): invalid FASTA file",
                    ));
                }
                sequence.extend_from_slice(line.as_bytes());
            }
        }

        if !name.is_empty() {
            if self.sequences.contains_key(&name) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("sequence name `{name}` occurred multiple times"),
                ));
            }
            self.sequences.insert(name.clone(), sequence);
            self.seqnames.push(name);
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
        for name in &self.seqnames {
            let sequence = self.sequences.get(name).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("OrderedStringSet::write_fasta(): missing sequence `{name}`"),
                )
            })?;
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

#[cfg(test)]
mod tests {

    use std::io::Cursor;

    use crate::orderedstringset::OrderedStringSet;

    #[test]
    fn read_fasta_preserves_order() {
        let mut sequences = OrderedStringSet::empty();
        sequences
            .read_fasta(Cursor::new(">chr2\nac\n>chr1\ngt\n"))
            .unwrap();

        assert_eq!(sequences.seqnames, vec!["chr2", "chr1"]);
        assert_eq!(sequences.sequences["chr1"], b"gt".to_vec());
    }

    #[test]
    fn duplicate_names_error() {
        let mut sequences = OrderedStringSet::empty();
        let error = sequences
            .read_fasta(Cursor::new(">chr1\nac\n>chr1\ngt\n"))
            .unwrap_err();
        assert!(error.to_string().contains("occurred multiple times"));
    }
}
