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

use std::io::{self, BufRead, Write};

use crate::bam::BamBlock;

/* -------------------------------------------------------------------------- */

/// Represents a single four-line FASTQ record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastqRecord {
    pub header: String,
    pub sequence: String,
    pub qualities: String,
}

/* -------------------------------------------------------------------------- */

impl FastqRecord {
    /// Creates a validated FASTQ record.
    pub fn new(header: String, sequence: String, qualities: String) -> io::Result<Self> {
        if header.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ header must not be empty",
            ));
        }
        if sequence.len() != qualities.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "FASTQ sequence and quality lengths differ ({} != {})",
                    sequence.len(),
                    qualities.len()
                ),
            ));
        }

        Ok(Self {
            header,
            sequence,
            qualities,
        })
    }

    /// Returns the first whitespace-delimited token from the FASTQ header.
    pub fn name(&self) -> &str {
        self.header.split_whitespace().next().unwrap_or("")
    }

    /// Creates a FASTQ record from a BAM alignment record.
    ///
    /// Reverse-strand BAM records are reverse-complemented and their qualities
    /// are reversed so the emitted FASTQ matches the original read orientation.
    pub fn from_bam_block(
        block: &BamBlock,
        append_pair_suffix: bool,
        fill_missing_quality: Option<char>,
    ) -> io::Result<Self> {
        if block.read_name.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "BAM record has an empty read name",
            ));
        }

        let mut header = block.read_name.clone();
        if append_pair_suffix {
            if block.flag.first_in_pair() {
                header.push_str("/1");
            } else if block.flag.second_in_pair() {
                header.push_str("/2");
            }
        }

        let mut sequence = block.seq.to_string();
        if block.flag.reverse_strand() {
            sequence = reverse_complement(&sequence);
        }

        let qualities = qualities_from_bam(block, sequence.len(), fill_missing_quality)?;

        FastqRecord::new(header, sequence, qualities)
    }
}

/* -------------------------------------------------------------------------- */

/// Reads FASTQ records from any buffered reader.
pub struct FastqReader<R: BufRead> {
    reader: R,
    line_number: usize,
}

/* -------------------------------------------------------------------------- */

impl<R: BufRead> FastqReader<R> {
    /// Creates a new FASTQ reader.
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_number: 0,
        }
    }

    fn read_line(&mut self, line: &mut String) -> io::Result<usize> {
        line.clear();
        let bytes = self.reader.read_line(line)?;
        if bytes > 0 {
            self.line_number += 1;
        }
        Ok(bytes)
    }

    /// Reads the next FASTQ record.
    pub fn read_record(&mut self) -> io::Result<Option<FastqRecord>> {
        let mut header = String::new();
        let mut sequence = String::new();
        let mut separator = String::new();
        let mut qualities = String::new();

        if self.read_line(&mut header)? == 0 {
            return Ok(None);
        }

        if !trim_line(&header).starts_with('@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid FASTQ header at line {}", self.line_number),
            ));
        }

        if self.read_line(&mut sequence)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "unexpected end of FASTQ after header at line {}",
                    self.line_number
                ),
            ));
        }

        if self.read_line(&mut separator)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "unexpected end of FASTQ after sequence at line {}",
                    self.line_number
                ),
            ));
        }
        if !trim_line(&separator).starts_with('+') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid FASTQ separator at line {}", self.line_number),
            ));
        }

        if self.read_line(&mut qualities)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "unexpected end of FASTQ after separator at line {}",
                    self.line_number
                ),
            ));
        }

        FastqRecord::new(
            trim_line(&header).trim_start_matches('@').to_string(),
            trim_line(&sequence).to_string(),
            trim_line(&qualities).to_string(),
        )
        .map(Some)
    }
}

/* -------------------------------------------------------------------------- */

/// Writes FASTQ records to any writer.
pub struct FastqWriter<W: Write> {
    writer: W,
}

/* -------------------------------------------------------------------------- */

impl<W: Write> FastqWriter<W> {
    /// Creates a new FASTQ writer.
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Writes a single FASTQ record.
    pub fn write_record(&mut self, record: &FastqRecord) -> io::Result<()> {
        writeln!(self.writer, "@{}", record.header)?;
        writeln!(self.writer, "{}", record.sequence)?;
        writeln!(self.writer, "+")?;
        writeln!(self.writer, "{}", record.qualities)?;
        Ok(())
    }

    /// Flushes the underlying writer.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

/* -------------------------------------------------------------------------- */

fn trim_line(line: &str) -> &str {
    line.trim_end_matches(&['\r', '\n'][..])
}

/* -------------------------------------------------------------------------- */

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(complement_base)
        .map(char::from)
        .collect()
}

/* -------------------------------------------------------------------------- */

fn complement_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'=' => b'=',
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'M' => b'K',
        b'R' => b'Y',
        b'S' => b'S',
        b'V' => b'B',
        b'W' => b'W',
        b'Y' => b'R',
        b'H' => b'D',
        b'K' => b'M',
        b'D' => b'H',
        b'B' => b'V',
        b'N' => b'N',
        _ => b'N',
    }
}

/* -------------------------------------------------------------------------- */

fn qualities_from_bam(
    block: &BamBlock,
    sequence_len: usize,
    fill_missing_quality: Option<char>,
) -> io::Result<String> {
    if block.qual.0.is_empty() && sequence_len == 0 {
        return Ok(String::new());
    }

    let all_missing = !block.qual.0.is_empty() && block.qual.0.iter().all(|&q| q == 0xff);
    if all_missing {
        if let Some(fill) = fill_missing_quality {
            return Ok(std::iter::repeat_n(fill, sequence_len).collect());
        }
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "BAM record `{}` has no quality scores; use --fill-missing-quality to emit FASTQ anyway",
                block.read_name
            ),
        ));
    }

    if block.qual.0.iter().any(|&q| q == 0xff) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "BAM record `{}` contains partially missing quality scores",
                block.read_name
            ),
        ));
    }

    if block.qual.0.len() != sequence_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "BAM record `{}` has inconsistent sequence and quality lengths ({} != {})",
                block.read_name,
                sequence_len,
                block.qual.0.len()
            ),
        ));
    }

    let mut qualities = block.qual.0.clone();
    if block.flag.reverse_strand() {
        qualities.reverse();
    }

    Ok(qualities.iter().map(|&q| (q + 33) as char).collect())
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::io::{BufReader, Cursor};

    use crate::bam::{BamBlock, BamFlag, BamQual, BamSeq};

    use super::{FastqReader, FastqRecord, FastqWriter};

    fn encode_base(base: u8) -> u8 {
        match base.to_ascii_uppercase() {
            b'=' => 0,
            b'A' => 1,
            b'C' => 2,
            b'M' => 3,
            b'G' => 4,
            b'R' => 5,
            b'S' => 6,
            b'V' => 7,
            b'T' => 8,
            b'W' => 9,
            b'Y' => 10,
            b'H' => 11,
            b'K' => 12,
            b'D' => 13,
            b'B' => 14,
            b'N' => 15,
            _ => 15,
        }
    }

    fn bam_seq_from_string(sequence: &str) -> BamSeq {
        let bytes = sequence.as_bytes();
        let mut encoded = Vec::with_capacity(bytes.len().div_ceil(2));
        let mut i = 0;

        while i < bytes.len() {
            let hi = encode_base(bytes[i]) << 4;
            let lo = if i + 1 < bytes.len() {
                encode_base(bytes[i + 1])
            } else {
                0
            };
            encoded.push(hi | lo);
            i += 2;
        }

        BamSeq(encoded)
    }

    fn bam_qual_from_string(qualities: &str) -> BamQual {
        BamQual(qualities.bytes().map(|q| q - 33).collect())
    }

    #[test]
    fn fastq_reader_reads_single_record() {
        let data = b"@read/1 lane:1\nACGT\n+\nIIII\n";
        let mut reader = FastqReader::new(BufReader::new(Cursor::new(&data[..])));
        let record = reader.read_record().unwrap().unwrap();

        assert_eq!(record.header, "read/1 lane:1");
        assert_eq!(record.name(), "read/1");
        assert_eq!(record.sequence, "ACGT");
        assert_eq!(record.qualities, "IIII");
        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn fastq_writer_roundtrips_record() {
        let mut buffer = Vec::new();
        let record =
            FastqRecord::new("read".to_string(), "ACGT".to_string(), "IIII".to_string()).unwrap();

        FastqWriter::new(&mut buffer).write_record(&record).unwrap();

        assert_eq!(String::from_utf8(buffer).unwrap(), "@read\nACGT\n+\nIIII\n");
    }

    #[test]
    fn bam_reverse_strand_converts_to_original_fastq_orientation() {
        let block = BamBlock {
            flag: BamFlag(0x001 | 0x080 | 0x010),
            read_name: "r001".to_string(),
            seq: bam_seq_from_string("CAGCGCCAT"),
            qual: bam_qual_from_string("ABCDEFGHI"),
            ..Default::default()
        };

        let record = FastqRecord::from_bam_block(&block, true, None).unwrap();

        assert_eq!(record.header, "r001/2");
        assert_eq!(record.sequence, "ATGGCGCTG");
        assert_eq!(record.qualities, "IHGFEDCBA");
    }

    #[test]
    fn bam_missing_quality_can_be_filled() {
        let block = BamBlock {
            read_name: "missing".to_string(),
            seq: bam_seq_from_string("ACGT"),
            qual: BamQual(vec![0xff; 4]),
            ..Default::default()
        };

        let record = FastqRecord::from_bam_block(&block, false, Some('!')).unwrap();

        assert_eq!(record.qualities, "!!!!");
    }
}
