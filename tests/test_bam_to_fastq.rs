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

#[cfg(test)]
mod tests {

    use std::env;
    use std::fs::{remove_file, File};
    use std::io::{BufReader, Cursor};
    use std::path::PathBuf;
    use std::process::Command;
    use std::time::{SystemTime, UNIX_EPOCH};

    use rustynetics::bam::{BamReader, BamReaderOptions};
    use rustynetics::fastq::{FastqReader, FastqRecord};

    #[derive(Default)]
    struct ExpectedCounts {
        total: usize,
        read1: usize,
        read2: usize,
        single: usize,
    }

    fn temp_path(filename: &str) -> PathBuf {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();

        env::temp_dir().join(format!("rustynetics-{}-{}", timestamp, filename))
    }

    fn read_fastq_records_from_bytes(bytes: &[u8]) -> Vec<FastqRecord> {
        let mut reader = FastqReader::new(BufReader::new(Cursor::new(bytes)));
        let mut records = Vec::new();

        while let Some(record) = reader.read_record().unwrap() {
            records.push(record);
        }

        records
    }

    fn read_fastq_records_from_file(filename: &PathBuf) -> Vec<FastqRecord> {
        let file = File::open(filename).unwrap();
        let mut reader = FastqReader::new(BufReader::new(file));
        let mut records = Vec::new();

        while let Some(record) = reader.read_record().unwrap() {
            records.push(record);
        }

        records
    }

    fn expected_counts(filename: &str) -> ExpectedCounts {
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let options = BamReaderOptions {
            read_name: true,
            read_cigar: false,
            read_sequence: false,
            read_auxiliary: false,
            read_qual: false,
        };
        let mut bam_reader = BamReader::new(reader, Some(options)).unwrap();
        let mut counts = ExpectedCounts::default();

        for result in bam_reader.read_single_end() {
            let block = result.unwrap().block;

            if block.flag.secondary_alignment()
                || block.flag.supplementary_alignment()
                || block.flag.not_passing_filters()
            {
                continue;
            }

            counts.total += 1;

            if block.flag.read_paired() && block.flag.first_in_pair() {
                counts.read1 += 1;
            } else if block.flag.read_paired() && block.flag.second_in_pair() {
                counts.read2 += 1;
            } else {
                counts.single += 1;
            }
        }

        counts
    }

    #[test]
    fn bam_to_fastq_writes_fastq_records_to_stdout() {
        let output = Command::new(env!("CARGO_BIN_EXE_bam-to-fastq"))
            .arg("tests/test_bam_1.bam")
            .arg("--pair-suffixes")
            .arg("--fill-missing-quality")
            .arg("I")
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "stderr: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        let records = read_fastq_records_from_bytes(&output.stdout);
        let counts = expected_counts("tests/test_bam_1.bam");

        assert_eq!(records.len(), counts.total);

        let mate1 = records
            .iter()
            .find(|record| record.header == "r001/1")
            .unwrap();
        assert_eq!(mate1.sequence, "ATGGCGCTG");
        assert_eq!(mate1.sequence.len(), mate1.qualities.len());
    }

    #[test]
    fn bam_to_fastq_can_split_paired_and_single_reads() {
        let read1 = temp_path("read1.fastq");
        let read2 = temp_path("read2.fastq");
        let single = temp_path("single.fastq");

        let output = Command::new(env!("CARGO_BIN_EXE_bam-to-fastq"))
            .arg("tests/test_bam_1.bam")
            .arg("--output1")
            .arg(&read1)
            .arg("--output2")
            .arg(&read2)
            .arg("--output-single")
            .arg(&single)
            .arg("--pair-suffixes")
            .arg("--fill-missing-quality")
            .arg("I")
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "stderr: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        let counts = expected_counts("tests/test_bam_1.bam");
        let read1_records = read_fastq_records_from_file(&read1);
        let read2_records = read_fastq_records_from_file(&read2);
        let single_records = read_fastq_records_from_file(&single);

        let _ = remove_file(&read1);
        let _ = remove_file(&read2);
        let _ = remove_file(&single);

        assert_eq!(read1_records.len(), counts.read1);
        assert_eq!(read2_records.len(), counts.read2);
        assert_eq!(single_records.len(), counts.single);
        assert!(read1_records
            .iter()
            .all(|record| record.header.ends_with("/1")));
        assert!(read2_records
            .iter()
            .all(|record| record.header.ends_with("/2")));
    }

    #[test]
    fn bam_to_fastq_reports_missing_qualities_without_fill_option() {
        let output = Command::new(env!("CARGO_BIN_EXE_bam-to-fastq"))
            .arg("tests/test_bam_1.bam")
            .output()
            .unwrap();

        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("has no quality scores"));
    }
}
