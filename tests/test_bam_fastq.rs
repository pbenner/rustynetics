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

    use std::collections::BTreeSet;
    use std::env;
    use std::fs::{remove_file, File};
    use std::io::{BufReader, BufWriter, Write};
    use std::path::{Path, PathBuf};
    use std::process::Command;
    use std::time::{SystemTime, UNIX_EPOCH};

    use rustynetics::bam::{BamReader, BamReaderOptions};

    fn temp_path(filename: &str) -> PathBuf {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();

        env::temp_dir().join(format!("rustynetics-{}-{}", timestamp, filename))
    }

    fn format_count(value: usize) -> String {
        let digits = value.to_string();
        let mut formatted = String::with_capacity(digits.len() + digits.len() / 3);

        for (index, ch) in digits.chars().enumerate() {
            if index > 0 && (digits.len() - index) % 3 == 0 {
                formatted.push(',');
            }
            formatted.push(ch);
        }

        formatted
    }

    fn collect_bam_read_names(filename: &str) -> Vec<String> {
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
        let mut names = BTreeSet::new();

        for result in bam_reader.read_single_end() {
            names.insert(result.unwrap().block.read_name);
        }

        names.into_iter().collect()
    }

    fn count_bam_unique_names(filename: &str) -> usize {
        collect_bam_read_names(filename).len()
    }

    fn count_bam_records(filename: &str) -> usize {
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
        let mut records = 0usize;

        for result in bam_reader.read_single_end() {
            result.unwrap();
            records += 1;
        }

        records
    }

    fn count_bam_records_with_name(filename: &str, expected_name: &str) -> usize {
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
        let mut records = 0usize;

        for result in bam_reader.read_single_end() {
            let block = result.unwrap().block;
            if block.read_name == expected_name {
                records += 1;
            }
        }

        records
    }

    fn write_fastq(filename: &Path, names: &[String]) {
        let file = File::create(filename).unwrap();
        let mut writer = BufWriter::new(file);

        for name in names {
            writeln!(writer, "@{}/1 lane:1", name).unwrap();
            writeln!(writer, "ACGT").unwrap();
            writeln!(writer, "+").unwrap();
            writeln!(writer, "IIII").unwrap();
        }
    }

    #[test]
    fn bam_check_fastq_accepts_matching_names() {
        let bam = "tests/test_bam_2.bam";
        let fastq = temp_path("matching.fastq");
        let names = collect_bam_read_names(bam);
        let bam_records = count_bam_records(bam);
        let bam_unique_names = count_bam_unique_names(bam);

        write_fastq(&fastq, &names);

        let output = Command::new(env!("CARGO_BIN_EXE_bam-check-fastq"))
            .arg(bam)
            .arg(&fastq)
            .output()
            .unwrap();

        let _ = remove_file(&fastq);

        assert!(output.status.success());
        let stdout = String::from_utf8_lossy(&output.stdout);

        assert!(stdout.contains("Result"));
        assert!(stdout.contains("BAM->FASTQ coverage"));
        assert!(stdout.contains("COMPLETE"));
        assert!(stdout.contains("Found BAM records"));
        assert!(stdout.contains("Missing BAM records"));
        assert!(stdout.contains("Found BAM unique names"));
        assert!(stdout.contains("Missing BAM unique names"));
        assert!(stdout.contains("BAM Record Coverage"));
        assert!(stdout.contains("BAM Unique Read Names"));
        assert!(stdout.contains("FASTQ Unique Read Names"));
        assert!(stdout.contains("BAM Record Classes"));
        assert!(stdout.contains("FASTQ File Summary"));
        assert!(stdout.contains(&format_count(bam_records)));
        assert!(stdout.contains(&format_count(bam_unique_names)));
        assert!(stdout.contains("100.00%"));
        assert!(stdout.contains("FASTQ-only unique read names"));
        assert!(stdout.contains(fastq.to_string_lossy().as_ref()));
    }

    #[test]
    fn bam_check_fastq_accepts_matching_names_across_multiple_fastqs() {
        let bam = "tests/test_bam_2.bam";
        let fastq1 = temp_path("matching-part1.fastq");
        let fastq2 = temp_path("matching-part2.fastq");
        let names = collect_bam_read_names(bam);
        let bam_records = count_bam_records(bam);
        let midpoint = names.len() / 2;

        write_fastq(&fastq1, &names[..midpoint]);
        write_fastq(&fastq2, &names[midpoint..]);

        let output = Command::new(env!("CARGO_BIN_EXE_bam-check-fastq"))
            .arg(bam)
            .arg(&fastq1)
            .arg(&fastq2)
            .output()
            .unwrap();

        let _ = remove_file(&fastq1);
        let _ = remove_file(&fastq2);

        assert!(output.status.success());
        let stdout = String::from_utf8_lossy(&output.stdout);

        assert!(stdout.contains("FASTQ File Summary"));
        assert!(stdout.contains(&format_count(bam_records)));
        assert!(stdout.contains(fastq1.to_string_lossy().as_ref()));
        assert!(stdout.contains(fastq2.to_string_lossy().as_ref()));
        assert!(stdout.contains("[1/2]"));
        assert!(stdout.contains("[2/2]"));
    }

    #[test]
    fn bam_check_fastq_reports_missing_names() {
        let bam = "tests/test_bam_2.bam";
        let fastq = temp_path("missing.fastq");
        let mut names = collect_bam_read_names(bam);
        let bam_records = count_bam_records(bam);
        let missing = names.remove(0);
        let missing_records = count_bam_records_with_name(bam, &missing);
        let matched_records = bam_records - missing_records;

        write_fastq(&fastq, &names);

        let output = Command::new(env!("CARGO_BIN_EXE_bam-check-fastq"))
            .arg(bam)
            .arg(&fastq)
            .arg("--report-limit")
            .arg("1")
            .output()
            .unwrap();

        let _ = remove_file(&fastq);

        assert!(!output.status.success());

        let stderr = String::from_utf8_lossy(&output.stderr);

        assert!(stderr.contains("Result"));
        assert!(stderr.contains("BAM->FASTQ coverage"));
        assert!(stderr.contains("INCOMPLETE"));
        assert!(stderr.contains("Missing BAM unique names"));
        assert!(stderr.contains("Missing BAM Read Names"));
        assert!(stderr.contains("Missing from FASTQ"));
        assert!(stderr.contains(&format_count(bam_records)));
        assert!(stderr.contains(&format!(
            "{} ({:.2}%)",
            format_count(matched_records),
            100.0 * matched_records as f64 / bam_records as f64
        )));
        assert!(stderr.contains(&format!(
            "{} ({:.2}%)",
            format_count(missing_records),
            100.0 * missing_records as f64 / bam_records as f64
        )));
        assert!(stderr.contains(&missing));
    }

    #[test]
    fn bam_check_fastq_can_keep_pair_suffixes() {
        let bam = "tests/test_bam_2.bam";
        let fastq = temp_path("strict.fastq");
        let names = collect_bam_read_names(bam);

        write_fastq(&fastq, &names);

        let output = Command::new(env!("CARGO_BIN_EXE_bam-check-fastq"))
            .arg(bam)
            .arg(&fastq)
            .arg("--keep-pair-suffixes")
            .output()
            .unwrap();

        let _ = remove_file(&fastq);

        assert!(!output.status.success());
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("keep trailing /1 and /2 suffixes")
        );
    }

    #[test]
    fn bam_check_fastq_reports_fastq_only_unique_names() {
        let bam = "tests/test_bam_2.bam";
        let fastq = temp_path("extra.fastq");
        let mut names = collect_bam_read_names(bam);
        names.push("read-only-in-fastq".to_string());

        write_fastq(&fastq, &names);

        let output = Command::new(env!("CARGO_BIN_EXE_bam-check-fastq"))
            .arg(bam)
            .arg(&fastq)
            .output()
            .unwrap();

        let _ = remove_file(&fastq);

        assert!(output.status.success());
        let stdout = String::from_utf8_lossy(&output.stdout);

        assert!(stdout.contains("FASTQ-only unique read names"));
        assert!(stdout.contains("1 ("));
        assert!(stdout.contains("Unique names also in BAM"));
    }

    #[test]
    fn bam_check_fastq_continues_after_truncated_fastq() {
        let bam = "tests/test_bam_2.bam";
        let fastq1 = temp_path("matching-part1.fastq");
        let fastq_truncated = temp_path("truncated.fastq");
        let fastq2 = temp_path("matching-part2.fastq");
        let names = collect_bam_read_names(bam);
        let midpoint = names.len() / 2;

        write_fastq(&fastq1, &names[..midpoint]);
        write_fastq(&fastq2, &names[midpoint..]);

        let file = File::create(&fastq_truncated).unwrap();
        let mut writer = BufWriter::new(file);
        writeln!(writer, "@truncated/1").unwrap();
        writeln!(writer, "ACGT").unwrap();
        writer.flush().unwrap();

        let output = Command::new(env!("CARGO_BIN_EXE_bam-check-fastq"))
            .arg(bam)
            .arg(&fastq1)
            .arg(&fastq_truncated)
            .arg(&fastq2)
            .output()
            .unwrap();

        let _ = remove_file(&fastq1);
        let _ = remove_file(&fastq_truncated);
        let _ = remove_file(&fastq2);

        assert!(output.status.success());

        let stdout = String::from_utf8_lossy(&output.stdout);

        assert!(stdout.contains("Result"));
        assert!(stdout.contains("BAM->FASTQ coverage"));
        assert!(stdout.contains("COMPLETE"));
        assert!(stdout.contains("FASTQ input status"));
        assert!(stdout.contains("INCOMPLETE (1 of 3 files truncated)"));
        assert!(stdout.contains("FASTQ Warnings"));
        assert!(stdout.contains("TRUNCATED"));
        assert!(stdout.contains(fastq_truncated.to_string_lossy().as_ref()));
        assert!(stdout.contains("unexpected end of FASTQ"));
        assert!(stdout.contains(fastq2.to_string_lossy().as_ref()));
    }

    #[test]
    fn bam_check_fastq_continues_after_malformed_fastq_tail() {
        let bam = "tests/test_bam_2.bam";
        let fastq_ok = temp_path("valid.fastq");
        let fastq_malformed_tail = temp_path("malformed-tail.fastq");
        let names = collect_bam_read_names(bam);

        write_fastq(&fastq_ok, &names[..1]);

        let file = File::create(&fastq_malformed_tail).unwrap();
        let mut writer = BufWriter::new(file);
        writeln!(writer, "@good/1").unwrap();
        writeln!(writer, "ACGT").unwrap();
        writeln!(writer, "+").unwrap();
        writeln!(writer, "IIII").unwrap();
        writeln!(writer, "@broken/1").unwrap();
        writeln!(writer, "ACGT").unwrap();
        writeln!(writer, "not-a-separator").unwrap();
        writeln!(writer, "IIII").unwrap();
        writer.flush().unwrap();

        let output = Command::new(env!("CARGO_BIN_EXE_bam-check-fastq"))
            .arg(bam)
            .arg(&fastq_ok)
            .arg(&fastq_malformed_tail)
            .output()
            .unwrap();

        let _ = remove_file(&fastq_ok);
        let _ = remove_file(&fastq_malformed_tail);

        assert!(!output.status.success());

        let stderr = String::from_utf8_lossy(&output.stderr);

        assert!(stderr.contains("Result"));
        assert!(stderr.contains("BAM->FASTQ coverage"));
        assert!(stderr.contains("INCOMPLETE"));
        assert!(stderr.contains("FASTQ input status"));
        assert!(stderr.contains("FASTQ Warnings"));
        assert!(stderr.contains("invalid FASTQ separator"));
        assert!(stderr.contains(fastq_malformed_tail.to_string_lossy().as_ref()));
    }

    #[test]
    fn bam_check_fastq_reports_the_broken_fastq_file() {
        let bam = "tests/test_bam_2.bam";
        let fastq_ok = temp_path("valid.fastq");
        let fastq_broken = temp_path("broken.fastq");
        let names = collect_bam_read_names(bam);

        write_fastq(&fastq_ok, &names[..1]);

        let file = File::create(&fastq_broken).unwrap();
        let mut writer = BufWriter::new(file);
        writeln!(writer, "@broken/1").unwrap();
        writeln!(writer, "ACGT").unwrap();
        writeln!(writer, "not-a-separator").unwrap();
        writeln!(writer, "IIII").unwrap();
        writer.flush().unwrap();

        let output = Command::new(env!("CARGO_BIN_EXE_bam-check-fastq"))
            .arg(bam)
            .arg(&fastq_ok)
            .arg(&fastq_broken)
            .output()
            .unwrap();

        let _ = remove_file(&fastq_ok);
        let _ = remove_file(&fastq_broken);

        assert!(!output.status.success());

        let stderr = String::from_utf8_lossy(&output.stderr);

        assert!(stderr.contains("failed to read FASTQ file"));
        assert!(stderr.contains(fastq_broken.to_string_lossy().as_ref()));
        assert!(stderr.contains("invalid FASTQ separator"));
    }
}
