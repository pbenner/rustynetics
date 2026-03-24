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

        write_fastq(&fastq, &names);

        let output = Command::new(env!("CARGO_BIN_EXE_bam-check-fastq"))
            .arg(bam)
            .arg(&fastq)
            .output()
            .unwrap();

        let _ = remove_file(&fastq);

        assert!(output.status.success());
        assert!(String::from_utf8_lossy(&output.stdout).contains("All"));
    }

    #[test]
    fn bam_check_fastq_reports_missing_names() {
        let bam = "tests/test_bam_2.bam";
        let fastq = temp_path("missing.fastq");
        let mut names = collect_bam_read_names(bam);
        let missing = names.remove(0);

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

        assert!(stderr.contains("missing BAM records"));
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
        assert!(String::from_utf8_lossy(&output.stderr).contains("missing BAM records"));
    }
}
