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

use flate2::read::GzDecoder;

use crate::alphabet::{Alphabet, ComplementableAlphabet, NucleotideAlphabet};
use crate::utility::is_gzip;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TFMatrix {
    pub values: Vec<Vec<f64>>,
}

/* -------------------------------------------------------------------------- */

impl TFMatrix {
    pub fn new(values: Vec<Vec<f64>>) -> Self {
        Self { values }
    }

    pub fn empty() -> Self {
        Self::default()
    }

    pub fn length(&self) -> usize {
        self.values.first().map_or(0, Vec::len)
    }

    pub fn get(&self, c: u8, j: usize) -> f64 {
        let i = NucleotideAlphabet.code(c).unwrap() as usize;
        self.values[i][j]
    }

    pub fn get_row(&self, c: u8) -> &[f64] {
        let i = NucleotideAlphabet.code(c).unwrap() as usize;
        &self.values[i]
    }

    pub fn revcomp(&self) -> Self {
        let alphabet = NucleotideAlphabet;
        let mut values = vec![Vec::new(); alphabet.length()];
        for i in 0..alphabet.length() {
            let j = alphabet.complement_coded(i as u8).unwrap() as usize;
            values[j] = self.values[i].iter().rev().copied().collect();
        }
        Self { values }
    }

    pub fn score<F>(&self, sequence: &[u8], revcomp: bool, x0: f64, f: F) -> Result<f64, String>
    where
        F: Fn(f64, f64) -> f64,
    {
        if sequence.len() != self.length() {
            return Err("TFMatrix::score(): sequence has invalid length".to_string());
        }

        let alphabet = NucleotideAlphabet;
        let mut x = x0;
        if revcomp {
            for j in 0..self.length() {
                let base = sequence[self.length() - j - 1];
                if base != b'N' && base != b'n' {
                    let complement = alphabet.complement(base).unwrap();
                    x = f(x, self.get(complement, j));
                }
            }
        } else {
            for (j, base) in sequence.iter().copied().enumerate() {
                if base != b'N' && base != b'n' {
                    x = f(x, self.get(base, j));
                }
            }
        }
        Ok(x)
    }

    pub fn read_matrix<R: io::Read>(&mut self, reader: R) -> io::Result<()> {
        let reader = BufReader::new(reader);
        let mut ncols = None;
        self.values = vec![Vec::new(); NucleotideAlphabet.length()];

        for line in reader.lines() {
            let line = line?;
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.is_empty() {
                continue;
            }
            if fields.len() <= 1 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "TFMatrix::read_matrix(): invalid TF matrix",
                ));
            }
            if ncols.is_none() {
                ncols = Some(fields.len() - 1);
            }
            if fields.len() != ncols.unwrap() + 1 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "TFMatrix::read_matrix(): invalid TF matrix",
                ));
            }
            let mut data = Vec::with_capacity(fields.len() - 1);
            for field in &fields[1..] {
                data.push(field.parse::<f64>().map_err(|error| {
                    io::Error::new(io::ErrorKind::InvalidData, error.to_string())
                })?);
            }
            let row = NucleotideAlphabet
                .code(fields[0].as_bytes()[0])
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?
                as usize;
            self.values[row] = data;
        }

        Ok(())
    }

    pub fn import_matrix(&mut self, filename: &str) -> io::Result<()> {
        let file = File::open(filename)?;
        if is_gzip(filename) {
            self.read_matrix(BufReader::new(GzDecoder::new(file)))
        } else {
            self.read_matrix(BufReader::new(file))
        }
    }

    pub fn write_matrix<W: Write>(&self, mut writer: W) -> io::Result<()> {
        for i in 0..self.values.len() {
            let c = NucleotideAlphabet
                .decode(i as u8)
                .unwrap()
                .to_ascii_uppercase();
            write!(writer, "{} ", c as char)?;
            for value in &self.values[i] {
                write!(writer, "{value:.6} ")?;
            }
            writeln!(writer)?;
        }
        Ok(())
    }

    pub fn write_jaspar<W: Write>(&self, mut writer: W) -> io::Result<()> {
        for i in 0..self.values.len() {
            let c = NucleotideAlphabet
                .decode(i as u8)
                .unwrap()
                .to_ascii_uppercase();
            write!(writer, "{} [ ", c as char)?;
            for value in &self.values[i] {
                write!(writer, "{value:.6} ")?;
            }
            writeln!(writer, "]")?;
        }
        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, PartialEq)]
pub struct PWM {
    pub matrix: TFMatrix,
}

/* -------------------------------------------------------------------------- */

impl PWM {
    pub fn new(matrix: TFMatrix) -> Self {
        Self { matrix }
    }

    pub fn scores(&self, sequence: &[u8], revcomp: bool) -> Vec<f64> {
        let motif_len = self.matrix.length();
        if motif_len == 0 || sequence.len() < motif_len {
            return Vec::new();
        }

        let mut result = Vec::with_capacity(sequence.len() - motif_len + 1);
        for i in 0..=sequence.len() - motif_len {
            result.push(
                self.matrix
                    .score(&sequence[i..i + motif_len], revcomp, 0.0, |a, b| a + b)
                    .unwrap(),
            );
        }
        result
    }

    pub fn max_score(&self, sequence: &[u8], revcomp: bool) -> f64 {
        self.scores(sequence, revcomp)
            .into_iter()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    pub fn mean_score(&self, sequence: &[u8], revcomp: bool) -> f64 {
        let motif_len = self.matrix.length();
        if motif_len == 0 || sequence.len() < motif_len {
            return f64::NEG_INFINITY;
        }

        let mut result = 0.0;
        for i in 0..=sequence.len() - motif_len {
            let score = self
                .matrix
                .score(&sequence[i..i + motif_len], revcomp, 0.0, |a, b| a + b)
                .unwrap();
            result = log_add(result, score);
        }
        result - ((sequence.len() - motif_len + 1) as f64).ln()
    }
}

/* -------------------------------------------------------------------------- */

fn log_add(mut a: f64, mut b: f64) -> f64 {
    if a > b {
        std::mem::swap(&mut a, &mut b);
    }
    if a.is_infinite() && a.is_sign_negative() {
        return b;
    }
    b + (a - b).exp().ln_1p()
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::path::PathBuf;

    use approx::assert_abs_diff_eq;

    use crate::tf::{TFMatrix, PWM};

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("gonetics")
            .join(name)
    }

    #[test]
    fn import_matrix_matches_go_fixture() {
        let mut matrix = TFMatrix::empty();
        matrix
            .import_matrix(fixture("tf_test.table").to_str().unwrap())
            .unwrap();

        assert_abs_diff_eq!(matrix.get(b'a', 0), -0.308122295362332, epsilon = 1e-8);
    }

    #[test]
    fn revcomp_scores_are_symmetric_for_fixture_motif() {
        let mut matrix = TFMatrix::empty();
        matrix
            .import_matrix(fixture("tf_test.table").to_str().unwrap())
            .unwrap();

        let pwm = PWM::new(matrix);
        let seq = b"cacgtg";

        assert_eq!(pwm.scores(seq, false)[0], pwm.scores(seq, true)[0]);
    }

    #[test]
    fn max_score_matches_reverse_complement() {
        let mut matrix = TFMatrix::empty();
        matrix
            .import_matrix(fixture("tf_test.table").to_str().unwrap())
            .unwrap();

        let pwm = PWM::new(matrix);
        let seq = b"cacgtgaaaccctttgg";

        assert_eq!(pwm.max_score(seq, false), pwm.max_score(seq, true));
    }
}
