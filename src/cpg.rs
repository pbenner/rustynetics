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

/* -------------------------------------------------------------------------- */

// Gardiner-Garden and Frommer, J. Mol. Biol. (1987) 196 (2), 261-282:
//   observed * length / (number of C * number of G)
pub fn observed_over_expected_cpg(sequence: &[u8]) -> f64 {
    let mut n_c   = 0;
    let mut n_g   = 0;
    let mut n_cpg = 0;

    for &base in sequence.iter() {
        match base {
            b'c' | b'C' => n_c += 1,
            b'g' | b'G' => n_g += 1,
            _ => {},
        }
    }

    for j in 0..sequence.len() - 1 {
        if (sequence[j] == b'c' || sequence[j] == b'C') && (sequence[j + 1] == b'g' || sequence[j + 1] == b'G') {
            n_cpg += 1;
        }
    }

    if n_cpg != 0 {
        (n_cpg as f64 * sequence.len() as f64) / (n_c as f64 * n_g as f64)
    } else {
        0.0
    }
}