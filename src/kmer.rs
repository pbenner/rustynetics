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

use std::fmt;

use crate::alphabet::ComplementableAlphabet;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Kmer(pub String);

/* -------------------------------------------------------------------------- */

impl Kmer {
    pub fn new<S: Into<String>>(value: S) -> Self {
        Self(value.into())
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub fn matches<A: ComplementableAlphabet>(&self, other: &Kmer, alphabet: &A) -> bool {
        let kmer1 = self.0.as_bytes();
        let kmer2 = other.0.as_bytes();
        assert!(
            kmer1.len() <= kmer2.len(),
            "Kmer::matches(): kmer1 must be smaller or equal in length than kmer2"
        );

        for i in 0..=kmer2.len() - kmer1.len() {
            if matches_window(kmer1, &kmer2[i..i + kmer1.len()], alphabet) {
                return true;
            }
        }
        false
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<&str> for Kmer {
    fn from(value: &str) -> Self {
        Self::new(value)
    }
}

/* -------------------------------------------------------------------------- */

fn matches_window<A: ComplementableAlphabet>(kmer1: &[u8], kmer2: &[u8], alphabet: &A) -> bool {
    if kmer1.len() != kmer2.len() {
        return false;
    }
    for i in 0..kmer1.len() {
        if kmer1[i].eq_ignore_ascii_case(&kmer2[i]) {
            continue;
        }
        let b1 = alphabet.bases(kmer1[i]).unwrap();
        let b2 = alphabet.bases(kmer2[i]).unwrap();
        if !byte_superset(&b1, &b2) {
            return false;
        }
    }
    true
}

fn byte_superset(a: &[u8], b: &[u8]) -> bool {
    b.iter().all(|item| a.contains(item))
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::alphabet::GappedNucleotideAlphabet;
    use crate::kmer::Kmer;

    #[test]
    fn kmer_matches_longer_kmer_with_ambiguity() {
        let a = Kmer::from("anctc");
        let b = Kmer::from("aagctc");
        assert!(a.matches(&b, &GappedNucleotideAlphabet));
    }

    #[test]
    fn kmer_does_not_match_when_ambiguity_is_insufficient() {
        let a = Kmer::from("anctc");
        let b = Kmer::from("aagntc");
        assert!(!a.matches(&b, &GappedNucleotideAlphabet));
    }
}
