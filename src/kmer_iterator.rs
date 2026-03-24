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

use crate::alphabet::ComplementableAlphabet;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct KmerIterator<A: ComplementableAlphabet + Clone> {
    c: Vec<u8>,
    alphabet: A,
    ok: bool,
    max_ambiguous: Option<usize>,
    num_ambiguous: usize,
}

/* -------------------------------------------------------------------------- */

impl<A: ComplementableAlphabet + Clone> KmerIterator<A> {
    pub fn new(k: usize, max_ambiguous: Option<usize>, alphabet: A) -> Self {
        let mut c = vec![0; k];
        for item in &mut c {
            *item = alphabet.decode(0).unwrap();
        }
        let num_ambiguous = count_ambiguous(&alphabet, &c);
        Self {
            c,
            alphabet,
            ok: true,
            max_ambiguous,
            num_ambiguous,
        }
    }

    pub fn get(&self) -> String {
        String::from_utf8(self.c.clone()).unwrap()
    }

    pub fn ok(&self) -> bool {
        self.ok
    }

    pub fn next(&mut self) {
        let k = self.c.len();
        let mut step = 0;
        while step < k {
            let idx = k - step - 1;
            let was_ambiguous = self.alphabet.is_ambiguous(self.c[idx]).unwrap();
            let mut ret = self.increment_position(idx);
            let is_ambiguous = self.alphabet.is_ambiguous(self.c[idx]).unwrap();
            if is_ambiguous && !was_ambiguous {
                self.num_ambiguous += 1;
            } else if !is_ambiguous && was_ambiguous {
                self.num_ambiguous -= 1;
            }
            if self
                .max_ambiguous
                .is_some_and(|limit| self.num_ambiguous > limit)
            {
                ret = false;
            } else {
                step += 1;
            }
            if ret {
                return;
            }
        }
        self.ok = false;
    }

    fn increment_position(&mut self, i: usize) -> bool {
        let code = self.alphabet.code(self.c[i]).unwrap();
        if (code as usize) + 1 < self.alphabet.length() {
            self.c[i] = self.alphabet.decode(code + 1).unwrap();
            true
        } else {
            self.c[i] = self.alphabet.decode(0).unwrap();
            false
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct KmerCylinderIterator<A: ComplementableAlphabet + Clone> {
    c: Vec<u8>,
    alphabet: A,
    ok: bool,
    max_ambiguous: Option<usize>,
    num_ambiguous: usize,
    j: usize,
    m: usize,
}

/* -------------------------------------------------------------------------- */

impl<A: ComplementableAlphabet + Clone> KmerCylinderIterator<A> {
    pub fn new(k: usize, max_ambiguous: Option<usize>, alphabet: A, j: usize, fixed: &str) -> Self {
        let bytes = fixed.as_bytes();
        let m = j + bytes.len();
        if m > k {
            panic!("KmerCylinderIterator::new(): invalid parameters");
        }

        let mut c = vec![0; k];
        for item in &mut c {
            *item = alphabet.decode(0).unwrap();
        }
        c[j..m].copy_from_slice(bytes);
        let num_ambiguous = count_ambiguous(&alphabet, &c);

        Self {
            c,
            alphabet,
            ok: true,
            max_ambiguous,
            num_ambiguous,
            j,
            m,
        }
    }

    pub fn get(&self) -> String {
        String::from_utf8(self.c.clone()).unwrap()
    }

    pub fn ok(&self) -> bool {
        self.ok
    }

    pub fn next(&mut self) {
        let k = self.c.len() as isize;
        let mut i = k - 1;
        if i >= self.j as isize && i < self.m as isize {
            i = self.j as isize - 1;
        }
        while i >= 0 {
            let idx = i as usize;
            let was_ambiguous = self.alphabet.is_ambiguous(self.c[idx]).unwrap();
            let mut ret = self.increment_position(idx);
            let is_ambiguous = self.alphabet.is_ambiguous(self.c[idx]).unwrap();
            if is_ambiguous && !was_ambiguous {
                self.num_ambiguous += 1;
            } else if !is_ambiguous && was_ambiguous {
                self.num_ambiguous -= 1;
            }
            if self
                .max_ambiguous
                .is_some_and(|limit| self.num_ambiguous > limit)
            {
                ret = false;
            } else {
                i -= 1;
                if i >= self.j as isize && i < self.m as isize {
                    i = self.j as isize - 1;
                }
            }
            if ret {
                return;
            }
        }
        self.ok = false;
    }

    fn increment_position(&mut self, i: usize) -> bool {
        let code = self.alphabet.code(self.c[i]).unwrap();
        if (code as usize) + 1 < self.alphabet.length() {
            self.c[i] = self.alphabet.decode(code + 1).unwrap();
            true
        } else {
            self.c[i] = self.alphabet.decode(0).unwrap();
            false
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct KmerInstantiationIterator<A: ComplementableAlphabet + Clone> {
    original: Vec<u8>,
    current: Vec<u8>,
    bases: Vec<Vec<u8>>,
    indices: Vec<usize>,
    ok: bool,
    partial: bool,
    alphabet: A,
}

/* -------------------------------------------------------------------------- */

impl<A: ComplementableAlphabet + Clone> KmerInstantiationIterator<A> {
    pub fn new(alphabet: A, pattern: &str, partial: bool) -> Self {
        let original = pattern.as_bytes().to_vec();
        let mut current = original.clone();
        let mut bases = vec![Vec::new(); original.len()];
        let mut indices = vec![0; original.len()];
        let mut ok = false;

        for (i, base) in original.iter().copied().enumerate() {
            if !alphabet.is_ambiguous(base).unwrap() {
                continue;
            }
            let choices = alphabet.bases(base).unwrap();
            if partial {
                bases[i] = choices;
            } else {
                bases[i] = choices;
                current[i] = bases[i][0];
                indices[i] = 1;
            }
            ok = true;
        }

        Self {
            original,
            current,
            bases,
            indices,
            ok,
            partial,
            alphabet,
        }
    }

    pub fn get(&self) -> String {
        String::from_utf8(self.current.clone()).unwrap()
    }

    pub fn ok(&self) -> bool {
        self.ok
    }

    pub fn next(&mut self) {
        let _ = &self.alphabet;
        let mut i = self.bases.len() as isize - 1;
        while i >= 0
            && (self.bases[i as usize].is_empty()
                || self.bases[i as usize].len() == self.indices[i as usize])
        {
            i -= 1;
        }
        if i < 0 {
            self.ok = false;
            self.current.clear();
            return;
        }

        let i = i as usize;
        self.current[i] = self.bases[i][self.indices[i]];
        for reset in (i + 1)..self.bases.len() {
            if self.bases[reset].is_empty() {
                continue;
            }
            if self.partial {
                self.indices[reset] = 0;
                self.current[reset] = self.original[reset];
            } else {
                self.indices[reset] = 1;
                self.current[reset] = self.bases[reset][0];
            }
        }
        self.indices[i] += 1;
    }
}

/* -------------------------------------------------------------------------- */

fn count_ambiguous<A: ComplementableAlphabet>(alphabet: &A, bytes: &[u8]) -> usize {
    bytes
        .iter()
        .copied()
        .filter(|base| alphabet.is_ambiguous(*base).unwrap_or(false))
        .count()
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::alphabet::{GappedNucleotideAlphabet, NucleotideAlphabet};
    use crate::kmer_iterator::{KmerCylinderIterator, KmerInstantiationIterator, KmerIterator};

    #[test]
    fn cylinder_iterator_matches_go_order() {
        let expected = vec![
            "acgtaa", "acgtac", "acgtag", "acgtat", "ccgtaa", "ccgtac", "ccgtag", "ccgtat",
            "gcgtaa", "gcgtac", "gcgtag", "gcgtat", "tcgtaa", "tcgtac", "tcgtag", "tcgtat",
        ];

        let mut it = KmerCylinderIterator::new(6, Some(2), NucleotideAlphabet, 1, "cgta");
        for expected_kmer in expected {
            assert!(it.ok());
            assert_eq!(it.get(), expected_kmer);
            it.next();
        }
        assert!(!it.ok());
    }

    #[test]
    fn instantiation_iterator_without_partial_matches_go_order() {
        let expected = vec![
            "ctaaa", "ctaca", "ctaga", "ctata", "ctcaa", "ctcca", "ctcga", "ctcta", "ctgaa",
            "ctgca", "ctgga", "ctgta", "cttaa", "cttca", "cttga", "cttta",
        ];

        let mut it = KmerInstantiationIterator::new(GappedNucleotideAlphabet, "ctnna", false);
        for expected_kmer in expected {
            assert!(it.ok());
            assert_eq!(it.get(), expected_kmer);
            it.next();
        }
        assert!(!it.ok());
    }

    #[test]
    fn instantiation_iterator_with_partial_matches_go_order() {
        let expected = vec![
            "ctnna", "ctnaa", "ctnca", "ctnga", "ctnta", "ctana", "ctaaa", "ctaca", "ctaga",
            "ctata", "ctcna", "ctcaa", "ctcca", "ctcga", "ctcta", "ctgna", "ctgaa", "ctgca",
            "ctgga", "ctgta", "cttna", "cttaa", "cttca", "cttga", "cttta",
        ];

        let mut it = KmerInstantiationIterator::new(GappedNucleotideAlphabet, "ctnna", true);
        for expected_kmer in expected {
            assert!(it.ok());
            assert_eq!(it.get(), expected_kmer);
            it.next();
        }
        assert!(!it.ok());
    }

    #[test]
    fn base_iterator_respects_ambiguity_limit() {
        let mut it = KmerIterator::new(2, Some(0), GappedNucleotideAlphabet);
        let mut seen = Vec::new();
        while it.ok() {
            seen.push(it.get());
            it.next();
        }
        assert_eq!(seen.len(), 16);
        assert!(seen.iter().all(|kmer| !kmer.contains('n')));
    }
}
