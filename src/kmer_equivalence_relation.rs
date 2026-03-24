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
use crate::kmer_class::KmerClass;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct KmerEquivalence<A: ComplementableAlphabet + Clone> {
    pub n: usize,
    pub m: usize,
    pub complement: bool,
    pub reverse: bool,
    pub revcomp: bool,
    pub max_ambiguous: Vec<Option<usize>>,
    pub alphabet: A,
}

/* -------------------------------------------------------------------------- */

impl<A: ComplementableAlphabet + Clone> KmerEquivalence<A> {
    pub fn new(
        n: usize,
        m: usize,
        complement: bool,
        reverse: bool,
        revcomp: bool,
        max_ambiguous: Option<Vec<usize>>,
        alphabet: A,
    ) -> Result<Self, String> {
        if n > m {
            return Err("KmerEquivalence::new(): n must be <= m".to_string());
        }

        let max_ambiguous = normalize_max_ambiguous(n, m, max_ambiguous)?;
        Ok(Self {
            n,
            m,
            complement,
            reverse,
            revcomp,
            max_ambiguous,
            alphabet,
        })
    }

    pub fn equals(&self, other: &Self) -> bool {
        self.n == other.n
            && self.m == other.m
            && self.complement == other.complement
            && self.reverse == other.reverse
            && self.revcomp == other.revcomp
            && self.alphabet.string() == other.alphabet.string()
            && self.max_ambiguous == other.max_ambiguous
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct KmerEquivalenceRelation<A: ComplementableAlphabet + Clone> {
    pub equivalence: KmerEquivalence<A>,
    powers: Vec<usize>,
}

/* -------------------------------------------------------------------------- */

impl<A: ComplementableAlphabet + Clone> KmerEquivalenceRelation<A> {
    pub fn new(
        n: usize,
        m: usize,
        complement: bool,
        reverse: bool,
        revcomp: bool,
        max_ambiguous: Option<Vec<usize>>,
        alphabet: A,
    ) -> Result<Self, String> {
        let equivalence =
            KmerEquivalence::new(n, m, complement, reverse, revcomp, max_ambiguous, alphabet)?;

        let mut powers = vec![1usize; m + 1];
        for i in 1..=m {
            powers[i] = powers[i - 1] * equivalence.alphabet.length();
        }

        Ok(Self {
            equivalence,
            powers,
        })
    }

    pub fn n(&self) -> usize {
        self.equivalence.n
    }

    pub fn m(&self) -> usize {
        self.equivalence.m
    }

    pub fn alphabet(&self) -> &A {
        &self.equivalence.alphabet
    }

    pub fn complement(&self) -> bool {
        self.equivalence.complement
    }

    pub fn reverse(&self) -> bool {
        self.equivalence.reverse
    }

    pub fn revcomp(&self) -> bool {
        self.equivalence.revcomp
    }

    pub fn max_ambiguous_for_k(&self, k: usize) -> Option<usize> {
        self.equivalence.max_ambiguous[k - self.n()]
    }

    pub fn equivalence_class(&self, kmer: &str) -> KmerClass {
        let mut c1: Vec<u8> = kmer.bytes().map(|base| base.to_ascii_lowercase()).collect();
        let k = c1.len();
        let mut c2 = vec![0; k];
        let mut c3 = vec![0; k];
        let mut c4 = vec![0; k];

        let mut i = 0usize;
        let mut i_c = 0usize;
        let mut i_r = 0usize;
        let mut i_rc = 0usize;
        for j in 0..k {
            let x1 = self.alphabet().code(c1[j]).unwrap() as usize;
            let x2 = self.alphabet().code(c1[k - j - 1]).unwrap() as usize;
            let y1 = self.alphabet().complement_coded(x1 as u8).unwrap() as usize;
            let y2 = self.alphabet().complement_coded(x2 as u8).unwrap() as usize;
            i += x2 * self.powers[j];
            i_c += y2 * self.powers[j];
            i_r += x1 * self.powers[j];
            i_rc += y1 * self.powers[j];
        }

        let mut selected_complement = Some(i_c);
        let mut selected_reverse = Some(i_r);
        let mut selected_revcomp = Some(i_rc);
        if self.complement() && i > i_c {
            i = i_c;
        } else {
            selected_complement = None;
        }
        if self.reverse() && i > i_r {
            i = i_r;
        } else {
            selected_reverse = None;
        }
        if self.revcomp() && i > i_rc {
            i = i_rc;
        } else {
            selected_revcomp = None;
        }

        let mut built_c2 = false;
        let mut built_c3 = false;
        let mut built_c4 = false;
        if Some(i) == selected_complement {
            self.comp(&mut c2, &c1);
            std::mem::swap(&mut c1, &mut c2);
            built_c2 = true;
        } else if Some(i) == selected_reverse {
            self.rev(&mut c3, &c1);
            std::mem::swap(&mut c1, &mut c3);
            built_c3 = true;
        } else if Some(i) == selected_revcomp {
            self.rev(&mut c2, &c1);
            self.comp(&mut c4, &c2);
            std::mem::swap(&mut c1, &mut c4);
            built_c2 = true;
            built_c4 = true;
        }

        if (!built_c2 && self.complement()) || self.revcomp() {
            self.comp(&mut c2, &c1);
        }
        if (!built_c3 && self.reverse()) || self.revcomp() {
            self.rev(&mut c3, &c1);
        }
        if !built_c4 && self.revcomp() {
            self.rev(&mut c4, &c2);
        }

        let mut elements = vec![String::from_utf8(c1).unwrap()];
        if self.complement() {
            elements.push(String::from_utf8(c2).unwrap());
        }
        if self.reverse() {
            elements.push(String::from_utf8(c3).unwrap());
        }
        if self.revcomp() {
            elements.push(String::from_utf8(c4).unwrap());
        }
        KmerClass::new(k, i, elements)
    }

    fn comp(&self, dest: &mut [u8], src: &[u8]) {
        for (d, s) in dest.iter_mut().zip(src.iter().copied()) {
            *d = self.alphabet().complement(s).unwrap();
        }
    }

    fn rev(&self, dest: &mut [u8], src: &[u8]) {
        for (target, source) in dest.iter_mut().zip(src.iter().rev().copied()) {
            *target = source;
        }
    }
}

/* -------------------------------------------------------------------------- */

fn normalize_max_ambiguous(
    n: usize,
    m: usize,
    max_ambiguous: Option<Vec<usize>>,
) -> Result<Vec<Option<usize>>, String> {
    let width = m - n + 1;
    match max_ambiguous {
        None => Ok(vec![None; width]),
        Some(values) if values.len() == 1 => Ok(vec![Some(values[0]); width]),
        Some(values) if values.len() == width => Ok(values.into_iter().map(Some).collect()),
        Some(_) => {
            Err("KmerEquivalence::new(): parameter `max_ambiguous` has invalid length".to_string())
        }
    }
}
