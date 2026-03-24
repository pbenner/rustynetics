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

use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt;

use crate::alphabet::ComplementableAlphabet;
use crate::kmer::Kmer;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct KmerClassId {
    pub k: usize,
    pub i: usize,
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct KmerClass {
    pub k: usize,
    pub i: usize,
    pub elements: Vec<String>,
}

/* -------------------------------------------------------------------------- */

impl KmerClass {
    pub fn new(k: usize, i: usize, elements: Vec<String>) -> Self {
        Self { k, i, elements }
    }

    pub fn id(&self) -> KmerClassId {
        KmerClassId {
            k: self.k,
            i: self.i,
        }
    }

    pub fn matches<A: ComplementableAlphabet>(&self, other: &KmerClass, alphabet: &A) -> bool {
        for left in &self.elements {
            for right in &other.elements {
                if Kmer::from(left.as_str()).matches(&Kmer::from(right.as_str()), alphabet) {
                    return true;
                }
            }
        }
        false
    }

    pub fn count_ambiguous<A: ComplementableAlphabet>(&self, alphabet: &A) -> usize {
        self.elements
            .first()
            .map(|sequence| {
                sequence
                    .bytes()
                    .filter(|base| alphabet.is_ambiguous(*base).unwrap_or(false))
                    .count()
            })
            .unwrap_or(0)
    }

    pub fn count_wildcard<A: ComplementableAlphabet>(&self, alphabet: &A) -> usize {
        self.elements
            .first()
            .map(|sequence| {
                sequence
                    .bytes()
                    .filter(|base| alphabet.is_wildcard(*base).unwrap_or(false))
                    .count()
            })
            .unwrap_or(0)
    }

    pub fn string(&self) -> String {
        self.to_string()
    }
}

/* -------------------------------------------------------------------------- */

impl Ord for KmerClass {
    fn cmp(&self, other: &Self) -> Ordering {
        self.id().cmp(&other.id())
    }
}

impl PartialOrd for KmerClass {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for KmerClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.elements.join("|"))
    }
}

/* -------------------------------------------------------------------------- */

pub type KmerClassSet = HashMap<KmerClassId, Vec<String>>;

pub(crate) fn union_kmer_classes(slices: &[&[KmerClass]]) -> Vec<KmerClass> {
    let mut set: KmerClassSet = HashMap::new();
    for slice in slices {
        for class in *slice {
            set.insert(class.id(), class.elements.clone());
        }
    }
    let mut kmers: Vec<_> = set
        .into_iter()
        .map(|(id, elements)| KmerClass::new(id.k, id.i, elements))
        .collect();
    kmers.sort();
    kmers
}
