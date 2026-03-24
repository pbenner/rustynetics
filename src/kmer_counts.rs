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

use std::collections::HashMap;

use crate::kmer_class::{union_kmer_classes, KmerClass, KmerClassId};

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct KmerCounts {
    pub kmers: Vec<KmerClass>,
    pub counts: HashMap<KmerClassId, usize>,
}

/* -------------------------------------------------------------------------- */

impl KmerCounts {
    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    pub fn is_empty(&self) -> bool {
        self.kmers.is_empty()
    }

    pub fn n(&self) -> usize {
        self.counts.len()
    }

    pub fn at(&self, i: usize) -> usize {
        self.counts.get(&self.kmers[i].id()).copied().unwrap_or(0)
    }

    pub fn get_count(&self, kmer: &KmerClass) -> usize {
        self.counts.get(&kmer.id()).copied().unwrap_or(0)
    }

    pub fn get_kmer(&self, i: usize) -> &KmerClass {
        &self.kmers[i]
    }

    pub fn iter(&self) -> KmerCountsIterator<'_> {
        KmerCountsIterator { counts: self, i: 0 }
    }

    pub fn set_kmers(&mut self, kmers: Vec<KmerClass>) {
        let mut counts = HashMap::new();
        for kmer in &kmers {
            if let Some(count) = self.counts.get(&kmer.id()) {
                counts.insert(kmer.id(), *count);
            }
        }
        self.kmers = kmers;
        self.counts = counts;
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct KmerCountsList {
    pub kmers: Vec<KmerClass>,
    pub counts: Vec<HashMap<KmerClassId, usize>>,
}

/* -------------------------------------------------------------------------- */

impl KmerCountsList {
    pub fn new(counts: Vec<KmerCounts>) -> Self {
        let result = Self::default();
        result.append(counts)
    }

    pub fn append(mut self, args: Vec<KmerCounts>) -> Self {
        if args.is_empty() {
            return self;
        }

        let mut slices: Vec<&[KmerClass]> = Vec::with_capacity(args.len() + 1);
        slices.push(&self.kmers);
        for counts in &args {
            slices.push(&counts.kmers);
        }
        self.kmers = union_kmer_classes(&slices);
        for counts in args {
            self.counts.push(counts.counts);
        }
        self
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    pub fn at(&self, i: usize) -> KmerCounts {
        KmerCounts {
            kmers: self.kmers.clone(),
            counts: self.counts[i].clone(),
        }
    }

    pub fn slice(&self, i: usize, j: usize) -> Self {
        Self {
            kmers: self.kmers.clone(),
            counts: self.counts[i..j].to_vec(),
        }
    }

    pub fn set_kmers(&mut self, kmers: Vec<KmerClass>) {
        for counts in &mut self.counts {
            let mut filtered = HashMap::new();
            for kmer in &kmers {
                if let Some(value) = counts.get(&kmer.id()) {
                    filtered.insert(kmer.id(), *value);
                }
            }
            *counts = filtered;
        }
        self.kmers = kmers;
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct KmerCountsIterator<'a> {
    counts: &'a KmerCounts,
    i: usize,
}

/* -------------------------------------------------------------------------- */

impl<'a> KmerCountsIterator<'a> {
    pub fn ok(&self) -> bool {
        self.i < self.counts.len()
    }

    pub fn get_kmer(&self) -> &KmerClass {
        &self.counts.kmers[self.i]
    }

    pub fn get_count(&self) -> usize {
        self.counts.at(self.i)
    }

    pub fn get_index(&self) -> usize {
        self.i
    }

    pub fn next(&mut self) {
        self.i += 1;
    }
}
