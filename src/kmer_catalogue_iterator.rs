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
use crate::kmer_catalogue::KmerCatalogue;
use crate::kmer_class::KmerClass;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct KmerCatalogueIterator {
    kmers: Vec<KmerClass>,
    i: usize,
}

/* -------------------------------------------------------------------------- */

impl KmerCatalogueIterator {
    pub fn new<A: ComplementableAlphabet + Clone>(catalogue: &KmerCatalogue<A>) -> Self {
        let mut kmers = Vec::new();
        for (offset, table) in catalogue.elements.iter().enumerate() {
            let k = catalogue.relation.n() + offset;
            for (id, elements) in table {
                kmers.push(KmerClass::new(k, *id, elements.clone()));
            }
        }
        kmers.sort();
        Self { kmers, i: 0 }
    }

    pub fn ok(&self) -> bool {
        self.i < self.kmers.len()
    }

    pub fn get(&self) -> &KmerClass {
        &self.kmers[self.i]
    }

    pub fn next(&mut self) {
        self.i += 1;
    }
}
