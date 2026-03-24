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

use std::collections::{HashMap, HashSet};

use crate::alphabet::ComplementableAlphabet;
use crate::kmer_catalogue::KmerCatalogue;
use crate::kmer_class::{KmerClass, KmerClassId};
use crate::kmer_counts::KmerCounts;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct KmerCounter<A: ComplementableAlphabet + Clone> {
    pub catalogue: KmerCatalogue<A>,
    kmap: Vec<HashMap<usize, Vec<usize>>>,
    frozen: bool,
}

/* -------------------------------------------------------------------------- */

impl<A: ComplementableAlphabet + Clone> KmerCounter<A> {
    pub fn new(
        n: usize,
        m: usize,
        complement: bool,
        reverse: bool,
        revcomp: bool,
        max_ambiguous: Option<Vec<usize>>,
        alphabet: A,
    ) -> Result<Self, String> {
        let catalogue =
            KmerCatalogue::new(n, m, complement, reverse, revcomp, max_ambiguous, alphabet)?;
        let mut kmap = Vec::with_capacity(m - n + 1);
        for _ in n..=m {
            kmap.push(HashMap::new());
        }
        Ok(Self {
            catalogue,
            kmap,
            frozen: false,
        })
    }

    pub fn new_with_kmers(
        n: usize,
        m: usize,
        complement: bool,
        reverse: bool,
        revcomp: bool,
        max_ambiguous: Option<Vec<usize>>,
        alphabet: A,
        kmers: Vec<KmerClass>,
    ) -> Result<Self, String> {
        let mut counter = Self::new(n, m, complement, reverse, revcomp, max_ambiguous, alphabet)?;
        for kmer in kmers {
            counter.add_kmer(kmer);
        }
        counter.freeze();
        Ok(counter)
    }

    pub fn freeze(&mut self) {
        self.frozen = true;
    }

    pub fn count_kmers(&mut self, sequence: &[u8]) -> KmerCounts {
        let mut kmers = Vec::new();
        let mut counts = HashMap::new();
        for k in self.catalogue.relation.n()..=self.catalogue.relation.m() {
            let result = self.count_kmers_fixed(sequence, k);
            kmers.extend(result.kmers);
            for (kmer_id, count) in result.counts {
                counts.insert(kmer_id, count);
            }
        }
        KmerCounts { kmers, counts }
    }

    pub fn identify_kmers(&mut self, sequence: &[u8]) -> KmerCounts {
        let mut kmers = Vec::new();
        let mut counts = HashMap::new();
        for k in self.catalogue.relation.n()..=self.catalogue.relation.m() {
            let result = self.identify_kmers_fixed(sequence, k);
            kmers.extend(result.kmers);
            for (kmer_id, count) in result.counts {
                counts.insert(kmer_id, count);
            }
        }
        KmerCounts { kmers, counts }
    }

    fn generalize_kmer_rec(
        &mut self,
        dest: &mut [u8],
        src: &[u8],
        result: &mut HashSet<usize>,
        i: usize,
    ) {
        if i == src.len() {
            let class = self
                .catalogue
                .get_kmer_class(std::str::from_utf8(dest).unwrap());
            result.insert(class.i);
            return;
        }
        let matching = self.catalogue.relation.alphabet().matching(src[i]).unwrap();
        for base in matching {
            if self
                .catalogue
                .relation
                .alphabet()
                .is_wildcard(base)
                .unwrap()
                && (i == 0 || i + 1 == src.len())
            {
                continue;
            }
            dest[i] = base;
            self.generalize_kmer_rec(dest, src, result, i + 1);
        }
    }

    fn generalize_kmer(&mut self, dest: &mut [u8], src: &[u8], result: &mut HashSet<usize>) {
        self.generalize_kmer_rec(dest, src, result, 0);
    }

    fn instantiate_kmer_rec(
        &self,
        dest: &mut [u8],
        src: &[u8],
        result: &mut HashSet<usize>,
        i: usize,
    ) {
        if i == src.len() {
            let class = self
                .catalogue
                .relation
                .equivalence_class(std::str::from_utf8(dest).unwrap());
            result.insert(class.i);
            return;
        }
        let bases = self.catalogue.relation.alphabet().bases(src[i]).unwrap();
        for base in bases {
            dest[i] = base;
            self.instantiate_kmer_rec(dest, src, result, i + 1);
        }
    }

    fn instantiate_kmer(&self, dest: &mut [u8], src: &[u8], result: &mut HashSet<usize>) {
        self.instantiate_kmer_rec(dest, src, result, 0);
    }

    fn add_observed_kmer(&mut self, kmer: KmerClass) -> Vec<usize> {
        let mut dest = vec![0; kmer.k];
        let mut result = HashSet::new();
        self.generalize_kmer(&mut dest, kmer.elements[0].as_bytes(), &mut result);
        let mut ids: Vec<_> = result.into_iter().collect();
        ids.sort_unstable();
        self.kmap[kmer.k - self.catalogue.relation.n()].insert(kmer.i, ids.clone());
        ids
    }

    fn add_kmer(&mut self, kmer: KmerClass) {
        let mut dest = vec![0; kmer.k];
        let mut result = HashSet::new();
        self.catalogue.add_kmer_class(kmer.clone());
        self.instantiate_kmer(&mut dest, kmer.elements[0].as_bytes(), &mut result);
        for id in result {
            self.kmap[kmer.k - self.catalogue.relation.n()]
                .entry(id)
                .or_default()
                .push(kmer.i);
        }
    }

    fn matching_kmers(&mut self, c: &[u8]) -> Vec<usize> {
        let class = self
            .catalogue
            .get_kmer_class(std::str::from_utf8(c).unwrap());
        if let Some(ids) = self.kmap[class.k - self.catalogue.relation.n()].get(&class.i) {
            ids.clone()
        } else if self.frozen {
            Vec::new()
        } else {
            self.add_observed_kmer(class)
        }
    }

    fn count_kmers_fixed(&mut self, sequence: &[u8], k: usize) -> KmerCounts {
        let sequence: Vec<u8> = sequence
            .iter()
            .map(|base| base.to_ascii_lowercase())
            .collect();
        let mut raw_counts: HashMap<usize, usize> = HashMap::new();
        for i in 0..sequence.len() {
            if i + k > sequence.len() {
                continue;
            }
            for id in self.matching_kmers(&sequence[i..i + k]) {
                *raw_counts.entry(id).or_insert(0) += 1;
            }
        }

        let mut kmers = Vec::with_capacity(raw_counts.len());
        let mut counts = HashMap::new();
        for (id, count) in raw_counts {
            let kmer = self.catalogue.get_kmer_class_from_id(k, id);
            counts.insert(kmer.id(), count);
            kmers.push(kmer);
        }
        kmers.sort();
        KmerCounts { kmers, counts }
    }

    fn identify_kmers_fixed(&mut self, sequence: &[u8], k: usize) -> KmerCounts {
        let sequence: Vec<u8> = sequence
            .iter()
            .map(|base| base.to_ascii_lowercase())
            .collect();
        let mut present = HashSet::new();
        for i in 0..sequence.len() {
            if i + k > sequence.len() {
                continue;
            }
            for id in self.matching_kmers(&sequence[i..i + k]) {
                present.insert(id);
            }
        }

        let mut kmers = Vec::with_capacity(present.len());
        let mut counts: HashMap<KmerClassId, usize> = HashMap::new();
        for id in present {
            let kmer = self.catalogue.get_kmer_class_from_id(k, id);
            counts.insert(kmer.id(), 1);
            kmers.push(kmer);
        }
        kmers.sort();
        KmerCounts { kmers, counts }
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::alphabet::GappedNucleotideAlphabet;
    use crate::kmer_counter::KmerCounter;

    #[test]
    fn counter_matches_go_reference_output() {
        let expected_kmers = vec![
            "acgt|tgca|tgca|acgt",
            "acnt|tgna|tnca|angt",
            "agcg|tcgc|gcga|cgct",
            "agng|tcnc|gnga|cnct",
            "ancg|tngc|gcna|cgnt",
            "anng|tnnc|gnna|cnnt",
            "annt|tnna|tnna|annt",
            "cagc|gtcg|cgac|gctg",
            "canc|gtng|cnac|gntg",
            "cgcg|gcgc|gcgc|cgcg",
            "cgtc|gcag|ctgc|gacg",
            "cgnc|gcng|cngc|gncg",
            "cgng|gcnc|gngc|cncg",
            "ctnc|gang|cntc|gnag",
            "cnnc|gnng|cnnc|gnng",
            "cnng|gnnc|gnnc|cnng",
            "acgtc|tgcag|ctgca|gacgt",
            "acgnc|tgcng|cngca|gncgt",
            "acntc|tgnag|ctnca|gangt",
            "acnnc|tgnng|cnnca|gnngt",
            "agcgc|tcgcg|cgcga|gcgct",
            "agcnc|tcgng|cncga|gngct",
            "agngc|tcncg|cgnga|gcnct",
            "agnnc|tcnng|cnnga|gnnct",
            "ancgc|tngcg|cgcna|gcgnt",
            "ancnc|tngng|cncna|gngnt",
            "angtc|tncag|ctgna|gacnt",
            "angnc|tncng|cngna|gncnt",
            "anngc|tnncg|cgnna|gcnnt",
            "anntc|tnnag|ctnna|gannt",
            "annnc|tnnng|cnnna|gnnnt",
            "cagcg|gtcgc|gcgac|cgctg",
            "cagng|gtcnc|gngac|cnctg",
            "cancg|gtngc|gcnac|cgntg",
            "canng|gtnnc|gnnac|cnntg",
            "cgacg|gctgc|gcagc|cgtcg",
            "cgang|gctnc|gnagc|cntcg",
            "cgcng|gcgnc|gncgc|cngcg",
            "cgtng|gcanc|gntgc|cnacg",
            "cgncg|gcngc|gcngc|cgncg",
            "cgnng|gcnnc|gnngc|cnncg",
            "cnang|gntnc|gnanc|cntng",
            "cncng|gngnc|gncnc|cngng",
            "cnnng|gnnnc|gnnnc|cnnng",
        ];
        let expected_counts = vec![
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2,
        ];

        let mut counter =
            KmerCounter::new(4, 5, true, true, true, None, GappedNucleotideAlphabet).unwrap();
        let counts = counter.count_kmers(b"acgtcgcg");

        let mut i = 0;
        let mut it = counts.iter();
        while it.ok() {
            assert_eq!(it.get_kmer().to_string(), expected_kmers[i]);
            assert_eq!(it.get_count(), expected_counts[i]);
            i += 1;
            it.next();
        }
    }
}
