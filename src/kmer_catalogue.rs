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

use crate::alphabet::ComplementableAlphabet;
use crate::kmer_class::KmerClass;
use crate::kmer_equivalence_relation::KmerEquivalenceRelation;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct KmerCatalogue<A: ComplementableAlphabet + Clone> {
    pub relation: KmerEquivalenceRelation<A>,
    pub elements: Vec<HashMap<usize, Vec<String>>>,
    pub idmap: Vec<HashMap<String, usize>>,
}

/* -------------------------------------------------------------------------- */

impl<A: ComplementableAlphabet + Clone> KmerCatalogue<A> {
    pub fn new(
        n: usize,
        m: usize,
        complement: bool,
        reverse: bool,
        revcomp: bool,
        max_ambiguous: Option<Vec<usize>>,
        alphabet: A,
    ) -> Result<Self, String> {
        let relation = KmerEquivalenceRelation::new(
            n,
            m,
            complement,
            reverse,
            revcomp,
            max_ambiguous,
            alphabet,
        )?;
        Ok(Self::from_relation(relation))
    }

    pub fn from_relation(relation: KmerEquivalenceRelation<A>) -> Self {
        let width = relation.m() - relation.n() + 1;
        let mut elements = Vec::with_capacity(width);
        let mut idmap = Vec::with_capacity(width);
        for _ in 0..width {
            elements.push(HashMap::new());
            idmap.push(HashMap::new());
        }
        Self {
            relation,
            elements,
            idmap,
        }
    }

    pub fn add_kmer_class(&mut self, kmer: KmerClass) {
        let offset = kmer.k - self.relation.n();
        for element in &kmer.elements {
            self.idmap[offset].insert(element.clone(), kmer.i);
        }
        self.elements[offset].insert(kmer.i, kmer.elements);
    }

    pub fn get_kmer_class_if_present(&self, kmer: &str) -> Option<KmerClass> {
        let kmer = normalize_kmer(kmer);
        let k = kmer.len();
        if k < self.relation.n() || k > self.relation.m() {
            panic!("KmerCatalogue::get_kmer_class_if_present(): k-mer has invalid length");
        }
        self.idmap[k - self.relation.n()]
            .get(&kmer)
            .map(|id| KmerClass::new(k, *id, self.elements[k - self.relation.n()][id].clone()))
    }

    pub fn get_kmer_class(&mut self, kmer: &str) -> KmerClass {
        let kmer = normalize_kmer(kmer);
        let k = kmer.len();
        if k < self.relation.n() || k > self.relation.m() {
            panic!("KmerCatalogue::get_kmer_class(): k-mer has invalid length");
        }
        if let Some(id) = self.idmap[k - self.relation.n()].get(&kmer).copied() {
            KmerClass::new(k, id, self.elements[k - self.relation.n()][&id].clone())
        } else {
            let class = self.relation.equivalence_class(&kmer);
            self.add_kmer_class(class.clone());
            class
        }
    }

    pub fn get_kmer_class_from_id(&self, k: usize, id: usize) -> KmerClass {
        if let Some(elements) = self.elements[k - self.relation.n()].get(&id) {
            KmerClass::new(k, id, elements.clone())
        } else {
            panic!("KmerCatalogue::get_kmer_class_from_id(): k-mer name not found");
        }
    }

    pub fn catalogue_size(&self) -> usize {
        self.elements.iter().map(HashMap::len).sum()
    }
}

/* -------------------------------------------------------------------------- */

fn normalize_kmer(kmer: &str) -> String {
    kmer.bytes()
        .map(|base| base.to_ascii_lowercase() as char)
        .collect()
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::alphabet::GappedNucleotideAlphabet;
    use crate::kmer_catalogue::KmerCatalogue;

    #[test]
    fn catalogue_matches_revcomp_equivalence_class() {
        let mut kmers =
            KmerCatalogue::new(4, 5, false, false, true, None, GappedNucleotideAlphabet).unwrap();
        let j = kmers.get_kmer_class("anntc").i;
        let i = kmers.get_kmer_class("anntc").i;
        let s = kmers.get_kmer_class("anntc").to_string();
        let t = kmers.get_kmer_class("gannt").to_string();
        assert_eq!(i, j);
        assert_eq!(s, "anntc|gannt");
        assert_eq!(s, t);
    }

    #[test]
    fn catalogue_reverse_lookup_matches_forward_lookup() {
        let mut kmers =
            KmerCatalogue::new(4, 5, false, false, true, None, GappedNucleotideAlphabet).unwrap();
        let j = kmers.get_kmer_class("gannt").i;
        let i = kmers.get_kmer_class("anntc").i;
        let s = kmers.get_kmer_class("anntc").to_string();
        let t = kmers.get_kmer_class("gannt").to_string();
        assert_eq!(i, j);
        assert_eq!(s, "anntc|gannt");
        assert_eq!(s, t);
    }

    #[test]
    fn catalogue_combines_complement_reverse_and_revcomp() {
        let mut kmers =
            KmerCatalogue::new(4, 5, true, true, true, None, GappedNucleotideAlphabet).unwrap();
        let j = kmers.get_kmer_class("tnnag").i;
        let i = kmers.get_kmer_class("anntc").i;
        let s = kmers.get_kmer_class("anntc").to_string();
        let t = kmers.get_kmer_class("gannt").to_string();
        assert_eq!(i, j);
        assert_eq!(s, "anntc|tnnag|ctnna|gannt");
        assert_eq!(s, t);
    }

    #[test]
    fn catalogue_lookup_by_reverse_sequence_matches_class() {
        let mut kmers =
            KmerCatalogue::new(4, 5, true, true, true, None, GappedNucleotideAlphabet).unwrap();
        let j = kmers.get_kmer_class("ctnna").i;
        let i = kmers.get_kmer_class("anntc").i;
        let s = kmers.get_kmer_class("anntc").to_string();
        let t = kmers.get_kmer_class("gannt").to_string();
        assert_eq!(i, j);
        assert_eq!(s, "anntc|tnnag|ctnna|gannt");
        assert_eq!(s, t);
    }
}
