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
use crate::kmer_catalogue::KmerCatalogue;
use crate::kmer_class::{KmerClass, KmerClassId};
use crate::kmer_equivalence_relation::KmerEquivalenceRelation;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct KmerGraph<A: ComplementableAlphabet + Clone> {
    pub catalogue: KmerCatalogue<A>,
    pub nodes: HashMap<KmerClassId, KmerGraphNode>,
}

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct KmerGraphNode {
    pub kmer: KmerClass,
    pub cluster_id: u64,
    pub rank: usize,
    pub infra: Vec<KmerClassId>,
    pub supra: Vec<KmerClassId>,
    pub intra: Vec<KmerClassId>,
}

/* -------------------------------------------------------------------------- */

impl<A: ComplementableAlphabet + Clone> KmerGraph<A> {
    pub fn new(kmers: Vec<KmerClass>, relation: KmerEquivalenceRelation<A>) -> Self {
        let mut catalogue = KmerCatalogue::from_relation(relation);
        for kmer in kmers {
            catalogue.add_kmer_class(kmer);
        }
        Self::from_catalogue(catalogue)
    }

    pub fn from_catalogue(catalogue: KmerCatalogue<A>) -> Self {
        let mut graph = Self {
            catalogue,
            nodes: HashMap::new(),
        };
        graph.construct_graph();
        graph
    }

    pub fn get_node(&self, kmer: &str) -> Option<&KmerGraphNode> {
        let class = self.catalogue.get_kmer_class_if_present(kmer)?;
        self.nodes.get(&class.id())
    }

    pub fn related_kmers(&self, kmer: &str) -> Vec<KmerClass> {
        let Some(node) = self.get_node(kmer) else {
            return Vec::new();
        };

        let mut result = Vec::new();
        for id in &node.infra {
            result.push(self.nodes[id].kmer.clone());
        }
        for id in &node.supra {
            result.push(self.nodes[id].kmer.clone());
        }
        result.sort();
        result
    }

    pub fn related_clusters(&self, cluster_id1: u64, cluster_id2: u64) -> bool {
        let mut r1 = 0u32;
        let mut r2 = 0u32;
        for k in 0..=self.catalogue.relation.m() {
            if cluster_id1 & (1 << k) != 0 {
                r1 = k as u32 + 1;
            }
            if cluster_id2 & (1 << k) != 0 {
                r2 = k as u32 + 1;
            }
        }
        if r2 > r1 {
            return false;
        }
        'outer: for offset in 0..=(r1 - r2) {
            for k in 0..r1 {
                let a = cluster_id1 & (1 << (k + offset)) != 0;
                let b = cluster_id2 & (1 << k) != 0;
                if !((a && b) || !b) {
                    continue 'outer;
                }
            }
            return true;
        }
        false
    }

    fn rank(&self, kmer: &KmerClass) -> usize {
        let mut count = 0usize;
        for base in kmer.elements[0].bytes() {
            if !self
                .catalogue
                .relation
                .alphabet()
                .is_wildcard(base)
                .unwrap()
            {
                count += 1;
            }
        }
        count.saturating_sub(1)
    }

    fn cluster_id(&self, kmer: &KmerClass) -> u64 {
        let mut id = 0u64;
        for (i, base) in kmer.elements[0].bytes().enumerate() {
            if !self
                .catalogue
                .relation
                .alphabet()
                .is_wildcard(base)
                .unwrap()
            {
                id += 1 << i;
            }
        }
        id
    }

    fn construct_graph(&mut self) {
        let n = self.catalogue.relation.n();
        let m = self.catalogue.relation.m();
        let mut layers: Vec<HashMap<u64, Vec<KmerClassId>>> = vec![HashMap::new(); m + 1];

        for k in n..=m {
            for (id, elements) in self.catalogue.elements[k - n].clone() {
                let kmer = KmerClass::new(k, id, elements);
                let node_id = self.new_node(kmer);
                let node = self.nodes.get(&node_id).unwrap();
                layers[node.rank]
                    .entry(node.cluster_id)
                    .or_default()
                    .push(node_id);
            }
        }

        let layers: Vec<_> = layers
            .into_iter()
            .filter(|layer| !layer.is_empty())
            .collect();
        for j in 1..layers.len() {
            for (cluster1_id, cluster1) in &layers[j] {
                for (cluster2_id, cluster2) in &layers[j - 1] {
                    if self.related_clusters(*cluster1_id, *cluster2_id) {
                        self.connect_clusters(cluster1, cluster2);
                    }
                }
            }
        }
    }

    fn new_node(&mut self, kmer: KmerClass) -> KmerClassId {
        let id = kmer.id();
        let node = KmerGraphNode {
            cluster_id: self.cluster_id(&kmer),
            rank: self.rank(&kmer),
            kmer,
            infra: Vec::new(),
            supra: Vec::new(),
            intra: Vec::new(),
        };
        self.nodes.insert(id, node);
        id
    }

    fn connect_clusters(&mut self, cluster1: &[KmerClassId], cluster2: &[KmerClassId]) {
        for left in cluster1 {
            for right in cluster2 {
                self.connect_nodes(*left, *right);
            }
        }
    }

    fn connect_nodes(&mut self, id1: KmerClassId, id2: KmerClassId) {
        let should_connect = {
            let node1 = self.nodes.get(&id1).unwrap();
            let node2 = self.nodes.get(&id2).unwrap();
            node2.kmer.elements[0].len() <= node1.kmer.elements[0].len()
                && node2
                    .kmer
                    .matches(&node1.kmer, self.catalogue.relation.alphabet())
        };
        if !should_connect {
            return;
        }
        self.nodes.get_mut(&id1).unwrap().supra.push(id2);
        self.nodes.get_mut(&id2).unwrap().infra.push(id1);
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::collections::HashMap;

    use crate::alphabet::GappedNucleotideAlphabet;
    use crate::kmer_catalogue::KmerCatalogue;
    use crate::kmer_graph::KmerGraph;

    #[test]
    fn graph_matches_go_reference_relationships() {
        let mut kmers =
            KmerCatalogue::new(2, 6, false, false, true, None, GappedNucleotideAlphabet).unwrap();
        for query in [
            "at", "tc", "gctc", "gcta", "annnc", "atnnc", "anntc", "anctc", "angtc", "agctc",
            "aagntc", "aagctc", "agctca",
        ] {
            kmers.get_kmer_class(query);
        }

        let graph = KmerGraph::from_catalogue(kmers);

        let expected: HashMap<&str, Vec<&str>> = HashMap::from([
            ("tc", vec!["anntc|gannt"]),
            ("at", vec!["atnnc|gnnat"]),
            ("annnc", vec!["atnnc|gnnat", "anntc|gannt"]),
            ("atnnc", vec!["at|at", "annnc|gnnnt"]),
            (
                "anntc",
                vec!["ga|tc", "anctc|gagnt", "angtc|gacnt", "annnc|gnnnt"],
            ),
            ("gctc", vec!["agctc|gagct"]),
            ("angtc", vec!["anntc|gannt"]),
            ("anctc", vec!["agctc|gagct", "anntc|gannt"]),
            (
                "agctc",
                vec!["gagc|gctc", "anctc|gagnt", "aagctc|gagctt", "agctca|tgagct"],
            ),
            ("aagntc", vec!["aagctc|gagctt"]),
            ("aagctc", vec!["agctc|gagct", "aagntc|ganctt"]),
            ("agctca", vec!["agctc|gagct"]),
        ]);

        for (query, expected_kmers) in expected {
            let related = graph.related_kmers(query);
            assert_eq!(related.len(), expected_kmers.len(), "query={query}");
            for (index, kmer) in related.iter().enumerate() {
                assert_eq!(kmer.to_string(), expected_kmers[index], "query={query}");
            }
        }

        let expected_cluster_ids = HashMap::from([
            ("at", 3u64),
            ("annnc", 17u64),
            ("atnnc", 19u64),
            ("anntc", 25u64),
            ("gctc", 15u64),
            ("gcta", 15u64),
            ("anctc", 29u64),
            ("aagntc", 55u64),
        ]);
        for (query, cluster_id) in expected_cluster_ids {
            let node = graph.get_node(query).unwrap();
            assert_eq!(node.cluster_id, cluster_id, "query={query}");
        }

        for (left, right, related) in [
            (19u64, 3u64, true),
            (19u64, 17u64, true),
            (25u64, 3u64, true),
            (25u64, 17u64, true),
            (15u64, 19u64, false),
            (15u64, 25u64, false),
            (31u64, 15u64, true),
            (31u64, 29u64, true),
            (55u64, 15u64, false),
            (55u64, 29u64, false),
            (63u64, 31u64, true),
            (63u64, 55u64, true),
        ] {
            assert_eq!(graph.related_clusters(left, right), related);
        }
    }
}
