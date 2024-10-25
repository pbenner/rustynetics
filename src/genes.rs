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
use std::collections::HashMap;
use std::error::Error;

use crate::range::Range;
use crate::genome::Genome;
use crate::granges::GRanges;
use crate::meta::MetaData;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct Genes {
    pub granges: GRanges,
    index  : HashMap<String, usize>,
}

/* -------------------------------------------------------------------------- */

impl Genes {
    fn new_impl(granges: GRanges) -> Genes {
        let names = granges.meta.get_column_str  ("names").unwrap();
        let mut index = HashMap::new();
        for i in 0..granges.num_rows() {
            // check if strand is valid
            if granges.strand[i] != '+' && granges.strand[i] != '-' {
                panic!("invalid strand");
            }
            index.insert(names[i].clone(), i);
        }
        let genes = Genes {
            granges,
            index,
        };
        genes
    }

    pub fn new(
        names   : Vec<String>,
        seqnames: Vec<String>,
        tx_from : Vec<usize>,
        tx_to   : Vec<usize>,
        cds_from: Vec<usize>,
        cds_to  : Vec<usize>,
        strand  : Vec<char>,
    ) -> Genes {
        assert!(names.len() == seqnames.len());
        assert!(names.len() == tx_from .len());
        assert!(names.len() == tx_to   .len());
        assert!(names.len() == cds_from.len());
        assert!(names.len() == cds_to  .len());
        assert!(names.len() == strand  .len());
        // construct cds ranges
        let mut cds = vec![];
        for i in 0..cds_from.len() {
            cds.push(Range::new(cds_from[i], cds_to[i]));
        }
        let mut granges = GRanges::new(
            seqnames,
            tx_from,
            tx_to,
            strand,
        );
        granges.meta.add("names", MetaData::StringArray(names)).unwrap();
        granges.meta.add("cds"  , MetaData::RangeArray(cds   )).unwrap();
        Genes::new_impl(granges)
    }

    pub fn names(&self) -> &Vec<String> {
        self.granges.meta.get_column_str("names").unwrap()
    }

    pub fn cds(&self) -> &Vec<Range> {
        self.granges.meta.get_column_range("cds").unwrap()
    }

    pub fn remove(&self, indices: &[usize]) -> Genes {
        let r = self.granges.remove(indices);
        Genes::new_impl(r)
    }

    pub fn remove_overlaps_with(&self, subject: &Genes) -> Genes {
        let r = self.granges.remove_overlaps_with(&subject.granges);
        Genes::new_impl(r)
    }

    pub fn keep_overlaps_with(&self, subject: &Genes) -> Genes {
        let r = self.granges.keep_overlaps_with(&subject.granges);
        Genes::new_impl(r)
    }

    pub fn subset(&self, indices: &[usize]) -> Genes {
        let r = self.granges.subset(indices);
        Genes::new_impl(r)
    }

    pub fn slice(&self, ifrom: usize, ito: usize) -> Genes {
        let r = self.granges.slice(ifrom, ito);
        Genes::new_impl(r)
    }

    pub fn intersection(&self, subject: &Genes) -> Genes {
        let r = self.granges.intersection(&subject.granges);
        Genes::new_impl(r)
    }

    pub fn sort(&self, name: &str, reverse: bool) -> Result<Genes, Box<dyn Error>> {
        let r = self.granges.sort(name, reverse)?;
        Ok(Genes::new_impl(r))
    }

    pub fn filter_genome(&self, genome: &Genome) -> Genes {
        let r = self.granges.filter_genome(genome);
        Genes::new_impl(r)
    }

    pub fn find_gene(&self, name: &str) -> Option<usize> {
        self.index.get(name).cloned()
    }

}

/* -------------------------------------------------------------------------- */

impl PartialEq for Genes {
    fn eq(&self, other: &Self) -> bool {
        self.granges == other.granges
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for Genes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad(&format!("{}", self.granges))
    }
}
