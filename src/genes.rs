/* Copyright (C) 2024 Philipp Benner
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

use std::collections::HashMap;

use crate::range::Range;
use crate::genome::Genome;
use crate::granges::GRanges;
use crate::meta::MetaData;
use crate::error::Error;

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
struct Genes<'a> {
    granges: GRanges,
    names  : &'a Vec<String>,
    cds    : &'a Vec<Range>,
    index  : HashMap<String, usize>,
}

/* -------------------------------------------------------------------------- */

impl<'a> Genes<'a> {
    fn new(granges: GRanges) -> Genes<'a> {
        let names = granges.meta.get_column_str  ("names").unwrap();
        let cds   = granges.meta.get_column_range("cds"  ).unwrap();

        let mut index = HashMap::new();
        for i in 0..granges.num_rows() {
            // check if strand is valid
            if granges.strand[i] != '+' && granges.strand[i] != '-' {
                panic!("NewGenes(): Invalid strand!");
            }
            index.insert(names[i].clone(), i);
        }
        Genes {
            granges,
            names,
            cds,
            index,
        }
    }

    fn new_genes(
        names   : Vec<&str>,
        seqnames: Vec<&str>,
        tx_from : Vec<usize>,
        tx_to   : Vec<usize>,
        cds_from: Vec<usize>,
        cds_to  : Vec<usize>,
        strand  : Vec<char>,
    ) -> Genes<'a> {
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
        granges.meta.add_meta("names", MetaData::StringArray(names.iter().map(|&x| x.into()).collect()));
        granges.meta.add_meta("cds"  , MetaData:: RangeArray(cds  .clone()));
        Genes::new(granges)
    }

    fn clone(&self) -> Genes {
        Genes::new(self.granges.clone())
    }

    fn remove(&self, indices: &[usize]) -> Genes {
        let r = self.granges.remove(indices);
        Genes::new(r)
    }

    fn remove_overlaps_with(&self, subject: &Genes) -> Genes {
        let r = self.granges.remove_overlaps_with(&subject.granges);
        Genes::new(r)
    }

    fn keep_overlaps_with(&self, subject: &Genes) -> Genes {
        let r = self.granges.keep_overlaps_with(&subject.granges);
        Genes::new(r)
    }

    fn subset(&self, indices: &[usize]) -> Genes {
        let r = self.granges.subset(indices);
        Genes::new(r)
    }

    fn slice(&self, ifrom: usize, ito: usize) -> Genes {
        let r = self.granges.slice(ifrom, ito);
        Genes::new(r)
    }

    fn intersection(&self, subject: &Genes) -> Genes {
        let r = self.granges.intersection(&subject.granges);
        Genes::new(r)
    }

    fn sort(&self, name: &str, reverse: bool) -> Result<Genes, Error> {
        let r = self.granges.sort(name, reverse)?;
        Ok(Genes::new(r))
    }

    fn filter_genome(&self, genome: &Genome) -> Genes {
        let r = self.granges.filter_genome(genome);
        Genes::new(r)
    }

    fn find_gene(&self, name: &str) -> Option<usize> {
        self.index.get(name).cloned()
    }

}
