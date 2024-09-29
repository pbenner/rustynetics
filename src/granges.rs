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

/* -------------------------------------------------------------------------- */

use std::fmt;
use std::collections::HashMap;
use std::error::Error;

use crate::range::Range;
use crate::genome::Genome;
use crate::granges_row::GRangesRow;
use crate::granges_find_overlaps::find_overlaps;
use crate::meta::Meta;
use crate::utility::remove_duplicates_int;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug)]
pub struct GRanges {
    pub seqnames: Vec<String>,
    pub ranges  : Vec<Range>,
    pub strand  : Vec<char>,
    pub meta    : Meta,
}

/* -------------------------------------------------------------------------- */

impl Default for GRanges {
    fn default() -> Self {
        let seqnames = vec!["".to_string(); 0];
        let ranges   = vec![Range::new(0, 0); 0];
        let strand   = vec!['*'; 0];
        GRanges {
            seqnames,
            ranges,
            strand,
            meta: Meta::default(),
        }
    }
}

/* -------------------------------------------------------------------------- */

impl GRanges {
    pub fn new(seqnames: Vec<String>, from: Vec<usize>, to: Vec<usize>, strand: Vec<char>) -> Self {
        let n = seqnames.len();
        if from.len() != n || to.len() != n || (strand.len() != 0 && strand.len() != n) {
            panic!("NewGRanges(): invalid arguments!");
        }
        let strand = if strand.len() == 0 {
            vec!['*'; n]
        } else {
            strand
        };
        let ranges = from
            .iter()
            .zip(to.iter())
            .map(|(&f, &t)| Range::new(f, t))
            .collect();
        GRanges {
            seqnames,
            ranges,
            strand,
            meta: Meta::default(),
        }
    }

    pub fn num_rows(&self) -> usize {
        self.ranges.len()
    }

    pub fn row(&self, i: usize) -> GRangesRow {
        GRangesRow::new(self, i)
    }

    pub fn append(&self, other: &GRanges) -> Result<Self, Box<dyn Error>> {
        let mut seqnames = self.seqnames.clone();
        seqnames.extend(other.seqnames.iter().cloned());
        let mut ranges = self.ranges.clone();
        ranges.extend(other.ranges.iter().cloned());
        let mut strand = self.strand.clone();
        strand.extend(other.strand.iter().cloned());
        let meta = self.meta.clone().append(&other.meta)?;
        Ok(GRanges {
            seqnames,
            ranges,
            strand,
            meta,
        })
    }

    pub fn remove(&self, indices: &[usize]) -> Self {
        if indices.is_empty() {
            return self.clone();
        }
        let indices = remove_duplicates_int(indices);
        let mut indices = indices.to_vec();
        indices.sort();
        let n = self.num_rows();
        let m = n - indices.len();
        let mut idx = Vec::with_capacity(m);
        let mut k = 0;
        for i in 0..n {
            while k < indices.len() - 1 && i > indices[k] {
                k += 1;
            }
            if i != indices[k] {
                idx.push(i);
            }
        }
        let result = self.subset(&idx);
        let meta = self.meta.subset(&idx);
        GRanges {
            seqnames: result.seqnames,
            ranges: result.ranges,
            strand: result.strand,
            meta,
        }
    }

    pub fn remove_overlaps_with(&self, subject: &GRanges) -> Self {
        let (query_hits, _) = find_overlaps(self, subject);
        self.remove(&query_hits)
    }

    pub fn keep_overlaps_with(&self, subject: &GRanges) -> Self {
        let (query_hits, _) = find_overlaps(self, subject);
        let query_hits = remove_duplicates_int(&query_hits);
        let mut query_hits = query_hits.to_vec();
        query_hits.sort();
        self.subset(&query_hits)
    }

    pub fn subset(&self, indices: &[usize]) -> Self {
        let seqnames = indices.iter().map(|&i| self.seqnames[i].clone()).collect();
        let from     = indices.iter().map(|&i| self.ranges  [i].from   ).collect();
        let to       = indices.iter().map(|&i| self.ranges  [i].to     ).collect();
        let strand   = indices.iter().map(|&i| self.strand  [i]        ).collect();
        let result   = GRanges::new(seqnames, from, to, strand);
        let meta     = self.meta.subset(indices);
        GRanges {
            seqnames: result.seqnames,
            ranges  : result.ranges,
            strand  : result.strand,
            meta,
        }
    }

    pub fn slice(&self, ifrom: usize, ito: usize) -> Self {
        let ifrom = ifrom.min(self.num_rows());
        let ito   = ito  .min(self.num_rows());
        let seqnames = (ifrom..ito).map(|i| self.seqnames[i].clone()).collect();
        let from     = (ifrom..ito).map(|i| self.ranges  [i].from   ).collect();
        let to       = (ifrom..ito).map(|i| self.ranges  [i].to     ).collect();
        let strand   = (ifrom..ito).map(|i| self.strand  [i]        ).collect();
        let result   = GRanges::new(seqnames, from, to, strand);
        let meta     = self.meta.slice(ifrom, ito);
        GRanges {
            seqnames: result.seqnames,
            ranges: result.ranges,
            strand: result.strand,
            meta,
        }
    }

    pub fn intersection(&self, s: &GRanges) -> Self {
        let (query_hits, subject_hits) = find_overlaps(self, s);
        let n = query_hits.len();

        let mut seqnames = Vec::with_capacity(n);
        let mut from     = Vec::with_capacity(n);
        let mut to       = Vec::with_capacity(n);
        let mut strand   = Vec::with_capacity(n);

        for i in 0..n {
            let i_q =   query_hits[i];
            let i_s = subject_hits[i];
            let gr  = self.ranges[i_q].intersection(&s.ranges[i_s]);

            seqnames.push(self.seqnames[i_q].clone());
            strand  .push(self.strand  [i_q]);
            from    .push(gr.from);
            to      .push(gr.to  );
        }
        let mut granges = GRanges::new(seqnames, from, to, strand);

        granges.meta = self.meta.subset(&query_hits);
        granges
    }

    pub fn filter_genome(&self, genome: &Genome) -> Self {
        let mut idx = Vec::new();
        let mut seqnames = HashMap::new();
        for i in 0..genome.len() {
            seqnames.insert(genome.seqnames[i].clone(), genome.lengths[i]);
        }
        for i in 0..self.num_rows() {
            let num_rows = seqnames.get(&self.seqnames[i]).cloned();
            if let Some(num_rows) = num_rows {
                if self.ranges[i].to <= num_rows {
                    continue;
                }
            }
            idx.push(i);
        }
        self.remove(&idx)
    }

    pub fn filter_strand(&self, s: char) -> Self {
        let mut idx = Vec::new();
        for i in 0..self.num_rows() {
            if self.strand[i] != s {
                idx.push(i);
            }
        }
        self.remove(&idx)
    }

    pub fn set_num_rowss(&self, n: usize) -> Self {
        let mut s = self.clone();

        for i in 0..s.num_rows() {
            if s.strand[i] == '+' {
                s.ranges[i].to = s.ranges[i].from + n;
            }
            if s.strand[i] == '-' {
                s.ranges[i].from = s.ranges[i].to - n;
            }
        }
        s
    }

}

/* -------------------------------------------------------------------------- */

impl PartialEq for GRanges {
    fn eq(&self, other: &Self) -> bool {
        self.seqnames == other.seqnames &&
        self.ranges   == other.ranges   &&
        self.strand   == other.strand   &&
        self.meta     == other.meta
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for GRanges {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad(&format!("{}", self.format_pretty(10).unwrap()))
    }
}
