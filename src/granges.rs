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
use std::cmp::Ordering;
use std::collections::HashMap;

use crate::range::Range;
use crate::genome::Genome;
use crate::granges_find_overlaps::find_overlaps;
use crate::meta::Meta;
use crate::error::Error;
use crate::utility::remove_duplicates_int;

/* -------------------------------------------------------------------------- */

pub struct GRangesRow<'a> {
    granges: &'a GRanges,
    row: usize,
}

impl<'a> GRangesRow<'a> {
    fn new(granges: &'a GRanges, row: usize) -> Self {
        GRangesRow { granges, row }
    }
}

impl<'a> fmt::Display for GRangesRow<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "GRangesRow(seqname={}, range=({}, {}), strand={})",
            self.granges.seqnames[self.row],
            self.granges.ranges[self.row].from,
            self.granges.ranges[self.row].to,
            self.granges.strand[self.row] as char
        )
    }
}

/* -------------------------------------------------------------------------- */

pub struct GRanges {
    pub seqnames: Vec<String>,
    pub ranges  : Vec<Range>,
    pub strand  : Vec<char>,
    pub meta    : Meta,
}

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
            meta: Meta::new_empty(n),
        }
    }

    pub fn new_empty(n: usize) -> Self {
        let seqnames = vec![String::new(); n];
        let ranges = vec![Range::new(0, 0); n];
        let strand = vec!['*'; n];
        GRanges {
            seqnames,
            ranges,
            strand,
            meta: Meta::new_empty(n),
        }
    }

    pub fn clone(&self) -> Self {
        let seqnames = self.seqnames.clone();
        let ranges = self.ranges.clone();
        let strand = self.strand.clone();
        let meta = self.meta.clone();
        GRanges {
            seqnames,
            ranges,
            strand,
            meta,
        }
    }

    pub fn length(&self) -> usize {
        self.ranges.len()
    }

    pub fn row(&self, i: usize) -> GRangesRow {
        GRangesRow::new(self, i)
    }

    pub fn append(&self, other: &GRanges) -> Result<Self, Error> {
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
        let n = self.length();
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
        let meta = self.meta.subset(indices);
        GRanges {
            seqnames: result.seqnames,
            ranges: result.ranges,
            strand: result.strand,
            meta,
        }
    }

    pub fn slice(&self, ifrom: usize, ito: usize) -> Self {
        let ifrom = ifrom.min(self.length());
        let ito   = ito  .min(self.length());
        let seqnames = (ifrom..ito)
            .map(|i| self.seqnames[i].clone())
            .collect();
        let from   = (ifrom..ito).map(|i| self.ranges[i].from).collect();
        let to     = (ifrom..ito).map(|i| self.ranges[i].to  ).collect();
        let strand = (ifrom..ito).map(|i| self.strand[i]     ).collect();
        let result = GRanges::new(seqnames, from, to, strand);
        let meta   = self.meta.slice(ifrom, ito);
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
        let seqnames = (0..n)
            .map(|i| {
                let i_q = query_hits[i];
                let i_s = subject_hits[i];
                let gr = self.ranges[i_q].intersection(&s.ranges[i_s]);
                self.seqnames[i_q].clone()
            })
            .collect();
        let strand = (0..n).map(|i| self.strand[query_hits[i]]     ).collect();
        let from   = (0..n).map(|i| self.ranges[query_hits[i]].from).collect();
        let to     = (0..n).map(|i| self.ranges[query_hits[i]].to  ).collect();
        let result = GRanges::new(seqnames, from, to, strand);
        let meta = self.meta.subset(&query_hits);
        GRanges {
            seqnames: result.seqnames,
            ranges: result.ranges,
            strand: result.strand,
            meta,
        }
    }

    pub fn sort(&self, name: &str, reverse: bool) -> Result<Self, Error> {
        if name.is_empty() {
            let mut l = GRangesSort::new(self);
            if reverse {
                l.indices.sort_by(|&i, &j| j.cmp(&i));
            } else {
                l.indices.sort();
            }
            let indices = l.indices.clone();
            Ok(self.subset(&indices))
        } else {
            let j = self.sorted_indices(name, reverse)?;
            Ok(self.subset(&j))
        }
    }

    pub fn filter_genome(&self, genome: &Genome) -> Self {
        let mut idx = Vec::new();
        let mut seqnames = HashMap::new();
        for i in 0..genome.length() {
            seqnames.insert(genome.seqnames[i].clone(), genome.lengths[i]);
        }
        for i in 0..self.length() {
            let length = seqnames.get(&self.seqnames[i]).cloned();
            if let Some(length) = length {
                if self.ranges[i].to <= length {
                    continue;
                }
            }
            idx.push(i);
        }
        self.remove(&idx)
    }

    pub fn filter_strand(&self, s: char) -> Self {
        let mut idx = Vec::new();
        for i in 0..self.length() {
            if self.strand[i] != s {
                idx.push(i);
            }
        }
        self.remove(&idx)
    }

    pub fn set_lengths(&self, n: usize) -> Self {
        let mut s = self.clone();

        for i in 0..s.length() {
            if s.strand[i] == '+' {
                s.ranges[i].to = s.ranges[i].from + n;
            }
            if s.strand[i] == '-' {
                s.ranges[i].from = s.ranges[i].to - n;
            }
        }
        s
    }

    pub fn sorted_indices(&self, name: &str, reverse: bool) -> Result<Vec<usize>, String> {
        let mut indices: Vec<usize> = (0..self.length()).collect();
        match name {
            "" => Err("Invalid sort name".to_string()),
            "seqnames" => {
                indices.sort_by(|&i, &j| {
                    let cmp = self.seqnames[i].cmp(&self.seqnames[j]);
                    if reverse {
                        cmp.reverse()
                    } else {
                        cmp
                    }
                });
                Ok(indices)
            }
            "ranges" => {
                indices.sort_by(|&i, &j| {
                    let cmp = self.ranges[i].from.cmp(&self.ranges[j].from);
                    if cmp == Ordering::Equal {
                        self.ranges[i].to.cmp(&self.ranges[j].to)
                    } else if reverse {
                        cmp.reverse()
                    } else {
                        cmp
                    }
                });
                Ok(indices)
            }
            "strand" => {
                indices.sort_by(|&i, &j| {
                    let cmp = self.strand[i].cmp(&self.strand[j]);
                    if reverse {
                        cmp.reverse()
                    } else {
                        cmp
                    }
                });
                Ok(indices)
            }
            _ => Err("Invalid sort name".to_string()),
        }
    }
}

/* -------------------------------------------------------------------------- */

struct GRangesSort<'a> {
    granges: &'a GRanges,
    indices: Vec<usize>,
}

impl<'a> GRangesSort<'a> {
    fn new(granges: &'a GRanges) -> Self {
        let indices: Vec<usize> = (0..granges.length()).collect();
        GRangesSort { granges, indices }
    }
}

impl<'a> fmt::Display for GRangesSort<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GRangesSort(indices={:?})", self.indices)
    }
}

impl<'a> PartialEq for GRangesSort<'a> {
    fn eq(&self, other: &Self) -> bool {
        self.indices == other.indices
    }
}

impl<'a> Eq for GRangesSort<'a> {}

impl<'a> PartialOrd for GRangesSort<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a> Ord for GRangesSort<'a> {
    fn cmp(&self, other: &Self) -> Ordering {
        let n = self.granges.length();
        for i in 0..n {
            let cmp = self.granges.seqnames[self.indices[i]].cmp(&self.granges.seqnames[other.indices[i]]);
            if cmp != Ordering::Equal {
                return cmp;
            }
        }
        Ordering::Equal
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for GRanges {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad(&format!("{}", self.pretty_string(10).unwrap()))
    }
}
