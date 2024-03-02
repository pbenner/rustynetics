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
 
use std::cmp::Ordering;
use std::collections::HashMap;

use crate::granges::GRanges;
use crate::granges_find_endpoint::{EndPoint, EndPointList};

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct FindNearestHits {
    pub query_hits  : Vec<usize>,
    pub subject_hits: Vec<usize>,
    pub distances   : Vec<i64>,
}

/* -------------------------------------------------------------------------- */

impl FindNearestHits {

    fn new(query_hits: Vec<usize>, subject_hits: Vec<usize>, distances: Vec<i64>) -> Self {
        FindNearestHits {
            query_hits,
            subject_hits,
            distances,
        }
    }

    fn len(&self) -> usize {
        self.query_hits.len()
    }

    fn sort(&mut self) {
        let mut r : Vec<FindNearestHitsItem> = (0..self.len()).map(
            |i| FindNearestHitsItem{
                query_hit  : self.  query_hits[i],
                subject_hit: self.subject_hits[i],
                distance   : self.   distances[i]}
            ).collect();

        r.sort();

        for i in 0..self.len() {
            self.  query_hits[i] = r[i].query_hit;
            self.subject_hits[i] = r[i].subject_hit;
            self.   distances[i] = r[i].distance;
        }
    }

}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct FindNearestHitsItem {
    pub query_hit  : usize,
    pub subject_hit: usize,
    pub distance   : i64,
}

/* -------------------------------------------------------------------------- */

impl PartialEq for FindNearestHitsItem {
    fn eq(&self, other: &Self) -> bool {
        self.query_hit   == other.query_hit   &&
        self.subject_hit == other.subject_hit &&
        self.distance    == other.distance
    }
}

impl Eq for FindNearestHitsItem {}

impl PartialOrd for FindNearestHitsItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for FindNearestHitsItem {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.query_hit == other.query_hit {
            if self.distance == other.distance {
                self.subject_hit.cmp(&other.subject_hit)
            } else {
                self.distance.cmp(&other.distance)
            }
        } else {
            self.query_hit.cmp(&other.query_hit)
        }
    }
}

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn find_nearest(query: &GRanges, subject: &GRanges, k: usize) -> FindNearestHits {
        let n = query  .num_rows();
        let m = subject.num_rows();

        let mut query_hits   = Vec::new();
        let mut subject_hits = Vec::new();
        let mut distances    = Vec::new();

        let mut rmap: HashMap<String, EndPointList> = HashMap::new();

        for i in 0..n {
            let start = EndPoint::new(query.ranges[i].from, i, true);
            let end   = EndPoint::new(query.ranges[i].to  , i, true);

            start.borrow_mut().end   = Some(end  .clone());
            end  .borrow_mut().start = Some(start.clone());

            let entry = rmap.entry(query.seqnames[i].clone()).or_insert(EndPointList::new());
            entry.push(start);
            entry.push(end);
        }

        for i in 0..m {
            let start = EndPoint::new(subject.ranges[i].from, i, false);
            let end   = EndPoint::new(subject.ranges[i].to  , i, false);

            start.borrow_mut().end   = Some(end  .clone());
            end  .borrow_mut().start = Some(start.clone());

            let entry = rmap.entry(subject.seqnames[i].clone()).or_insert(EndPointList::new());
            entry.push(start);
            entry.push(end);
        }

        for (_, entry) in rmap.iter_mut() {
            entry.sort();
        }

        for (_, mut entry) in rmap.iter_mut() {
            EndPointList::find_overlaps_entry(&mut query_hits, &mut subject_hits, &mut entry);

            for _ in 0..query_hits.len() {
                distances.push(0);
            }

            for i in 0..entry.len() {
                let r = &entry[i];
                // loop over query start regions
                if r.borrow().is_query && !r.borrow().is_end() {
                    let mut i1 = i as i64 - 1;
                    let mut i2 = i as i64 + 1;
                    // find k nearest neighbors
                    for _ in 0..k {
                        if i1 < 0 && (i2 as usize) >= entry.len() {
                            break;
                        }
                        // find next subject end to the left
                        for _ in 0..entry.len() {
                            if !(i1 >= 0) {
                                break;
                            }
                            if !entry[i1 as usize].borrow().is_query && entry[i1 as usize].borrow().is_end() {
                                break;
                            }
                            i1 -= 1;
                        }
                        // find next subject start to the right (and drop overlaps)
                        for _ in 0..entry.len() {
                            if !((i2 as usize) < entry.len()) {
                                break;
                            }
                            if !entry[i2 as usize].borrow().is_query && entry[i2 as usize].borrow().is_start() && entry[i2 as usize].borrow().position > r.borrow().get_end() {
                                break;
                            }
                            i2 += 1;
                        }
                        // are there two elements to compare?
                        let mut ir = -1;
                        let mut dr = -1;

                        if i1 >= 0 && i2 < entry.len() as i64 {
                            let (d1, s1) = EndPoint::distance(&r.borrow(), &entry[i1 as usize].borrow());
                            let (d2, s2) = EndPoint::distance(&r.borrow(), &entry[i2 as usize].borrow());

                            if d1 <= d2 {
                                dr = d1*s1; ir = i1; i1 -= 1;
                            } else {
                                dr = d2*s2; ir = i2; i2 += 1;
                            }
                        } else {
                            if i1 >= 0 {
                                let (d1, s1) = EndPoint::distance(&r.borrow(), &entry[i1 as usize].borrow());
                                dr = d1*s1; ir = i1; i1 -= 1;
                            }
                            if (i2 as usize) < entry.len() {
                                let (d2, s2) = EndPoint::distance(&r.borrow(), &entry[i2 as usize].borrow());
                                dr = d2*s2; ir = i2; i2 += 1;
                            }
                        }
                        if ir != -1 {
                            query_hits  .push(entry[i  as usize].borrow().src_idx);
                            subject_hits.push(entry[ir as usize].borrow().src_idx);
                            distances   .push(dr);
                        }
                    }
                }
            }
        }

        let mut find_nearest_hits = FindNearestHits::new(query_hits, subject_hits, distances);

        find_nearest_hits.sort();
        find_nearest_hits
    }

}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::granges::GRanges;

    #[test]
    fn test_nearest() {

        let r_query = GRanges::new(
            vec!["chr4", "chr4"],
            vec![600, 850],
            vec![950, 950],
            vec![]
        );
        let r_subjects = GRanges::new(
            vec!["chr4", "chr4", "chr4", "chr4"],
            vec![100, 200, 300, 400],
            vec![900, 300, 700, 600],
            vec![]
        );

        let rq = vec![0, 0, 0,   0, 1,   1,   1];
        let rs = vec![0, 2, 3,   1, 0,   2,   3];
        let rd = vec![0, 0, 0, 300, 0, 150, 250];

        let r = GRanges::find_nearest(&r_query, &r_subjects, 2);

        assert!(r.query_hits  .len() == rq.len());
        assert!(r.subject_hits.len() == rs.len());
        assert!(r.distances   .len() == rd.len());

        for i in 0..rd.len() {
            assert!(rq[i] == r.query_hits  [i]);
            assert!(rs[i] == r.subject_hits[i]);
            assert!(rd[i] == r.distances   [i]);
        }

    }
}
