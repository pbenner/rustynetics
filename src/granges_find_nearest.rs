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
use std::rc::Rc;

use crate::granges::GRanges;
use crate::granges_find_endpoint::{EndPoint, EndPointList};

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct FindNearestHits {
    query_hits  : Vec<i32>,
    subject_hits: Vec<i32>,
    distances   : Vec<i64>,
}

/* -------------------------------------------------------------------------- */

impl FindNearestHits {

    fn new(query_hits: Vec<i32>, subject_hits: Vec<i32>, distances: Vec<i64>) -> Self {
        FindNearestHits {
            query_hits,
            subject_hits,
            distances,
        }
    }

}

/* -------------------------------------------------------------------------- */

fn distance(r1: &EndPoint, r2: &EndPoint) -> (usize, i8) {
    let mut sign = -1;

    let (r1, r2) = if r1.position > r2.position {
        sign = 1;
        (r2, r1)
    } else {
        (r1, r2)
    };

    if r1.start.is_none() || r2.position <= r1.start.as_ref().unwrap().position {
        return (0, sign);
    }

    let d1 = r2.position - r1.start.as_ref().unwrap().position;
    let d2 = r2.position - r1.end  .as_ref().unwrap().position;

    if d1 < d2 {
        (d1, sign)
    } else {
        (d2, sign)
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
            let start = Rc::new(EndPoint {
                position: query.ranges[i].from,
                start   : None,
                end     : None,
                src_idx : i,
                is_query: true,
            });
            let end = Rc::new(EndPoint {
                position: query.ranges[i].to - 1,
                start   : Some(Rc::clone(&start)),
                end     : None,
                src_idx : i,
                is_query: true,
            });

            let entry = rmap.entry(query.seqnames[i].clone()).or_insert(EndPointList::new());
            entry.push(start);
            entry.push(end);
        }

        for i in 0..m {
            let start = Rc::new(EndPoint {
                position: query.ranges[i].from,
                start   : None,
                end     : None,
                src_idx : i,
                is_query: false,
            });
            let end = Rc::new(EndPoint {
                position: query.ranges[i].to - 1,
                start   : Some(Rc::clone(&start)),
                end     : None,
                src_idx : i,
                is_query: false,
            });

            let entry = rmap.entry(subject.seqnames[i].clone()).or_insert(EndPointList::new());

            entry.push(start);
            entry.push(end);
        }

        for (_, entry) in rmap.iter_mut() {
            entry.sort();
        }

        for (_, mut entry) in rmap.iter_mut() {
            EndPointList::find_overlaps_entry(&mut query_hits, &mut subject_hits, &mut entry);

            for i in 0..entry.len() {
                let r = &entry[i];
                // loop over query start regions
                if r.is_query && r.end.is_some() {
                    let mut i1 = i as i64 - 1;
                    let mut i2 = i as i64 + 1;
                    // find k nearest neighbors
                    for _ in 0..k {

                        if i1 < 0 && (i2 as usize) >= entry.len() {
                            break;
                        }
                        // find next subject end to the left
                        for _ in 0..entry.len() {
                            if i1 >= 0 && !entry[i1 as usize].is_query && entry[i1 as usize].end.is_some() {
                                break;
                            }
                            i1 -= 1;
                        }
                        // find next subject start to the right (and drop overlaps)
                        for _ in 0..entry.len() {
                            if (i2 as usize) < entry.len() && !entry[i2 as usize].is_query && entry[i2 as usize].start.is_some() && entry[i2 as usize].position > r.end.as_ref().unwrap().position {
                                break;
                            }
                            i2 += 1;
                        }
                        // are there two elements to compare?
                        let mut ir = -1;
                        let mut dr = -1;

                        if i1 >= 0 && i2 < entry.len() as i64 {
                            let (d1, s1) = distance(r, &entry[i1 as usize]);
                            let (d2, s2) = distance(r, &entry[i2 as usize]);

                            if d1 <= d2 {
                                dr = d1*s1; ir = i1; i1 -= 1;
                            } else {
                                dr = d2*s2; ir = i2; i2 += 1;
                            }
                        } else {
                            if i1 >= 0 {
                                let (d1, s1) = distance(r, &entry[i1 as usize]);
                                dr = d1*s1; ir = i1; i1 -= 1;
                            }
                            if (i2 as usize) < entry.len() {
                                let (d2, s2) = distance(r, &entry[i2 as usize]);
                                dr = d2*s2; ir = i2; i2 += 1;
                            }
                        }
                        if ir != -1 {
                            query_hits  .push(entry[i  as usize].index as i32);
                            subject_hits.push(entry[ir as usize].index as i32);
                            distances   .push(dr);
                        }
                    }
                }
            }
        }

        let find_nearest_hits = FindNearestHits::new(query_hits, subject_hits, distances);

        find_nearest_hits
    }

}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::granges::GRanges;

    #[test]
    fn test_nearest() {

        let mut rQuery = GRanges::new(
            vec!["chr4", "chr4"],
            vec![600, 850],
            vec![950, 950],
            vec![]
        );
        let mut rSubjects = GRanges::new(
            vec!["chr4", "chr4", "chr4", "chr4"],
            vec![100, 200, 300, 400],
            vec![900, 300, 700, 600],
            vec![]
        );

        let r = GRanges::find_nearest(&rQuery, &rSubjects, 2);

        println!("{:?}", r.query_hits);
        println!("{:?}", r.subject_hits);
        println!("{:?}", r.distances);

    }
}
