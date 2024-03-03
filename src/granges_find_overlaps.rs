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

use crate::granges::GRanges;
use crate::granges_find_endpoint::{EndPoint, EndPointList};

/* -------------------------------------------------------------------------- */

pub fn find_overlaps(query: &GRanges, subject: &GRanges) -> (Vec<usize>, Vec<usize>) {
    let n = query  .num_rows();
    let m = subject.num_rows();

    let mut query_hits   = Vec::new();
    let mut subject_hits = Vec::new();

    let mut rmap: HashMap<String, EndPointList> = HashMap::new();

    for i in 0..n {
        let start = EndPoint::new(query.ranges[i].from, i, true);
        let end   = EndPoint::new(query.ranges[i].to-1, i, true);

        end.borrow_mut().start = Some(start.clone());

        let entry = rmap.entry(query.seqnames[i].clone()).or_insert_with(EndPointList::new);
        entry.push(start);
        entry.push(end);
    }

    for i in 0..m {
        let start = EndPoint::new(subject.ranges[i].from, i, false);
        let end   = EndPoint::new(subject.ranges[i].to-1, i, false);

        end.borrow_mut().start = Some(start.clone());

        let entry = rmap.entry(subject.seqnames[i].clone()).or_insert_with(EndPointList::new);
        entry.push(start);
        entry.push(end);
    }

    for (_, entry) in rmap.iter_mut() {
        entry.sort();
    }

    for (_, mut entry) in rmap.into_iter() {
        EndPointList::find_overlaps_entry(&mut query_hits, &mut subject_hits, &mut entry);
    }

    (query_hits, subject_hits)
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use super::*;
    use crate::granges::GRanges;

    #[test]
    fn test_overlaps_list() {

        let r1 = EndPoint::new(100, 1, true);
        let r2 = EndPoint::new(200, 1, true);
        let r3 = EndPoint::new(300, 1, true);
        let r4 = EndPoint::new(300, 1, true);
        let r5 = r2.clone();

        let mut s = EndPointList::new();

        s.push(r1.clone());
        s.push(r2.clone());
        s.push(r3.clone());
        s.push(r4.clone());

        s.remove(&r5);

        assert!(r1.clone() == s[0]);
        assert!(r3.clone() == s[1]);
        assert!(r3.clone() == s[2]);
        assert!(r4.clone() == s[2]);

    }

    #[test]
    fn test_overlaps() {

        let granges1 = GRanges::new(
            vec!["chr4", "chr4"],
            vec![600, 850],
            vec![950, 950],
            vec![]
        );
        let granges2 = GRanges::new(
            vec!["chr4", "chr4", "chr4", "chr4"],
            vec![100, 200, 300, 400],
            vec![900, 800, 700, 600],
            vec![]
        );

        let (query_hits, subject_hits) = find_overlaps(&granges1, &granges2);

        assert!(query_hits[0] == 0);
        assert!(query_hits[1] == 0);
        assert!(query_hits[2] == 0);
        assert!(query_hits[3] == 1);

        assert!(subject_hits[0] == 0);
        assert!(subject_hits[1] == 1);
        assert!(subject_hits[2] == 2);
        assert!(subject_hits[3] == 0);
    }

}
