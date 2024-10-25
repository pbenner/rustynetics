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
            vec!["chr4", "chr4"].iter().map(|&x| x.into()).collect(),
            vec![600, 850],
            vec![950, 950],
            vec![]
        );
        let granges2 = GRanges::new(
            vec!["chr4", "chr4", "chr4", "chr4"].iter().map(|&x| x.into()).collect(),
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
