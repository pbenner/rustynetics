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

impl GRanges {

    fn merge_impl(r: &GRanges, seqname: String, entry: &EndPointList) -> GRanges {

        let mut seqnames = vec![];
        let mut from     = vec![];
        let mut to       = vec![];

        let mut i = 0;
        while i < entry.len() {
            // next item must be a start position
            assert!(entry[i].is_start());

            let     r_from = entry[i].get_position();
            let mut r_to   = entry[i].get_position() + 1;

            // k: number open intervals
            let mut k = 1;
            while k > 0 {
                // go to next item
                i += 1;
                if entry[i].is_start() {
                    // start of an interval
                    k += 1;
                    r_to = entry[i].get_position();
                } else {
                    // end of an interval
                    k -= 1;
                    r_to = entry[i].get_position() + 1;
                }
            }
            seqnames.push(seqname.clone());
            from    .push(r_from);
            to      .push(r_to);
            // go to next item
            i += 1;
        }
        r.append(&GRanges::new(seqnames, from, to, vec![])).unwrap()
    }

    pub fn merge(granges: &[GRanges]) -> GRanges {

        let mut r = GRanges::default();
        let mut rmap: HashMap<String, EndPointList> = HashMap::new();

        for g in granges.iter() {
            for i in 0..g.num_rows() {

                let start = EndPoint::new(g.ranges[i].from, i, true);
                let end   = EndPoint::new(g.ranges[i].to-1, i, true);
    
                end.borrow_mut().start = Some(start.clone());
    
                let entry = rmap.entry(g.seqnames[i].clone()).or_insert_with(EndPointList::new);
                entry.push(start);
                entry.push(end);
            }
        }

        let mut seqnames: Vec<String> = rmap.keys().cloned().collect();
        seqnames.sort();

        for entry in rmap.values_mut() {
            entry.sort();
        }

        for seqname in seqnames {
            r = GRanges::merge_impl(&r, seqname.clone(), &rmap[&seqname]);
        }
        r
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::granges::GRanges;

    #[test]
    fn test_merge() {

        let seqnames = vec!["chr1", "chr1", "chr1", "chr2", "chr2"].iter().map(|&x| x.into()).collect();
        let from     = vec![ 6, 10, 24,  6, 10];
        let to       = vec![21, 31, 81, 21, 31];
        let strand   = vec![];

        let granges  = GRanges::new(seqnames, from, to, strand);

        let r = GRanges::merge(&[granges]);

        assert!(r.num_rows() == 2);

        assert!(r.ranges[0].from == 6);
        assert!(r.ranges[0].to   == 81);

        assert!(r.ranges[1].from == 6);
        assert!(r.ranges[1].to   == 31);

    }

}
