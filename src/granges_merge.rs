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

impl GRanges {

    fn merge_impl(r: &GRanges, seqname: &str, entry: &EndPointList) -> GRanges {

        let mut seqnames = vec![];
        let mut from     = vec![];
        let mut to       = vec![];

        let mut i = 0;
        while i < entry.len() {
            // next item must be a start position
            assert!(entry[i].is_start());

            let r_from = entry[i].get_position();
            let r_to   = entry[i].get_position() + 1;

            let mut k = 1;
            while k > 0 {
                i += 1;
                if entry[i].is_start() {
                    k += 1;
                    to[i] = entry[i].get_position();
                } else {
                    k -= 1;
                    to[i] = entry[i].get_position() + 1;
                }
            }

            seqnames.push(seqname);
            from    .push(r_from);
            to      .push(r_to);
        }
        r.append(&GRanges::new(seqnames, from, to, vec![])).unwrap()
    }

    pub fn merge(granges: &[GRanges]) -> GRanges {

        let mut r = GRanges::new_empty();
        let mut rmap: HashMap<String, EndPointList> = HashMap::new();

        for g in granges.iter() {
            for i in 0..g.num_rows() {

                let start = EndPoint::new(g.ranges[i].from, i, true);
                let end   = EndPoint::new(g.ranges[i].to  , i, true);
    
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
            r = GRanges::merge_impl(&r, &seqname, &rmap[&seqname]);
        }
        r
    }
}
