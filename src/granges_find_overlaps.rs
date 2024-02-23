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
use std::rc::Rc;

use crate::granges::{GRanges};

/* -------------------------------------------------------------------------- */

#[derive(Debug, Clone)]
struct EndPoint {
    position: usize,
    start   : Option<Rc<EndPoint>>,
    end     : Option<Rc<EndPoint>>,
    src_idx : usize,
    is_query: bool,
}

impl EndPoint {
    fn is_start(&self) -> bool {
        self.start.is_none()
    }

    fn is_end(&self) -> bool {
        self.end.is_none()
    }

    fn get_start(&self) -> usize {
        if let Some(start) = &self.start {
            start.position
        } else {
            self.position
        }
    }

    fn get_end(&self) -> usize {
        if let Some(end) = &self.end {
            end.position
        } else {
            self.position
        }
    }
}

impl PartialEq for EndPoint {
    fn eq(&self, other: &Self) -> bool {
        self.position == other.position
    }
}

impl Eq for EndPoint {}

impl PartialOrd for EndPoint {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for EndPoint {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.position != other.position {
            self.position.cmp(&other.position)
        } else if self.start.is_none() && other.start.is_some() {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    }
}

#[derive(Debug)]
struct EndPointList(Vec<Rc<EndPoint>>);

impl EndPointList {
    fn new() -> Self {
        EndPointList(Vec::new())
    }

    fn append(&mut self, endpoint: Rc<EndPoint>) {
        self.0.push(endpoint);
    }

    fn remove(&mut self, endpoint: &Rc<EndPoint>) {
        if let Some(index) = self.0.iter().position(|e| **e == **endpoint) {
            self.0.remove(index);
        }
    }
}

impl PartialEq for EndPointList {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl Eq for EndPointList {}

impl PartialOrd for EndPointList {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for EndPointList {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

fn find_overlaps_entry(
    query_hits: &mut Vec<usize>,
    subject_hits: &mut Vec<usize>,
    entry: &mut EndPointList,
) {
    let mut query_list = EndPointList::new();
    let mut subject_list = EndPointList::new();

    for endpoint in &entry.0 {
        if endpoint.is_query {
            if endpoint.start.is_none() {
                query_list.append(endpoint.clone());
                for subject_endpoint in &subject_list.0 {
                    query_hits.push(endpoint.src_idx);
                    subject_hits.push(subject_endpoint.src_idx);
                }
            } else {
                query_list.remove(endpoint.start.as_ref().unwrap());
            }
        } else {
            if endpoint.start.is_none() {
                subject_list.append(endpoint.clone());
                for query_endpoint in &query_list.0 {
                    query_hits.push(query_endpoint.src_idx);
                    subject_hits.push(endpoint.src_idx);
                }
            } else {
                subject_list.remove(endpoint.start.as_ref().unwrap());
            }
        }
    }
}

pub fn find_overlaps(query: &GRanges, subject: &GRanges) -> (Vec<usize>, Vec<usize>) {
    let n = query.length();
    let m = subject.length();

    let mut query_hits = Vec::new();
    let mut subject_hits = Vec::new();

    let mut rmap: HashMap<String, EndPointList> = HashMap::new();

    for i in 0..n {
        let start = Rc::new(EndPoint {
            position: query.ranges[i].from,
            start: None,
            end: None,
            src_idx: i,
            is_query: true,
        });
        let end = Rc::new(EndPoint {
            position: query.ranges[i].to - 1,
            start: Some(Rc::clone(&start)),
            end: None,
            src_idx: i,
            is_query: true,
        });
        let entry = rmap.entry(query.seqnames[i].clone()).or_insert_with(EndPointList::new);
        entry.append(start);
        entry.append(end);
    }

    for i in 0..m {
        let start = Rc::new(EndPoint {
            position: subject.ranges[i].from,
            start: None,
            end: None,
            src_idx: i,
            is_query: false,
        });
        let end = Rc::new(EndPoint {
            position: subject.ranges[i].to - 1,
            start: Some(Rc::clone(&start)),
            end: None,
            src_idx: i,
            is_query: false,
        });
        let entry = rmap.entry(subject.seqnames[i].clone()).or_insert_with(EndPointList::new);
        entry.append(start);
        entry.append(end);
    }

    for (_, entry) in rmap.iter_mut() {
        entry.0.sort();
    }

    for (_, mut entry) in rmap.into_iter() {
        find_overlaps_entry(&mut query_hits, &mut subject_hits, &mut entry);
    }

    (query_hits, subject_hits)
}
