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

use std::fmt;
use std::cmp::Ordering;
use std::rc::Rc;
use std::cell::RefCell;
use std::clone::Clone;

/* -------------------------------------------------------------------------- */
 
#[derive(Clone)]
pub struct EndPointNode {
    pub position: usize,
    pub start   : Option<EndPoint>,
    pub end     : Option<EndPoint>,
    pub src_idx : usize,
    pub is_query: bool,
}

/* -------------------------------------------------------------------------- */

impl fmt::Debug for EndPointNode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[position: {}, src_idx: {}, is_query: {}]", self.position, self.src_idx, self.is_query
        )
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct EndPoint(Rc<RefCell<EndPointNode>>);

impl EndPoint {

    pub fn new(position: usize, src_idx: usize, is_query: bool) -> EndPoint {
        EndPoint(Rc::new(RefCell::new(EndPointNode {
            position: position,
            start   : None,
            end     : None,
            src_idx : src_idx,
            is_query: is_query,
        })))
    }

    pub fn is_start(&self) -> bool {
        self.0.borrow().start.is_none()
    }

    pub fn is_end(&self) -> bool {
        self.0.borrow().end.is_none()
    }

    pub fn is_query(&self) -> bool {
        self.0.borrow().is_query
    }

    pub fn src_idx(&self) -> usize {
        self.0.borrow().src_idx
    }

    pub fn get_start(&self) -> usize {
        if let Some(start) = &self.0.borrow().start {
            start.0.borrow().position
        } else {
            self.0.borrow().position
        }
    }

    pub fn get_end(&self) -> usize {
        if let Some(end) = &self.0.borrow().end {
            end.0.borrow().position
        } else {
            self.0.borrow().position
        }
    }

    pub fn get_position(&self) -> usize {
        self.0.borrow().position
    }

    pub fn distance(r1: &EndPoint, r2: &EndPoint) -> (i64, i64) {
        let mut sign = -1;

        let (r1, r2) = if r1.get_position() > r2.get_position() {
            sign = 1;
            (r2, r1)
        } else {
            (r1, r2)
        };

        if r1.0.borrow().start.is_none() || r2.get_position() <= r1.0.borrow().start.as_ref().unwrap().0.borrow().position {
            return (0, sign);
        }

        let d1 = r2.get_start() as i64 - r1.get_end  () as i64;
        let d2 = r2.get_end  () as i64 - r1.get_start() as i64;

        if d1 < d2 {
            (d1, sign)
        } else {
            (d2, sign)
        }
    }
}

impl PartialEq for EndPoint {
    fn eq(&self, other: &Self) -> bool {
        self.get_position() == other.get_position()
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
        if self.get_position() != other.get_position() {
            self.get_position().cmp(&other.get_position())
        } else if self.is_start() && other.is_end() {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    }
}

impl Clone for EndPoint {
    fn clone(&self) -> Self {
        EndPoint(Rc::clone(&self.0))
    }
}

impl std::ops::Deref for EndPoint {
    type Target = Rc<RefCell<EndPointNode>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for EndPoint {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
 
/* -------------------------------------------------------------------------- */
 
#[derive(Debug)]
pub struct EndPointList(Vec<EndPoint>);

impl EndPointList {
    pub fn new() -> Self {
            EndPointList(Vec::new())
        }

    pub fn push(&mut self, endpoint: EndPoint) {
        self.0.push(endpoint);
    }

    pub fn remove(&mut self, endpoint: &EndPoint) {
        if let Some(index) = self.0.iter().position(|e| *e == *endpoint) {
            self.0.remove(index);
        }
    }

    pub fn sort(&mut self) {
        self.0.sort();
    }
}

impl std::ops::Deref for EndPointList {
    type Target = Vec<EndPoint>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for EndPointList {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
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

/* -------------------------------------------------------------------------- */

impl EndPointList {

    pub fn find_overlaps_entry(
        query_hits  : &mut Vec<usize>,
        subject_hits: &mut Vec<usize>,
        entry       : &mut EndPointList,
    ) {
        let mut query_list   = EndPointList::new();
        let mut subject_list = EndPointList::new();

        for endpoint in &entry.0 {
            if endpoint.borrow().is_query {
                if endpoint.is_start() {
                    query_list.push(endpoint.clone());
                    for subject_endpoint in &subject_list.0 {
                        query_hits  .push(        endpoint.src_idx());
                        subject_hits.push(subject_endpoint.src_idx());
                    }
                } else {
                    query_list.remove(endpoint.borrow().start.as_ref().unwrap());
                }
            } else {
                if endpoint.is_start() {
                    subject_list.push(endpoint.clone());
                    for query_endpoint in &query_list.0 {
                        query_hits  .push(query_endpoint.src_idx());
                        subject_hits.push(      endpoint.src_idx());
                    }
                } else {
                    subject_list.remove(endpoint.borrow().start.as_ref().unwrap());
                }
            }
        }
    }

}
