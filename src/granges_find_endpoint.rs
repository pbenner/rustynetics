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

#[derive(Debug)]
pub struct Link(Rc<RefCell<EndPoint>>);

impl Clone for Link {
    fn clone(&self) -> Self {
        Link(Rc::clone(self))
    }
}

impl std::ops::Deref for Link {
    type Target = Rc<RefCell<EndPoint>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for Link {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl PartialEq for Link {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}
impl Eq for Link {}

impl PartialOrd for Link {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl Ord for Link {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

/* -------------------------------------------------------------------------- */
 
 #[derive(Clone)]
pub struct EndPoint {
    pub position: usize,
    pub start   : Option<Link>,
    pub end     : Option<Link>,
    pub src_idx : usize,
    pub is_query: bool,
}
 
impl EndPoint {

    pub fn new(position: usize, src_idx: usize, is_query: bool) -> Link {
        Link(Rc::new(RefCell::new(EndPoint {
            position: position,
            start   : None,
            end     : None,
            src_idx : src_idx,
            is_query: is_query,
        })))
    }

    pub fn is_start(&self) -> bool {
        self.start.is_none()
    }

    pub fn is_end(&self) -> bool {
        self.end.is_none()
    }

    pub fn get_start(&self) -> usize {
        if let Some(start) = &self.start {
            start.borrow().position
        } else {
            self.position
        }
    }

    pub fn get_end(&self) -> usize {
        if let Some(end) = &self.end {
            end.borrow().position
        } else {
            self.position
        }
    }

    pub fn distance(r1: &EndPoint, r2: &EndPoint) -> (i64, i64) {
        let mut sign = -1;

        let (r1, r2) = if r1.position > r2.position {
            sign = 1;
            (r2, r1)
        } else {
            (r1, r2)
        };

        if r1.start.is_none() || r2.position <= r1.start.as_ref().unwrap().borrow().position {
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
        } else if self.is_start() && other.is_end() {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    }
}

impl fmt::Debug for EndPoint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[position: {}, src_idx: {}, is_query: {}]", self.position, self.src_idx, self.is_query
        )
    }
}
 
/* -------------------------------------------------------------------------- */
 
#[derive(Debug)]
pub struct EndPointList(Vec<Link>);

impl EndPointList {
    pub fn new() -> Self {
            EndPointList(Vec::new())
        }

    pub fn push(&mut self, endpoint: Link) {
        self.0.push(endpoint);
    }

    pub fn remove(&mut self, endpoint: &Link) {
        if let Some(index) = self.0.iter().position(|e| **e == **endpoint) {
            self.0.remove(index);
        }
    }

    pub fn sort(&mut self) {
        self.0.sort();
    }
}

impl std::ops::Deref for EndPointList {
    type Target = Vec<Link>;

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
                if endpoint.borrow().is_start() {
                    query_list.push(endpoint.clone());
                    for subject_endpoint in &subject_list.0 {
                        query_hits  .push(        endpoint.borrow().src_idx);
                        subject_hits.push(subject_endpoint.borrow().src_idx);
                    }
                } else {
                    query_list.remove(endpoint.borrow().start.as_ref().unwrap());
                }
            } else {
                if endpoint.borrow().is_start() {
                    subject_list.push(endpoint.clone());
                    for query_endpoint in &query_list.0 {
                        query_hits  .push(query_endpoint.borrow().src_idx);
                        subject_hits.push(      endpoint.borrow().src_idx);
                    }
                } else {
                    subject_list.remove(endpoint.borrow().start.as_ref().unwrap());
                }
            }
        }
    }

}
