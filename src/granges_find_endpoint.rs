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
use std::collections::HashMap;
use std::rc::Rc;

use crate::granges::GRanges;
 
/* -------------------------------------------------------------------------- */
 
 #[derive(Clone)]
pub struct EndPoint {
    pub position: usize,
    pub start   : Option<Rc<EndPoint>>,
    pub end     : Option<Rc<EndPoint>>,
    pub src_idx : usize,
    pub is_query: bool,
}
 
impl EndPoint {
    pub fn is_start(&self) -> bool {
        self.start.is_none()
    }

    pub fn is_end(&self) -> bool {
        self.end.is_none()
    }

    pub fn get_start(&self) -> usize {
        if let Some(start) = &self.start {
            start.position
        } else {
            self.position
        }
    }

    pub fn get_end(&self) -> usize {
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
            "<{},{}>", self.get_start(), self.get_end()
        )
    }
}
 
 /* -------------------------------------------------------------------------- */
 
#[derive(Debug)]
pub struct EndPointList(pub Vec<Rc<EndPoint>>);

impl EndPointList {
pub fn new() -> Self {
        EndPointList(Vec::new())
    }

pub fn append(&mut self, endpoint: Rc<EndPoint>) {
    self.0.push(endpoint);
}

pub fn remove(&mut self, endpoint: &Rc<EndPoint>) {
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
