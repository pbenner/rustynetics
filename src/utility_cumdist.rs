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

use std::hash::Hash;
use std::hash::Hasher;
use std::collections::HashMap;
use std::cmp::Ordering;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Copy)]
pub struct OrderedFloat(f64);

/* -------------------------------------------------------------------------- */

impl PartialEq for OrderedFloat {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl Eq for OrderedFloat {}

impl PartialOrd for OrderedFloat {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl Ord for OrderedFloat {
    fn cmp(&self, other: &Self) -> Ordering {
        // Handle NaN values, treating them as less than any other number
        if self.0.is_nan() && other.0.is_nan() {
            Ordering::Equal
        } else if self.0.is_nan() {
            Ordering::Less
        } else if other.0.is_nan() {
            Ordering::Greater
        } else {
            self.0.partial_cmp(&other.0).unwrap()
        }
    }
}

impl Hash for OrderedFloat {
    fn hash<H: Hasher>(&self, state: &mut H) {
        if self.0.is_nan() {
            // All NaNs should hash to the same value
            0x7ff8_0000_0000_0000_u64.hash(state);
        } else {
            // Treat +0.0 and -0.0 as equal
            let bits = if self.0 == 0.0 {
                0_u64
            } else {
                self.0.to_bits()
            };
            bits.hash(state);
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct CumDist {
    pub x: Vec<f64>,
    pub y: Vec<usize>,
}

/* -------------------------------------------------------------------------- */

// Implement methods for CumDist
impl CumDist {

    pub fn new_from_counts(x: Vec<f64>, y: Vec<usize>) -> Self {
        // Sort the arrays based on x values and keep y aligned with x
        let mut indices: Vec<usize> = (0..x.len()).collect();

        indices.sort_by(|&i, &j| x[i].partial_cmp(&x[j]).unwrap_or(Ordering::Equal));

        let     sorted_x: Vec<f64>   = indices.iter().map(|&i| x[i]).collect();
        let mut sorted_y: Vec<usize> = indices.iter().map(|&i| y[i]).collect();

        // Compute the cumulative distribution
        let mut n = 0;
        for val in &mut sorted_y {
            n += *val;
            *val = n;
        }

        CumDist { x: sorted_x, y: sorted_y }
    }

    pub fn new(m: HashMap<OrderedFloat, usize>) -> Self {

        let mut x : Vec<f64>   = vec![];
        let mut y : Vec<usize> = vec![];

        for (key, val) in m.iter() {
            x.push(key.0);
            y.push(*val);
        }

        Self::new_from_counts(x, y)
    }
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::collections::HashMap;

    use crate::utility_cumdist::CumDist;
    use crate::utility_cumdist::OrderedFloat;

    #[test]
    fn test_utility_cumdist_1() {
        let mut data = HashMap::new();
        data.insert(OrderedFloat(3.0), 1);
        data.insert(OrderedFloat(1.0), 2);
        data.insert(OrderedFloat(2.0), 3);

        let cum_dist = CumDist::new(data);

        assert_eq!(cum_dist.x, vec![1.0, 2.0, 3.0]);
        assert_eq!(cum_dist.y, vec![2, 5, 6]);

        println!("{:?}", cum_dist);
    }
}
