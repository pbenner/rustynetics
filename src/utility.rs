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

/* -------------------------------------------------------------------------- */

use num::traits::PrimInt;

use std::collections::HashSet;

/* -------------------------------------------------------------------------- */

pub fn remove_duplicates_int(s: &[usize]) -> Vec<usize> {
    let mut m: HashSet<usize> = HashSet::new();
    let mut r: Vec<usize> = Vec::new();
    for &v in s {
        if m.insert(v) {
            r.push(v);
        }
    }
    r
}

/* -------------------------------------------------------------------------- */

// Helper function for integer division rounding up
pub fn div_int_up<T : PrimInt>(a: T, b: T) -> T {
    (a + b - T::one()) / b
}

// Helper function for integer division rounding down
pub fn div_int_down<T : PrimInt>(n: T, d: T) -> T {
    n / d
}
