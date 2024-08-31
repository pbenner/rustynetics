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

// Gardiner-Garden and Frommer, J. Mol. Biol. (1987) 196 (2), 261-282:
//   observed * length / (number of C * number of G)
pub fn observed_over_expected_cpg(sequence: &[u8]) -> f64 {
    let mut n_c   = 0;
    let mut n_g   = 0;
    let mut n_cpg = 0;

    for &base in sequence.iter() {
        match base {
            b'c' | b'C' => n_c += 1,
            b'g' | b'G' => n_g += 1,
            _ => {},
        }
    }

    for j in 0..sequence.len() - 1 {
        if (sequence[j] == b'c' || sequence[j] == b'C') && (sequence[j + 1] == b'g' || sequence[j + 1] == b'G') {
            n_cpg += 1;
        }
    }

    if n_cpg != 0 {
        (n_cpg as f64 * sequence.len() as f64) / (n_c as f64 * n_g as f64)
    } else {
        0.0
    }
}