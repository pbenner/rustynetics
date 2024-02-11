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

#[cfg(test)]
mod tests {

    use rustynetics::genome::{Genome};

    #[test]
    fn test_genome() {
        let mut genome = Genome::new(vec![], vec![]);

        let r = genome.import("data/hg19.genome");

        assert_eq!(r.is_ok(), true);

        assert_eq!(genome.length(), 25)
    }
}
