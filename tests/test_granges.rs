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

    use std::fs::remove_file;

    use rustynetics::granges::{GRanges};

    #[test]
    fn test_granges_bed3() {

        let mut granges = GRanges::new_empty();
        
        assert!(
            granges.import_bed("tests/test_granges.bed", 3, false).is_ok());

        println!("{}", granges);

        assert!(
            granges.export_bed3("tests/test_granges.bed.tmp", false).is_ok());
        assert!(
            remove_file("tests/test_granges.bed.tmp").is_ok());

        }
}
