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

    use rand::{thread_rng, Rng};
    use std::fs::remove_file;

    use rustynetics::meta::{Meta, MetaData};
    use rustynetics::range::Range;
    use rustynetics::granges::GRanges;

    #[test]
    fn test_meta() {

        let mut rng = thread_rng();

        let n = 20;

        let names = vec!["name", "score", "float", "range", "matrix"];
        let data  = vec![
            MetaData::StringArray(
                (0..n).map(|i: i32| format!("Entry {}", i)).collect()),
            MetaData::IntArray(
                (0..n).map(|_| rng.gen_range(-100..100)).collect()),
            MetaData::FloatArray(
                (0..n).map(|_| rng.gen_range(-100.0..100.0)).collect()),
            MetaData::RangeArray(
                (0..n).map(|_| rng.gen_range(0..10000)).map(|x| Range::new(x, x+500)).collect()),
            MetaData::FloatMatrix(
                (0..n).map(|_| (0..5).map(|_| rng.gen_range(0.0..1000.0)).collect()).collect())
        ];

        let mut granges1 = GRanges::new_empty();
        
        assert!(
            granges1.import_bed("tests/test_meta.bed", 3, false).is_ok());

        granges1.meta = Meta::new(names, data).unwrap();

        let granges2 = granges1.clone();

        assert!(
            granges1 == granges2);

        match granges1.meta.get_column_str_mut("name") {
            Some(v) => v[1] = String::from("Test"),
            _ => ()
        };

        println!("{}", granges1);
        println!("{}", granges2);

        assert!(
            granges1 != granges2);

    }

    #[test]
    fn test_meta_bed6() {

        let mut granges1 = GRanges::new_empty();
        let mut granges2 = GRanges::new_empty();

        assert!(
            granges1.import_bed("tests/test_meta.bed", 6, false).is_ok());
        assert!(
            granges1.export_bed6("tests/test_meta.bed.tmp", false).is_ok());
        assert!(
            granges2.import_bed6("tests/test_meta.bed.tmp", false).is_ok());
        assert!(
            granges1 == granges2);
        assert!(
            remove_file("tests/test_meta.bed.tmp").is_ok());

    }

}
