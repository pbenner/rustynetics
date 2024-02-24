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

    use rustynetics::meta::{Meta, MetaData};
    use rustynetics::range::Range;
    use rustynetics::granges::GRanges;

    #[test]
    fn test_meta() {

        let mut rng = thread_rng();

        let n = 20;

        let names = vec!["name", "hello", "genomics", "score", "yeehaa"];
        let data  = vec![
            MetaData::StringArray(
                (0..n).map(|i: i32| format!("Entry {}", i)).collect()),
            MetaData::FloatArray(
                (0..n).map(|_| rng.gen_range(-10.0..100.0)).collect()),
            MetaData::RangeArray(vec![Range::new(0, 1000), Range::new(300, 400), Range::new(20,10000), Range::new(0, 1000), Range::new(300, 400), Range::new(20,10000), Range::new(0, 1000), Range::new(300, 400), Range::new(20,10000), Range::new(0, 1000), Range::new(300, 400), Range::new(20,10000), Range::new(0, 1000), Range::new(300, 400), Range::new(20,10000), Range::new(0, 1000), Range::new(300, 400), Range::new(20,10000), Range::new(0, 1000), Range::new(300, 400)]),
            MetaData::IntArray(vec![3, 2, 5, 3, 2, 5, 3, 2, 5, 3, 2, 5, 3, 2, 5, 3, 2, 5, 3, 2]),
            MetaData::FloatMatrix(vec![
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],
                vec![3.2, 2.2, 5.6],                
                vec![3.2, 2.2, 5.6]]
            )
        ];

        let meta = Meta::new(names, data).unwrap();

        let mut granges = GRanges::new_empty(0);
        
        match granges.import_bed("tests/test_granges.bed", 3, false) {
            Err(r) => panic!("{}", r),
            Ok (_) => (),
        };
        granges.meta = meta;

        match granges.meta.get_column_str_mut("name") {
            Some(v) => v[1] = String::from("hello"),
            _ => ()
        };

        println!("{}", granges);


        match granges.export_bed6("test.bed", false) {
            Ok(_) => (),
            Err(e) => panic!("Error message: {}", e)
        }
    }
}
