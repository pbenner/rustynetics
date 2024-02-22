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

    use rustynetics::meta::{Meta, MetaData};
    use rustynetics::range::Range;
    use rustynetics::granges::GRanges;

    #[test]
    fn test_meta() {

        let names = vec![String::from("hello"), String::from("genomics"), String::from("world"), String::from("yeehaa")];
        let data  = vec![
            MetaData::FloatArray(vec![3.2, 2.2, 5.6, 3.2, 2.2, 5.6, 3.2, 2.2, 5.6, 3.4, 3.2, 2.2, 5.6, 3.2, 2.2, 5.6, 3.2, 2.2, 5.6, 3.4]),
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

        let mut meta = Meta::new(names, data).unwrap();

        let mut granges = GRanges::new_empty(0);
        
        match granges.import_bed("tests/test_granges.bed", 3, false) {
            Err(r) => panic!("{}", r),
            Ok (_) => (),
        };
        granges.meta = meta;

        println!("{}", granges)

    }
}
