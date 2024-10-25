// Copyright (C) 2024 Philipp Benner
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the “Software”), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

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
                (0..n).map(|i: i32| format!("Entry_{}", i)).collect()),
            MetaData::IntArray(
                (0..n).map(|_| rng.gen_range(-100..100)).collect()),
            MetaData::FloatArray(
                (0..n).map(|_| rng.gen_range(-100.0..100.0)).collect()),
            MetaData::RangeArray(
                (0..n).map(|_| rng.gen_range(0..10000)).map(|x| Range::new(x, x+500)).collect()),
            MetaData::FloatMatrix(
                (0..n).map(|_| (0..5).map(|_| rng.gen_range(0.0..1000.0)).collect()).collect())
        ];

        let mut granges1 = GRanges::default();
        
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

        //println!("{}", granges1);
        //println!("{}", granges2);

        assert!(
            granges1 != granges2);

    }

    #[test]
    fn test_meta_bed6() {

        let mut granges1 = GRanges::default();
        let mut granges2 = GRanges::default();

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
