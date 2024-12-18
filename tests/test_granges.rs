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

    use std::fs::remove_file;
    use rand::{thread_rng, Rng};

    use rustynetics::range::Range;
    use rustynetics::granges::GRanges;
    use rustynetics::meta::{Meta, MetaData};

    #[test]
    fn test_granges_bed3() {

        let mut granges1 = GRanges::default();
        let mut granges2 = GRanges::default();

        // Import given granges
        assert!(
            granges1.import_bed("tests/test_granges_1.bed", 3, false).is_ok());

        assert_eq!(
            granges1.num_rows(), 20);

        // Export to new file and import again
        assert!(
            granges1.export_bed3("tests/test_granges_1.bed.tmp", false).is_ok());
        assert!(
            granges2.import_bed3("tests/test_granges_1.bed.tmp", false).is_ok());
        assert!(
            granges1 == granges2);

        // Modify granges and test for inequality
        granges2.seqnames[2] = String::from("test");
        assert!(
            granges1 != granges2);

        assert!(
            remove_file("tests/test_granges_1.bed.tmp").is_ok());

    }

    #[test]
    fn test_granges_bed9() {

        let mut granges1 = GRanges::default();
        let mut granges2 = GRanges::default();

        // Import given granges
        assert!(
            granges1.import_bed("tests/test_granges_2.bed", 9, false).is_ok());

        assert_eq!(
            granges1.num_rows(), 9);

        // Export to new file and import again
        assert!(
            granges1.export_bed9("tests/test_granges_2.bed.tmp", false).is_ok());
        assert!(
            granges2.import_bed9("tests/test_granges_2.bed.tmp", false).is_ok());
        assert!(
            granges1 == granges2);

        // Modify granges and test for inequality
        granges2.seqnames[2] = String::from("test");
        assert!(
            granges1 != granges2);

        assert!(
            remove_file("tests/test_granges_2.bed.tmp").is_ok());

    }

    #[test]
    fn test_granges_intersection() {

        let mut granges1 = GRanges::default();
        let mut granges2 = GRanges::default();

        // Import given granges
        assert!(
            granges1.import_bed("tests/test_granges_1.bed", 3, false).is_ok()
        );
        assert!(
            granges2.import_bed("tests/test_granges_1.bed", 3, false).is_ok()
        );
        granges1 = granges1.subset(&[0,1]);
        granges2 = granges1.subset(&[0,1]);

        let granges3 = granges1.intersection(&granges2);

        assert!(
            granges3.ranges[0] == Range::new(100000266, 100000291)
        );
        assert!(
            granges3.ranges[1] == Range::new(100000271, 100000291)
        );
        assert!(
            granges3.ranges[2] == Range::new(100000271, 100000291)
        );
        assert!(
            granges3.ranges[3] == Range::new(100000271, 100000296)
        );

    }

    #[test]
    fn test_granges_table() {

        let r = |x: f64| (x * 100.0).round() / 100.0;

        let mut rng = thread_rng();

        let n = 20;

        let names = vec!["name", "score", "float", "matrix"];
        let data  = vec![
            MetaData::StringArray(
                (0..n).map(|i: i32| format!("Entry_{}", i)).collect()),
            MetaData::IntArray(
                (0..n).map(|_| rng.gen_range(-100..100)).collect()),
            MetaData::FloatArray(
                (0..n).map(|_| r(rng.gen_range(-100.0..100.0))).collect()),
            MetaData::FloatMatrix(
                (0..n).map(|_| (0..5).map(|_| r(rng.gen_range(0.0..1000.0))).collect()).collect())
        ];

        let mut granges1 = GRanges::default();
        let mut granges2 = GRanges::default();

        assert!(
            granges1.import_bed("tests/test_meta.bed", 3, false).is_ok());

        granges1.meta = Meta::new(names, data).unwrap();

        // Export to new file
        assert!(
            granges1.export_table("tests/test_granges.table.tmp", false, &[]).is_ok());
        // Import table again
        if let Err(v) = granges2.import_table(
                "tests/test_granges.table.tmp",
                &["name", "score", "float", "matrix"],
                &["String", "Int", "Float", "Vec<Float>"],
                false) {
            println!("{}", v);
        }

        assert!(
            remove_file("tests/test_granges.table.tmp").is_ok());

            assert!(
            granges1 == granges2);

    }

}
