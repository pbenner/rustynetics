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

    use approx::assert_relative_eq;

    use rustynetics::bigwig::{BigWigFile, BigWigQueryType};

    #[test]
    fn bigwig_test_1() {

        let result =  BigWigFile::new_reader("tests/test_bigwig_1.bw");

        assert!(result.is_ok());

        if let Ok(mut bw) = result {

            assert_eq!(bw.genome().len(), 2);

            assert_eq!(bw.genome().seqnames[0], "test1");
            assert_eq!(bw.genome().seqnames[1], "test2");

            let mut sum_from = 0;
            let mut sum_to   = 0;
            let mut sum_min  = 0.0;
            let mut sum_max  = 0.0;

            for result in bw.query("test1", 0, 100, 10) {
                if let Ok(item) = result {
                    assert_eq!(item.data.chrom, "test1");
                    sum_from += item.data.from;
                    sum_to   += item.data.to;
                    sum_min  += item.data.statistics.min;
                    sum_max  += item.data.statistics.max;
                }
            }

            assert_eq!(sum_from, 450);
            assert_eq!(sum_to  , 550);

            assert_relative_eq!(sum_min, 49.5, epsilon = 1e-6);
            assert_relative_eq!(sum_max, 49.5, epsilon = 1e-6);
        }
    }

    #[test]
    fn bigwig_test_2() {

        let result = BigWigFile::new_reader("tests/test_bigwig_2.bw");

        assert!(result.is_ok());

        if let Ok(mut bw) = result {

            // Query raw data
            let r1 : Vec<BigWigQueryType> = bw.query("chrY", 1838100, 1838600, 100)
                .filter_map(|r| r.ok())
                .collect();

            // Query zoom data
            let r2 : Vec<BigWigQueryType> = bw.query("chrY", 1838100, 1838600, 400)
                .filter_map(|r| r.ok())
                .collect();

            assert_eq!(r1.len(), 5);
            assert_eq!(r2.len(), 2);

            assert_eq!(r1[0].data.from, 1838100);
            assert_eq!(r1[0].data.to  , 1838200);
            assert_eq!(r1[0].data.statistics.valid, 1.0);
            assert_eq!(r1[0].data.statistics.min  , 1.0);
            assert_eq!(r1[0].data.statistics.max  , 1.0);
            assert_eq!(r1[0].data.statistics.sum  , 1.0);

            assert_eq!(r1[1].data.from, 1838200);
            assert_eq!(r1[1].data.to  , 1838300);
            assert_eq!(r1[1].data.statistics.valid, 1.0);
            assert_eq!(r1[1].data.statistics.min  , 1.0);
            assert_eq!(r1[1].data.statistics.max  , 1.0);
            assert_eq!(r1[1].data.statistics.sum  , 1.0);

            assert_eq!(r1[2].data.from, 1838300);
            assert_eq!(r1[2].data.to  , 1838400);
            assert_eq!(r1[2].data.statistics.valid, 1.0);
            assert_eq!(r1[2].data.statistics.min  , 0.0);
            assert_eq!(r1[2].data.statistics.max  , 0.0);
            assert_eq!(r1[2].data.statistics.sum  , 0.0);

            assert_eq!(r1[3].data.from, 1838400);
            assert_eq!(r1[3].data.to  , 1838500);
            assert_eq!(r1[3].data.statistics.valid, 1.0);
            assert_eq!(r1[3].data.statistics.min  , 0.0);
            assert_eq!(r1[3].data.statistics.max  , 0.0);
            assert_eq!(r1[3].data.statistics.sum  , 0.0);

            assert_eq!(r1[4].data.from, 1838500);
            assert_eq!(r1[4].data.to  , 1838600);
            assert_eq!(r1[4].data.statistics.valid, 1.0);
            assert_eq!(r1[4].data.statistics.min  , 0.0);
            assert_eq!(r1[4].data.statistics.max  , 0.0);
            assert_eq!(r1[4].data.statistics.sum  , 0.0);

            assert_eq!(r2[0].data.from, 1838000);
            assert_eq!(r2[0].data.to  , 1838400);
            assert_eq!(r2[0].data.statistics.valid, 4.0);
            assert_eq!(r2[0].data.statistics.min  , 0.0);
            assert_eq!(r2[0].data.statistics.max  , 1.0);
            assert_eq!(r2[0].data.statistics.sum  , 2.0);

            assert_eq!(r2[1].data.from, 1838400);
            assert_eq!(r2[1].data.to  , 1838800);
            assert_eq!(r2[1].data.statistics.valid, 4.0);
            assert_eq!(r2[1].data.statistics.min  , 0.0);
            assert_eq!(r2[1].data.statistics.max  , 0.0);
            assert_eq!(r2[1].data.statistics.sum  , 0.0);

        }
    }

}
