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

    use approx::assert_relative_eq;

    use rustynetics::bbi::BbiQueryType;
    use rustynetics::bigwig::BigWigFile;

    #[test]
    fn bigwig_test_1() {

        let result =  BigWigFile::new_reader("tests/test_bigwig_1.bw");

        assert!(result.is_ok());

        if let Ok(mut bw) = result {

            assert_eq!(bw.genome().len(), 2);

            assert_eq!(bw.genome().seqnames[0], "test1");
            assert_eq!(bw.genome().seqnames[1], "test2");

            let mut sum_id   = 0;
            let mut sum_from = 0;
            let mut sum_to   = 0;
            let mut sum_min  = 0.0;
            let mut sum_max  = 0.0;

            for result in bw.query("test1", 0, 100, 10) {
                if let Ok(item) = result {
                    sum_id   += item.data.chrom_id;
                    sum_from += item.data.from;
                    sum_to   += item.data.to;
                    sum_min  += item.data.statistics.min;
                    sum_max  += item.data.statistics.max;
                }
            }

            assert_eq!(sum_id  ,   0);
            assert_eq!(sum_from, 450);
            assert_eq!(sum_to  , 550);

            assert_relative_eq!(sum_min, 49.5, epsilon = 1e-6);
            assert_relative_eq!(sum_max, 49.5, epsilon = 1e-6);
        }
    }

    #[test]
    fn bigwig_test_2() {

        let result =  BigWigFile::new_reader("tests/test_bigwig_2.bw");

        assert!(result.is_ok());

        if let Ok(mut bw) = result {

            // Query raw data
            let r1 : Vec<BbiQueryType> = bw.query("chrY", 1838100, 1838600, 100)
                .filter_map(|r| r.ok())
                .collect();

            // Query zoom data
            let r2 : Vec<BbiQueryType> = bw.query("chrY", 1838100, 1838600, 400)
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
