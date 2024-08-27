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

    use byteorder::LittleEndian;

    use crate::bigwig::{BigWigFile, BigWigOrder};
    use crate::netfile::NetFile;

    #[test]
    fn test_bigwig_1() {

        let result =  BigWigFile::open("tests/test_bigwig_1.bw");

        assert!(result.is_ok());

        if let Ok(mut bw) = result {

            assert_eq!(bw.genome().length(), 2);

            assert_eq!(bw.genome().seqnames[0], "test1");
            assert_eq!(bw.genome().seqnames[1], "test2");

            let mut sum_id   = 0;
            let mut sum_from = 0;
            let mut sum_to   = 0;
            let mut sum_min  = 0.0;
            let mut sum_max  = 0.0;

            for result in bw.bwf.query::<LittleEndian, NetFile>(&mut bw.reader, 0, 0, 100, 1) {
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

            assert_relative_eq!(sum_min, 7.0, epsilon = f32::EPSILON as f64);
            assert_relative_eq!(sum_max, 7.0, epsilon = f32::EPSILON as f64);

        }
    }
}
