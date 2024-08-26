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

    use rustynetics::bigwig::BigWigFile;

    #[test]
    fn test_bigwig_1() {

        match BigWigFile::open("tests/test_bigwig_1.bw") {
        Err(err) => println!("Err: {}", err),
        Ok (mut bw) => {
                println!("Genome: {}", bw.genome());

                for result in bw.query_iterator("test1", 0, 100, 1) {
                    if let Ok(item) = result {
                        println!("result {}", item);
                    }
                }

            }
        }

    }


}
