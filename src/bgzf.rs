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

use std::io::{self, Read};
use flate2::read::MultiGzDecoder;
use byteorder::{LittleEndian, ReadBytesExt};

/* -------------------------------------------------------------------------- */

// Structure to hold the BGZF extra fields
#[derive(Clone, Debug, PartialEq)]
pub struct BgzfExtra {
    pub si1  : u8,
    pub si2  : u8,
    pub slen : u16,
    pub bsize: u16,
}

/* -------------------------------------------------------------------------- */

// BGZF reader structure that wraps GzDecoder
#[derive(Debug)]
pub struct BgzfReader<R: Read> {
    decoder: MultiGzDecoder<R>,
}

/* -------------------------------------------------------------------------- */

impl<R: Read> BgzfReader<R> {
    // Function to create a new BgzfReader
    pub fn new(reader: R) -> io::Result<BgzfReader<R>> {
        let decoder = MultiGzDecoder::new(reader);
        Ok(BgzfReader { decoder })
    }

    // Function to extract the BGZF extra fields
    pub fn get_extra(&mut self) -> io::Result<BgzfExtra> {
        let header = self.decoder.header().ok_or_else(|| {
            io::Error::new(io::ErrorKind::Other, "Failed to read gzip header")
        })?;
        
        if let Some(extra) = header.extra() {
            if extra.len() != 6 {
                return Err(io::Error::new(io::ErrorKind::Other, "No extra information available"));
            }

            let mut cursor = io::Cursor::new(extra);
            let si1   = cursor.read_u8()?;
            let si2   = cursor.read_u8()?;
            let slen  = cursor.read_u16::<LittleEndian>()?;
            let bsize = cursor.read_u16::<LittleEndian>()?;

            Ok(BgzfExtra { si1, si2, slen, bsize })
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "No extra information available"))
        }
    }
}

/* -------------------------------------------------------------------------- */

impl<R: Read> Read for BgzfReader<R> {

    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.decoder.read(buf)
    }

}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::fs::File;
    use byteorder::ReadBytesExt;
    use byteorder::LittleEndian;

    use crate::bgzf::{BgzfExtra, BgzfReader};
    use crate::netfile::NetFile;

    #[test]
    fn test_bgzf() {

        let file = File::open("src/bam_test.1.bam");
        let bgzf = BgzfReader::new(file.unwrap());

        let mut reader = bgzf.unwrap();

        assert_eq!(
            reader.read_i64::<LittleEndian>().unwrap(),
            150345695554
        );

        assert_eq!(
            reader.get_extra().unwrap(),
            BgzfExtra{si1: 66, si2: 67, slen: 2, bsize: 77}
        );

    }

    #[test]
    fn test_bgzf_netfile() {

        let file = NetFile::open("src/bam_test.1.bam");
        let bgzf = BgzfReader::new(file.unwrap());

        let mut reader = bgzf.unwrap();

        assert_eq!(
            reader.read_i64::<LittleEndian>().unwrap(),
            150345695554
        );

        assert_eq!(
            reader.get_extra().unwrap(),
            BgzfExtra{si1: 66, si2: 67, slen: 2, bsize: 77}
        );

    }

}
