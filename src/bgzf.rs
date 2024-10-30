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

use std::io::{self, Read};
use flate2::read::MultiGzDecoder;
use byteorder::{LittleEndian, ReadBytesExt};

/* -------------------------------------------------------------------------- */

/// Structure to hold the BGZF extra fields.
///
/// This struct represents additional metadata fields that are
/// specific to the BGZF format. These fields include:
/// - `si1`: Subfield identifier 1.
/// - `si2`: Subfield identifier 2.
/// - `slen`: Length of the subfield data.
/// - `bsize`: Block size of the BGZF compressed block.
#[derive(Clone, Debug, PartialEq)]
pub struct BgzfExtra {
    pub si1  : u8,
    pub si2  : u8,
    pub slen : u16,
    pub bsize: u16,
}

/* -------------------------------------------------------------------------- */

/// BGZF reader structure that wraps around a GzDecoder.
///
/// `BgzfReader` is a wrapper around `MultiGzDecoder` that provides
/// additional functionality for reading BGZF compressed files, specifically
/// the ability to extract BGZF-specific extra fields and metadata.
#[derive(Debug)]
pub struct BgzfReader<R: Read> {
    decoder: MultiGzDecoder<R>,
}

/* -------------------------------------------------------------------------- */

impl<R: Read> BgzfReader<R> {
    /// Creates a new `BgzfReader` from a given reader.
    ///
    /// # Parameters
    ///
    /// - `reader`: The reader from which BGZF-compressed data is read.
    ///
    /// # Returns
    ///
    /// A `Result` containing a new `BgzfReader` instance if successful,
    /// or an `io::Error` if initialization fails.
    pub fn new(reader: R) -> io::Result<BgzfReader<R>> {
        let decoder = MultiGzDecoder::new(reader);
        Ok(BgzfReader { decoder })
    }

    /// Extracts the BGZF-specific extra fields from the compressed data.
    ///
    /// This function reads the gzip header of the BGZF-compressed data and
    /// retrieves the extra fields if they are present. These fields contain
    /// metadata specific to the BGZF format.
    ///
    /// # Returns
    ///
    /// A `Result` containing a `BgzfExtra` struct with BGZF metadata,
    /// or an `io::Error` if the header or extra fields cannot be read.
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

    /// Reads data from the BGZF compressed stream into the provided buffer.
    ///
    /// This implementation of the `Read` trait allows reading decompressed
    /// data from the BGZF stream directly, making it compatible with standard
    /// Rust I/O operations.
    ///
    /// # Parameters
    ///
    /// - `buf`: A mutable buffer where the decompressed data will be stored.
    ///
    /// # Returns
    ///
    /// The number of bytes read into `buf`, or an `io::Error` if reading fails.
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

        let file = File::open("tests/test_bam_1.bam");
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

        let file = NetFile::open("tests/test_bam_1.bam");
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
