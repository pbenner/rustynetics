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

use std::fmt;
use std::io::{self, Read};

/* -------------------------------------------------------------------------- */

pub fn indent_fmt(f: &mut fmt::Formatter<'_>, indent: usize, content: &str) -> fmt::Result {
    writeln!(f, "{:indent$}{}", "", content, indent = indent)
}

/* -------------------------------------------------------------------------- */

pub fn read_until_null<R: Read>(reader: &mut R) -> io::Result<Vec<u8>> {
    let mut buffer = Vec::new();
    let mut byte = [0; 1]; // Buffer to hold a single byte

    loop {
        let bytes_read = reader.read(&mut byte)?;
        if bytes_read == 0 || byte[0] == 0 {
            break;
        }
        buffer.push(byte[0]);
    }

    Ok(buffer)
}

/* -------------------------------------------------------------------------- */

pub fn skip_n_bytes<R: Read>(reader: &mut R, n: usize) -> io::Result<()> {
    const BUFFER_SIZE: usize = 8;  // A small fixed-size buffer
    let mut buffer = [0u8; BUFFER_SIZE];  // No dynamic allocation
    let mut bytes_to_skip = n;

    while bytes_to_skip > 0 {
        // Determine how many bytes to read in this iteration
        let bytes_to_read = BUFFER_SIZE.min(bytes_to_skip);
        // Read and discard the bytes
        let bytes_read = reader.read(&mut buffer[..bytes_to_read])?;

        // If the reader has no more data, break early
        if bytes_read == 0 {
            break;
        }

        // Decrease the remaining number of bytes to skip
        bytes_to_skip -= bytes_read;
    }

    Ok(())
}
