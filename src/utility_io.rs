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

/// Formats and writes a string with a specified indentation level.
///
/// This function indents the provided content by a specified number of spaces
/// and writes it to the provided formatter. The indentation is applied by
/// using the `write!` macro to create a new line with leading spaces.
///
/// # Parameters
///
/// - `f`: A mutable reference to a `fmt::Formatter` where the formatted output
///   will be written.
/// - `indent`: The number of spaces to indent the content.
/// - `content`: The string content to be indented and written.
///
/// # Returns
///
/// A `fmt::Result` indicating the success or failure of the formatting operation.
pub fn indent_fmt(f: &mut fmt::Formatter<'_>, indent: usize, content: &str) -> fmt::Result {
    writeln!(f, "{:indent$}{}", "", content, indent = indent)
}

/* -------------------------------------------------------------------------- */

/// Reads bytes from the given reader until a null byte (0) is encountered.
///
/// This function reads data from a type that implements the `Read` trait until
/// it encounters a null byte or the end of the stream. It collects the read
/// bytes into a vector and returns it.
///
/// # Type Parameters
///
/// - `R`: A type that implements the `Read` trait.
///
/// # Parameters
///
/// - `reader`: A mutable reference to the reader from which bytes are to be read.
///
/// # Returns
///
/// A `Result` containing a `Vec<u8>` with the read bytes, or an `io::Error` if
/// an error occurs during reading.
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

/// Skips a specified number of bytes in the given reader.
///
/// This function reads and discards bytes from a type that implements the `Read`
/// trait without storing them. It uses a fixed-size buffer to read bytes in chunks,
/// minimizing dynamic allocation. The operation will stop either when the specified
/// number of bytes have been skipped or when the end of the stream is reached.
///
/// # Type Parameters
///
/// - `R`: A type that implements the `Read` trait.
///
/// # Parameters
///
/// - `reader`: A mutable reference to the reader from which bytes are to be skipped.
/// - `n`: The number of bytes to skip.
///
/// # Returns
///
/// A `Result` indicating the success or failure of the skip operation.
/// The function returns `Ok(())` if successful, or an `io::Error` if an error occurs.
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
