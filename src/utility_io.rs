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
