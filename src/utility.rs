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

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
use flate2::write::GzEncoder;
use flate2::Compression;
use regex::Regex;

/* -------------------------------------------------------------------------- */

pub fn remove_duplicates_int(s: Vec<i32>) -> Vec<i32> {
    let mut m: HashSet<i32> = HashSet::new();
    let mut r: Vec<i32> = Vec::new();
    for v in s {
        if m.insert(v) {
            r.push(v);
        }
    }
    r
}

pub fn i_min(a: i32, b: i32) -> i32 {
    std::cmp::min(a, b)
}

pub fn i_max(a: i32, b: i32) -> i32 {
    std::cmp::max(a, b)
}

pub fn i_pow(x: i32, k: i32) -> i32 {
    x.pow(k as u32)
}

pub fn div_int_down(a: i32, b: i32) -> i32 {
    a / b
}

pub fn div_int_up(a: i32, b: i32) -> i32 {
    (a + b - 1) / b
}

pub fn write_file(filename: &str, r: &mut dyn Read, compress: bool) -> io::Result<()> {
    let mut buffer = Vec::new();
    if compress {
        let mut w = GzEncoder::new(&mut buffer, Compression::default());
        io::copy(r, &mut w)?;
        w.finish()?;
    } else {
        io::copy(r, &mut buffer)?;
    }
    std::fs::write(filename, buffer)
}

pub fn is_gzip(filename: &str) -> bool {
    let path = Path::new(filename);
    if let Ok(mut f) = File::open(&path) {
        let mut b = [0; 2];
        if f.read(&mut b).is_ok() {
            return b[0] == 0x1f && b[1] == 0x8b;
        }
    }
    false
}

pub fn fields_quoted(line: &str) -> Vec<String> {
    let mut q = false;
    let mut fields = Vec::new();
    let mut current = String::new();
    for c in line.chars() {
        if c == '"' {
            q = !q;
        } else if c.is_whitespace() && !q {
            if !current.is_empty() {
                fields.push(current.clone());
                current.clear();
            }
        } else {
            current.push(c);
        }
    }
    if !current.is_empty() {
        fields.push(current);
    }
    fields
}

pub fn remove_quotes(s: &str) -> String {
    let re = Regex::new(r#""([^"]*)""#).unwrap();
    re.replace_all(s, "$1").to_string()
}

pub fn reverse_float64(x: Vec<f64>) -> Vec<f64> {
    let mut y = Vec::with_capacity(x.len());
    for i in x.iter().rev() {
        y.push(*i);
    }
    y
}

pub fn bufio_read_line(reader: &mut BufReader<File>) -> io::Result<String> {
    let mut line = String::new();
    match reader.read_line(&mut line) {
        Ok(0) => Err(io::Error::new(io::ErrorKind::UnexpectedEof, "EOF reached")),
        Ok(_) => {
            if line.ends_with('\n') {
                line.pop();
                if line.ends_with('\r') {
                    line.pop();
                }
            }
            Ok(line)
        }
        Err(e) => Err(e),
    }
}

pub fn byte_superset(a: &[u8], b: &[u8]) -> bool {
    let m: HashSet<u8> = a.iter().cloned().collect();
    b.iter().all(|&c| m.contains(&c))
}
