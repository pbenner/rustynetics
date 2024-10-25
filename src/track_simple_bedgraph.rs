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

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::str::FromStr;
use std::{cell::RefCell, rc::Rc};

use flate2::read::GzDecoder;

use crate::track_simple::SimpleTrack;
use crate::utility::is_gzip;

/* -------------------------------------------------------------------------- */

impl SimpleTrack {

    pub fn read_bedgraph<R: BufRead>(&mut self, reader: R) -> io::Result<()> {

        let mut cur_seq  = &mut Rc::new(RefCell::new(vec![]));
        let mut cur_name = String::new();
        let bin_size     = self.bin_size;

        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split_whitespace().collect();

            if fields.is_empty() {
                break;
            }
            if fields.len() != 4 {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, "BedGraph file must have four columns!"));
            }

            let from  = usize::from_str(fields[1]).map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid integer in column 2"))?;
            let to    = usize::from_str(fields[2]).map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid integer in column 3"))?;
            let value =   f64::from_str(fields[3]).map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid float in column 4"))?;
            let name  = fields[0].to_string();

            if name != cur_name {

                cur_seq = self.data.entry(name.clone()).or_insert_with(
                    || Rc::new(RefCell::new(vec![0.0; (to / bin_size) + 1]))
                );

                cur_name = name;
            }

            let mut s = cur_seq.borrow_mut();

            if from >= s.len() * bin_size || to >= s.len() * bin_size {
                continue;
            }

            for i in (from / bin_size)..(to / bin_size) {
                s[i] = value;
            }
        }

        Ok(())
    }

    pub fn import_bedgraph(&mut self, filename: &str) -> io::Result<()> {
        let file = File::open(&filename)?;

        let reader: Box<dyn BufRead> = if is_gzip(&filename) {
            let decoder = GzDecoder::new(file);
            Box::new(BufReader::new(decoder))
        } else {
            Box::new(BufReader::new(file))
        };

        self.read_bedgraph(reader)
    }
}
