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

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;
use std::{cell::RefCell, rc::Rc};

use flate2::read::GzDecoder;

use crate::track_simple::SimpleTrack;

/* -------------------------------------------------------------------------- */

impl SimpleTrack {

    pub fn read_bedgraph<R: BufRead>(&mut self, reader: R) -> io::Result<()> {

        let mut cur_seq  = Vec::<f64>::new();
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
                ).clone();

                cur_name = name;
            }

            if from >= cur_seq.len() * bin_size || to >= cur_seq.len() * bin_size {
                continue;
            }

            for i in (from / bin_size)..(to / bin_size) {
                cur_seq[i] = value;
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

fn is_gzip<P: AsRef<Path>>(filename: P) -> bool {
    filename.as_ref().extension().map_or(false, |ext| ext == "gz")
}
