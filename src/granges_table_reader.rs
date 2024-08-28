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

 use std::io;
 use std::str::FromStr;
 
 use crate::range::Range;
 use crate::granges::GRanges;
 
 /* -------------------------------------------------------------------------- */
 
pub struct GRangesTableReader {
    col_seqname: i32,
    col_from   : i32,
    col_to     : i32,
    col_strand : i32,
    granges    : GRanges,
}

/* -------------------------------------------------------------------------- */

impl GRangesTableReader {

    pub fn new() -> Self {
        GRangesTableReader{
            col_seqname: -1,
            col_from   : -1,
            col_to     : -1,
            col_strand : -1,
            granges    : GRanges::default(),
        }
    }

    pub fn read_header(&mut self, line: &String) -> io::Result<()> {

        let fields: Vec<&str> = line.trim().split_whitespace().collect();
        for (i, field) in fields.iter().enumerate() {
            match *field {
                "seqnames" => self.col_seqname = i as i32,
                "from"     => self.col_from    = i as i32,
                "to"       => self.col_to      = i as i32,
                "strand"   => self.col_strand  = i as i32,
                "start" if self.col_from == -1 => self.col_from = i as i32,
                "end"   if self.col_to   == -1 => self.col_to   = i as i32,
                _ => (),
            }
        }

        if self.col_seqname == -1 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "is missing a seqnames column"));
        }
        if self.col_from == -1 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "is missing a from column"));
        }
        if self.col_to == -1 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "is missing a to column"));
        }
        Ok(())
    }

    pub fn read_line(&mut self, line: &String, i: i32) -> io::Result<()> {

        let fields: Vec<&str> = line.trim().split_whitespace().collect();
        if fields.len() < self.col_seqname as usize || fields.len() < self.col_from as usize || fields.len() < self.col_to as usize {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid table"));
        }

        let seqname = fields[self.col_seqname as usize].to_string();
        let from    = usize::from_str(fields[self.col_from as usize]).map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("parsing `from' column `{}` failed at line `{}`", self.col_from + 1, i + 1)))?;
        let to      = usize::from_str(fields[self.col_to   as usize]).map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("parsing `to' column `{}` failed at line `{}`"  , self.col_to   + 1, i + 1)))?;

        self.granges.seqnames.push(seqname);
        self.granges.ranges  .push(Range{ from, to });

        if self.col_strand != -1 {
            if fields.len() < self.col_strand as usize {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid table"));
            }
            let strand = fields[self.col_strand as usize].chars().next().unwrap();
            self.granges.strand.push(strand);
        } else {
            self.granges.strand.push('*');
        }
        Ok(())
    }

    pub fn push(&mut self, granges: &mut GRanges) {
        granges.seqnames = self.granges.seqnames.clone();
        granges.ranges   = self.granges.ranges  .clone();
        granges.strand   = self.granges.strand  .clone();
    }

}