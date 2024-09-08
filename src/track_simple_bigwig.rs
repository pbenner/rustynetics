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
use std::error::Error;
use std::io::{self, Read, Seek, Write};

use crate::bigwig::BigWigReader;
use crate::track_generic::GenericTrack;
use crate::track_simple::SimpleTrack;
use crate::track_statistics::BinSummaryStatistics;

/* -------------------------------------------------------------------------- */

impl SimpleTrack {

    fn read_bigwig<R: Read + Seek>(
        &mut self,
        reader: &mut R,
        name  : &str,
        f     : BinSummaryStatistics,
        mut bin_size: usize,
        bin_overlap: usize,
        init: f64,
    ) -> Result<(), Box<dyn Error>> {
        let mut bwr = BigWigReader::new(reader)?;

        let mut sequences: Vec<Vec<f64>> = Vec::new();

        for seqname in &bwr.genome.seqnames {
            let (s, b) = bwr.query_sequence(seqname, &f, bin_size, bin_overlap, init)?;
            if bin_size == 0 {
                bin_size = b;
            }
            sequences.push(s);
        }

        let tmp = SimpleTrack::new(name.to_string(), sequences, bwr.genome, bin_size)?;
        *self = tmp;

        Ok(())
    }

    fn import_bigwig(
        &mut self,
        filename   : &str,
        name       : &str,
        s          : BinSummaryStatistics,
        bin_size   : usize,
        bin_overlap: usize,
        init       : f64,
    ) -> Result<(), Box<dyn Error>> {

        let mut f = File::open(filename)?;
        self.read_big_wig(&mut f, name, s, bin_size, bin_overlap, init)
            .map_err(|e| {
                io::Error::new(io::ErrorKind::Other, format!("importing bigWig file from `{}` failed: {}", filename.display(), e))
            }
        )

    }

    fn write_bigwig<W: Write + Seek>(
        &self,
        writer: &mut W,
        args: &[impl ToString]
    ) -> Result<(), Box<dyn Error>> {

        GenericTrack { track: self }.write_big_wig(writer, args)

    }

    fn export_bigwig(
        &self,
        filename: &str,
        args: &[impl ToString]
    ) -> Result<(), Box<dyn Error>> {

        let mut f = File::create(filename)?;
        self.write_bigwig(&mut f, args).map_err(|e| {
            io::Error::new(io::ErrorKind::Other, format!("exporting bigWig file to `{}` failed: {}", filename.display(), e))
        })

    }

}