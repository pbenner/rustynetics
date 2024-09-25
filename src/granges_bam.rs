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

use std::io::Read;
use std::error::Error;

use crate::bam::{BamReader, BamReaderOptions};
use crate::granges::GRanges;
use crate::meta::MetaData;
use crate::netfile::NetFile;

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn read_bam_single_end<R: Read>(&mut self, reader: R, options_arg: Option<BamReaderOptions>) -> Result<(), Box<dyn Error>> {

        let options = options_arg.unwrap_or_default();

        let mut bam_reader = BamReader::new(reader, Some(options))?;

        let mut seqnames = Vec::new();
        let mut from     = Vec::new();
        let mut to       = Vec::new();
        let mut strand   = Vec::new();
        let mut sequence = Vec::new();
        let mut mapq     = Vec::new();
        let mut cigar    = Vec::new();
        let mut flag     = Vec::new();
        let mut qual     = Vec::new();

        let genome = bam_reader.get_genome().clone();

        for item in bam_reader.read_single_end() {

            let block = item?.block;

            if block.ref_id == -1 || block.flag.unmapped() || block.ref_id < 0 || block.ref_id as usize >= genome.seqnames.len() {
                continue;
            }

            let pos = block.position as usize;
            let len = block.cigar.alignment_length();

            seqnames.push(genome.seqnames[block.ref_id as usize].clone());
            from    .push(pos);
            to      .push(pos+len);
            strand  .push(if block.flag.reverse_strand() { '-' } else { '+' });

            flag    .push(block.flag.0 as i64);
            mapq    .push(block.mapq   as i64);

            if options.read_sequence {
                sequence.push(block.seq.to_string());
            }
            if options.read_cigar {
                cigar.push(block.cigar.to_string());
            }
            if options.read_qual {
                qual.push(block.qual.to_string());
            }
        }

        *self = GRanges::new(seqnames, from, to, strand);
        self.meta.add("flag", MetaData::IntArray(flag))?;
        self.meta.add("mapq", MetaData::IntArray(mapq))?;
        if options.read_sequence {
            self.meta.add("sequence", MetaData::StringArray(sequence))?;
        }
        if options.read_cigar {
            self.meta.add("cigar", MetaData::StringArray(cigar))?;
        }
        if options.read_qual {
            self.meta.add("qual", MetaData::StringArray(qual))?;
        }

        Ok(())
    }

    pub fn import_bam_single_end(&mut self, filename: &str, options: Option<BamReaderOptions>) -> Result<(), Box<dyn Error>> {
        let file = NetFile::open(filename)?;
        self.read_bam_single_end(file, options)
    }
}

