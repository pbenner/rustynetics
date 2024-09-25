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
use std::io::{self, Read};

use crate::bam::{BamReader, BamReaderOptions};
use crate::granges::GRanges;
use crate::meta::MetaData;

/* -------------------------------------------------------------------------- */

impl GRanges {

    fn read_bam_single_end<R: Read>(&mut self, reader: R, options: BamReaderOptions) -> io::Result<()> {

        let mut bam_reader = BamReader::new(reader, options)?;

        let mut seqnames = Vec::new();
        let mut from     = Vec::new();
        let mut to       = Vec::new();
        let mut strand   = Vec::new();
        let mut sequence = Vec::new();
        let mut mapq     = Vec::new();
        let mut cigar    = Vec::new();
        let mut flag     = Vec::new();
        let mut qual     = Vec::new();

        for item in bam_reader.read_single_end() {

            let block = item?.block;

            if let Some(err) = block.error {
                return Err(err);
            }
            if block.ref_id == -1 || block.flag.unmapped() || block.ref_id < 0 || block.ref_id as usize >= bam_reader.genome.seqnames.len() {
                continue;
            }

            seqnames.push(bam_reader.genome.seqnames[block.ref_id as usize].clone());
            from    .push(block.position as i32);
            to      .push((block.position + block.alignment_length() as u32) as i32);
            strand  .push(if block.flag.reverse_strand() { '-' } else { '+' });
            flag    .push(block.flag.0 as i32);
            mapq    .push(block.mapq as i32);

            if options.read_sequence {
                sequence.push(block.seq);
            }
            if options.read_cigar {
                cigar.push(block.cigar);
            }
            if options.read_qual {
                qual.push(block.qual);
            }
        }

        *self = GRanges::new(seqnames, from, to, strand);
        self.meta.add("flag", flag);
        self.meta.add("mapq", mapq);
        if options.read_sequence {
            self.meta.add("sequence", sequence);
        }
        if options.read_cigar {
            self.meta.add("cigar", cigar);
        }
        if options.read_qual {
            self.meta.add("qual", qual);
        }

        Ok(())
    }

    fn import_bam_single_end(&mut self, filename: &str, options: BamReaderOptions) -> io::Result<()> {
        let file = File::open(filename)?;
        self.read_bam_single_end(file, options)
    }
}

