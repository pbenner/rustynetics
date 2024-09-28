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

const BUFSIZE : usize = 10000;

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

            if seqnames.capacity() == seqnames.len() {
                seqnames.reserve(BUFSIZE);
                from    .reserve(BUFSIZE);
                to      .reserve(BUFSIZE);
                strand  .reserve(BUFSIZE);
                sequence.reserve(BUFSIZE);
                mapq    .reserve(BUFSIZE);
                cigar   .reserve(BUFSIZE);
                flag    .reserve(BUFSIZE);
                qual    .reserve(BUFSIZE);
            }

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
        seqnames.shrink_to_fit();
        from    .shrink_to_fit();
        to      .shrink_to_fit();
        strand  .shrink_to_fit();
        sequence.shrink_to_fit();
        mapq    .shrink_to_fit();
        cigar   .shrink_to_fit();
        flag    .shrink_to_fit();
        qual    .shrink_to_fit();

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

    pub fn read_bam_paired_end<R: Read>(&mut self, reader: R, options_arg: Option<BamReaderOptions>) -> Result<(), Box<dyn Error>> {

        let mut options = options_arg.unwrap_or_default();

        // Ensure that CIGAR strings are always read
        options.read_cigar = true;

        let mut bam_reader = BamReader::new(reader, Some(options))?;

        let mut seqnames  = Vec::new();
        let mut from      = Vec::new();
        let mut to        = Vec::new();
        let mut strand    = Vec::new();
        let mut sequence1 = Vec::new();
        let mut sequence2 = Vec::new();
        let mut mapq1     = Vec::new();
        let mut mapq2     = Vec::new();
        let mut cigar1    = Vec::new();
        let mut cigar2    = Vec::new();
        let mut flag1     = Vec::new();
        let mut flag2     = Vec::new();
        let mut qual1     = Vec::new();
        let mut qual2     = Vec::new();

        let genome = bam_reader.get_genome().clone();

        for item in bam_reader.read_paired_end() {

            if seqnames.capacity() == seqnames.len() {
                seqnames .reserve(BUFSIZE);
                from     .reserve(BUFSIZE);
                to       .reserve(BUFSIZE);
                strand   .reserve(BUFSIZE);
                sequence1.reserve(BUFSIZE);
                sequence2.reserve(BUFSIZE);
                mapq1    .reserve(BUFSIZE);
                mapq2    .reserve(BUFSIZE);
                cigar1   .reserve(BUFSIZE);
                cigar2   .reserve(BUFSIZE);
                flag1    .reserve(BUFSIZE);
                flag2    .reserve(BUFSIZE);
                qual1    .reserve(BUFSIZE);
                qual2    .reserve(BUFSIZE);
            }

            let r = item?;

            let block1 = &r.block1;
            let block2 = &r.block2;

            // Skip unmapped reads or reads not mapped properly as pairs
            if block1.flag.unmapped() || !block1.flag.read_mapped_proper_paired() {
                continue;
            }
            if block2.flag.unmapped() || !block2.flag.read_mapped_proper_paired() {
                continue;
            }

            // Check for valid reference ID
            if block1.ref_id == -1 || block1.ref_id != block2.ref_id {
                continue;
            }

            if block1.ref_id < 0 || block1.ref_id as usize >= genome.seqnames.len() {
                return Err(format!("bam file contains invalid RefID '{}'", block1.ref_id).into());
            }

            let pos1 = block1.position as usize;
            let pos2 = block2.position as usize;
            let len2 = block2.cigar.alignment_length();

            seqnames.push(genome.seqnames[block1.ref_id as usize].clone());
            from    .push(pos1);
            to      .push(pos2 + len2);
            strand  .push('*');  // Paired-end, strand is ambiguous, so it's set to '*'
            // Capture metadata
            flag1   .push(block1.flag.0 as i64);
            flag2   .push(block2.flag.0 as i64);
            mapq1   .push(block1.mapq   as i64);
            mapq2   .push(block2.mapq   as i64);

            if options.read_sequence {
                sequence1.push(block1.seq.to_string());
                sequence2.push(block2.seq.to_string());
            }
            if options.read_cigar {
                cigar1.push(block1.cigar.to_string());
                cigar2.push(block2.cigar.to_string());
            }
            if options.read_qual {
                qual1.push(block1.qual.to_string());
                qual2.push(block2.qual.to_string());
            }
        }

        seqnames .shrink_to_fit();
        from     .shrink_to_fit();
        to       .shrink_to_fit();
        strand   .shrink_to_fit();
        sequence1.shrink_to_fit();
        sequence2.shrink_to_fit();
        mapq1    .shrink_to_fit();
        mapq2    .shrink_to_fit();
        cigar1   .shrink_to_fit();
        cigar2   .shrink_to_fit();
        flag1    .shrink_to_fit();
        flag2    .shrink_to_fit();
        qual1    .shrink_to_fit();
        qual2    .shrink_to_fit();

        *self = GRanges::new(seqnames, from, to, strand);
        self.meta.add("flag1", MetaData::IntArray(flag1))?;
        self.meta.add("flag2", MetaData::IntArray(flag2))?;
        self.meta.add("mapq1", MetaData::IntArray(mapq1))?;
        self.meta.add("mapq2", MetaData::IntArray(mapq2))?;

        if options.read_sequence {
            self.meta.add("sequence1", MetaData::StringArray(sequence1))?;
            self.meta.add("sequence2", MetaData::StringArray(sequence2))?;
        }
        if options.read_cigar {
            self.meta.add("cigar1", MetaData::StringArray(cigar1))?;
            self.meta.add("cigar2", MetaData::StringArray(cigar2))?;
        }
        if options.read_qual {
            self.meta.add("qual1", MetaData::StringArray(qual1))?;
            self.meta.add("qual2", MetaData::StringArray(qual2))?;
        }

        Ok(())
    }

}

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn import_bam_single_end(&mut self, filename: &str, options: Option<BamReaderOptions>) -> Result<(), Box<dyn Error>> {
        let file = NetFile::open(filename)?;
        self.read_bam_single_end(file, options)
    }

    pub fn import_bam_paired_end(&mut self, filename: &str, options: Option<BamReaderOptions>) -> Result<(), Box<dyn Error>> {
        let file = NetFile::open(filename)?;
        self.read_bam_paired_end(file, options)
    }
}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::bam::BamReaderOptions;
    use crate::granges::GRanges;

    #[test]
    fn test_granges_bam_read_single_end() {

        let n = 4891;

        let mut granges = GRanges::default();
        let mut options = BamReaderOptions::default();

        options.read_cigar = true;
        options.read_qual  = true;

        assert!(
            granges.import_bam_single_end("src/bam_test.2.bam", Some(options)).is_ok()
        );
        assert_eq!(
            granges.num_rows(), n
        );

        let flag  = granges.meta.get_column_int("flag" ).unwrap();
        let mapq  = granges.meta.get_column_int("mapq" ).unwrap();
        let cigar = granges.meta.get_column_str("cigar").unwrap();
        let qual  = granges.meta.get_column_str("qual" ).unwrap();

        // Check values of last row, which should be ok if everything before was
        // read without error
        assert_eq!(
            granges.seqnames[n-1], "chr11"
        );
        assert_eq!(
            granges.ranges[n-1].from, 4737786
        );
        assert_eq!(
            granges.ranges[n-1].to, 4737837
        );
        assert_eq!(
            flag[n-1], 163
        );
        assert_eq!(
            mapq[n-1], 60
        );
        assert_eq!(
            cigar[n-1], "51M"
        );
        assert_eq!(
            qual[n-1], "@@<DDBDDFD+C?A:1CFDHBFHC<?F9+CGGI:49CCGFACE99?DC990"
        );

    }

    #[test]
    fn test_granges_bam_read_paired_end() {

        let n = 2335;

        let mut granges = GRanges::default();
        let mut options = BamReaderOptions::default();

        options.read_cigar    = true;
        options.read_qual     = true;
        options.read_sequence = true;

        assert!(
            granges.import_bam_paired_end("src/bam_test.2.bam", Some(options)).is_ok()
        );
        assert_eq!(
            granges.num_rows(), n
        );

        let flag1  = granges.meta.get_column_int("flag1" ).unwrap();
        let flag2  = granges.meta.get_column_int("flag2" ).unwrap();
        let mapq1  = granges.meta.get_column_int("mapq1" ).unwrap();
        let mapq2  = granges.meta.get_column_int("mapq2" ).unwrap();
        let cigar1 = granges.meta.get_column_str("cigar1").unwrap();
        let cigar2 = granges.meta.get_column_str("cigar2").unwrap();
        let seq1   = granges.meta.get_column_str("sequence1").unwrap();
        let seq2   = granges.meta.get_column_str("sequence2").unwrap();
        let qual1  = granges.meta.get_column_str("qual1" ).unwrap();
        let qual2  = granges.meta.get_column_str("qual2" ).unwrap();

        // Check values of last row, which should be ok if everything before was
        // read without error

        assert_eq!(
            granges.seqnames[n-1], "chr11"
        );
        assert_eq!(
            granges.ranges[n-1].from, 4737786
        );
        assert_eq!(
            granges.ranges[n-1].to, 4737889
        );
        assert_eq!(
            flag1[n-1], 163
        );
        assert_eq!(
            flag2[n-1], 83
        );
        assert_eq!(
            mapq1[n-1], 60
        );
        assert_eq!(
            mapq2[n-1], 60
        );
        assert_eq!(
            cigar1[n-1], "51M"
        );
        assert_eq!(
            cigar2[n-1], "51M"
        );
        assert_eq!(
            seq1[n-1], "AGGGACAGATCCGTACTATGTTGCTCAGGATGGTCCCAAATTCCTGATTCA"
        );
        assert_eq!(
            seq2[n-1], "ACAATCCTCTTGCCTCAGCCTCTTAACTAAGTAGCTGTCACACTGCATGGT"
        );
        assert_eq!(
            qual1[n-1], "@@<DDBDDFD+C?A:1CFDHBFHC<?F9+CGGI:49CCGFACE99?DC990"
        );
        assert_eq!(
            qual2[n-1], "DB9;HCD?D??:?:):)CCA<C2:@HFAHEEHF@<?<?:ACADB;:BB1@?"
        );

    }

}
