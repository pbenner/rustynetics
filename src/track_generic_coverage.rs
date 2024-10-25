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

use std::error::Error;

use futures::executor::block_on_stream;

use crate::coverage::CoverageConfig;
use crate::genome::Genome;
use crate::read_stream::ReadStream;
use crate::track_generic::GenericMutableTrack;
use crate::bam::BamFile;

/* -------------------------------------------------------------------------- */

impl<'a> GenericMutableTrack<'a> {

    pub fn coverage_from_bam(
        mut config         : CoverageConfig,
        track1             : GenericMutableTrack,
        track2_arg         : Option<GenericMutableTrack>,
        filenames_treatment: &Vec<&str>,
        filenames_control  : &Vec<&str>,
        fraglen_treatment  : &Vec<usize>,
        fraglen_control    : &Vec<usize>,
        genome             : Genome,
    ) -> Result<(), Box<dyn Error>> {
        // Treatment data
        let mut n_treatment = 0;
        let mut n_control   = 0;

        for (i, filename) in filenames_treatment.iter().enumerate() {
            let mut err_opt = None;
            let fraglen = fraglen_treatment[i];

            log!(config.logger, "Reading treatment tags from `{}`", filename);
            let mut bam = BamFile::open(filename, None)?;

            let treatment = Box::pin(bam.reader.read_simple_stream(!config.paired_as_single_end, config.paired_end_strand_specific));
            // First round of filtering
            let treatment = ReadStream::filter_paired_end   (treatment, Some(&config.logger),  config.filter_paired_end);
            let treatment = ReadStream::filter_single_end   (treatment, Some(&config.logger),  config.filter_single_end);
            let treatment = ReadStream::paired_as_single_end(treatment, Some(&config.logger),  config.paired_as_single_end);
            let treatment = ReadStream::filter_read_length  (treatment, Some(&config.logger), &config.filter_read_lengths);
            let treatment = ReadStream::filter_duplicates   (treatment, Some(&config.logger),  config.filter_duplicates);
            let treatment = ReadStream::filter_mapq         (treatment, Some(&config.logger),  config.filter_mapq);
            // Second round of filtering
            let treatment = ReadStream::filter_strand       (treatment, Some(&config.logger),  config.filter_strand);
            let treatment = ReadStream::shift_reads         (treatment, Some(&config.logger), &config.shift_reads);

            let treatment_iter = block_on_stream(treatment).map_while(|item| {
                match item {
                    Ok(read) => Some(read),
                    Err(err) => {
                        err_opt = Some(err);
                        None
                    }
                }
            });
        
            n_treatment += track1.add_reads(treatment_iter, fraglen, &config.binning_method);

            if let Some(err) = err_opt {
                return Err(Box::new(err));
            }
        }

        // Normalization for treatment
        if config.normalize_track == "rpkm" {
            log!(config.logger, "Normalizing treatment track (rpkm)");
            let c = 1_000_000.0 / (n_treatment as f64 * config.bin_size as f64);
            track1.map(|_name, _i, x| c * x)?;
            config.pseudocounts[0] *= c;
        }

        if config.normalize_track == "cpm" {
            log!(config.logger, "Normalizing treatment track (cpm)");
            let c = 1_000_000.0 / n_treatment as f64;
            track1.map(|_name, _i, x| c * x)?;
            config.pseudocounts[0] *= c;
        }

        if let Some(track2) = track2_arg {
            // Control data

            for (i, filename) in filenames_control.iter().enumerate() {
                let mut err_opt = None;
                let fraglen = fraglen_control[i];

                log!(config.logger, "Reading control tags from `{}`", filename);
                let mut bam = BamFile::open(filename, None)?;
                let control = Box::pin(bam.reader.read_simple_stream(!config.paired_as_single_end, config.paired_end_strand_specific));

                // First round of filtering
                let control = ReadStream::filter_paired_end   (control, Some(&config.logger),  config.filter_paired_end);
                let control = ReadStream::filter_single_end   (control, Some(&config.logger),  config.filter_single_end);
                let control = ReadStream::paired_as_single_end(control, Some(&config.logger),  config.paired_as_single_end);
                let control = ReadStream::filter_read_length  (control, Some(&config.logger), &config.filter_read_lengths);
                let control = ReadStream::filter_duplicates   (control, Some(&config.logger),  config.filter_duplicates);
                let control = ReadStream::filter_mapq         (control, Some(&config.logger),  config.filter_mapq);
                // Second round of filtering
                let control = ReadStream::filter_strand       (control, Some(&config.logger),  config.filter_strand);
                let control = ReadStream::shift_reads         (control, Some(&config.logger), &config.shift_reads);

                let control_iter = block_on_stream(control).map_while(|item| {
                    match item {
                        Ok(read) => Some(read),
                        Err(err) => {
                            err_opt = Some(err);
                            None
                        }
                    }
                });
        
                n_control += track2.add_reads(control_iter, fraglen, &config.binning_method);

                if let Some(err) = err_opt {
                    return Err(Box::new(err));
                }
            }

            // Normalization for control
            if config.normalize_track == "rpkm" {
                log!(config.logger, "Normalizing control track (rpkm)");
                let c = 1_000_000.0 / (n_control as f64 * config.bin_size as f64);
                track2.map(|_name, _i, x| c * x)?;
                config.pseudocounts[1] *= c;
            }

            if config.normalize_track == "cpm" {
                log!(config.logger, "Normalizing control track (cpm)");
                let c = 1_000_000.0 / n_control as f64;
                track2.map(|_name, _i, x| c * x)?;
                config.pseudocounts[1] *= c;
            }

            if config.smoothen_control {
                track2.smoothen(config.smoothen_min, config.smoothen_sizes.clone())?;
            }

            log!(config.logger, "Combining treatment and control tracks...");
            track1.normalize(
                track2.track.as_track(), config.pseudocounts[0], config.pseudocounts[1], config.log_scale)?;
        } else {
            // No control data
            if config.pseudocounts[0] != 0.0 {
                log!(config.logger, "Adding pseudocount `{}`", config.pseudocounts[0]);
                track1.map(|_name, _i, x| x + config.pseudocounts[0])?;
            }
            if config.log_scale {
                log!(config.logger, "Log-transforming data");
                track1.map(|_name, _i, x| x.ln())?;
            }
        }

        // Filtering chromosomes
        if config.remove_filtered_chroms {
            if !config.filter_chroms.is_empty() {
                log!(config.logger, "Removing chromosomes `{}`", config.filter_chroms.join(", "));
                track1.filter_genome(|name, _length| {
                    !config.filter_chroms.contains(&name.to_string())
                });
            }
        } else {
            if !config.filter_chroms.is_empty() {
                log!(config.logger, "Removing all reads from `{}`", config.filter_chroms.join(", "));
                for chr in &config.filter_chroms {
                    if let Ok(mut s) = track1.track.get_sequence_mut(chr) {
                        for i in 0..s.n_bins() {
                            s.set_bin(i, 0.0);
                        }
                    }
                }
            }
        }

        Ok(())

    }

}

