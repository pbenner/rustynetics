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

use crate::bam::{BamFile, bam_import_genome};
use crate::genome::Genome;
use crate::error::ArgumentError;
use crate::log;

use crate::coverage::{CoverageError, CoverageConfig, OptionCoverage, FraglenEstimate};
use crate::read_stream::ReadStream;
use crate::track_simple::SimpleTrack;
use crate::track_statistics::estimate_fragment_length;

/* -------------------------------------------------------------------------- */

/// Estimates the mean fragment length from BAM data.
///
/// This function processes a BAM file to read, filter, and analyze genomic reads according to
/// parameters defined in `CoverageConfig`. It returns an estimate of the mean fragment length,
/// using various filters specified in the configuration, such as mapping quality, read length,
/// and duplicate read removal.
///
/// # Parameters
///
/// - `config`: Reference to `CoverageConfig`, defining filtering, logging, and configuration
///             for fragment length estimation. Key fields used include:
///   - `logger`: A `Logger` for logging progress and details during processing.
///   - `fraglen_range`: A `(i32, i32)` tuple specifying the range of fragment lengths to consider.
///   - `fraglen_bin_size`: Size of bins for fragment length estimation.
///   - `filter_mapq`: Minimum mapping quality threshold for filtering reads.
///   - `filter_read_lengths`: `[usize; 2]` array specifying minimum and maximum read lengths.
///   - `filter_duplicates`: Boolean flag to exclude duplicate reads.
///
/// - `filename`: Path to the BAM file containing reads for fragment length estimation.
/// - `genome`: A `Genome` reference.
///
/// # Returns
///
/// - `Ok(FraglenEstimate)`: Contains estimated mean fragment length and associated data if
///   successful.
/// - `Err`: Returns an error if there is an issue reading the BAM file or estimating the
///          fragment length.
///
/// # Errors
///
/// Returns an error in case of:
/// - BAM file reading issues.
/// - Any filtering or read stream errors encountered during processing.
pub fn estimate_fraglen(config: &CoverageConfig, filename: &str, genome: &Genome) -> Result<FraglenEstimate, Box<dyn Error>> {

    log!(config.logger, "Reading tags from `{}`", filename);

    let bam_result = BamFile::open(filename, None);
    let mut bam = match bam_result {
        Ok (b)   => b,
        Err(err) => return Err(err),
    };

    // Read the reads
    let reads = Box::pin(bam.reader.read_simple_stream(false, false));

    // First round of filtering
    let reads = ReadStream::filter_single_end (reads, Some(&config.logger), true);
    let reads = ReadStream::filter_read_length(reads, Some(&config.logger), &config.filter_read_lengths);
    let reads = ReadStream::filter_duplicates (reads, Some(&config.logger),  config.filter_duplicates);
    let reads = ReadStream::filter_mapq       (reads, Some(&config.logger),  config.filter_mapq);

    let mut err_opt = None;
    // Convert stream to iterator and catch errors
    let reads_iter = block_on_stream(reads).map_while(|item| {
        match item {
            Ok(read) => Some(read),
            Err(err) => {
                err_opt = Some(err);
                None
            }
        }
    });

    // Estimate fragment length
    log!(config.logger, "Estimating mean fragment length");

    let r = match estimate_fragment_length(reads_iter, genome, 2000, config.fraglen_bin_size, config.fraglen_range) {
        Ok((fraglen, x, y, _n)) => {
            log!(config.logger, "Estimated mean fragment length: {}", fraglen);
            Ok(FraglenEstimate { fraglen: fraglen as usize, x, y })
        },
        Err(err) => Err(err),
    };

    // Check if reading the BAM file yielded an error
    if let Some(err) = err_opt {
        return Err(Box::new(err));
    }
    r
}

/* -------------------------------------------------------------------------- */

/// Calculates coverage tracks from treatment and control BAM files, with optional fragment length estimation.
///
/// This function reads and processes BAM files for both treatment and control groups to calculate
/// coverage tracks, optionally estimating fragment lengths if none are provided. It supports
/// configurable filtering and normalization through `CoverageConfig` options, including settings
/// for bin size, read shifting, and duplicate filtering.
///
/// # Parameters
///
/// - `filenames_treatment`: Vector of file paths to treatment BAM files.
/// - `filenames_control`: Vector of file paths to control BAM files.
/// - `fraglen_treatment`: Vector of optional fragment length estimates for each treatment BAM file.
/// - `fraglen_control`: Vector of optional fragment length estimates for each control BAM file.
/// - `options`: A vector of `OptionCoverage` enums specifying configuration options such as:
///   - `BinningMethod`, `BinSize`, `BinOverlap`: Define binning parameters for coverage.
///   - `NormalizeTrack`: Specifies normalization strategy for the output track.
///   - `ShiftReads`: Allows for shifting reads, e.g., for strand-specific analysis.
///   - `EstimateFraglen`: Boolean flag to control whether fragment lengths should be estimated.
///
/// # Returns
///
/// - `Ok((SimpleTrack, Vec<FraglenEstimate>, Vec<FraglenEstimate>))`: Tuple containing the coverage
///   track, treatment fragment length estimates, and control fragment length estimates.
/// - `Err(CoverageError)`: Indicates an error due to input inconsistency, genome mismatch, or
///   calculation issues.
///
/// # Errors
///
/// - Returns `CoverageError` if:
///   - Provided fragment length vectors do not match the number of input files.
///   - Treatment and control BAM files have differing genome references.
///   - An error occurs during coverage calculation or file processing.
pub fn bam_coverage(
    filenames_treatment: &Vec<&str>,
    filenames_control  : &Vec<&str>,
    fraglen_treatment  : &Vec<Option<usize>>,
    fraglen_control    : &Vec<Option<usize>>,
    options            : Vec<OptionCoverage>,
) -> Result<(SimpleTrack, Vec<FraglenEstimate>, Vec<FraglenEstimate>), CoverageError> {

    let mut config = CoverageConfig::default();

    // Parse options
    for option in options {
        config.insert_option(option);
    }

    let mut fraglen_treatment = fraglen_treatment.clone();
    let mut fraglen_control   = fraglen_control  .clone();

    // Check fraglen arguments
    if fraglen_treatment.is_empty() {
        fraglen_treatment = vec![None; filenames_treatment.len()];
    }
    if fraglen_control.is_empty() {
        fraglen_control   = vec![None; filenames_control.len()];
    }
    if fraglen_treatment.len() != filenames_treatment.len() {
        let e = ArgumentError(
            format!("Number of provided treatment fragment lengths `{}` does not match number of treatment files `{}`",
                fraglen_treatment.len(), filenames_treatment.len())
        );
        return Err(CoverageError::new_empty(Box::new(e)));
    }
    if fraglen_control.len() != filenames_control.len() {
        let e = ArgumentError(
            format!("Number of provided control fragment lengths `{}` does not match number of control files `{}`",
                fraglen_control.len(), filenames_control.len())
        );
        return Err(CoverageError::new_empty(Box::new(e)));
    }

    // Read genome
    let mut genome: Genome = Genome::default();
    for filename in filenames_treatment.iter().chain(filenames_control.iter()) {
        let g = bam_import_genome(filename).map_err(
            |e| CoverageError::new_empty(e)
        )?;

        if genome.len() == 0 {
            genome = g;
        } else if genome != g {
            let e = ArgumentError("Treatment and control tracks have different genomes".to_string());
            return Err(CoverageError::new_empty(Box::new(e)));
        }
    }

    let mut treatment_fraglen_estimates = vec![FraglenEstimate::default(); filenames_treatment.len()];
    let mut   control_fraglen_estimates = vec![FraglenEstimate::default(); filenames_control  .len()];

    // No fragment length estimation, copy fragment lengths from arguments
    if !config.estimate_fraglen {
        for (i, fraglen) in fraglen_treatment.iter().enumerate() {
            if let Some(k) = *fraglen {
                treatment_fraglen_estimates[i].fraglen = k;
            } else {
                // If no fragment length is provided and also estimation is
                // switched off, then do not extend reads
                treatment_fraglen_estimates[i].fraglen = 0;
            }
        }
        for (i, fraglen) in fraglen_control.iter().enumerate() {
            if let Some(k) = *fraglen {
                control_fraglen_estimates[i].fraglen = k;
            } else {
                // If no fragment length is provided and also estimation is
                // switched off, then do not extend reads
                control_fraglen_estimates[i].fraglen = 0;
            }
        }
    }
    // Fragment length estimation
    else {

        for (i, filename) in filenames_treatment.iter().enumerate() {

            match fraglen_treatment[i] {
                Some(k) => treatment_fraglen_estimates[i].fraglen = k,
                None    => {
                    let estimate = estimate_fraglen(&config, filename, &genome).map_err(
                        |e| CoverageError::new_empty(e)
                    )?;

                    treatment_fraglen_estimates[i] = estimate;
                }
            }
        }

        for (i, filename) in filenames_control.iter().enumerate() {

            match fraglen_control[i] {
                Some(k) => control_fraglen_estimates[i].fraglen = k,
                None    => {
                    let estimate = estimate_fraglen(&config, filename, &genome).map_err(
                        |e| CoverageError::new_empty(e)
                    )?;

                    control_fraglen_estimates[i] = estimate;
                }
            }
        }
    }
    let fraglen_treatment_arg : Vec<usize> = treatment_fraglen_estimates.iter().map(|x| x.fraglen as usize).collect();
    let fraglen_control_arg   : Vec<usize> =   control_fraglen_estimates.iter().map(|x| x.fraglen as usize).collect();

    let track = SimpleTrack::coverage_from_bam(config, filenames_treatment, filenames_control, &fraglen_treatment_arg, &fraglen_control_arg, genome).map_err(
        |e| CoverageError::new(e, treatment_fraglen_estimates.clone(), control_fraglen_estimates.clone())
    )?;

    Ok((track, treatment_fraglen_estimates, control_fraglen_estimates))
}
