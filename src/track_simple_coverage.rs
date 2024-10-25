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

use crate::coverage::CoverageConfig;
use crate::genome::Genome;
use crate::track_generic::GenericMutableTrack;
use crate::track_simple::SimpleTrack;

/* -------------------------------------------------------------------------- */

impl SimpleTrack {

    pub fn coverage_from_bam(
        config             : CoverageConfig,
        filenames_treatment: &Vec<&str>,
        filenames_control  : &Vec<&str>,
        fraglen_treatment  : &Vec<usize>,
        fraglen_control    : &Vec<usize>,
        genome             : Genome,
    ) -> Result<SimpleTrack, Box<dyn Error>> {

        let mut track1 = SimpleTrack::alloc("treatment".to_string(), genome.clone(), config.initial_value, config.bin_size);
        if !filenames_control.is_empty() {
            // Control data
            let mut track2 = SimpleTrack::alloc("control".to_string(), genome.clone(), config.initial_value, config.bin_size);

            GenericMutableTrack::coverage_from_bam(
                config,
                GenericMutableTrack::wrap(&mut track1),
                Some(GenericMutableTrack::wrap(&mut track2)),
                filenames_treatment,
                filenames_control,
                fraglen_treatment,
                fraglen_control)?;

        } else {

            GenericMutableTrack::coverage_from_bam(
                config,
                GenericMutableTrack::wrap(&mut track1),
                None,
                filenames_treatment,
                filenames_control,
                fraglen_treatment,
                fraglen_control)?;

        };

        Ok(track1)

    }

}
