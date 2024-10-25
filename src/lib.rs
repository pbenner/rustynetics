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

/* -------------------------------------------------------------------------- */

// Public crates
pub mod alphabet;
pub mod bam;
pub mod bbi;
pub mod bgzf;
pub mod bigwig;
pub mod cpg;
pub mod genes;
pub mod genes_ucsc;
pub mod genome;
pub mod granges;
pub mod granges_bam;
pub mod granges_bed;
pub mod granges_error;
pub mod granges_gtf;
pub mod granges_row;
pub mod granges_table;
pub mod granges_pretty;
pub mod granges_random;
pub mod granges_merge;
pub mod granges_sort;
#[macro_use]
pub mod infologger;
pub mod meta;
pub mod range;
pub mod read;
pub mod read_stream;
pub mod error;
pub mod track;
pub mod track_bed;
pub mod track_bigwig;
pub mod track_coverage;
pub mod track_generic;
pub mod track_granges;
pub mod track_simple;
pub mod track_simple_bedgraph;
pub mod track_simple_bigwig;
pub mod track_simple_wig;
pub mod track_statistics;

// Private crates
mod granges_find_endpoint;
mod granges_find_nearest;
mod granges_find_overlaps;
mod granges_table_reader;
mod granges_table_writer;
mod meta_pretty;
mod meta_table;
mod meta_table_reader;
mod meta_table_writer;
mod netfile;
mod utility;
mod utility_cumdist;
mod utility_io;

// Macros
extern crate approx;
