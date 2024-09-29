
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
pub mod meta;
pub mod range;
pub mod reads;
pub mod error;
pub mod track;
pub mod track_bed;
pub mod track_bigwig;
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
