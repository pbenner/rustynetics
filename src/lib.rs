
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

pub mod alphabet;
pub mod bbi;
pub mod bigwig;
pub mod genes;
pub mod genes_ucsc;
pub mod genome;
pub mod granges;
pub mod granges_row;
pub mod meta;
pub mod range;
pub mod error;
pub mod options_print;
pub mod track;
pub mod track_simple;

// Private crates
mod granges_find_endpoint;
mod granges_find_nearest;
mod granges_find_overlaps;
mod granges_bed;
mod granges_pretty;
mod granges_random;
mod granges_merge;
mod granges_sort;
mod granges_table;
mod granges_table_reader;
mod granges_table_writer;
mod meta_pretty;
mod meta_table;
mod meta_table_reader;
mod meta_table_writer;
mod netfile;
mod utility;

// Tests
mod bigwig_test;

// Macros
extern crate approx;
