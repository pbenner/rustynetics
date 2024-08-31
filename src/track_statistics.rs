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

use std::f64;

use crate::genome::Genome;
use crate::reads::Read;
use crate::track::Track;
use crate::track_simple::SimpleTrack;
use crate::track_generic::GenericMutableTrack;

/* -------------------------------------------------------------------------- */

// Type alias for BinSummaryStatistics function
pub type BinSummaryStatistics = fn(f64, f64, f64, f64, f64) -> f64;

/* -------------------------------------------------------------------------- */

// Function implementations for BinSummaryStatistics
fn bin_mean(sum: f64, _sum_squares: f64, _min: f64, _max: f64, n: f64) -> f64 {
    sum / n
}

fn bin_max(_sum: f64, _sum_squares: f64, _min: f64, max: f64, _n: f64) -> f64 {
    max
}

fn bin_min(_sum: f64, _sum_squares: f64, min: f64, _max: f64, _n: f64) -> f64 {
    min
}

fn bin_discrete_mean(sum: f64, _sum_squares: f64, _min: f64, _max: f64, n: f64) -> f64 {
    (sum / n).floor() + 0.5
}

fn bin_discrete_max(_sum: f64, _sum_squares: f64, _min: f64, max: f64, _n: f64) -> f64 {
    max.floor()
}

fn bin_discrete_min(_sum: f64, _sum_squares: f64, min: f64, _max: f64, _n: f64) -> f64 {
    min.floor()
}

fn bin_variance(sum: f64, sum_squares: f64, _min: f64, _max: f64, n: f64) -> f64 {
    sum_squares / n - (sum / n) * (sum / n)
}

/* -------------------------------------------------------------------------- */

pub fn bin_summary_statistics_from_string(s: &str) -> Option<BinSummaryStatistics> {
    match s {
        "mean"          => Some(bin_mean),
        "max"           => Some(bin_max),
        "min"           => Some(bin_min),
        "discrete mean" => Some(bin_discrete_mean),
        "discrete max"  => Some(bin_discrete_max),
        "discrete min"  => Some(bin_discrete_min),
        "variance"      => Some(bin_variance),
        _               => None,
    }
}

/* -------------------------------------------------------------------------- */

// Helper function for integer division rounding up
fn div_int_up(a: i32, b: i32) -> i32 {
    (a + b - 1) / b
}

/* -------------------------------------------------------------------------- */

// Function to compute the cross-correlation between two tracks
pub fn track_crosscorrelation(
    track1: &dyn Track,
    track2: &dyn Track,
    from: i32,
    to: i32,
    normalize: bool,
) -> Result<(Vec<i32>, Vec<f64>), String> {
    if from < 0 || to < from {
        return Err("Crosscorrelation(): invalid parameters".to_string());
    }

    if track1.get_bin_size() != track2.get_bin_size() {
        return Err("Crosscorrelation(): track binSizes do not match".to_string());
    }

    let mut x = Vec::new();
    let mut y = Vec::new();

    for name in track1.get_seq_names() {
        let sequence1 = track1.get_sequence(&name)?;
        if let Ok(sequence2) = track2.get_sequence(&name) {
            if sequence1.n_bins() != sequence2.n_bins() {
                return Err("Crosscorrelation(): track sequence lengths do not match".to_string());
            }
        } else {
            continue;
        }
    }

    let b = track1.get_bin_size();
    let n = div_int_up(to - from, b as i32);
    let mut m = 0.0;
    let mut mean1 = 0.0;
    let mut mean2 = 0.0;
    let mut variance1 = 1.0;
    let mut variance2 = 1.0;

    x.resize(n as usize, 0);
    y.resize(n as usize, 0.0);

    for (j, l) in (from..to).step_by(b as usize).enumerate() {
        x[j] = l / (b as i32);
    }

    if normalize {
        for name in track1.get_seq_names() {
            let sequence1 = track1.get_sequence(&name)?;
            if let Ok(sequence2) = track2.get_sequence(&name) {
                let (mut s1, mut s2, mut t1, mut t2) = (0.0, 0.0, 0.0, 0.0);

                for i in 0..sequence1.n_bins() {
                    s1 += sequence1.at_bin(i);
                    s2 += sequence2.at_bin(i);
                    t1 += sequence1.at_bin(i) * sequence1.at_bin(i);
                    t2 += sequence2.at_bin(i) * sequence2.at_bin(i);
                }

                let k = sequence1.n_bins() as f64;
                mean1 = m / (m + k) * mean1 + 1.0 / (m + k) * s1;
                mean2 = m / (m + k) * mean2 + 1.0 / (m + k) * s2;
                variance1 = m / (m + k) * variance1 + 1.0 / (m + k) * t1;
                variance2 = m / (m + k) * variance2 + 1.0 / (m + k) * t2;
                m += k;
            } else {
                continue;
            }
        }
        variance1 -= mean1 * mean1;
        variance2 -= mean2 * mean2;
    }

    m = 0.0;
    for name in track1.get_seq_names() {
        let sequence1 = track1.get_sequence(&name)?;
        if let Ok(sequence2) = track2.get_sequence(&name) {
            let mut s = vec![0.0; n as usize];

            for i in 0..sequence1.n_bins() {
                for j in 0..n as usize {
                    if i + (x[j] as usize) < sequence1.n_bins() {
                        s[j] += (sequence1.at_bin(i) - mean1) * (sequence2.at_bin(i + x[j] as usize) - mean2);
                    }
                }
            }

            let k = sequence1.n_bins() as f64;
            for j in 0..n as usize {
                y[j] = m / (m + k) * y[j] + 1.0 / (m + k) * s[j];
            }
            m += k;
        } else {
            continue;
        }
    }

    for j in 0..n as usize {
        x[j] *= b as i32;
        y[j] /= (variance1 * variance2).sqrt();
    }

    Ok((x, y))
}

/* -------------------------------------------------------------------------- */

// Function to cross-correlate reads on forward and reverse strands
pub fn crosscorrelate_reads(
    reads    : impl Iterator<Item = Read>,
    genome   : &Genome,
    max_delay: i32,
    bin_size : usize,
) -> Result<(Vec<i32>, Vec<f64>, i32, u64), String> {

    let mut track1 = SimpleTrack::alloc("forward".to_string(), genome.clone(), bin_size);
    let mut track2 = SimpleTrack::alloc("reverse".to_string(), genome.clone(), bin_size);
    let mut n = 0_u64;
    let mut read_length = 0_u64;

    for read in reads {
        if read.strand == '+' {
            let mut r = read.clone();
            r.range.to = r.range.from + 1;
            if GenericMutableTrack::wrap(&mut track1).add_read(&r, 0).is_ok() {
                read_length += (r.range.to - r.range.from) as u64;
                n += 1;
            }
        } else if read.strand == '-' {
            let mut r = read.clone();
            r.range.from = r.range.to - 1;
            r.range.from += 1;
            r.range.to += 1;
            if GenericMutableTrack::wrap(&mut track2).add_read(&r, 0).is_ok() {
                read_length += (r.range.to - r.range.from) as u64;
                n += 1;
            }
        }
    }

    if n == 0 {
        return Err("computing cross-correlation failed: no reads available".to_string());
    }

    read_length /= n;
    let (x, y) = track_crosscorrelation(&track1, &track2, 0, max_delay, true)?;

    Ok((x, y, read_length as i32, n))
}

/* -------------------------------------------------------------------------- */

// Function to estimate fragment length
pub fn estimate_fragment_length(
    reads        : impl Iterator<Item = Read>,
    genome       : &Genome,
    max_delay    : i32,
    bin_size     : usize,
    fraglen_range: (i32, i32),
) -> Result<(i32, Vec<i32>, Vec<f64>, u64), String> {
    let (x, y, read_length, n) = crosscorrelate_reads(reads, genome, max_delay, bin_size)?;

    let mut from = (read_length + read_length / 2) as i32;
    let mut to   = max_delay;

    if fraglen_range.0 != -1 {
        from = fraglen_range.0;
    }

    if fraglen_range.1 != -1 {
        to = fraglen_range.1;
    }

    let mut i        = from / (bin_size as i32);
    let mut max      = y[i as usize];
    let mut frag_len = x[i as usize];

    while i < to / (bin_size as i32) {
        if y[i as usize] > max {
            max = y[i as usize];
            frag_len = x[i as usize];
        }
        i += 1;
    }

    Ok((frag_len, x, y, n))
}
