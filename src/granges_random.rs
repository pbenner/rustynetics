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

use rand::Rng;
use rand::prelude::SliceRandom;

use crate::genome::Genome;
use crate::granges::GRanges;

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
struct GenomeRng {
    weights: Vec<f64>,
    genome : Genome,
    max_len: usize,
}

/* -------------------------------------------------------------------------- */

impl GenomeRng {
    fn new(genome: &Genome) -> GenomeRng {
        let max_len = genome.lengths.iter().max().unwrap_or(&0);
        let weights: Vec<f64> = genome
            .lengths
            .iter()
            .map(|&length| length as f64)
            .collect();
        let sum    : f64      = weights.iter().sum();
        let weights: Vec<f64> = weights.iter().map(|&weight| weight / sum).collect();
        let mut cumulative_probabilities = vec![weights[0]];
        for i in 1..weights.len() {
            cumulative_probabilities.push(cumulative_probabilities[i - 1] + weights[i]);
        }
        GenomeRng {
            weights: cumulative_probabilities,
            genome : genome.clone(),
            max_len: *max_len,
        }
    }

    fn draw(&self, wsize: usize) -> Option<(usize, usize)> {
        if wsize > self.max_len {
            return None;
        }
        let mut rng = rand::thread_rng();
        let mut k = 0;
        let mut t = 0.0;
        loop {
            let p = rng.gen::<f64>();
            for i in 0..self.weights.len() {
                if t <= p && p < self.weights[i] {
                    k = i;
                    break;
                }
                t = self.weights[i];
            }
            if self.genome.lengths[k] >= wsize {
                break;
            }
        }
        let i = rng.gen_range(0..=self.genome.lengths[k] - wsize);
        Some((k, i))
    }
}

/* -------------------------------------------------------------------------- */

impl GRanges {

    pub fn random(n: usize, wsize: usize, genome: &Genome, use_strand: bool) -> GRanges {
        let gnome_rng = GenomeRng::new(genome);
        let mut seqnames = Vec::with_capacity(n);
        let mut from     = Vec::with_capacity(n);
        let mut to       = Vec::with_capacity(n);
        let mut strand   = Vec::with_capacity(n);
        let mut rng      = rand::thread_rng();
        for _ in 0..n {
            let (j, position) = match gnome_rng.draw(wsize) {
                Some(v) => v,
                None    => return GRanges::new_empty()
            };
            seqnames.push(genome.seqnames[j]);
            from    .push(position);
            to      .push(position + wsize);
            if use_strand {
                let k = rng.gen_range(0..2);
                strand.push(['+', '-'][k]);
            }
        }
        GRanges::new(seqnames, from, to, strand)
    }

    pub fn random_permutation(&self) -> GRanges {
        let mut rng = rand::thread_rng();
        let mut idx: Vec<usize> = (0..self.num_rows()).collect();
        idx.shuffle(&mut rng);
        self.subset(&idx)
    }
}

