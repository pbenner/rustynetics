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
                None    => return GRanges::default()
            };
            seqnames.push(genome.seqnames[j].clone());
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

