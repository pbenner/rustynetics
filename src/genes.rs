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

use std::fmt;
use std::collections::HashMap;
use std::error::Error;

use crate::range::Range;
use crate::genome::Genome;
use crate::granges::GRanges;
use crate::meta::MetaData;

/* -------------------------------------------------------------------------- */

/// The `Genes` struct represents a collection of genes with associated genomic ranges.
/// Each gene is characterized by its name, sequence name, transcript start and end positions,
/// coding sequence (CDS) start and end positions, and strand orientation. Additionally,
/// each gene is indexed for efficient retrieval by name.
///
/// # Fields
///
/// - `granges`: A `GRanges` object that stores the genomic ranges and associated metadata.
/// - `index`: A `HashMap` that maps gene names to their indices in `granges` for fast lookup.
///
/// # Examples
///
/// ```
/// use rustynetics::genes::Genes;
///
/// let names = vec!["gene1".to_string(), "gene2".to_string()];
/// let seqnames = vec!["chr1".to_string(), "chr2".to_string()];
/// let tx_from = vec![100, 200];
/// let tx_to = vec![150, 250];
/// let cds_from = vec![120, 220];
/// let cds_to = vec![140, 240];
/// let strand = vec!['+', '-'];
///
/// let genes = Genes::new(names, seqnames, tx_from, tx_to, cds_from, cds_to, strand);
/// ```
#[derive(Clone, Debug)]
pub struct Genes {
    pub granges: GRanges,
    index      : HashMap<String, usize>,
}

/* -------------------------------------------------------------------------- */

impl Genes {
    fn new_impl(granges: GRanges) -> Genes {
        let names = granges.meta.get_column_str  ("names").unwrap();
        let mut index = HashMap::new();
        for i in 0..granges.num_rows() {
            // check if strand is valid
            if granges.strand[i] != '+' && granges.strand[i] != '-' {
                panic!("invalid strand");
            }
            index.insert(names[i].clone(), i);
        }
        let genes = Genes {
            granges,
            index,
        };
        genes
    }

    /// Constructs a new `Genes` object.
    ///
    /// # Parameters
    ///
    /// - `names`: Names of the genes.
    /// - `seqnames`: Chromosome or sequence names for each gene.
    /// - `tx_from`: Transcript start positions for each gene.
    /// - `tx_to`: Transcript end positions for each gene.
    /// - `cds_from`: Coding sequence (CDS) start positions for each gene.
    /// - `cds_to`: CDS end positions for each gene.
    /// - `strand`: Strand orientation ('+' or '-') for each gene.
    ///
    /// # Panics
    ///
    /// Panics if any strand value is not '+' or '-', or if the input vectors have mismatched lengths.
    pub fn new(
        names   : Vec<String>,
        seqnames: Vec<String>,
        tx_from : Vec<usize>,
        tx_to   : Vec<usize>,
        cds_from: Vec<usize>,
        cds_to  : Vec<usize>,
        strand  : Vec<char>,
    ) -> Genes {
        assert!(names.len() == seqnames.len());
        assert!(names.len() == tx_from .len());
        assert!(names.len() == tx_to   .len());
        assert!(names.len() == cds_from.len());
        assert!(names.len() == cds_to  .len());
        assert!(names.len() == strand  .len());
        // construct cds ranges
        let mut cds = vec![];
        for i in 0..cds_from.len() {
            cds.push(Range::new(cds_from[i], cds_to[i]));
        }
        let mut granges = GRanges::new(
            seqnames,
            tx_from,
            tx_to,
            strand,
        );
        granges.meta.add("names", MetaData::StringArray(names)).unwrap();
        granges.meta.add("cds"  , MetaData::RangeArray(cds   )).unwrap();
        Genes::new_impl(granges)
    }

    /// Retrieves the names of all genes.
    ///
    /// # Returns
    ///
    /// A reference to a vector of gene names.
    pub fn names(&self) -> &Vec<String> {
        self.granges.meta.get_column_str("names").unwrap()
    }

    /// Retrieves the coding sequences (CDS) ranges for each gene.
    ///
    /// # Returns
    ///
    /// A reference to a vector of `Range` objects representing the CDS.
    pub fn cds(&self) -> &Vec<Range> {
        self.granges.meta.get_column_range("cds").unwrap()
    }

    /// Removes genes by their indices.
    ///
    /// # Parameters
    ///
    /// - `indices`: The indices of genes to remove.
    ///
    /// # Returns
    ///
    /// A new `Genes` object with the specified genes removed.
    pub fn remove(&self, indices: &[usize]) -> Genes {
        let r = self.granges.remove(indices);
        Genes::new_impl(r)
    }

    /// Removes genes that overlap with any gene in the specified `subject`.
    ///
    /// # Parameters
    ///
    /// - `subject`: A `Genes` object containing genes to check for overlaps.
    ///
    /// # Returns
    ///
    /// A new `Genes` object without overlaps with `subject`.
    pub fn remove_overlaps_with(&self, subject: &Genes) -> Genes {
        let r = self.granges.remove_overlaps_with(&subject.granges);
        Genes::new_impl(r)
    }

    /// Retains only genes that overlap with any gene in the specified `subject`.
    ///
    /// # Parameters
    ///
    /// - `subject`: A `Genes` object containing genes to check for overlaps.
    ///
    /// # Returns
    ///
    /// A new `Genes` object with only the overlapping genes.
    pub fn keep_overlaps_with(&self, subject: &Genes) -> Genes {
        let r = self.granges.keep_overlaps_with(&subject.granges);
        Genes::new_impl(r)
    }

    /// Creates a subset of genes by the specified indices.
    ///
    /// # Parameters
    ///
    /// - `indices`: Indices of genes to include in the subset.
    ///
    /// # Returns
    ///
    /// A new `Genes` object containing only the specified genes.
    pub fn subset(&self, indices: &[usize]) -> Genes {
        let r = self.granges.subset(indices);
        Genes::new_impl(r)
    }

    /// Slices a range of genes between specified indices.
    ///
    /// # Parameters
    ///
    /// - `ifrom`: The starting index of the slice.
    /// - `ito`: The ending index of the slice (exclusive).
    ///
    /// # Returns
    ///
    /// A new `Genes` object containing the genes within the specified range.
    pub fn slice(&self, ifrom: usize, ito: usize) -> Genes {
        let r = self.granges.slice(ifrom, ito);
        Genes::new_impl(r)
    }

    /// Finds the intersection of genes with another `Genes` object.
    ///
    /// # Parameters
    ///
    /// - `subject`: The `Genes` object to intersect with.
    ///
    /// # Returns
    ///
    /// A new `Genes` object containing the intersection.
    pub fn intersection(&self, subject: &Genes) -> Genes {
        let r = self.granges.intersection(&subject.granges);
        Genes::new_impl(r)
    }

    /// Sorts genes based on the specified metadata column name.
    ///
    /// # Parameters
    ///
    /// - `name`: The name of the metadata column to sort by.
    /// - `reverse`: A boolean indicating whether to sort in descending order.
    ///
    /// # Returns
    ///
    /// A sorted `Genes` object or an error if sorting fails.
    pub fn sort(&self, name: &str, reverse: bool) -> Result<Genes, Box<dyn Error>> {
        let r = self.granges.sort(name, reverse)?;
        Ok(Genes::new_impl(r))
    }

    /// Filters genes based on a specified `Genome`.
    ///
    /// # Parameters
    ///
    /// - `genome`: The `Genome` object used to filter genes.
    ///
    /// # Returns
    ///
    /// A new `Genes` object filtered by the genome.
    pub fn filter_genome(&self, genome: &Genome) -> Genes {
        let r = self.granges.filter_genome(genome);
        Genes::new_impl(r)
    }

    /// Finds the index of a gene by name.
    ///
    /// # Parameters
    ///
    /// - `name`: The name of the gene to find.
    ///
    /// # Returns
    ///
    /// An `Option` containing the gene index, or `None` if the gene is not found.
    pub fn find_gene(&self, name: &str) -> Option<usize> {
        self.index.get(name).cloned()
    }

}

/* -------------------------------------------------------------------------- */

impl PartialEq for Genes {
    fn eq(&self, other: &Self) -> bool {
        self.granges == other.granges
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for Genes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad(&format!("{}", self.granges))
    }
}
