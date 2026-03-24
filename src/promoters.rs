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

use std::collections::HashMap;
use std::error::Error;
use std::io;

use crate::genes::Genes;
use crate::granges::GRanges;
use crate::meta::{Meta, MetaData};

/* -------------------------------------------------------------------------- */

pub fn promoters(genes: &Genes, offset1: isize, offset2: isize) -> Result<GRanges, Box<dyn Error>> {
    let mut seqnames = Vec::with_capacity(genes.granges.num_rows());
    let mut from = Vec::with_capacity(genes.granges.num_rows());
    let mut to = Vec::with_capacity(genes.granges.num_rows());
    let mut strand = Vec::with_capacity(genes.granges.num_rows());
    let mut names = Vec::with_capacity(genes.granges.num_rows());

    for i in 0..genes.granges.num_rows() {
        let tx_from = genes.granges.ranges[i].from as isize;
        let tx_to = genes.granges.ranges[i].to as isize;

        let (promoter_from, promoter_to) = match genes.granges.strand[i] {
            '+' => (tx_from - offset1, tx_from + offset2 + 1),
            '-' => (tx_to - 1 - offset2, tx_to - 1 + offset1 + 1),
            _ => {
                return Err(Box::new(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("gene `{i}` has no strand information"),
                )))
            }
        };

        if promoter_to < promoter_from {
            return Err(Box::new(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("gene `{i}` produces an invalid promoter window"),
            )));
        }

        seqnames.push(genes.granges.seqnames[i].clone());
        from.push(promoter_from.max(0) as usize);
        to.push(promoter_to.max(0) as usize);
        strand.push(genes.granges.strand[i]);
        names.push(genes.names()[i].clone());
    }

    let mut result = GRanges::new(seqnames, from, to, strand);
    result.meta = genes.granges.meta.clone();
    result.meta.delete_meta("cds");
    replace_meta_column(&mut result.meta, "names", MetaData::StringArray(names))?;

    merge_duplicate_promoters(result)
}

/* -------------------------------------------------------------------------- */

impl Genes {
    pub fn promoters(&self, offset1: isize, offset2: isize) -> Result<GRanges, Box<dyn Error>> {
        promoters(self, offset1, offset2)
    }
}

/* -------------------------------------------------------------------------- */

fn replace_meta_column(meta: &mut Meta, name: &str, data: MetaData) -> Result<(), Box<dyn Error>> {
    meta.delete_meta(name);
    meta.add(name, data)
}

fn merge_duplicate_promoters(granges: GRanges) -> Result<GRanges, Box<dyn Error>> {
    let mut groups = Vec::<Vec<usize>>::new();
    let mut index = HashMap::<(String, usize, usize), usize>::new();

    for i in 0..granges.num_rows() {
        let key = (
            granges.seqnames[i].clone(),
            granges.ranges[i].from,
            granges.ranges[i].to,
        );
        if let Some(group) = index.get(&key).copied() {
            groups[group].push(i);
        } else {
            index.insert(key, groups.len());
            groups.push(vec![i]);
        }
    }

    let first_indices: Vec<usize> = groups.iter().map(|group| group[0]).collect();
    let mut result = granges.subset(&first_indices);
    result.meta = merge_meta_groups(&granges.meta, &groups)?;
    Ok(result)
}

fn merge_meta_groups(meta: &Meta, groups: &[Vec<usize>]) -> Result<Meta, Box<dyn Error>> {
    let mut merged = Meta::default();

    for (name, data) in meta.iter() {
        let data = match data {
            MetaData::StringArray(values) => MetaData::StringMatrix(
                groups
                    .iter()
                    .map(|group| group.iter().map(|&i| values[i].clone()).collect())
                    .collect(),
            ),
            MetaData::FloatArray(values) => MetaData::FloatMatrix(
                groups
                    .iter()
                    .map(|group| group.iter().map(|&i| values[i]).collect())
                    .collect(),
            ),
            MetaData::IntArray(values) => MetaData::IntMatrix(
                groups
                    .iter()
                    .map(|group| group.iter().map(|&i| values[i]).collect())
                    .collect(),
            ),
            MetaData::RangeArray(values) => {
                let mut merged_ranges = Vec::with_capacity(groups.len());
                for group in groups {
                    let first = values[group[0]];
                    if group.iter().any(|&i| values[i] != first) {
                        return Err(Box::new(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!(
                                "cannot merge range metadata column `{name}` with differing values"
                            ),
                        )));
                    }
                    merged_ranges.push(first);
                }
                MetaData::RangeArray(merged_ranges)
            }
            MetaData::StringMatrix(_) | MetaData::FloatMatrix(_) | MetaData::IntMatrix(_) => {
                return Err(Box::new(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("cannot merge matrix metadata column `{name}`"),
                )));
            }
        };
        merged.add(name, data)?;
    }

    Ok(merged)
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::PathBuf;

    use crate::genes::Genes;
    use crate::meta::MetaData;
    use crate::promoters::promoters;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("gonetics")
            .join(name)
    }

    fn import_ucsc_genes_fixture(name: &str) -> Genes {
        let path = fixture(name);
        let file = File::open(path).unwrap();
        let decoder = flate2::read::GzDecoder::new(file);
        let reader = BufReader::new(decoder);

        let mut names = Vec::new();
        let mut seqnames = Vec::new();
        let mut tx_from = Vec::new();
        let mut tx_to = Vec::new();
        let mut cds_from = Vec::new();
        let mut cds_to = Vec::new();
        let mut strand = Vec::new();

        for line in reader.lines() {
            let line = line.unwrap();
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.is_empty() {
                continue;
            }
            names.push(fields[0].to_string());
            seqnames.push(fields[1].to_string());
            strand.push(fields[2].chars().next().unwrap());
            tx_from.push(fields[3].parse().unwrap());
            tx_to.push(fields[4].parse().unwrap());
            cds_from.push(fields[5].parse().unwrap());
            cds_to.push(fields[6].parse().unwrap());
        }

        Genes::new(names, seqnames, tx_from, tx_to, cds_from, cds_to, strand)
    }

    #[test]
    fn promoters_merge_duplicate_windows_and_metadata() {
        let mut genes = Genes::new(
            vec!["g1".into(), "g2".into(), "g3".into()],
            vec!["chr1".into(), "chr1".into(), "chr2".into()],
            vec![100, 100, 200],
            vec![150, 150, 260],
            vec![110, 110, 210],
            vec![140, 140, 250],
            vec!['+', '+', '-'],
        );
        genes
            .granges
            .meta
            .add("expr", MetaData::FloatArray(vec![1.0, 2.0, 3.0]))
            .unwrap();

        let promoters = promoters(&genes, 10, 5).unwrap();

        assert_eq!(promoters.num_rows(), 2);
        assert_eq!(
            promoters.seqnames,
            vec!["chr1".to_string(), "chr2".to_string()]
        );
        assert_eq!(promoters.ranges[0].from, 90);
        assert_eq!(promoters.ranges[0].to, 106);
        assert_eq!(promoters.ranges[1].from, 254);
        assert_eq!(promoters.ranges[1].to, 270);
        assert!(promoters.meta.get_column_range("cds").is_none());
        assert_eq!(
            promoters.meta.get_column("names"),
            Some(&MetaData::StringMatrix(vec![
                vec!["g1".into(), "g2".into()],
                vec!["g3".into()]
            ]))
        );
        assert_eq!(
            promoters.meta.get_column("expr"),
            Some(&MetaData::FloatMatrix(vec![vec![1.0, 2.0], vec![3.0]]))
        );
    }

    #[test]
    fn promoters_match_go_fixture_count() {
        let genes = import_ucsc_genes_fixture("Data/hg19.knownGene.txt.gz");
        let promoters = genes.promoters(500, 500).unwrap();

        assert_eq!(promoters.num_rows(), 51384);
    }
}
