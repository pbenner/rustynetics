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

use std::fs::File;
use std::io::{BufRead, BufReader};

use flate2::read::GzDecoder;

use crate::genes::Genes;
use crate::granges::GRanges;
use crate::meta::MetaData;
use crate::utility::is_gzip;

/* -------------------------------------------------------------------------- */

impl Genes {
    fn set_expr_column(&mut self, expr: Vec<f64>) -> Result<(), Box<dyn std::error::Error>> {
        self.granges.meta.delete_meta("expr");
        self.granges.meta.add("expr", MetaData::FloatArray(expr))?;
        Ok(())
    }

    /// Parses expression data from a GTF file and stores it as the `expr` metadata column.
    ///
    /// Expression values are summed by matching the GTF `gene_id_name` field against the
    /// gene names stored in this `Genes` object.
    pub fn read_gtf_expr<R: BufRead>(
        &mut self,
        reader: R,
        gene_id_name: &str,
        expr_id_name: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let granges = GRanges::read_gtf(
            reader,
            vec![gene_id_name, expr_id_name],
            vec!["str", "float"],
            vec![],
        )?;

        let gene_ids = granges
            .meta
            .get_column_str(gene_id_name)
            .ok_or_else(|| format!("invalid gene_id_name `{}`", gene_id_name))?;
        let expr_values = granges
            .meta
            .get_column_float(expr_id_name)
            .ok_or_else(|| format!("invalid expr_id_name `{}`", expr_id_name))?;

        let mut expr = vec![0.0; self.granges.num_rows()];

        for i in 0..granges.num_rows() {
            if let Some(j) = self.find_gene(&gene_ids[i]) {
                expr[j] += expr_values[i];
            }
        }

        self.set_expr_column(expr)
    }

    /// Imports expression data from a GTF file and stores it as the `expr` metadata column.
    pub fn import_gtf_expr(
        &mut self,
        filename: &str,
        gene_id_name: &str,
        expr_id_name: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::open(filename)?;
        let reader: Box<dyn BufRead> = if is_gzip(filename) {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        self.read_gtf_expr(reader, gene_id_name, expr_id_name)
    }

    /// Parses a Cufflinks `genes.fpkm_tracking` file and stores the result as the `expr` metadata column.
    pub fn read_cufflinks_fpkm_tracking<R: BufRead>(
        &mut self,
        reader: R,
        verbose: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut expr = vec![0.0; self.granges.num_rows()];
        let mut lines = reader.lines();

        let header = lines.next().ok_or("invalid header")??;
        let header_fields: Vec<&str> = header.split_whitespace().collect();

        if header_fields.len() != 13
            || header_fields[0] != "tracking_id"
            || (header_fields[9] != "FPKM" && header_fields[9] != "RPKM")
        {
            return Err("invalid header".into());
        }

        for line in lines {
            let line = line?;
            let fields: Vec<&str> = line.split_whitespace().collect();

            if fields.is_empty() {
                continue;
            }
            if fields.len() != 13 {
                return Err("file must have 13 columns".into());
            }
            if fields[12] != "OK" {
                continue;
            }

            let gene_id = fields[0];
            let expr_value = fields[9].parse::<f64>()?;

            if let Some(i) = self.find_gene(gene_id) {
                expr[i] += expr_value;
            } else if verbose {
                eprintln!("`{}` not present in gene list!", gene_id);
            }
        }

        self.set_expr_column(expr)
    }

    /// Imports a Cufflinks `genes.fpkm_tracking` file and stores the result as the `expr` metadata column.
    pub fn import_cufflinks_fpkm_tracking(
        &mut self,
        filename: &str,
        verbose: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::open(filename)?;
        let reader: Box<dyn BufRead> = if is_gzip(filename) {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        self.read_cufflinks_fpkm_tracking(reader, verbose)
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::io::{BufReader, Cursor};

    use crate::genes::Genes;

    fn test_genes() -> Genes {
        Genes::new(
            vec!["geneA".into(), "geneB".into(), "geneC".into()],
            vec!["chr1".into(), "chr1".into(), "chr2".into()],
            vec![100, 300, 500],
            vec![200, 400, 600],
            vec![120, 320, 520],
            vec![180, 380, 580],
            vec!['+', '-', '+'],
        )
    }

    #[test]
    fn read_gtf_expr_sums_expression_by_gene() {
        let data = br#"chr1	src	transcript	100	200	.	+	.	transcript_id "geneA"; FPKM "1.5";
chr1	src	exon	120	180	.	+	.	transcript_id "geneA"; FPKM "2.5";
chr1	src	transcript	300	400	.	-	.	transcript_id "geneB"; FPKM "4.0";
"#;
        let mut genes = test_genes();

        genes
            .read_gtf_expr(
                BufReader::new(Cursor::new(&data[..])),
                "transcript_id",
                "FPKM",
            )
            .unwrap();

        assert_eq!(
            genes.granges.meta.get_column_float("expr").unwrap(),
            &vec![4.0, 4.0, 0.0]
        );
    }

    #[test]
    fn read_gtf_expr_replaces_existing_expr_column() {
        let data1 = br#"chr1	src	transcript	100	200	.	+	.	transcript_id "geneA"; FPKM "1.0";
"#;
        let data2 = br#"chr1	src	transcript	300	400	.	-	.	transcript_id "geneB"; FPKM "3.0";
"#;
        let mut genes = test_genes();

        genes
            .read_gtf_expr(
                BufReader::new(Cursor::new(&data1[..])),
                "transcript_id",
                "FPKM",
            )
            .unwrap();
        genes
            .read_gtf_expr(
                BufReader::new(Cursor::new(&data2[..])),
                "transcript_id",
                "FPKM",
            )
            .unwrap();

        assert_eq!(
            genes.granges.meta.get_column_float("expr").unwrap(),
            &vec![0.0, 3.0, 0.0]
        );
        assert_eq!(
            genes
                .granges
                .meta
                .meta_name
                .iter()
                .filter(|name| name.as_str() == "expr")
                .count(),
            1
        );
    }

    #[test]
    fn read_cufflinks_fpkm_tracking_imports_ok_rows() {
        let data = br#"tracking_id class_code nearest_ref_id gene_id gene_short_name tss_id locus length coverage FPKM FPKM_conf_lo FPKM_conf_hi status
geneA - - - - - chr1:100-200 100 10.0 1.5 0.0 2.0 OK
geneB - - - - - chr1:300-400 100 8.0 3.0 0.0 4.0 OK
geneX - - - - - chr1:500-600 100 5.0 7.0 0.0 8.0 OK
geneA - - - - - chr1:100-200 100 10.0 2.5 0.0 3.0 FAIL
"#;
        let mut genes = test_genes();

        genes
            .read_cufflinks_fpkm_tracking(BufReader::new(Cursor::new(&data[..])), false)
            .unwrap();

        assert_eq!(
            genes.granges.meta.get_column_float("expr").unwrap(),
            &vec![1.5, 3.0, 0.0]
        );
    }

    #[test]
    fn read_cufflinks_fpkm_tracking_rejects_invalid_header() {
        let data =
            br#"tracking_id class_code nearest_ref_id gene_id gene_short_name tss_id locus length coverage TPM FPKM_conf_lo FPKM_conf_hi status
"#;
        let mut genes = test_genes();

        assert!(genes
            .read_cufflinks_fpkm_tracking(BufReader::new(Cursor::new(&data[..])), false)
            .is_err());
    }
}
