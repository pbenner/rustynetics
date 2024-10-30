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
use mysql::prelude::*;
use mysql::from_row;
use mysql::Pool;

use crate::genes::Genes;
use crate::utility::is_gzip;

/* -------------------------------------------------------------------------- */

impl Genes {

    /// Imports gene information from a file.
    ///
    /// # Parameters
    ///
    /// - `filename`: The path to the file containing gene data. The file must contain
    ///   seven whitespace-separated columns: `name`, `seqname`, `tx_from`, `tx_to`,
    ///   `cds_from`, `cds_to`, and `strand`. The file can be plain text or gzip-compressed.
    ///
    /// # Returns
    ///
    /// A `Result` with a `Genes` object on success, or an error if the file format
    /// is incorrect or any field parsing fails.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The file cannot be opened.
    /// - The file contains rows with fewer than seven columns.
    /// - Any numeric field fails to parse.
    /// - Strand column contains an invalid character.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use rustynetics::genes::Genes;
    ///
    /// let genes = Genes::import_genes("path/to/genes.txt")?;
    /// ```
    pub fn import_genes(filename: &str) -> Result<Genes, Box<dyn std::error::Error>> {
        let mut names    = vec![];
        let mut seqnames = vec![];
        let mut tx_from  = vec![];
        let mut tx_to    = vec![];
        let mut cds_from = vec![];
        let mut cds_to   = vec![];
        let mut strand   = vec![];

        let file = File::open(filename)?;
        let reader: Box<dyn BufRead> = if is_gzip(filename) {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        for line_ in reader.lines() {
            let line = line_?;
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.is_empty() {
                continue;
            }
            if fields.len() != 7 {
                return Err("file must have seven columns".into());
            }
            names   .push(String::from(fields[0]));
            seqnames.push(String::from(fields[1]));
            tx_from .push(fields[2].parse()?);
            tx_to   .push(fields[3].parse()?);
            cds_from.push(fields[4].parse()?);
            cds_to  .push(fields[5].parse()?);
            strand  .push(fields[6].chars().next().unwrap());
        }
        Ok(Genes::new(names, seqnames, tx_from, tx_to, cds_from, cds_to, strand))
    }

    /// Imports gene data from the UCSC genome database.
    ///
    /// # Parameters
    ///
    /// - `genome`: The UCSC genome assembly name (e.g., `hg19` or `hg38`).
    /// - `table`: The UCSC database table name (e.g., `refGene` or `knownGene`).
    ///
    /// # Returns
    ///
    /// A `Result` with a `Genes` object on success, or an error if the connection or query fails.
    ///
    /// # Details
    ///
    /// This function connects to the UCSC MySQL database for the specified genome and retrieves
    /// gene information from the specified table. It queries for `name`, `chrom`, `strand`,
    /// `txStart`, `txEnd`, `cdsStart`, and `cdsEnd`, converting these fields into the required
    /// format to create a `Genes` object.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The MySQL connection fails.
    /// - The table is not accessible.
    /// - There are issues in parsing or processing the retrieved rows.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use rustynetics::genes::Genes;
    ///
    /// let genes = Genes::import_genes_from_ucsc("hg38", "refGene")?;
    /// ```
    pub fn import_genes_from_ucsc(genome: &str, table: &str) -> Result<Genes, Box<dyn std::error::Error>> {
        let mut names    = vec![];
        let mut seqnames = vec![];
        let mut tx_from  = vec![];
        let mut tx_to    = vec![];
        let mut cds_from = vec![];
        let mut cds_to   = vec![];
        let mut strand   = vec![];

        let url      = format!("mysql://genome@genome-mysql.cse.ucsc.edu:3306/{}", genome);
        let pool     = Pool::new(url.as_str())?;
        let mut conn = pool.get_conn()?;
        let query    = format!("SELECT name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd FROM {}", table);

        let mut result = conn.query_iter(query)?;

        while let Some(result_set) = result.iter() {

            for row in result_set {
                let r: (String, String, String, i32, i32, i32, i32) = from_row(row.unwrap());

                names   .push(r.0.clone());
                seqnames.push(r.1.clone());
                strand  .push(r.2.chars().next().unwrap());
                tx_from .push(r.3 as usize);
                tx_to   .push(r.4 as usize);
                cds_from.push(r.5 as usize);
                cds_to  .push(r.6 as usize);
            }
        }

        Ok(Genes::new(names, seqnames, tx_from, tx_to, cds_from, cds_to, strand))
    }

}
