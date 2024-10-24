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
