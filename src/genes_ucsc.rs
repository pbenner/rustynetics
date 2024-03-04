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
use mysql::Pool;

use crate::genes::Genes;

/* -------------------------------------------------------------------------- */

impl Genes {

    fn read_ucsc_genes(filename: &str) -> Result<Genes, Box<dyn std::error::Error>> {
        let names    = vec![];
        let seqnames = vec![];
        let tx_from  = vec![];
        let tx_to    = vec![];
        let cds_from = vec![];
        let cds_to   = vec![];
        let strand   = vec![];

        let file = File::open(filename)?;
        let reader: Box<dyn BufRead> = if filename.ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.is_empty() {
                continue;
            }
            if fields.len() != 7 {
                return Err("file must have seven columns".into());
            }
            names   .push(fields[0]);
            seqnames.push(fields[1]);
            tx_from .push(fields[2].parse()?);
            tx_to   .push(fields[3].parse()?);
            cds_from.push(fields[4].parse()?);
            cds_to  .push(fields[5].parse()?);
            strand  .push(fields[6].chars().next().unwrap());
        }
        Ok(Genes::new(names, seqnames, tx_from, tx_to, cds_from, cds_to, strand))
    }

    fn import_genes_from_ucsc(genome: &str, table: &str) -> Result<Genes, Box<dyn std::error::Error>> {
        let names    = vec![];
        let seqnames = vec![];
        let tx_from  = vec![];
        let tx_to    = vec![];
        let cds_from = vec![];
        let cds_to   = vec![];
        let strand   = vec![];

        let url      = format!("genome@tcp(genome-mysql.cse.ucsc.edu:3306)/{}", genome).as_str();
        let pool     = Pool::new(url)?;
        let mut conn = pool.get_conn()?;
        let query    = format!("SELECT name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd FROM {}", table);

        conn.query_map(query, |(name_, seqname_, strand_, tx_from_, tx_to_, cds_from_, cds_to_)| {
            names   .push(name_    );
            seqnames.push(seqname_ );
            strand  .push(strand_  );
            tx_from .push(tx_from_ );
            tx_to   .push(tx_to_   );
            cds_from.push(cds_from_);
            cds_to  .push(cds_to_  );
        })?;

        Ok(Genes::new(names, seqnames, tx_from, tx_to, cds_from, cds_to, strand))
    }

}
