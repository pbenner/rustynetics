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

use std::process;

use clap::{Arg, Command};

use rustynetics::bigwig::BigWigFile;

/* -------------------------------------------------------------------------- */

fn query(filename_in: &str, chrom: &str, from: usize, to: usize, bin_size: usize, verbose: bool) {

    if verbose {
        eprintln!("Opening bigWig file {}", filename_in);
    }

    // Open the BigWig file
    let mut reader = BigWigFile::new_reader(filename_in).unwrap_or_else(|err| {
        eprintln!("Error opening file: {}", err);
        process::exit(1);
    });
    let genome = reader.genome().clone();

    // Query the BigWig file
    for result in reader.query(chrom, from, to, bin_size) {
        match result {
            Ok(record) => {
                println!("{}:[{}, {})={}",
                    genome.seqnames[record.data.chrom_id as usize],
                    record.data.from, record.data.to, record.data.statistics);
            }
            Err(err) => {
                eprintln!("Error querying BigWig file: {}", err);
                process::exit(1);
            }
        }
    }
}

/* -------------------------------------------------------------------------- */

fn main() {
    let matches = Command::new("BigWig Query")
        .version("1.0")
        .author("Philipp Benner [https://github.com/pbenner]")
        .about("Query BigWig files")
        .arg(
            Arg::new("input")
                .help("The input BigWig file")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("chrom")
                .help("The chromosome to query")
                .required(true)
                .index(2),
        )
        .arg(
            Arg::new("from")
                .help("The start position")
                .required(true)
                .index(3),
        )
        .arg(
            Arg::new("to")
                .help("The end position")
                .required(true)
                .index(4),
        )
        .arg(
            Arg::new("binsize")
                .help("The bin size for the query")
                .required(true)
                .index(5),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .help("Be verbose"))
        .get_matches();

    let filename_in = matches.get_one::<String>("input").expect("Input file is required");
    let chrom = matches.get_one::<String>("chrom").expect("Chromosome is required");
    let from: usize = matches
        .get_one::<String>("from")
        .expect("Start position is required")
        .parse()
        .unwrap_or_else(|_| {
            eprintln!("Invalid start position");
            process::exit(1);
        });
    let to: usize = matches
        .get_one::<String>("to")
        .expect("End position is required")
        .parse()
        .unwrap_or_else(|_| {
            eprintln!("Invalid end position");
            process::exit(1);
        });
    let bin_size: usize = matches
        .get_one::<String>("binsize")
        .expect("Bin size is required")
        .parse()
        .unwrap_or_else(|_| {
            eprintln!("Invalid bin size");
            process::exit(1);
        });
    let verbose = matches.get_flag("verbose");

    query(filename_in, chrom, from, to, bin_size, verbose);

}
