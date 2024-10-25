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

    // Query the BigWig file
    for result in reader.query(chrom, from, to, bin_size) {
        match result {
            Ok(record) => {
                println!("{}:[{}, {})={}",
                    record.data.chrom,
                    record.data.from,
                    record.data.to,
                    record.data.statistics);
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
                .action(clap::ArgAction::SetTrue)
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
