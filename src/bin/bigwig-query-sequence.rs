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
use std::error::Error;

use clap::{Arg, Command};

use rustynetics::bigwig::BigWigFile;
use rustynetics::track_statistics::{BinSummaryStatistics, bin_summary_statistics_from_string};

/* -------------------------------------------------------------------------- */

fn query(
    filename_in: &str,
    chrom      : &str,
    from       : usize,
    to         : usize,
    bin_size   : usize,
    bin_overlap: usize,
    summary    : BinSummaryStatistics,
    verbose    : bool,
) -> Result<(), Box<dyn Error>> {

    if verbose {
        eprintln!("Opening bigWig file {}", filename_in);
    }

    // Open BigWig file
    let mut reader = BigWigFile::new_reader(filename_in).unwrap_or_else(|err| {
        eprintln!("Error opening file: {}", err);
        process::exit(1);
    });

    // Query BigWig file
    let result = reader.query_slice(
        chrom,
        from,
        to,
        summary,
        bin_size,
        bin_overlap,
        f64::NAN,
    )?;

    println!("{:?}", result.0);

    Ok(())
}

/* -------------------------------------------------------------------------- */

fn main() {
    let matches = Command::new("BigWig Query")
        .version("1.0")
        .about("Query BigWig files")
        .arg(Arg::new("input")
            .help("The input BigWig file")
            .required(true))
        .arg(Arg::new("chrom")
            .help("The chromosome to query")
            .required(true))
        .arg(Arg::new("from")
            .help("The start position")
            .required(true))
        .arg(Arg::new("to")
            .help("The end position")
            .required(true))
        .arg(Arg::new("binsize")
            .help("The bin size for the query")
            .required(true))
        .arg(Arg::new("bin-overlap")
            .long("bin-overlap")
            .default_value("0")
            .help("Number of overlapping bins when computing the summary"))
        .arg(Arg::new("bin-summary")
            .long("bin-summary")
            .default_value("mean")
            .help("Bin summary statistic [mean, max, min]"))
        .arg(Arg::new("verbose")
            .short('v')
            .long("verbose")
            .action(clap::ArgAction::SetTrue)
            .help("Be verbose"))
        .get_matches();

    let filename_in = matches.get_one::<String>("input").unwrap();
    let chrom = matches.get_one::<String>("chrom").unwrap();
    let from = matches.get_one::<String>("from").unwrap().parse::<usize>().unwrap_or_else(|_| {
        eprintln!("Invalid start position");
        process::exit(1);
    });
    let to = matches.get_one::<String>("to").unwrap().parse::<usize>().unwrap_or_else(|_| {
        eprintln!("Invalid end position");
        process::exit(1);
    });
    let bin_size = matches.get_one::<String>("binsize").unwrap().parse::<usize>().unwrap_or_else(|_| {
        eprintln!("Invalid bin size");
        process::exit(1);
    });
    let bin_overlap = matches.get_one::<String>("bin-overlap").unwrap().parse::<usize>().unwrap();
    let summary_str = matches.get_one::<String>("bin-summary").unwrap();
    let verbose = matches.get_flag("verbose");

    let summary = bin_summary_statistics_from_string(summary_str).unwrap_or_else(|| {
        eprintln!("Invalid bin summary statistic. Choose from: mean, max, min");
        process::exit(1);
    });

    if let Err(e) = query(filename_in, chrom, from, to, bin_size, bin_overlap, summary, verbose) {
        eprintln!("{}", e);
        process::exit(1);
    }
}
