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

    println!("{:?}", result);

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
            .help("Be verbose"))
        .arg(Arg::new("help")
            .short('h')
            .long("help")
            .help("Print help"))
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
