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
use rustynetics::bam::bam_import_genome;

/* -------------------------------------------------------------------------- */

fn get_genome(filename_in: &str, verbose: bool) -> Result<(), Box<dyn Error>> {
    if verbose {
        println!("Opening BAM file `{}`...", filename_in)
    }
    // Import the genome from the BAM file
    match bam_import_genome(filename_in) {
        Ok(genome) => {
            eprint!("{}", genome);
        }
        Err(err) => {
            eprintln!("Error importing genome: {}", err);
            process::exit(1);
        }
    }
    Ok(())
}

/* -------------------------------------------------------------------------- */

fn main() {
    let matches = Command::new("Get BAM Genome")
        .version("1.0")
        .author("Philipp Benner [https://github.com/pbenner]")
        .about("Extract and print the genome information from a BAM file")
        .arg(
            Arg::new("input")
                .help("The input BAM file")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(clap::ArgAction::SetTrue)
                .help("Enable verbose output")
        )
        .get_matches();

    let filename_in = matches.get_one::<String>("input").expect("Input file is required");
    let verbose     = matches.get_flag("verbose");

    if let Err(e) = get_genome(filename_in, verbose) {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}
