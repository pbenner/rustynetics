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

use rustynetics::bigwig::bigwig_import_genome;

/* -------------------------------------------------------------------------- */

fn get_genome(filename_in: &str, verbose: bool) {
    if verbose {
        eprintln!("Reading genome from file: {}", filename_in);
    }

    // Placeholder for reading BigWig genome; Replace with actual library call
    match bigwig_import_genome(filename_in) {
        Ok(genome) => print!("{}", genome),
        Err(e) => {
            eprintln!("Error reading genome: {}", e);
            process::exit(1);
        }
    }
}

/* -------------------------------------------------------------------------- */

fn main() {
    let matches = Command::new("BigWig Genome")
        .version("1.0")
        .author("Philipp Benner [https://github.com/pbenner]")
        .about("Print bigWig genome entries")
        .arg(Arg::new("input")
            .required(true)
            .index(1)
            .help("Input BigWig file"))
        .arg(Arg::new("verbose")
            .short('v')
            .long("verbose")
            .action(clap::ArgAction::SetTrue)
            .help("Enable verbose output"))
        .get_matches();

    let filename_in = matches.get_one::<String>("input").unwrap();
    let verbose     = matches.get_flag("verbose");

    get_genome(filename_in, verbose);
}
