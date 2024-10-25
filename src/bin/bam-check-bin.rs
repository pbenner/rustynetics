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
use std::io::BufReader;
use std::process;
use std::error::Error;

use clap::{Arg, Command};

use rustynetics::bam::{BamReader, BamReaderOptions};

/* -------------------------------------------------------------------------- */

fn reg2bin(beg: i32, mut end: i32) -> i32 {
    end -= 1;
    if beg >> 14 == end >> 14 {
        return ((1 << 15) - 1) / 7 + (beg >> 14);
    }
    if beg >> 17 == end >> 17 {
        return ((1 << 12) - 1) / 7 + (beg >> 17);
    }
    if beg >> 20 == end >> 20 {
        return ((1 << 9) - 1) / 7 + (beg >> 20);
    }
    if beg >> 23 == end >> 23 {
        return ((1 << 6) - 1) / 7 + (beg >> 23);
    }
    if beg >> 26 == end >> 26 {
        return ((1 << 3) - 1) / 7 + (beg >> 26);
    }
    0
}

/* -------------------------------------------------------------------------- */

fn check_bin(filename_in: &str, verbose: bool) -> Result<(), Box<dyn Error>> {
    let file = File::open(filename_in)?;
    let reader = BufReader::new(file);

    let options = BamReaderOptions {
        read_name     : true,
        read_cigar    : true,
        read_sequence : true,
        read_auxiliary: false,
        read_qual     : false,
    };

    let mut bam_reader = BamReader::new(reader, Some(options))?;
    let mut i = 0;

    for result in bam_reader.read_single_end() {
        let block = result?.block;
        let bin;

        if block.flag.unmapped() {
            let from = block.position as i32;
            let to = from + 1;
            bin = reg2bin(from, to);
        } else {
            let from = block.position as i32;
            let to = from + block.cigar.alignment_length() as i32;
            bin = reg2bin(from, to);
        }

        if bin != block.bin as i32 {
            println!(
                "record `{}` has invalid bin value (expected `{}` but bin value is `{}`)",
                i + 1,
                bin,
                block.bin
            );
            println!(" -> {:?}\n", block);
        } else if verbose {
            println!(
                "record `{}` has correct bin value `{}`",
                i + 1,
                bin
            );
            println!(" -> {:?}\n", block);
        }

        i += 1;
    }

    Ok(())
}

/* -------------------------------------------------------------------------- */

fn main() {
    let matches = Command::new("Check BAM Bin Values")
        .version("1.0")
        .author("Philipp Benner [https://github.com/pbenner]")
        .about("Check the bin values of BAM records")
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

    if let Err(e) = check_bin(filename_in, verbose) {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}
