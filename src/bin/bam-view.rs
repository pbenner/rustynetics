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

struct Config {
    print_read_name: bool,
    print_cigar    : bool,
    print_sequence : bool,
    print_auxiliary: bool,
}

/* -------------------------------------------------------------------------- */

fn bam_view(config: Config, filename_in: &str) -> Result<(), Box<dyn Error>> {

    let file   = File::open(filename_in)?;
    let reader = BufReader::new(file);
    let star   = "*".to_string();
    let mut pos_str;

    // Options for the BAM reader
    let options = BamReaderOptions {
        read_name     : config.print_read_name,
        read_cigar    : config.print_cigar,
        read_sequence : config.print_sequence,
        read_auxiliary: config.print_auxiliary,
        read_qual     : false,
    };

    let mut bam_reader = BamReader::new(reader, Some(options))?;

    let genome = bam_reader.get_genome().clone();

    // Print header
    print!("{:>10} {:>15} {:>17} {:>4}", "Seqname", "Position", "Flag", "MapQ");

    if options.read_cigar {
        print!(" {:>20}", "Cigar");
    }
    if options.read_name {
        print!(" {:>40}", "ReadName");
    }
    if options.read_sequence {
        print!(" {}", "Sequence");
    }
    if options.read_auxiliary {
        print!(" {}", "Auxiliary");
    }
    println!();

    // Iterate over BAM records
    for result in bam_reader.read_single_end() {

        let block = result?.block;

        let seqname = if block.ref_id < 0 {
            &star
        } else {
            &genome.seqnames[block.ref_id as usize]
        };
        let position = if block.position < 0 {
            &star
        } else {
            pos_str = block.position.to_string();
            &pos_str
        };

        print!("{:>10} {:>15} {:>5}:{:011b} {:>4}",
            seqname,
            position,
            block.flag.0,
            block.flag.0,
            block.mapq
        );

        if options.read_cigar {
            if block.cigar.0.len() > 0 {
                print!(" {:>20}", block.cigar.to_string());
            } else {
                print!(" {:>20}", "-");
            }
        }
        if options.read_name {
            if block.read_name.len() > 0 {
                print!(" {:>40}", block.read_name);
            } else {
                print!(" {:>40}", "-");
            }
        }
        if options.read_sequence {
            if block.seq.0.len() > 0 {
                print!(" {}", block.seq);
            } else {
                print!(" -");
            }
        }
        if options.read_auxiliary {
            if block.auxiliary.is_empty() {
                print!("-");
            } else {
                for aux in block.auxiliary.iter() {
                    print!(" {}", aux.to_string());
                }
            }
        }
        println!();
    }

    Ok(())
}

/* -------------------------------------------------------------------------- */

fn main() {
    let matches = Command::new("BAM Viewer")
        .version("1.0")
        .author("Philipp Benner [https://github.com/pbenner]")
        .about("Prints various fields from a BAM file")
        .arg(
            Arg::new("no-read-name")
                .long("no-read-name")
                .action(clap::ArgAction::SetTrue)
                .help("Do not print read names"),
        )
        .arg(
            Arg::new("no-cigar")
                .long("no-cigar")
                .action(clap::ArgAction::SetTrue)
                .help("Do not print cigar strings"),
        )
        .arg(
            Arg::new("sequence")
                .long("sequence")
                .action(clap::ArgAction::SetTrue)
                .help("Print read sequences"),
        )
        .arg(
            Arg::new("auxiliary")
                .long("auxiliary")
                .action(clap::ArgAction::SetTrue)
                .help("Print auxiliary information"),
        )
        .arg(
            Arg::new("input")
                .help("The input BAM file")
                .required(true)
                .index(1),
        )
        .get_matches();

    let config = Config {
        print_read_name: !matches.get_flag("no-read-name"),
        print_cigar    : !matches.get_flag("no-cigar"),
        print_sequence :  matches.get_flag("sequence"),
        print_auxiliary:  matches.get_flag("auxiliary"),
    };

    let filename_in = matches.get_one::<String>("input").unwrap();

    if let Err(e) = bam_view(config, filename_in) {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}
