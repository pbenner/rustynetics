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
use std::io::BufReader;
use std::process;
use std::error::Error;

use clap::{Arg, Command};
use rustynetics::bam::{BamReader, BamReaderOptions};

/* -------------------------------------------------------------------------- */

struct Config {
    print_read_name: bool,
    print_cigar: bool,
    print_sequence: bool,
    print_auxiliary: bool,
}

/* -------------------------------------------------------------------------- */

fn bam_view(config: Config, filename_in: &str) -> Result<(), Box<dyn Error>> {
    let file = File::open(filename_in)?;
    let reader = BufReader::new(file);

    // Options for the BAM reader
    let options = BamReaderOptions {
        read_name: config.print_read_name,
        read_cigar: config.print_cigar,
        read_sequence: config.print_sequence,
        read_auxiliary: config.print_auxiliary,
        read_qual: false,
    };

    let mut bam_reader = BamReader::new(reader, Some(options))?;

    let genome = bam_reader.get_genome().clone();

    // Print header
    print!("{:<10} {:<15} {:<17} {:<4}", "Seqname", "Position", "Flag", "MapQ");

    if options.read_cigar {
        print!("{:<20}", "Cigar");
    }
    if options.read_name {
        print!("{:<40}", "ReadName");
    }
    if options.read_sequence {
        print!("{}", "Sequence");
    }
    if options.read_auxiliary {
        print!("{}", "Auxiliary");
    }
    println!();

    // Iterate over BAM records
    for result in bam_reader.read_single_end() {
        let block = result?.block;

        print!("{:<10} {:<15} {:<5}:{:011b} {:<4}", 
            genome.seqnames[block.ref_id as usize],
            block.position,
            block.flag.0,
            block.flag.0,
            block.mapq
        );

        if options.read_cigar {
            if block.cigar.0.len() > 0 {
                print!("{:<20}", block.cigar);
            } else {
                print!("{:<20}", "-");
            }
        }
        if options.read_name {
            if block.read_name.len() > 0 {
                print!("{:<40}", block.read_name);
            } else {
                print!("{:<40}", "-");
            }
        }
        if options.read_sequence {
            if block.seq.0.len() > 0 {
                print!("{}", block.seq);
            } else {
                print!("-");
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
