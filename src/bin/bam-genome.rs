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
    let verbose = matches.get_flag("verbose");

    if let Err(e) = get_genome(filename_in, verbose) {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}
