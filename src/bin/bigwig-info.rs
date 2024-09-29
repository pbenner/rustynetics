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

fn print_info(filename_in: &str, verbose: bool) {
    if verbose {
        eprintln!("Opening bigWig file {}", filename_in);
    }

    // Open the BigWig file
    let reader = BigWigFile::new_reader(filename_in).unwrap_or_else(|err| {
        eprintln!("Error opening file: {}", err);
        process::exit(1);
    });

    println!("BigWig Header:\n{}", reader.header());
}

/* -------------------------------------------------------------------------- */

fn main() {
    let matches = Command::new("BigWig Genome")
        .version("1.0")
        .author("Philipp Benner [https://github.com/pbenner]")
        .about("Print bigWig information")
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

    print_info(filename_in, verbose);
}
