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

use serde::Serialize;
use serde_json;

use rustynetics::bbi::{BbiHeader, BbiHeaderZoom};
use rustynetics::bbi::{BBI_TYPE_FIXED, BBI_TYPE_VARIABLE, BBI_TYPE_BED_GRAPH};
use rustynetics::bigwig::BigWigFile;

/* -------------------------------------------------------------------------- */

// Wrapper struct for BbiZoomHeader to implement Serialize
#[derive(Serialize)]
pub struct SerializableBbiHeaderZoom {
    pub reduction_level : u32,
    pub reserved        : u32,
    pub data_offset     : u64,
    pub index_offset    : u64,
    pub n_blocks        : u32,
}

/* -------------------------------------------------------------------------- */

impl From<&BbiHeaderZoom> for SerializableBbiHeaderZoom {
    fn from(header: &BbiHeaderZoom) -> Self {
        SerializableBbiHeaderZoom {
            reduction_level : header.reduction_level,
            reserved        : header.reserved,
            data_offset     : header.data_offset,
            index_offset    : header.index_offset,
            n_blocks        : header.n_blocks,
        }
    }
}

/* -------------------------------------------------------------------------- */

// Wrapper struct for BbiHeader to implement Serialize
#[derive(Serialize)]
struct SerializableBbiHeader {
    magic              : u32,
    version            : u16,
    zoom_levels        : u16,
    ct_offset          : u64,
    data_offset        : u64,
    index_offset       : u64,
    field_count        : u16,
    defined_field_count: u16,
    sql_offset         : u64,
    summary_offset     : u64,
    uncompress_buf_size: u32,
    extension_offset   : u64,
    n_bases_covered    : u64,
    min_val            : f64,
    max_val            : f64,
    sum_data           : f64,
    sum_squares        : f64,
    zoom_headers       : Vec<SerializableBbiHeaderZoom>,
    n_blocks           : u64,
}

/* -------------------------------------------------------------------------- */

impl From<&BbiHeader> for SerializableBbiHeader {
    fn from(header: &BbiHeader) -> Self {
        SerializableBbiHeader {
            magic              : header.magic,
            version            : header.version,
            zoom_levels        : header.zoom_levels,
            ct_offset          : header.ct_offset,
            data_offset        : header.data_offset,
            index_offset       : header.index_offset,
            field_count        : header.field_count,
            defined_field_count: header.defined_field_count,
            sql_offset         : header.sql_offset,
            summary_offset     : header.summary_offset,
            uncompress_buf_size: header.uncompress_buf_size,
            extension_offset   : header.extension_offset,
            n_bases_covered    : header.n_bases_covered,
            min_val            : header.min_val,
            max_val            : header.max_val,
            sum_data           : header.sum_data,
            sum_squares        : header.sum_squares,
            zoom_headers       : header.zoom_headers.iter().map(|h| SerializableBbiHeaderZoom::from(h)).collect(),
            n_blocks           : header.n_blocks,
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(Serialize)]
struct TrackInfo {
    seqname   : String,
    track_type: String,
    length    : usize,
    binsize   : Option<usize>,
}

#[derive(Serialize)]
struct FileInfo {
    filename: String,
    header  : SerializableBbiHeader,  // Including header info
    tracks  : Vec<TrackInfo>,
}

/* -------------------------------------------------------------------------- */

fn print_json(filename_in: &str, verbose: bool) {
    if verbose {
        eprintln!("Opening bigWig file {}", filename_in);
    }

    // Open the BigWig file
    let mut reader = BigWigFile::new_reader(filename_in).unwrap_or_else(|err| {
        eprintln!("Error opening file: {}", err);
        process::exit(1);
    });

    // Get header and genome from the reader
    let header = reader.header().clone();
    let genome = reader.genome().clone();

    // Collect tracks information
    let mut tracks: Vec<TrackInfo> = Vec::new();

    for (seqname, length) in genome.iter() {
        for r in reader.query(&seqname, 0, *length, 0) {
            if let Ok(record) = r {
                let track_type = if record.data_type == BBI_TYPE_BED_GRAPH {
                    "BedGraph".to_string()
                } else if record.data_type == BBI_TYPE_FIXED {
                    "Fixed".to_string()
                } else if record.data_type == BBI_TYPE_VARIABLE {
                    "Variable".to_string()
                } else {
                    eprintln!("Invalid track data type");
                    process::exit(1);
                };

                let binsize = if record.data_type == BBI_TYPE_BED_GRAPH {
                    None
                } else {
                    Some((record.data.to - record.data.from) as usize)
                };

                tracks.push(TrackInfo {
                    seqname: seqname.clone(),
                    track_type,
                    length: *length,
                    binsize,
                });
                break; // Exit the loop after processing one record for each seqname
            }
        }
    }

    // Prepare file information, now including the header
    let file_info = FileInfo {
        filename: filename_in.to_string(),
        header  : SerializableBbiHeader::from(&header),
        tracks,
    };

    // Print the file information as JSON
    let json_output = serde_json::to_string_pretty(&file_info).unwrap_or_else(|err| {
        eprintln!("Error serializing to JSON: {}", err);
        process::exit(1);
    });

    println!("{}", json_output);
}

/* -------------------------------------------------------------------------- */

fn print_info(filename_in: &str, verbose: bool) {
    if verbose {
        eprintln!("Opening bigWig file {}", filename_in);
    }

    // Open the BigWig file
    let mut reader = BigWigFile::new_reader(filename_in).unwrap_or_else(|err| {
        eprintln!("Error opening file: {}", err);
        process::exit(1);
    });
    let genome = reader.genome().clone();

    print!("{}", reader.header());

    println!("Tracks:");
    for (seqname, length) in genome.iter() {
        for r in reader.query(&seqname, 0, *length, 0) {
            if let Ok(record) = r {
                if record.data_type == BBI_TYPE_BED_GRAPH {
                    println!("  {}:", seqname);
                    println!("    type   : bedgraph");
                    println!("    length : {}", length);
                    println!("    binsize: *");
                } else {
                    let binsize = (record.data.to - record.data.from) as usize;
                    println!("  {}:", seqname);
                    if record.data_type == BBI_TYPE_FIXED {
                        println!("    type   : fixed");
                    } else
                    if record.data_type == BBI_TYPE_VARIABLE {
                        println!("    type   : variable");
                    } else {
                        eprintln!("Invalid track data type");
                        process::exit(1);                
                    }
                    println!("    length : {}", length);
                    println!("    binsize: {}", binsize);
                }
                break;
            }
        }
    }
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
        .arg(Arg::new("json")
            .short('j')
            .long("json")
            .action(clap::ArgAction::SetTrue)
            .help("Print output in json format"))
        .get_matches();

    let filename_in = matches.get_one::<String>("input").unwrap();
    let verbose     = matches.get_flag("verbose");
    let json        = matches.get_flag("json");

    if json {
        print_json(filename_in, verbose);
    } else {
        print_info(filename_in, verbose);
    }
}
