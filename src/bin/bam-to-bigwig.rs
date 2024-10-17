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
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::process;
use std::error::Error;

use clap::{Arg, ArgAction, Command};
use plotters::prelude::*;

use rustynetics::bam::bam_import_genome;
use rustynetics::track_coverage::estimate_fraglen;
use rustynetics::track_coverage::FraglenEstimate;
use rustynetics::track_coverage::bam_coverage;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, Default)]
struct Config {
    bw_zoom_levels: Vec<i32>,
    save_fraglen: bool,
    save_cross_corr: bool,
    save_cross_corr_plot: bool,
    verbose: u8,
}

/* -------------------------------------------------------------------------- */

macro_rules! print_stderr {
    ($config:expr, $level:expr, $($arg:tt)*) => {
        if $config.verbose >= $level {
            eprintln!($($arg)*);
        }
    };
}


/* -------------------------------------------------------------------------- */

fn parse_filename(filename: &str) -> (String, i32) {
    let parts: Vec<&str> = filename.split(':').collect();
    if parts.len() == 2 {
        let t = parts[1].parse::<i32>().unwrap_or_else(|err| {
            eprintln!("Error parsing filename: {}", err);
            process::exit(1);
        });
        return (parts[0].to_string(), t);
    } else if parts.len() >= 2 {
        eprintln!("Invalid input file description `{}`", filename);
        process::exit(1);
    }
    (filename.to_string(), -1)
}

/* -------------------------------------------------------------------------- */

fn save_fraglen(config: &Config, filename: &str, fraglen: i32) -> io::Result<()> {
    let basename = filename.trim_end_matches(Path::new(filename).extension().unwrap_or_default().to_str().unwrap());
    let out_filename = format!("{}.fraglen.txt", basename);

    let mut file = File::create(&out_filename)?;
    writeln!(file, "{}", fraglen)?;

    print_stderr!(config, 1, "Wrote fragment length estimate to `{}`\n", out_filename);
    Ok(())
}

/* -------------------------------------------------------------------------- */

fn save_cross_corr(config: &Config, filename: &str, x: &[i32], y: &[f64]) -> io::Result<()> {
    let basename = filename.trim_end_matches(Path::new(filename).extension().unwrap_or_default().to_str().unwrap());
    let out_filename = format!("{}.fraglen.table", basename);

    let mut file = File::create(&out_filename)?;
    for (xi, yi) in x.iter().zip(y.iter()) {
        writeln!(file, "{} {}", xi, yi)?;
    }

    print_stderr!(config, 1, "Wrote cross-correlation table to `{}`\n", out_filename);
    Ok(())
}

/* -------------------------------------------------------------------------- */

fn save_cross_corr_plot(
    filename: &str,
    cross_corr: &[f64],
    fraglen_est: FraglenEstimate,
    outliers: &[usize],
    _scale_x: bool, // You can implement scaling if needed
) -> Result<(), Box<dyn Error>> {
    // Create a 600x400 image for the plot
    let root = BitMapBackend::new(filename, (600, 400)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_y = cross_corr.iter().cloned().fold(f64::MIN, f64::max);
    let min_y = cross_corr.iter().cloned().fold(f64::MAX, f64::min);

    let mut chart = ChartBuilder::on(&root)
        .caption("Cross-correlation", ("sans-serif", 20))
        .x_label_area_size(30)
        .y_label_area_size(30)
        .margin(5)
        .build_cartesian_2d(0..cross_corr.len(), min_y..max_y)?;

    chart.configure_mesh()
        .x_desc("Fragment length")
        .y_desc("Cross-correlation")
        .draw()?;

    // Plot the cross-correlation line
    chart.draw_series(LineSeries::new(
        cross_corr.iter().enumerate().map(|(i, &y)| (i, y)),
        &RED,
    ))?;

    // Mark the estimated fragment length with a vertical line
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(fraglen_est.fraglen as usize, min_y as f64), (fraglen_est.fraglen as usize, max_y as f64)],
        BLUE.stroke_width(1),
    )))?;

    // Mark outliers with blue lines
    for &outlier in outliers {
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(outlier, min_y), (outlier, max_y)],
            GREEN.stroke_width(1),
        )))?;
    }

    // Ensure the plot is saved
    root.present()?;

    Ok(())
}

/* -------------------------------------------------------------------------- */

fn import_fraglen(config: &Config, filename: &str) -> i32 {
    // Try reading the fragment length from file
    let basename = Path::new(filename).with_extension(""); // Remove file extension
    let fraglen_filename = format!("{}.fraglen.txt", basename.display());

    match File::open(&fraglen_filename) {
        Ok(file) => {
            print_stderr!(config, 1, "Reading fragment length from `{}`... ", fraglen_filename);
            let reader = io::BufReader::new(file);
            let mut lines = reader.lines();

            if let Some(Ok(line)) = lines.next() {
                if let Ok(fraglen) = line.parse::<i64>() {
                    print_stderr!(config, 1, "done\n");
                    return fraglen as i32;
                }
            }

            print_stderr!(config, 1, "failed\n");
            -1
        }
        Err(_) => -1,
    }
}

/* -------------------------------------------------------------------------- */

fn main() {
    let mut app = Command::new("bam-to-bigwig")
        .version("1.0")
        .author("Philipp Benner [https://github.com/pbenner]")
        .about("Convert BAM to BigWig files")
        // bigWig options
        .arg(Arg::new("bigwig-zoom-levels")
            .long("bigwig-zoom-levels")
            .num_args(1)
            .help("Comma separated list of BigWig zoom levels"))
        // read options
        .arg(Arg::new("shift-reads")
            .long("shift-reads")
            .num_args(1)
            .help("Shift reads on the positive strand by `x` bps and those on the negative strand by `y` bps [format: x,y]"))
        .arg(Arg::new("paired-as-single-end")
            .long("paired-as-single-end")
            .action(ArgAction::SetTrue)
            .help("Treat paired reads as single end reads"))
        .arg(Arg::new("paired-end-strand-specific")
            .long("paired-end-strand-specific")
            .action(ArgAction::SetTrue)
            .help("Strand specific paired-end sequencing"))
        // options for filtering reads
        .arg(Arg::new("filter-strand")
            .long("filter-strand")
            .num_args(1)
            .help("Use reads on either the forward `+` or reverse `-` strand"))
        .arg(Arg::new("filter-read-lengths")
            .long("filter-read-lengths")
            .num_args(1)
            .help("Feasible range of read lengths [format: min:max]"))
        .arg(Arg::new("filter-mapq")
            .long("filter-mapq")
            .num_args(1)
            .default_value("0")
            .help("Filter reads for minimum mapping quality [default: 0]"))
        .arg(Arg::new("filter-duplicates")
            .long("filter-duplicates")
            .action(ArgAction::SetTrue)
            .help("Remove reads marked as duplicates"))
        .arg(Arg::new("filter-paired-end")
            .long("filter-paired-end")
            .action(ArgAction::SetTrue)
            .help("Remove all single end reads"))
        .arg(Arg::new("filter-single-end")
            .long("filter-single-end")
            .action(ArgAction::SetTrue)
            .help("Remove all paired end reads"))
        .arg(Arg::new("filter-chromosomes")
            .long("filter-chromosomes")
            .num_args(1)
            .help("Remove all reads on the given chromosomes [comma separated list]"))
        // track options
        .arg(Arg::new("binning-method")
            .long("binning-method")
            .num_args(1)
            .help("Binning method"))
        .arg(Arg::new("bin-size")
            .long("bin-size")
            .num_args(1)
            .default_value("10")
            .help("Track bin size [default: 10]"))
        .arg(Arg::new("normalize-track")
            .long("normalize-track")
            .num_args(1)
            .help("Normalize track with the specified method"))
        .arg(Arg::new("pseudocounts")
            .long("pseudocounts")
            .num_args(1)
            .help("Pseudocounts added to treatment and control signal [default: `0.0,0.0']"))
        .arg(Arg::new("smoothen-control")
            .long("smoothen-control")
            .action(ArgAction::SetTrue)
            .help("Smoothen control with an adaptive window method"))
        .arg(Arg::new("smoothen-window-sizes")
            .long("smoothen-window-sizes")
            .num_args(1)
            .help("Feasible window sizes for the smoothening method [format: s1,s2,...]"))
        .arg(Arg::new("smoothen-min-counts")
            .long("smoothen-min-counts")
            .num_args(1)
            .help("Minimum number of counts for the smoothening method"))
        .arg(Arg::new("log-scale")
            .long("log-scale")
            .action(ArgAction::SetTrue)
            .help("Log-transform data"))
        // options for estimating and setting fragment lengths
        .arg(Arg::new("fragment-length")
            .long("fragment-length")
            .num_args(1)
            .help("Fragment length for all input files [reads are extended to the given length]"))
        .arg(Arg::new("fragment-length-range")
            .long("fragment-length-range")
            .num_args(1)
            .help("Feasible range of fragment lengths [format from:to]"))
        .arg(Arg::new("fragment-length-bin-size")
            .long("fragment-length-bin-size")
            .num_args(1)
            .default_value("10")
            .help("Bin size used when estimating the fragment length [default: 10]"))
        .arg(Arg::new("estimate-fragment-length")
            .long("estimate-fragment-length")
            .action(ArgAction::SetTrue)
            .help("Use crosscorrelation to estimate the fragment length"))
        .arg(Arg::new("save-fraglen")
            .long("save-fraglen")
            .action(ArgAction::SetTrue)
            .help("Save estimated fragment length in a file"))
        .arg(Arg::new("save-crosscorrelation")
            .long("save-crosscorrelation")
            .action(ArgAction::SetTrue)
            .help("Save crosscorrelation between forward and reverse strands"))
        .arg(Arg::new("save-crosscorrelation-plot")
            .long("save-crosscorrelation-plot")
            .action(ArgAction::SetTrue)
            .help("Save crosscorrelation plot"))
        // generic options
        .arg(Arg::new("verbose")
            .short('v')
            .action(ArgAction::Count)
            .help("Verbose level [-v or -vv]"))
        .arg(Arg::new("files")
            .num_args(1..)
            .required(true)
            .help("<TREATMENT1.bam[:FRAGLEN],TREATMENT2.bam[:FRAGLEN],...> [<CONTROL1.bam[:FRAGLEN],CONTROL2.bam[:FRAGLEN],...>] <RESULT.bw>"));

    let matches = app.clone().get_matches();
    let mut config = Config::default();

    // handle verbosity
    config.verbose = matches.get_count("verbose") as u8;

    // handle file arguments
    let files: Vec<_> = matches.get_many::<String>("files").unwrap().collect();
    if files.len() != 2 && files.len() != 3 {
        eprintln!("Error: Invalid number of arguments.");
        eprintln!("{}", app.render_usage());
        process::exit(1);
    }

    print_stderr!(config, 1, "Verbose mode enabled with level: {}", config.verbose);

}
