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

use plotters::prelude::*;

use rustynetics::bam::bam_import_genome;
use rustynetics::track_coverage::estimate_fraglen;
use rustynetics::track_coverage::FraglenEstimate;
use rustynetics::track_coverage::bam_coverage;

/* -------------------------------------------------------------------------- */

struct Config {
    bw_zoom_levels: Vec<i32>,
    save_fraglen: bool,
    save_cross_corr: bool,
    save_cross_corr_plot: bool,
    verbose: usize,
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

}