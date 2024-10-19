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

use rustynetics::bigwig::BigWigParameters;
use rustynetics::track_generic::GenericTrack;
use rustynetics::track_coverage::OptionBamCoverage;
use rustynetics::track_coverage::bam_coverage;

/* -------------------------------------------------------------------------- */

#[derive(Clone, Debug, Default)]
struct Config {
    bw_zoom_levels      : Vec<i32>,
    save_fraglen        : bool,
    save_cross_corr     : bool,
    save_cross_corr_plot: bool,
    verbose             : u8,
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

fn parse_filename(filename: &str) -> (String, Option<usize>) {
    let parts: Vec<&str> = filename.split(':').collect();
    if parts.len() == 2 {
        let t = parts[1].parse::<usize>().unwrap_or_else(|err| {
            eprintln!("Error parsing filename: {}", err);
            process::exit(1);
        });
        return (parts[0].to_string(), Some(t));
    } else if parts.len() >= 2 {
        eprintln!("Invalid input file description `{}`", filename);
        process::exit(1);
    }
    (filename.to_string(), None)
}

/* -------------------------------------------------------------------------- */

fn save_fraglen(config: &Config, filename: &str, fraglen: usize) -> io::Result<()> {
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
    config     : &Config,
    filename   : &str,
    fraglen    : usize,
    x          : &Vec<i32>,
    y          : &Vec<f64>,
) -> Result<(), Box<dyn Error>> {

    // Create a 600x400 image for the plot
    let root = BitMapBackend::new(filename, (600, 400)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_y = y.iter().cloned().fold(f64::MIN, f64::max);
    let min_y = y.iter().cloned().fold(f64::MAX, f64::min);

    let z : Vec<(usize, f64)> = x.iter().map(|&x| x as usize).zip(y.iter().map(|&x| x)).collect();

    let mut chart = ChartBuilder::on(&root)
        .caption("Cross-correlation", ("sans-serif", 20))
        .x_label_area_size(30)
        .y_label_area_size(30)
        .margin(5)
        .build_cartesian_2d(0..y.len(), min_y..max_y)?;

    chart.configure_mesh()
        .x_desc("Fragment length")
        .y_desc("Cross-correlation")
        .draw()?;

    // Plot the cross-correlation line
    chart.draw_series(LineSeries::new(
        z,
        &BLACK,
    ))?;

    // Mark the estimated fragment length with a vertical line
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(fraglen, min_y as f64), (fraglen, max_y as f64)],
        RED.stroke_width(1),
    )))?;

    // Ensure the plot is saved
    root.present()?;

    print_stderr!(config, 1, "Wrote cross-correlation plot to `{}`\n", filename);

    Ok(())
}

/* -------------------------------------------------------------------------- */

fn import_fraglen(config: &Config, filename: &str) -> Option<usize> {
    // Try reading the fragment length from file
    let basename = Path::new(filename).with_extension(""); // Remove file extension
    let fraglen_filename = format!("{}.fraglen.txt", basename.display());

    match File::open(&fraglen_filename) {
        Ok(file) => {
            print_stderr!(config, 1, "Reading fragment length from `{}`... ", fraglen_filename);
            let reader = io::BufReader::new(file);
            let mut lines = reader.lines();

            if let Some(Ok(line)) = lines.next() {
                if let Ok(fraglen) = line.parse::<usize>() {
                    print_stderr!(config, 1, "done\n");
                    return Some(fraglen);
                }
            }

            print_stderr!(config, 1, "failed\n");
            None
        }
        Err(_) => None,
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
        .arg(Arg::new("fraglen-range")
            .long("fraglen-range")
            .num_args(1)
            .help("Feasible range of fragment lengths [format from:to]"))
        .arg(Arg::new("fraglen-bin-size")
            .long("fraglen-bin-size")
            .num_args(1)
            .default_value("10")
            .help("Bin size used when estimating the fragment length [default: 10]"))
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

    let mut options_list: Vec<OptionBamCoverage> = Vec::new();

    if let Some(opt_bin_size) = matches.get_one::<String>("bin-size") {
        let bin_size : usize = opt_bin_size.parse().unwrap_or_else(|_| {
            eprintln!("Invalid bin size");
            process::exit(1);
        });
        if bin_size < 1 {
            eprintln!("{}", app.render_usage());
            process::exit(1);
        } else {
            options_list.push(OptionBamCoverage::BinSize(bin_size));
        }
    }

    if let Some(opt_binning_method) = matches.get_one::<String>("binning-method") {
        match opt_binning_method.as_str() {
            "simple" | "default" | "overlap" | "mean overlap" => {}
            _ => {
                eprintln!("invalid binning method `{}`", opt_binning_method);
                process::exit(1);
            }
        }
        options_list.push(OptionBamCoverage::BinningMethod(opt_binning_method.clone()));
    }

    if let Some(opt_pseudocounts) = matches.get_one::<String>("pseudocounts") {
        let tmp: Vec<&str> = opt_pseudocounts.split(',').collect();
        if tmp.len() != 2 {
            eprintln!("{}", app.render_usage());
            process::exit(1);
        }
        let t1 = tmp[0].parse::<f64>().expect("Invalid pseudocount value");
        let t2 = tmp[1].parse::<f64>().expect("Invalid pseudocount value");
        options_list.push(OptionBamCoverage::Pseudocounts([t1, t2]));
    }

    if let Some(opt_read_length) = matches.get_one::<String>("filter-read-lengths") {
        let tmp: Vec<&str> = opt_read_length.split(':').collect();
        if tmp.len() != 2 {
            eprintln!("{}", app.render_usage());
            process::exit(1);
        }
        let t1 = tmp[0].parse::<usize>().expect("Invalid read length");
        let t2 = tmp[1].parse::<usize>().expect("Invalid read length");
        if t1 > t2 {
            eprintln!("{}", app.render_usage());
            process::exit(1);
        }
        options_list.push(OptionBamCoverage::FilterReadLengths([t1, t2]));
    }

    if let Some(opt_filter_mapq) = matches.get_one::<String>("filter-mapq") {
        let filter_mapq : i64 = opt_filter_mapq.parse().unwrap_or_else(|_| {
            eprintln!("Invalid bin size");
            process::exit(1);
        });
        if filter_mapq < 0 {
            eprintln!("{}", app.render_usage());
            process::exit(1);
        }
        options_list.push(OptionBamCoverage::FilterMapQ(filter_mapq));
    }

    if let Some(opt_filter_strand) = matches.get_one::<String>("filter-strand") {
        match opt_filter_strand.as_str() {
            "+" => options_list.push(OptionBamCoverage::FilterStrand('+')),
            "-" => options_list.push(OptionBamCoverage::FilterStrand('-')),
            _ => {
                eprintln!("{}", app.render_usage());
                process::exit(1);
            }
        }
    }

    if let Some(opt_shift_reads) = matches.get_one::<String>("shift-reads") {
        let tmp: Vec<&str> = opt_shift_reads.split(',').collect();
        if tmp.len() != 2 {
            eprintln!("{}", app.render_usage());
            process::exit(1);
        }
        let t1 = tmp[0].parse::<usize>().expect("Invalid shift value");
        let t2 = tmp[1].parse::<usize>().expect("Invalid shift value");
        options_list.push(OptionBamCoverage::ShiftReads([t1, t2]));
    }

    if matches.get_flag("smoothen-control") {
        options_list.push(OptionBamCoverage::SmoothenControl(true));
    }

    if let Some(opt_smoothen_sizes) = matches.get_one::<String>("smoothen-window-sizes") {
        let tmp: Vec<&str> = opt_smoothen_sizes.split(',').collect();
        let smoothen_sizes: Vec<usize> = tmp.iter()
            .map(|&s| s.parse::<usize>().expect("Invalid smoothen size"))
            .collect();
        options_list.push(OptionBamCoverage::SmoothenSizes(smoothen_sizes));
    }

    if let Some(opt_smoothen_min) = matches.get_one::<String>("smoothen-min-counts") {
        let t = opt_smoothen_min.parse::<f64>().expect("Invalid minimum counts");
        options_list.push(OptionBamCoverage::SmoothenMin(t));
    }

    if let Some(opt_bw_zoom_levels) = matches.get_one::<String>("bigwig-zoom-levels") {
        let tmp: Vec<&str> = opt_bw_zoom_levels.split(',').collect();
        let bw_zoom_levels: Vec<i32> = tmp.iter()
            .map(|&s| s.parse::<i32>().expect("Invalid zoom level"))
            .collect();
        config.bw_zoom_levels = bw_zoom_levels;
    }

    if let Some(opt_normalize_track) = matches.get_one::<String>("normalize-track") {
        match opt_normalize_track.to_lowercase().as_str() {
            "rpkm" | "cpm" => {},
            _ => {
                eprintln!("invalid normalization method `{}`", opt_normalize_track);
                process::exit(1);
            }
        }
        options_list.push(OptionBamCoverage::NormalizeTrack(opt_normalize_track.to_lowercase()));
    }
    
    if let Some(opt_fraglen_range) = matches.get_one::<String>("fraglen-range") {
        let tmp: Vec<&str> = opt_fraglen_range.split(':').collect();
        if tmp.len() != 2 {
            eprintln!("{}", app.render_usage());
            process::exit(1);
        }
        let t1 = tmp[0].parse::<i32>().expect("parsing fragment length range failed");
        let t2 = tmp[1].parse::<i32>().expect("parsing fragment length range failed");
        options_list.push(OptionBamCoverage::FraglenRange((t1, t2)));
    }
    
    if let Some(opt_fraglen_bin_size) = matches.get_one::<String>("fraglen-bin-size") {
        let fraglen_bin_size : usize = opt_fraglen_bin_size.parse().unwrap_or_else(|_| {
            eprintln!("Invalid bin size");
            process::exit(1);
        });
        if fraglen_bin_size > 0 {
            options_list.push(OptionBamCoverage::FraglenBinSize(fraglen_bin_size));
        }
    }
    
    if matches.get_flag("filter-paired-end") && matches.get_flag("estimate-fraglen") {
        eprintln!("cannot estimate fragment length for paired end reads");
        process::exit(1);
    }
    
    if matches.get_flag("filter-paired-end") && matches.get_flag("filter-single-end") {
        eprintln!("cannot filter for paired and single end reads");
        process::exit(1);
    }
    
    if let Some(opt_filter_chroms) = matches.get_one::<String>("filter-chromosomes") {
        options_list.push(OptionBamCoverage::FilterChroms(
            opt_filter_chroms.split(',').map(String::from).collect(),
        ));
    }
    
    options_list.push(OptionBamCoverage::EstimateFraglen(matches.get_flag("estimate-fragment-length")));
    options_list.push(OptionBamCoverage::LogScale(matches.get_flag("log-scale")));
    options_list.push(OptionBamCoverage::PairedAsSingleEnd(matches.get_flag("paired-as-single-end")));
    options_list.push(OptionBamCoverage::PairedEndStrandSpecific(matches.get_flag("paired-end-strand-specific")));
    options_list.push(OptionBamCoverage::FilterDuplicates(matches.get_flag("filter-duplicates")));
    options_list.push(OptionBamCoverage::FilterPairedEnd(matches.get_flag("filter-paired-end")));
    options_list.push(OptionBamCoverage::FilterSingleEnd(matches.get_flag("filter-single-end")));
    
    config.save_fraglen         = matches.get_flag("save-fraglen");
    config.save_cross_corr      = matches.get_flag("save-crosscorrelation");
    config.save_cross_corr_plot = matches.get_flag("save-crosscorrelation-plot");
    
    // Parse arguments
    //////////////////////////////////////////////////////////////////////////////
    
    let mut filenames_treatment: Vec<String> = files[0]
        .split(',')
        .map(String::from)
        .collect();
    
    let mut filenames_control: Vec<String> = Vec::new();
    let filename_track : String;
    
    if files.len() == 3 {
        filenames_control = files[1]
            .split(',')
            .map(String::from)
            .collect();
        filename_track = String::from(files[2]);
    } else {
        filename_track = String::from(files[1]);
    }
    
    let mut fraglen_treatment = vec![None; filenames_treatment.len()];
    let mut fraglen_control   = vec![None; filenames_control  .len()];
    
    // Split filename:fraglen for treatment filenames
    for (filename, fraglen) in filenames_treatment.iter_mut().zip(fraglen_treatment.iter_mut()) {
        let (new_filename, new_fraglen) = parse_filename(filename);
        *filename = new_filename;
        *fraglen  = new_fraglen;
    }

    // Split filename:fraglen for control filenames
    for (filename, fraglen) in filenames_control.iter_mut().zip(fraglen_control.iter_mut()) {
        let (new_filename, new_fraglen) = parse_filename(filename);
        *filename = new_filename;
        *fraglen  = new_fraglen;
    }

    
    if let Some(&opt_fraglen) = matches.get_one::<usize>("fragment-length") {
        if opt_fraglen != 0 {
            for fraglen in &mut fraglen_treatment {
                *fraglen = Some(opt_fraglen);
            }
            for fraglen in &mut fraglen_control {
                *fraglen = Some(opt_fraglen);
            }
        }
    }
    
    // Import fragment length
    //////////////////////////////////////////////////////////////////////////////
    
    for (i, filename) in filenames_treatment.iter().enumerate() {
        if fraglen_treatment[i] == None {
            fraglen_treatment[i] = import_fraglen(&config, filename);
        }
    }
    
    for (i, filename) in filenames_control.iter().enumerate() {
        if fraglen_control[i] == None {
            fraglen_control[i] = import_fraglen(&config, filename);
        }
    }

    let result = bam_coverage(
        &filenames_treatment,
        &filenames_control,
        &fraglen_treatment,
        &fraglen_control,
        options_list
    );

    let fraglen_treatment_estimate;
    let fraglen_control_estimate;

    // Get fraglen estimates
    match result {
        Ok((ref _track, ref fraglen_treatment_estimate_, ref fraglen_control_estimate_)) => {
            fraglen_treatment_estimate = fraglen_treatment_estimate_.clone();
            fraglen_control_estimate   = fraglen_control_estimate_.clone();
        },
        Err(ref err) => {
            fraglen_treatment_estimate = err.treatment_fraglen_estimates.clone();
            fraglen_control_estimate   = err.control_fraglen_estimates.clone();
        }
    }

    // Save fragment length estimates if the option is set
    if matches.get_flag("estimate-fragment-length") {
        for (i, estimate) in fraglen_treatment_estimate.iter().enumerate() {
            let filename = &filenames_treatment[i];

            if config.save_fraglen {
                if let Err(e) = save_fraglen(&config, filename, estimate.fraglen) {
                    eprintln!("{}", e);
                    std::process::exit(1);
                }
            }
            if config.save_cross_corr && estimate.x.len() > 0 {
                if let Err(e) = save_cross_corr(&config, filename, &estimate.x, &estimate.y) {
                    eprintln!("{}", e);
                    std::process::exit(1);
                }
            }
            if config.save_cross_corr_plot && estimate.x.len() > 0 {
                if let Err(e) = save_cross_corr_plot(&config, filename, estimate.fraglen, &estimate.x, &estimate.y) {
                    eprintln!("{}", e);
                    std::process::exit(1);
                }
            }
        }

        for (i, estimate) in fraglen_control_estimate.iter().enumerate() {
            let filename = &filenames_control[i];

            if config.save_fraglen {
                if let Err(e) = save_fraglen(&config, filename, estimate.fraglen) {
                    eprintln!("{}", e);
                    std::process::exit(1);
                }
            }
            if config.save_cross_corr && estimate.x.len() > 0 {
                if let Err(e) = save_cross_corr(&config, filename, &estimate.x, &estimate.y) {
                    eprintln!("{}", e);
                    std::process::exit(1);
                }
            }
            if config.save_cross_corr_plot && estimate.x.len() > 0 {
                if let Err(e) = save_cross_corr_plot(&config, filename, estimate.fraglen, &estimate.x, &estimate.y) {
                    eprintln!("{}", e);
                    std::process::exit(1);
                }
            }
        }

        // Exit on error
        if let Err(err) = result {
            eprintln!("Error: {}", err);
            std::process::exit(1);
        }
        // Process the result
        if let Ok((track, _, _)) = result {

            eprint!("Writing track `{}`... ", filename_track);

            let mut parameters = BigWigParameters::default();
            parameters.reduction_levels = config.bw_zoom_levels.clone();

            if let Err(err) = GenericTrack::wrap(&track).export_bigwig(&filename_track, Some(parameters)) {
                eprintln!("failed");
                eprintln!("{}", err);
                std::process::exit(1);
            } else {
                eprintln!("done");
            }

        }

    }

}
