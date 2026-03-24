use std::cmp::Ordering;
use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::granges::GRanges;
use rustynetics::granges_table::OptionPrintScientific;
use rustynetics::meta::MetaData;
use rustynetics::track::Track;

mod common;

fn parse_input_spec(spec: &str, default_threshold: f64) -> (String, f64) {
    let parts: Vec<_> = spec.split(':').collect();
    if parts.len() == 1 {
        return (parts[0].to_string(), default_threshold);
    }
    if parts.len() != 2 {
        eprintln!("invalid input specification `{spec}`");
        process::exit(1);
    }
    let threshold = parts[1].parse().unwrap_or_else(|error| {
        eprintln!("invalid threshold in `{spec}`: {error}");
        process::exit(1);
    });
    (parts[0].to_string(), threshold)
}

fn get_joint_peaks(
    tracks: &[rustynetics::track_simple::SimpleTrack],
    thresholds: &[f64],
) -> GRanges {
    if tracks.is_empty() {
        return GRanges::default();
    }

    let bin_size = tracks[0].get_bin_size();
    let mut seqnames = Vec::new();
    let mut from = Vec::new();
    let mut to = Vec::new();
    let mut counts = Vec::new();

    for seqname in tracks[0].get_seq_names() {
        let sequences: Vec<_> = tracks
            .iter()
            .map(|track| {
                track.get_sequence(&seqname).unwrap_or_else(|error| {
                    eprintln!("reading sequence `{seqname}` failed: {error}");
                    process::exit(1);
                })
            })
            .collect();

        let n_bins = sequences[0].n_bins();
        for sequence in &sequences[1..] {
            if sequence.n_bins() != n_bins || sequence.get_bin_size() != bin_size {
                eprintln!("track sequences are not compatible for `{seqname}`");
                process::exit(1);
            }
        }

        let mut i = 0usize;
        while i < n_bins {
            let is_positive = |idx: usize| -> bool {
                sequences
                    .iter()
                    .zip(thresholds.iter().copied())
                    .all(|(sequence, threshold)| {
                        let value = sequence.at_bin(idx);
                        !value.is_nan() && value > threshold
                    })
            };

            if !is_positive(i) {
                i += 1;
                continue;
            }

            let start = i;
            let mut max_idx = i;
            let mut max_sum = sequences
                .iter()
                .map(|sequence| sequence.at_bin(i))
                .sum::<f64>();
            i += 1;
            while i < n_bins && is_positive(i) {
                let sum = sequences
                    .iter()
                    .map(|sequence| sequence.at_bin(i))
                    .sum::<f64>();
                if sum > max_sum {
                    max_sum = sum;
                    max_idx = i;
                }
                i += 1;
            }

            seqnames.push(seqname.clone());
            from.push(start * bin_size);
            to.push(i * bin_size);
            counts.push(
                sequences
                    .iter()
                    .map(|sequence| sequence.at_bin(max_idx))
                    .collect::<Vec<_>>(),
            );
        }
    }

    let mut peaks = GRanges::new(seqnames, from, to, Vec::new());
    peaks
        .meta
        .add("counts", MetaData::FloatMatrix(counts))
        .unwrap_or_else(|error| {
            eprintln!("building peak table failed: {error}");
            process::exit(1);
        });

    let sums: Vec<f64> = peaks
        .meta
        .get_column("counts")
        .and_then(|column| match column {
            MetaData::FloatMatrix(values) => Some(
                values
                    .iter()
                    .map(|row| row.iter().sum::<f64>())
                    .collect::<Vec<_>>(),
            ),
            _ => None,
        })
        .unwrap();
    let mut indices: Vec<_> = (0..peaks.num_rows()).collect();
    indices.sort_by(|a, b| sums[*b].partial_cmp(&sums[*a]).unwrap_or(Ordering::Equal));
    peaks.subset(&indices)
}

fn import_regions_table(path: &str, extra_meta: &[String]) -> GRanges {
    let mut reader = common::open_reader(Some(path)).unwrap_or_else(|error| {
        eprintln!("opening regions failed: {error}");
        process::exit(1);
    });
    let mut granges = GRanges::default();
    let mut names = vec!["name"];
    let mut types = vec!["String"];
    for column in extra_meta {
        names.push(column.as_str());
        types.push("String");
    }
    granges
        .read_table(&mut reader, &names, &types)
        .unwrap_or_else(|error| {
            eprintln!("reading regions table failed: {error}");
            process::exit(1);
        });
    granges
}

fn find_nearest_regions(
    mut peaks: GRanges,
    regions_path: Option<&str>,
    k_nearest: usize,
    region_meta: &[String],
) -> GRanges {
    let Some(regions_path) = regions_path else {
        return peaks;
    };

    let regions = if regions_path.ends_with(".bed") || regions_path.ends_with(".bed.gz") {
        let mut reader = common::open_reader(Some(regions_path)).unwrap_or_else(|error| {
            eprintln!("opening regions BED failed: {error}");
            process::exit(1);
        });
        let mut regions = GRanges::default();
        regions.read_bed6(&mut reader).unwrap_or_else(|error| {
            eprintln!("reading regions BED failed: {error}");
            process::exit(1);
        });
        regions
    } else {
        import_regions_table(regions_path, region_meta)
    };

    let names = regions.meta.get_column_str("name").unwrap_or_else(|| {
        eprintln!("regions file has no `name` column");
        process::exit(1);
    });
    let hits = GRanges::find_nearest(&peaks, &regions, k_nearest);

    let mut nearest_names = vec![Vec::new(); peaks.num_rows()];
    let mut distances = vec![Vec::new(); peaks.num_rows()];
    let mut extra_values: Vec<Vec<Vec<String>>> = region_meta
        .iter()
        .map(|_| vec![Vec::new(); peaks.num_rows()])
        .collect();

    let region_columns: Vec<_> = region_meta
        .iter()
        .map(|name| {
            regions.meta.get_column_str(name).unwrap_or_else(|| {
                eprintln!("regions file has no `{name}` column");
                process::exit(1);
            })
        })
        .collect();

    for i in 0..hits.query_hits.len() {
        let q = hits.query_hits[i];
        let s = hits.subject_hits[i];
        nearest_names[q].push(names[s].clone());
        distances[q].push(hits.distances[i]);
        for (j, column) in region_columns.iter().enumerate() {
            extra_values[j][q].push(column[s].clone());
        }
    }

    peaks
        .meta
        .add("names", MetaData::StringMatrix(nearest_names))
        .unwrap();
    peaks
        .meta
        .add("distances", MetaData::IntMatrix(distances))
        .unwrap();
    for (name, values) in region_meta.iter().zip(extra_values.into_iter()) {
        peaks
            .meta
            .add(name, MetaData::StringMatrix(values))
            .unwrap();
    }

    peaks
}

fn main() {
    let matches = Command::new("bigwig-positive")
        .about("Call joint positive regions across one or more BigWig tracks")
        .arg(Arg::new("bin-size").long("bin-size").default_value("0"))
        .arg(
            Arg::new("bin-summary")
                .long("bin-summary")
                .default_value("mean")
                .value_parser([
                    "mean",
                    "max",
                    "min",
                    "discrete mean",
                    "discrete max",
                    "discrete min",
                    "variance",
                ]),
        )
        .arg(
            Arg::new("bin-overlap")
                .long("bin-overlap")
                .default_value("0"),
        )
        .arg(Arg::new("exclude").long("exclude").value_name("BED"))
        .arg(
            Arg::new("regions")
                .long("regions")
                .value_name("BED_OR_TABLE"),
        )
        .arg(Arg::new("k-nearest").long("k-nearest").default_value("1"))
        .arg(
            Arg::new("regions-meta")
                .long("regions-meta")
                .default_value(""),
        )
        .arg(Arg::new("threshold").long("threshold").default_value("1.0"))
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("output").required(true).index(1))
        .arg(Arg::new("inputs").required(true).index(2).num_args(1..))
        .get_matches();

    let output = matches.get_one::<String>("output").unwrap();
    let bin_size: usize = matches
        .get_one::<String>("bin-size")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid bin size: {error}");
            process::exit(1);
        });
    let bin_overlap: usize = matches
        .get_one::<String>("bin-overlap")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid bin overlap: {error}");
            process::exit(1);
        });
    let default_threshold: f64 = matches
        .get_one::<String>("threshold")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid default threshold: {error}");
            process::exit(1);
        });
    let k_nearest: usize = matches
        .get_one::<String>("k-nearest")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid k-nearest value: {error}");
            process::exit(1);
        });
    let region_meta: Vec<String> = matches
        .get_one::<String>("regions-meta")
        .unwrap()
        .split(',')
        .filter(|value| !value.is_empty())
        .map(str::to_string)
        .collect();
    let summary = common::parse_bin_summary(matches.get_one::<String>("bin-summary").unwrap())
        .unwrap_or_else(|error| {
            eprintln!("{error}");
            process::exit(1);
        });

    let inputs: Vec<(String, f64)> = matches
        .get_many::<String>("inputs")
        .unwrap()
        .map(|spec| parse_input_spec(spec, default_threshold))
        .collect();

    let tracks: Vec<_> = inputs
        .iter()
        .map(|(path, _)| {
            common::import_simple_track(path, "", summary, bin_size, bin_overlap, f64::NAN)
        })
        .collect::<Result<_, _>>()
        .unwrap_or_else(|error| {
            eprintln!("importing BigWig tracks failed: {error}");
            process::exit(1);
        });
    let thresholds: Vec<f64> = inputs.iter().map(|(_, threshold)| *threshold).collect();

    let mut peaks = get_joint_peaks(&tracks, &thresholds);
    if let Some(exclude_path) = matches.get_one::<String>("exclude") {
        let mut reader = common::open_reader(Some(exclude_path)).unwrap_or_else(|error| {
            eprintln!("opening exclude BED failed: {error}");
            process::exit(1);
        });
        let mut exclude = GRanges::default();
        exclude.read_bed3(&mut reader).unwrap_or_else(|error| {
            eprintln!("reading exclude BED failed: {error}");
            process::exit(1);
        });
        peaks = peaks.remove_overlaps_with(&exclude);
    }

    peaks = find_nearest_regions(
        peaks,
        matches.get_one::<String>("regions").map(String::as_str),
        k_nearest,
        &region_meta,
    );

    let mut writer = common::open_writer(Some(output)).unwrap_or_else(|error| {
        eprintln!("opening output failed: {error}");
        process::exit(1);
    });
    if let Err(error) = peaks.write_table(&mut writer, &[&OptionPrintScientific(true)]) {
        eprintln!("writing result table failed: {error}");
        process::exit(1);
    }
}
