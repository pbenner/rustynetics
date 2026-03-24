use std::cmp::Ordering;
use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::bigwig::BigWigFile;
use rustynetics::granges::GRanges;
use rustynetics::meta::MetaData;

mod common;

fn import_bed6(path: &str) -> GRanges {
    let mut reader = common::open_reader(Some(path)).unwrap_or_else(|error| {
        eprintln!("opening BED failed: {error}");
        process::exit(1);
    });
    let mut granges = GRanges::default();
    granges.read_bed6(&mut reader).unwrap_or_else(|error| {
        eprintln!("reading BED6 failed: {error}");
        process::exit(1);
    });
    granges
}

fn filter_states(granges: &GRanges, states: &[String]) -> GRanges {
    let names = granges.meta.get_column_str("name").unwrap_or_else(|| {
        eprintln!("segmentation BED is missing the `name` column");
        process::exit(1);
    });
    let indices: Vec<_> = names
        .iter()
        .enumerate()
        .filter_map(|(i, name)| states.contains(name).then_some(i))
        .collect();
    granges.subset(&indices)
}

fn region_score(path: &str, seqname: &str, from: usize, to: usize, summary: &str) -> f64 {
    let summary = common::parse_bin_summary(summary).unwrap();
    let mut reader = BigWigFile::new_reader(path).unwrap_or_else(|error| {
        eprintln!("opening BigWig `{path}` failed: {error}");
        process::exit(1);
    });
    let values = reader
        .query_slice(seqname, from, to, summary, to - from, 0, f64::NAN)
        .unwrap_or_else(|error| {
            eprintln!("querying BigWig `{path}` failed: {error}");
            process::exit(1);
        })
        .0;
    values.first().copied().unwrap_or(f64::NAN)
}

fn apply_bigwig_filter(granges: &GRanges, path: &str, threshold: f64, keep_above: bool) -> GRanges {
    let mut keep = Vec::new();
    for i in 0..granges.num_rows() {
        let value = region_score(
            path,
            &granges.seqnames[i],
            granges.ranges[i].from,
            granges.ranges[i].to,
            "mean",
        );
        let passes = if keep_above {
            value.is_nan() || value >= threshold
        } else {
            value.is_nan() || value <= threshold
        };
        if passes {
            keep.push(i);
        }
    }
    granges.subset(&keep)
}

fn overlaps(left: &GRanges, li: usize, right: &GRanges, ri: usize) -> bool {
    left.seqnames[li] == right.seqnames[ri]
        && left.ranges[li].from < right.ranges[ri].to
        && right.ranges[ri].from < left.ranges[li].to
}

fn compute_labels(merged: &GRanges, regions: &[GRanges]) -> Vec<Vec<i64>> {
    let mut labels = vec![Vec::new(); merged.num_rows()];
    for (j, region_set) in regions.iter().enumerate() {
        for i in 0..merged.num_rows() {
            if (0..region_set.num_rows()).any(|k| overlaps(merged, i, region_set, k)) {
                labels[i].push(j as i64);
            }
        }
    }
    labels
}

fn compute_scores(granges: &GRanges, labels: &[Vec<i64>], score_files: &[String]) -> Vec<f64> {
    let mut scores = vec![0.0; granges.num_rows()];
    for i in 0..granges.num_rows() {
        let n_pos = labels[i].len();
        if n_pos == 0 {
            continue;
        }
        if n_pos == score_files.len() {
            scores[i] = f64::NEG_INFINITY;
            continue;
        }
        let mut pos_sum = 0.0;
        let mut neg_sum = 0.0;
        for (j, path) in score_files.iter().enumerate() {
            let value = region_score(
                path,
                &granges.seqnames[i],
                granges.ranges[i].from,
                granges.ranges[i].to,
                "max",
            );
            if labels[i].contains(&(j as i64)) {
                pos_sum += value;
            } else {
                neg_sum += value;
            }
        }
        scores[i] = ((pos_sum + 1.0) / n_pos as f64).ln()
            - ((neg_sum + 1.0) / (score_files.len() - n_pos) as f64).ln();
    }
    scores
}

fn parse_threshold_pairs(value: &str) -> (Vec<String>, Vec<f64>) {
    let mut files = Vec::new();
    let mut thresholds = Vec::new();
    if value.is_empty() {
        return (files, thresholds);
    }
    for item in value.split(',') {
        let parts: Vec<_> = item.split(':').collect();
        if parts.len() != 2 {
            eprintln!("invalid threshold spec `{item}`");
            process::exit(1);
        }
        files.push(parts[0].to_string());
        thresholds.push(parts[1].parse().unwrap_or_else(|error| {
            eprintln!("invalid threshold in `{item}`: {error}");
            process::exit(1);
        }));
    }
    (files, thresholds)
}

fn main() {
    let matches = Command::new("segmentation-differential")
        .about("Merge and score differential chromatin states across segmentations")
        .arg(
            Arg::new("filter-bigwig")
                .long("filter-bigwig")
                .default_value(""),
        )
        .arg(
            Arg::new("filter-out-bigwig")
                .long("filter-out-bigwig")
                .default_value(""),
        )
        .arg(
            Arg::new("filter-max-size")
                .long("filter-max-size")
                .default_value("0"),
        )
        .arg(Arg::new("scores").long("scores").default_value(""))
        .arg(
            Arg::new("region-size")
                .long("region-size")
                .default_value("0"),
        )
        .arg(
            Arg::new("print-statistics")
                .long("print-statistics")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("states").required(true).index(1))
        .arg(
            Arg::new("segmentations")
                .required(true)
                .index(2)
                .num_args(2..),
        )
        .get_matches();

    let states: Vec<String> = matches
        .get_one::<String>("states")
        .unwrap()
        .split(',')
        .map(str::to_string)
        .collect();
    let segmentation_files: Vec<String> = matches
        .get_many::<String>("segmentations")
        .unwrap()
        .map(|s| s.to_string())
        .collect();
    let (filter_files, filter_thr) =
        parse_threshold_pairs(matches.get_one::<String>("filter-bigwig").unwrap());
    let (filter_out_files, filter_out_thr) =
        parse_threshold_pairs(matches.get_one::<String>("filter-out-bigwig").unwrap());
    let score_files: Vec<String> = matches
        .get_one::<String>("scores")
        .unwrap()
        .split(',')
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect();
    let filter_max_size: usize = matches
        .get_one::<String>("filter-max-size")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid filter-max-size: {error}");
            process::exit(1);
        });
    let region_size: usize = matches
        .get_one::<String>("region-size")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid region-size: {error}");
            process::exit(1);
        });

    if !filter_files.is_empty() && filter_files.len() != segmentation_files.len() {
        eprintln!("--filter-bigwig must provide one file per segmentation");
        process::exit(1);
    }
    if !filter_out_files.is_empty() && filter_out_files.len() != segmentation_files.len() {
        eprintln!("--filter-out-bigwig must provide one file per segmentation");
        process::exit(1);
    }
    if !score_files.is_empty() && score_files.len() != segmentation_files.len() {
        eprintln!("--scores must provide one file per segmentation");
        process::exit(1);
    }

    let mut regions = Vec::with_capacity(segmentation_files.len());
    for (i, path) in segmentation_files.iter().enumerate() {
        let mut granges = filter_states(&import_bed6(path), &states);
        if !filter_files.is_empty() {
            granges = apply_bigwig_filter(&granges, &filter_files[i], filter_thr[i], true);
        }
        if !filter_out_files.is_empty() {
            granges = apply_bigwig_filter(&granges, &filter_out_files[i], filter_out_thr[i], false);
        }
        regions.push(granges);
    }

    let region_refs: Vec<_> = regions.iter().collect();
    let mut merged = GRanges::merge(&region_refs);
    let labels = compute_labels(&merged, &regions);

    merged
        .meta
        .add(
            "name",
            MetaData::StringArray(vec![states.join(","); merged.num_rows()]),
        )
        .unwrap();
    merged
        .meta
        .add("labels", MetaData::IntMatrix(labels.clone()))
        .unwrap();

    if !score_files.is_empty() {
        let scores = compute_scores(&merged, &labels, &score_files);
        merged
            .meta
            .add("score", MetaData::FloatArray(scores.clone()))
            .unwrap();
        let mut indices: Vec<_> = (0..merged.num_rows()).collect();
        indices.sort_by(|a, b| {
            scores[*b]
                .partial_cmp(&scores[*a])
                .unwrap_or(Ordering::Equal)
        });
        merged = merged.subset(&indices);
    }

    if filter_max_size > 0 {
        let keep: Vec<_> = (0..merged.num_rows())
            .filter(|&i| merged.ranges[i].to - merged.ranges[i].from <= filter_max_size)
            .collect();
        merged = merged.subset(&keep);
    }
    if region_size > 0 {
        for i in 0..merged.num_rows() {
            let center = (merged.ranges[i].from + merged.ranges[i].to) / 2;
            merged.ranges[i].from = center.saturating_sub(region_size / 2);
            merged.ranges[i].to = merged.ranges[i].from + region_size;
        }
    }

    let mut writer = common::open_writer(None).unwrap_or_else(|error| {
        eprintln!("opening stdout failed: {error}");
        process::exit(1);
    });
    if let Err(error) = merged.write_table(&mut writer, &[]) {
        eprintln!("writing result failed: {error}");
        process::exit(1);
    }

    if matches.get_flag("print-statistics") {
        let labels = merged
            .meta
            .get_column("labels")
            .and_then(|column| match column {
                MetaData::IntMatrix(values) => Some(values),
                _ => None,
            })
            .unwrap();
        eprintln!("Total number of regions: {}", merged.num_rows());
        for tissues in 1..=segmentation_files.len() {
            let count = labels.iter().filter(|row| row.len() == tissues).count();
            let fraction = if merged.num_rows() == 0 {
                0.0
            } else {
                100.0 * count as f64 / merged.num_rows() as f64
            };
            eprintln!(
                "Fraction of regions present in {} tissues: {:7.3}%",
                tissues, fraction
            );
        }
    }
}
