use std::collections::HashMap;
use std::process;

use clap::{Arg, ArgAction, Command};

use rustynetics::granges::GRanges;
use rustynetics::meta::MetaData;

mod common;

#[derive(Clone, Copy)]
struct GtfRecord {
    from: usize,
    to: usize,
    strand: char,
}

impl GtfRecord {
    fn merge(&mut self, other: GtfRecord) {
        self.from = self.from.min(other.from);
        self.to = self.to.max(other.to);
        if self.strand != other.strand {
            self.strand = '*';
        }
    }
}

fn merge_rows(granges: &GRanges, merge_by: &str) -> Result<GRanges, String> {
    let values = granges
        .meta
        .get_column_str(merge_by)
        .ok_or_else(|| format!("attribute `{merge_by}` missing"))?;

    let mut entries: HashMap<String, HashMap<String, GtfRecord>> = HashMap::new();

    for (i, value) in values.iter().enumerate() {
        if value.is_empty() {
            continue;
        }
        let record = GtfRecord {
            from: granges.ranges[i].from,
            to: granges.ranges[i].to,
            strand: granges.strand[i],
        };
        let per_seq = entries.entry(value.clone()).or_default();
        if let Some(existing) = per_seq.get_mut(&granges.seqnames[i]) {
            existing.merge(record);
        } else {
            per_seq.insert(granges.seqnames[i].clone(), record);
        }
    }

    let mut names = Vec::new();
    let mut seqnames = Vec::new();
    let mut from = Vec::new();
    let mut to = Vec::new();
    let mut strand = Vec::new();

    let mut keys: Vec<_> = entries.into_iter().collect();
    keys.sort_by(|a, b| a.0.cmp(&b.0));

    for (name, per_seq) in keys {
        let mut per_seq: Vec<_> = per_seq.into_iter().collect();
        per_seq.sort_by(|a, b| a.0.cmp(&b.0));
        for (seqname, record) in per_seq {
            names.push(name.clone());
            seqnames.push(seqname);
            from.push(record.from);
            to.push(record.to);
            strand.push(record.strand);
        }
    }

    let mut merged = GRanges::new(seqnames, from, to, strand);
    merged
        .meta
        .add("name", MetaData::StringArray(names))
        .map_err(|error| error.to_string())?;
    merged
        .meta
        .add("score", MetaData::IntArray(vec![0; merged.num_rows()]))
        .map_err(|error| error.to_string())?;
    Ok(merged)
}

fn main() {
    let matches = Command::new("gtf-to-bed")
        .about("Convert GTF records to BED6")
        .arg(Arg::new("input").long("input").value_name("FILE"))
        .arg(Arg::new("output").long("output").value_name("FILE"))
        .arg(
            Arg::new("merge-by")
                .long("merge-by")
                .value_name("ATTRIBUTE"),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::SetTrue),
        )
        .get_matches();

    let input_path = matches.get_one::<String>("input").map(String::as_str);
    let output_path = matches.get_one::<String>("output").map(String::as_str);
    let merge_by = matches
        .get_one::<String>("merge-by")
        .map(String::as_str)
        .unwrap_or("");
    let verbose = matches.get_flag("verbose");

    let opt_names = if merge_by.is_empty() {
        Vec::new()
    } else {
        vec![merge_by]
    };
    let opt_types = if merge_by.is_empty() {
        Vec::new()
    } else {
        vec!["str"]
    };
    let defaults = if merge_by.is_empty() {
        Vec::new()
    } else {
        vec![Some("")]
    };

    let mut reader = common::open_reader(input_path).unwrap_or_else(|error| {
        eprintln!("opening GTF failed: {error}");
        process::exit(1);
    });
    if verbose {
        eprintln!("Reading GTF...");
    }
    let mut granges = GRanges::read_gtf(&mut reader, opt_names, opt_types, defaults)
        .unwrap_or_else(|error| {
            eprintln!("reading GTF failed: {error}");
            process::exit(1);
        });

    granges.meta.rename_meta("score", "gtfScore");

    let mut bed = if merge_by.is_empty() {
        granges.meta.rename_meta("feature", "name");
        granges
    } else {
        merge_rows(&granges, merge_by).unwrap_or_else(|error| {
            eprintln!("merging rows failed: {error}");
            process::exit(1);
        })
    };

    if bed.meta.get_column_int("score").is_none() {
        if let Err(error) = bed
            .meta
            .add("score", MetaData::IntArray(vec![0; bed.num_rows()]))
        {
            eprintln!("adding BED score column failed: {error}");
            process::exit(1);
        }
    }

    let mut writer = common::open_writer(output_path).unwrap_or_else(|error| {
        eprintln!("opening output failed: {error}");
        process::exit(1);
    });
    if let Err(error) = bed.write_bed6(&mut writer) {
        eprintln!("writing BED failed: {error}");
        process::exit(1);
    }
}
