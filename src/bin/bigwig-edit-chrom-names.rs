use std::fs::{File, OpenOptions};
use std::io::{Seek, SeekFrom, Write};
use std::process;

use byteorder::LittleEndian;
use clap::{Arg, ArgAction, Command};
use regex::Regex;

use rustynetics::bbi::BbiFile;
use rustynetics::genome::Genome;
use rustynetics::track::Track;
use rustynetics::track_simple::SimpleTrack;

mod common;

const BIGWIG_MAGIC: u32 = 0x888f_fc26;

fn decode_chrom_key(key: &[u8]) -> String {
    let end = key.iter().position(|byte| *byte == 0).unwrap_or(key.len());
    String::from_utf8_lossy(&key[..end]).to_string()
}

fn encode_chrom_key(name: &str, key_size: usize) -> Vec<u8> {
    let mut key = vec![0; key_size];
    let copy_len = name.as_bytes().len().min(key_size);
    key[..copy_len].copy_from_slice(&name.as_bytes()[..copy_len]);
    key
}

fn collect_edits(
    input: &str,
    regex: &Regex,
    replacement: &str,
) -> Result<Vec<(String, String, u64, usize)>, String> {
    let mut file =
        File::open(input).map_err(|error| format!("opening `{input}` failed: {error}"))?;
    let mut bwf = BbiFile::default();
    bwf.open::<LittleEndian, _>(&mut file, BIGWIG_MAGIC)
        .map_err(|error| format!("reading `{input}` failed: {error}"))?;

    Ok(bwf
        .chrom_data
        .keys
        .iter()
        .zip(bwf.chrom_data.key_offsets().iter())
        .map(|(key, offset)| {
            let old_name = decode_chrom_key(key);
            let new_name = regex.replace_all(&old_name, replacement).to_string();
            (old_name, new_name, *offset as u64, key.len())
        })
        .collect())
}

fn rewrite_in_place(input: &str, edits: &[(String, String, u64, usize)]) -> Result<(), String> {
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .open(input)
        .map_err(|error| format!("opening `{input}` for update failed: {error}"))?;

    for (_, new_name, offset, key_size) in edits {
        let key = encode_chrom_key(new_name, *key_size);
        file.seek(SeekFrom::Start(*offset))
            .map_err(|error| format!("seeking in `{input}` failed: {error}"))?;
        file.write_all(&key)
            .map_err(|error| format!("writing chromosome names in `{input}` failed: {error}"))?;
    }

    Ok(())
}

fn main() {
    let matches = Command::new("bigwig-edit-chrom-names")
        .about("Edit chromosome names in a BigWig with a regex transform")
        .arg(Arg::new("bin-size").long("bin-size").default_value("0"))
        .arg(Arg::new("output").long("output").value_name("FILE"))
        .arg(
            Arg::new("dry-run")
                .long("dry-run")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("input").required(true).index(1))
        .arg(Arg::new("regex").required(true).index(2))
        .arg(Arg::new("replacement").required(true).index(3))
        .get_matches();

    let input = matches.get_one::<String>("input").unwrap();
    let output = matches.get_one::<String>("output").map(String::as_str);
    let bin_size: usize = matches
        .get_one::<String>("bin-size")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid bin size: {error}");
            process::exit(1);
        });
    let regex = Regex::new(matches.get_one::<String>("regex").unwrap()).unwrap_or_else(|error| {
        eprintln!("invalid regular expression: {error}");
        process::exit(1);
    });
    let replacement = matches.get_one::<String>("replacement").unwrap();
    let dry_run = matches.get_flag("dry-run");
    let verbose = matches.get_flag("verbose");

    let edits = collect_edits(input, &regex, replacement).unwrap_or_else(|error| {
        eprintln!("{error}");
        process::exit(1);
    });

    let track = common::import_simple_track(
        input,
        "",
        common::parse_bin_summary("mean").unwrap(),
        bin_size,
        0,
        f64::NAN,
    );

    if dry_run {
        for (old_name, new_name, _, _) in &edits {
            println!("`{old_name}` -> `{new_name}`");
        }
        return;
    }

    if let Some(output) = output {
        let track = track.unwrap_or_else(|error| {
            eprintln!("importing BigWig failed: {error}");
            process::exit(1);
        });
        let new_seqnames = edits
            .iter()
            .map(|(_, new_name, _, _)| new_name.clone())
            .collect::<Vec<_>>();
        let lengths = track.get_genome().lengths.clone();
        let sequences: Vec<_> = track
            .get_genome()
            .seqnames
            .iter()
            .map(|name| track.get_sequence(name).unwrap().clone_as_vec())
            .collect();
        let genome = Genome::new(new_seqnames, lengths);
        let renamed = SimpleTrack::new(String::new(), sequences, genome, track.get_bin_size())
            .unwrap_or_else(|error| {
                eprintln!("building renamed track failed: {error}");
                process::exit(1);
            });
        if let Err(error) = common::export_simple_track(&renamed, output) {
            eprintln!("writing output failed: {error}");
            process::exit(1);
        }
        return;
    }

    if verbose {
        eprintln!("Editing chromosome names in `{input}`...");
    }
    if let Err(error) = rewrite_in_place(input, &edits) {
        eprintln!("{error}");
        process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    use rustynetics::bigwig::BigWigReader;
    use rustynetics::genome::Genome;
    use rustynetics::track::MutableTrack;
    use rustynetics::track_simple::SimpleTrack;

    fn temp_bigwig_path(name: &str) -> String {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir()
            .join(format!("rustynetics-{name}-{nanos}.bw"))
            .to_string_lossy()
            .to_string()
    }

    #[test]
    fn rewrite_in_place_updates_bigwig_sequence_names() {
        let path = temp_bigwig_path("chrom-rename");
        let genome = Genome::new(vec!["chr1".into(), "chr2".into()], vec![4, 4]);
        let mut track = SimpleTrack::alloc(String::new(), genome, 0.0, 1);
        track.get_sequence_mut("chr1").unwrap().set(0, 1.0);
        track.get_sequence_mut("chr2").unwrap().set(0, 2.0);
        track.export_bigwig(&path, vec![]).unwrap();

        let regex = Regex::new("^chr").unwrap();
        let edits = collect_edits(&path, &regex, "CHR").unwrap();
        rewrite_in_place(&path, &edits).unwrap();

        let reader = BigWigReader::new(File::open(&path).unwrap()).unwrap();
        assert_eq!(
            reader.genome().seqnames,
            vec!["CHR1".to_string(), "CHR2".to_string()]
        );

        let _ = fs::remove_file(path);
    }
}
