use std::collections::HashMap;
use std::ffi::CString;
use std::os::raw::c_char;
use std::process;

use clap::{Arg, ArgAction, Command};
use libloading::{Library, Symbol};

use rustynetics::track::Track;
use rustynetics::track_bigwig::LazyTrackFile;
use rustynetics::track_generic::{GenericMutableTrack, GenericTrack};
use rustynetics::track_simple::SimpleTrack;

mod common;

type MapFn = unsafe extern "C" fn(*const c_char, usize, *const f64, usize) -> f64;

fn load_tracks(
    inputs: &[String],
    summary_name: &str,
    bin_size: usize,
    bin_overlap: usize,
    verbose: u8,
) -> Vec<LazyTrackFile> {
    let summary = common::parse_bin_summary(summary_name).unwrap_or_else(|error| {
        eprintln!("{error}");
        process::exit(1);
    });

    inputs
        .iter()
        .map(|path| {
            if verbose > 0 {
                eprintln!("Opening BigWig `{path}`...");
            }
            common::import_lazy_track(path, path, summary, bin_size, bin_overlap, f64::NAN)
                .unwrap_or_else(|error| {
                    eprintln!("opening `{path}` failed: {error}");
                    process::exit(1);
                })
        })
        .collect()
}

fn build_output_track(tracks: &[LazyTrackFile], requested_bin_size: usize) -> SimpleTrack {
    let bin_size = if requested_bin_size == 0 {
        tracks.first().map(Track::get_bin_size).unwrap_or(0)
    } else {
        requested_bin_size
    };

    if bin_size == 0 {
        eprintln!("could not determine track bin size");
        process::exit(1);
    }

    SimpleTrack::alloc(
        String::new(),
        tracks[0].get_genome().clone(),
        f64::NAN,
        bin_size,
    )
}

fn build_seqname_cache(track: &dyn Track) -> HashMap<String, CString> {
    track
        .get_seq_names()
        .into_iter()
        .map(|seqname| {
            let c_string = CString::new(seqname.clone()).unwrap_or_else(|error| {
                eprintln!("invalid sequence name `{seqname}` for plugin interface: {error}");
                process::exit(1);
            });
            (seqname, c_string)
        })
        .collect()
}

fn write_output(track: &SimpleTrack, output: &str, verbose: u8) {
    if verbose > 0 {
        eprintln!("Writing BigWig `{output}`...");
    }
    GenericTrack::wrap(track)
        .export_bigwig(output, vec![])
        .unwrap_or_else(|error| {
            eprintln!("writing `{output}` failed: {error}");
            process::exit(1);
        });
}

fn main() {
    let matches = Command::new("bigwig-map")
        .about("Apply a shared-library mapping function to one or more BigWig tracks")
        .after_help(
            "The shared library must export a symbol with the signature:\n\
             `extern \"C\" fn(*const c_char, usize, *const f64, usize) -> f64`\n\
             The arguments are `(seqname, genomic_position, values_ptr, values_len)`.\n\
             Use `--symbol` to select a symbol name other than `bigwig_map`.",
        )
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
        .arg(
            Arg::new("symbol")
                .long("symbol")
                .default_value("bigwig_map"),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("plugin").required(true).index(1))
        .arg(Arg::new("output").required(true).index(2))
        .arg(Arg::new("inputs").required(true).num_args(1..).index(3))
        .get_matches();

    let plugin_path = matches.get_one::<String>("plugin").unwrap();
    let output_path = matches.get_one::<String>("output").unwrap();
    let input_paths: Vec<String> = matches
        .get_many::<String>("inputs")
        .unwrap()
        .map(|s| s.to_string())
        .collect();
    let summary_name = matches.get_one::<String>("bin-summary").unwrap();
    let symbol_name = matches.get_one::<String>("symbol").unwrap();
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
    let verbose = matches.get_count("verbose");

    let tracks = load_tracks(&input_paths, summary_name, bin_size, bin_overlap, verbose);
    let mut output_track = build_output_track(&tracks, bin_size);
    let seqname_cache = build_seqname_cache(&output_track);
    let track_refs: Vec<&dyn Track> = tracks.iter().map(|track| track as &dyn Track).collect();

    let library = unsafe { Library::new(plugin_path) }.unwrap_or_else(|error| {
        eprintln!("loading plugin `{plugin_path}` failed: {error}");
        process::exit(1);
    });
    let mapper: Symbol<MapFn> =
        unsafe { library.get(symbol_name.as_bytes()) }.unwrap_or_else(|error| {
            eprintln!("loading symbol `{symbol_name}` from `{plugin_path}` failed: {error}");
            process::exit(1);
        });

    if verbose > 0 {
        eprintln!(
            "Mapping {} track(s) with symbol `{}`...",
            track_refs.len(),
            symbol_name
        );
    }

    GenericMutableTrack::wrap(&mut output_track)
        .map_list(&track_refs, |seqname, position, values| {
            let c_seqname = seqname_cache.get(seqname).unwrap_or_else(|| {
                eprintln!("internal error: missing cached sequence name `{seqname}`");
                process::exit(1);
            });
            unsafe { mapper(c_seqname.as_ptr(), position, values.as_ptr(), values.len()) }
        })
        .unwrap_or_else(|error| {
            eprintln!("mapping tracks failed: {error}");
            process::exit(1);
        });

    write_output(&output_track, output_path, verbose);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::CStr;

    #[test]
    fn seqname_cache_contains_c_strings() {
        let track = SimpleTrack::alloc(
            String::new(),
            rustynetics::genome::Genome::new(vec!["chr1".into(), "chr2".into()], vec![10, 10]),
            f64::NAN,
            1,
        );
        let cache = build_seqname_cache(&track);

        let chr1 = cache.get("chr1").unwrap();
        let value = unsafe { CStr::from_ptr(chr1.as_ptr()) }.to_str().unwrap();
        assert_eq!(value, "chr1");
    }
}
