use std::collections::HashMap;
use std::ffi::CString;
use std::process;

use clap::{parser::ValueSource, Arg, ArgAction, Command};
use libloading::{Library, Symbol};

use rustynetics::bigwig_map_plugin::{
    BigWigMapInput, BigWigMapLegacyFn, BigWigMapRustFn, BIGWIG_MAP_LEGACY_SYMBOL,
    BIGWIG_MAP_RUST_SYMBOL,
};
use rustynetics::track::Track;
use rustynetics::track_bigwig::LazyTrackFile;
use rustynetics::track_generic::{GenericMutableTrack, GenericTrack};
use rustynetics::track_simple::SimpleTrack;

mod common;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum PluginAbi {
    Auto,
    Legacy,
    Rust,
}

enum Mapper<'a> {
    Legacy(Symbol<'a, BigWigMapLegacyFn>),
    Rust(Symbol<'a, BigWigMapRustFn>),
}

struct LoadedMapper<'a> {
    abi: PluginAbi,
    symbol_name: String,
    mapper: Mapper<'a>,
}

impl Mapper<'_> {
    unsafe fn call(&self, input: &BigWigMapInput) -> f64 {
        match self {
            Mapper::Legacy(mapper) => mapper(
                input.seqname_ptr,
                input.position,
                input.values_ptr,
                input.values_len,
            ),
            Mapper::Rust(mapper) => mapper(input as *const BigWigMapInput),
        }
    }
}

fn parse_plugin_abi(value: &str) -> PluginAbi {
    match value {
        "auto" => PluginAbi::Auto,
        "legacy" => PluginAbi::Legacy,
        "rust" => PluginAbi::Rust,
        _ => {
            eprintln!("invalid plugin ABI `{value}`");
            process::exit(1);
        }
    }
}

fn try_load_legacy_mapper<'a>(
    library: &'a Library,
    symbol_name: &str,
) -> Result<LoadedMapper<'a>, libloading::Error> {
    let mapper = unsafe { library.get::<BigWigMapLegacyFn>(symbol_name.as_bytes()) }?;
    Ok(LoadedMapper {
        abi: PluginAbi::Legacy,
        symbol_name: symbol_name.to_string(),
        mapper: Mapper::Legacy(mapper),
    })
}

fn try_load_rust_mapper<'a>(
    library: &'a Library,
    symbol_name: &str,
) -> Result<LoadedMapper<'a>, libloading::Error> {
    let mapper = unsafe { library.get::<BigWigMapRustFn>(symbol_name.as_bytes()) }?;
    Ok(LoadedMapper {
        abi: PluginAbi::Rust,
        symbol_name: symbol_name.to_string(),
        mapper: Mapper::Rust(mapper),
    })
}

fn load_default_legacy_mapper<'a>(library: &'a Library) -> Result<LoadedMapper<'a>, String> {
    let mut errors = Vec::new();
    for symbol_name in ["F", BIGWIG_MAP_LEGACY_SYMBOL] {
        match try_load_legacy_mapper(library, symbol_name) {
            Ok(mapper) => return Ok(mapper),
            Err(error) => errors.push(format!("`{symbol_name}`: {error}")),
        }
    }
    Err(format!(
        "loading legacy mapper failed: {}",
        errors.join("; ")
    ))
}

fn load_mapper<'a>(
    library: &'a Library,
    abi: PluginAbi,
    symbol_name: Option<&str>,
    symbol_user_provided: bool,
) -> Result<LoadedMapper<'a>, String> {
    match abi {
        PluginAbi::Legacy => match symbol_name {
            Some(symbol_name) => try_load_legacy_mapper(library, symbol_name)
                .map_err(|error| format!("loading legacy symbol `{symbol_name}` failed: {error}")),
            None => load_default_legacy_mapper(library),
        },
        PluginAbi::Rust => {
            let symbol_name = symbol_name.unwrap_or(BIGWIG_MAP_RUST_SYMBOL);
            try_load_rust_mapper(library, symbol_name)
                .map_err(|error| format!("loading Rust symbol `{symbol_name}` failed: {error}"))
        }
        PluginAbi::Auto => {
            if !symbol_user_provided {
                if let Ok(mapper) = try_load_rust_mapper(library, BIGWIG_MAP_RUST_SYMBOL) {
                    return Ok(mapper);
                }
                return load_default_legacy_mapper(library);
            }

            match symbol_name {
                Some(BIGWIG_MAP_RUST_SYMBOL) => {
                    try_load_rust_mapper(library, BIGWIG_MAP_RUST_SYMBOL).map_err(|error| {
                        format!(
                            "loading Rust symbol `{}` failed: {}",
                            BIGWIG_MAP_RUST_SYMBOL, error
                        )
                    })
                }
                Some("F") | Some(BIGWIG_MAP_LEGACY_SYMBOL) => {
                    let symbol_name = symbol_name.unwrap();
                    try_load_legacy_mapper(library, symbol_name).map_err(|error| {
                        format!("loading legacy symbol `{symbol_name}` failed: {error}")
                    })
                }
                Some(symbol_name) => Err(format!(
                    "cannot infer plugin ABI for custom symbol `{symbol_name}`; use `--abi rust` or `--abi legacy`"
                )),
                None => unreachable!(),
            }
        }
    }
}

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
            "Supported plugin ABIs:\n\
             - `--abi rust`: export `rustynetics_bigwig_map` with signature\n\
               `extern \"C\" fn(*const BigWigMapInput) -> f64`\n\
               using `rustynetics::bigwig_map_plugin::BigWigMapInput`.\n\
             - `--abi legacy`: export `F` (or `--symbol`) with signature\n\
               `extern \"C\" fn(*const c_char, usize, *const f64, usize) -> f64`.\n\
             Without `--abi`, the tool first looks for the Rust symbol\n\
             `rustynetics_bigwig_map`, then falls back to legacy `F` and `bigwig_map`.",
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
            Arg::new("abi")
                .long("abi")
                .default_value("auto")
                .value_parser(["auto", "legacy", "rust"]),
        )
        .arg(Arg::new("symbol").long("symbol"))
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
    let abi = parse_plugin_abi(matches.get_one::<String>("abi").unwrap());
    let summary_name = matches.get_one::<String>("bin-summary").unwrap();
    let symbol_name = matches.get_one::<String>("symbol").map(String::as_str);
    let symbol_user_provided = matches.value_source("symbol") == Some(ValueSource::CommandLine);
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
    let mapper =
        load_mapper(&library, abi, symbol_name, symbol_user_provided).unwrap_or_else(|error| {
            eprintln!("{error}");
            process::exit(1);
        });

    if verbose > 0 {
        let abi_name = match mapper.abi {
            PluginAbi::Auto => "auto",
            PluginAbi::Legacy => "legacy",
            PluginAbi::Rust => "rust",
        };
        eprintln!(
            "Loaded `{}` plugin symbol `{}` from `{}`...",
            abi_name, mapper.symbol_name, plugin_path
        );
    }

    if verbose > 0 {
        eprintln!(
            "Mapping {} track(s) with symbol `{}`...",
            track_refs.len(),
            mapper.symbol_name
        );
    }

    GenericMutableTrack::wrap(&mut output_track)
        .map_list(&track_refs, |seqname, position, values| {
            let c_seqname = seqname_cache.get(seqname).unwrap_or_else(|| {
                eprintln!("internal error: missing cached sequence name `{seqname}`");
                process::exit(1);
            });
            let input = BigWigMapInput {
                seqname_ptr: c_seqname.as_ptr(),
                seqname_len: seqname.len(),
                position,
                values_ptr: values.as_ptr(),
                values_len: values.len(),
            };
            unsafe { mapper.mapper.call(&input) }
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
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::process::Command as ProcessCommand;
    use std::time::{SystemTime, UNIX_EPOCH};

    const RUST_PLUGIN_SOURCE: &str = r#"
#[repr(C)]
pub struct BigWigMapInput {
    pub seqname_ptr: *const std::os::raw::c_char,
    pub seqname_len: usize,
    pub position: usize,
    pub values_ptr: *const f64,
    pub values_len: usize,
}

#[no_mangle]
pub unsafe extern "C" fn rustynetics_bigwig_map(input: *const BigWigMapInput) -> f64 {
    let input = &*input;
    let values = std::slice::from_raw_parts(input.values_ptr, input.values_len);
    values.iter().copied().sum::<f64>() + input.position as f64
}
"#;

    const LEGACY_PLUGIN_SOURCE: &str = r#"
#[no_mangle]
pub unsafe extern "C" fn F(
    _seqname_ptr: *const std::os::raw::c_char,
    position: usize,
    values_ptr: *const f64,
    values_len: usize,
) -> f64 {
    let values = std::slice::from_raw_parts(values_ptr, values_len);
    values.iter().copied().product::<f64>() + position as f64
}
"#;

    fn temp_plugin_dir(name: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("rustynetics-{name}-{nanos}"))
    }

    fn compile_plugin(source: &str, name: &str) -> (PathBuf, PathBuf) {
        let dir = temp_plugin_dir(name);
        fs::create_dir_all(&dir).unwrap();
        let source_path = dir.join("plugin.rs");
        let output_path = dir.join(format!("plugin.{}", std::env::consts::DLL_EXTENSION));
        fs::write(&source_path, source).unwrap();

        let output = ProcessCommand::new("rustc")
            .arg("--crate-type")
            .arg("cdylib")
            .arg(&source_path)
            .arg("-O")
            .arg("-o")
            .arg(&output_path)
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "compiling test plugin failed:\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );

        (dir, output_path)
    }

    fn sample_input(seqname: &CString, values: &[f64], position: usize) -> BigWigMapInput {
        BigWigMapInput {
            seqname_ptr: seqname.as_ptr(),
            seqname_len: seqname.as_bytes().len(),
            position,
            values_ptr: values.as_ptr(),
            values_len: values.len(),
        }
    }

    fn cleanup(dir: &Path) {
        let _ = fs::remove_dir_all(dir);
    }

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

    #[test]
    fn auto_loads_rust_plugin_symbol() {
        let (dir, plugin_path) = compile_plugin(RUST_PLUGIN_SOURCE, "bigwig-map-rust");
        let library = unsafe { Library::new(&plugin_path) }.unwrap();
        let mapper = load_mapper(&library, PluginAbi::Auto, None, false).unwrap();
        let seqname = CString::new("chr1").unwrap();
        let values = [1.0, 2.0, 3.0];
        let input = sample_input(&seqname, &values, 5);

        assert_eq!(mapper.abi, PluginAbi::Rust);
        assert_eq!(mapper.symbol_name, BIGWIG_MAP_RUST_SYMBOL);
        let result = unsafe { mapper.mapper.call(&input) };
        assert_eq!(result, 11.0);

        drop(library);
        cleanup(&dir);
    }

    #[test]
    fn auto_falls_back_to_legacy_plugin_symbol() {
        let (dir, plugin_path) = compile_plugin(LEGACY_PLUGIN_SOURCE, "bigwig-map-legacy");
        let library = unsafe { Library::new(&plugin_path) }.unwrap();
        let mapper = load_mapper(&library, PluginAbi::Auto, None, false).unwrap();
        let seqname = CString::new("chr2").unwrap();
        let values = [2.0, 4.0];
        let input = sample_input(&seqname, &values, 3);

        assert_eq!(mapper.abi, PluginAbi::Legacy);
        assert_eq!(mapper.symbol_name, "F");
        let result = unsafe { mapper.mapper.call(&input) };
        assert_eq!(result, 11.0);

        drop(library);
        cleanup(&dir);
    }
}
