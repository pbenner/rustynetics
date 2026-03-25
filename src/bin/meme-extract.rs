use std::collections::HashMap;
use std::fs;
use std::io;
use std::process;

use clap::{Arg, ArgAction, Command};
use regex::Regex;

use rustynetics::alphabet::{Alphabet, NucleotideAlphabet};
use rustynetics::tf::TFMatrix;

mod common;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum InputFormat {
    Meme,
    Dreme,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum OutputFormat {
    Table,
    Jaspar,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum OutputType {
    Pwm,
    Ppm,
}

#[derive(Clone, Debug, PartialEq)]
struct ParsedMotif {
    e_value: Option<f64>,
    columns: Vec<[f64; 4]>,
}

#[derive(Clone, Debug, PartialEq)]
struct ParsedDocument {
    background: [f64; 4],
    motifs: Vec<ParsedMotif>,
}

#[derive(Clone, Copy, Debug)]
struct Config {
    alpha: f64,
    filter_e_value: f64,
    input_format: InputFormat,
    output_format: OutputFormat,
    output_type: OutputType,
    verbose: u8,
}

fn log_verbose(config: Config, level: u8, message: &str) {
    if config.verbose >= level {
        eprintln!("{message}");
    }
}

fn parse_input_format(value: &str) -> InputFormat {
    match value.to_ascii_lowercase().as_str() {
        "meme" => InputFormat::Meme,
        "dreme" => InputFormat::Dreme,
        _ => {
            eprintln!("invalid input format `{value}`");
            process::exit(1);
        }
    }
}

fn parse_output_format(value: &str) -> OutputFormat {
    match value.to_ascii_lowercase().as_str() {
        "table" => OutputFormat::Table,
        "jaspar" => OutputFormat::Jaspar,
        _ => {
            eprintln!("invalid output format `{value}`");
            process::exit(1);
        }
    }
}

fn parse_output_type(value: &str) -> OutputType {
    match value.to_ascii_lowercase().as_str() {
        "pwm" => OutputType::Pwm,
        "ppm" => OutputType::Ppm,
        _ => {
            eprintln!("invalid output type `{value}`");
            process::exit(1);
        }
    }
}

fn parse_attributes(tag: &str) -> HashMap<String, String> {
    let attr_re = Regex::new(r#"([A-Za-z_:][A-Za-z0-9_.:-]*)\s*=\s*"([^"]*)""#).unwrap();
    attr_re
        .captures_iter(tag)
        .map(|caps| (caps[1].to_string(), xml_unescape(&caps[2])))
        .collect()
}

fn xml_unescape(text: &str) -> String {
    text.replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&apos;", "'")
        .replace("&amp;", "&")
        .trim()
        .to_string()
}

fn extract_single_capture(xml: &str, pattern: &str, what: &str) -> Result<String, String> {
    let re = Regex::new(pattern).unwrap();
    let caps = re.captures(xml).ok_or_else(|| format!("missing {what}"))?;
    Ok(xml_unescape(caps.get(1).unwrap().as_str()))
}

fn parse_background_array(xml: &str, background_tag: &str) -> Result<[f64; 4], String> {
    let alphabet_array_re = Regex::new(background_tag).unwrap();
    let array_xml = alphabet_array_re
        .captures(xml)
        .and_then(|caps| caps.get(1).map(|m| m.as_str().to_string()))
        .ok_or_else(|| "missing background frequencies".to_string())?;

    let value_re = Regex::new(r#"(?s)<value\b([^>]*)>(.*?)</value>"#).unwrap();
    let alphabet = NucleotideAlphabet;
    let mut result = [0.0; 4];

    for caps in value_re.captures_iter(&array_xml) {
        let attrs = parse_attributes(caps.get(1).unwrap().as_str());
        let letter = attrs
            .get("letter_id")
            .ok_or_else(|| "background has invalid letter".to_string())?;
        if letter.len() != 1 {
            return Err("background has invalid letter".to_string());
        }
        let value: f64 = xml_unescape(caps.get(2).unwrap().as_str())
            .parse()
            .map_err(|error| format!("background has invalid value: {error}"))?;
        let code = alphabet
            .code(letter.as_bytes()[0])
            .map_err(|_| "background has invalid letter".to_string())?;
        result[code as usize] = value;
    }

    Ok(result)
}

fn parse_meme_background(xml: &str) -> Result<[f64; 4], String> {
    parse_background_array(
        xml,
        r#"(?s)<background_frequencies>\s*<alphabet_array\b[^>]*>(.*?)</alphabet_array>\s*</background_frequencies>"#,
    )
}

fn parse_meme_probabilities(xml: &str) -> Result<Vec<[f64; 4]>, String> {
    let matrix_re = Regex::new(
        r#"(?s)<probabilities>\s*<alphabet_matrix\b[^>]*>(.*?)</alphabet_matrix>\s*</probabilities>"#,
    )
    .unwrap();
    let matrix_xml = matrix_re
        .captures(xml)
        .and_then(|caps| caps.get(1).map(|m| m.as_str().to_string()))
        .ok_or_else(|| "missing motif probability matrix".to_string())?;

    parse_meme_alphabet_matrix(&matrix_xml)
}

fn parse_meme_alphabet_matrix(xml: &str) -> Result<Vec<[f64; 4]>, String> {
    let array_re = Regex::new(r#"(?s)<alphabet_array\b[^>]*>(.*?)</alphabet_array>"#).unwrap();
    let value_re = Regex::new(r#"(?s)<value\b([^>]*)>(.*?)</value>"#).unwrap();
    let alphabet = NucleotideAlphabet;
    let mut columns = Vec::new();

    for array_caps in array_re.captures_iter(xml) {
        let array_xml = array_caps.get(1).unwrap().as_str();
        let mut column = [0.0; 4];

        for value_caps in value_re.captures_iter(array_xml) {
            let attrs = parse_attributes(value_caps.get(1).unwrap().as_str());
            let letter = attrs
                .get("letter_id")
                .ok_or_else(|| "matrix has invalid letter".to_string())?;
            if letter.len() != 1 {
                return Err("matrix has invalid letter".to_string());
            }
            let value: f64 = xml_unescape(value_caps.get(2).unwrap().as_str())
                .parse()
                .map_err(|error| format!("matrix has invalid value: {error}"))?;
            let code = alphabet
                .code(letter.as_bytes()[0])
                .map_err(|_| "matrix has invalid letter".to_string())?;
            column[code as usize] = value;
        }

        columns.push(column);
    }

    if columns.is_empty() {
        return Err("missing matrix values".to_string());
    }

    Ok(columns)
}

fn parse_meme(xml: &str) -> Result<ParsedDocument, String> {
    let alphabet_name =
        extract_single_capture(xml, r#"(?s)<alphabet\b[^>]*name="([^"]+)""#, "alphabet")?;
    if !alphabet_name.eq_ignore_ascii_case("DNA") {
        return Err(format!("invalid alphabet `{alphabet_name}`"));
    }

    let background = parse_meme_background(xml)?;
    let motif_re = Regex::new(r#"(?s)<motif\b([^>]*)>(.*?)</motif>"#).unwrap();
    let mut motifs = Vec::new();

    for caps in motif_re.captures_iter(xml) {
        let attrs = parse_attributes(caps.get(1).unwrap().as_str());
        let e_value = attrs
            .get("e_value")
            .map(|value| value.parse::<f64>())
            .transpose()
            .map_err(|error| format!("invalid motif e-value: {error}"))?;
        let columns = parse_meme_probabilities(caps.get(2).unwrap().as_str())?;
        motifs.push(ParsedMotif { e_value, columns });
    }

    Ok(ParsedDocument { background, motifs })
}

fn parse_dreme_background(xml: &str) -> Result<[f64; 4], String> {
    let background_re = Regex::new(r#"(?s)<background\b([^>]*)/?>"#).unwrap();
    let attrs = background_re
        .captures(xml)
        .map(|caps| parse_attributes(caps.get(1).unwrap().as_str()))
        .ok_or_else(|| "missing background frequencies".to_string())?;
    let alphabet = NucleotideAlphabet;
    let mut result = [0.0; 4];

    for (key, value) in attrs {
        if key.len() != 1 {
            continue;
        }
        let code = alphabet
            .code(key.as_bytes()[0])
            .map_err(|_| "background has invalid letter".to_string())?;
        result[code as usize] = value
            .parse()
            .map_err(|error| format!("background has invalid value: {error}"))?;
    }

    Ok(result)
}

fn parse_dreme(xml: &str) -> Result<ParsedDocument, String> {
    let alphabet_name =
        extract_single_capture(xml, r#"(?s)<alphabet\b[^>]*name="([^"]+)""#, "alphabet")?;
    if !alphabet_name.eq_ignore_ascii_case("DNA") {
        return Err(format!("invalid alphabet `{alphabet_name}`"));
    }

    let background = parse_dreme_background(xml)?;
    let motif_re = Regex::new(r#"(?s)<motif\b[^>]*>(.*?)</motif>"#).unwrap();
    let pos_re = Regex::new(r#"(?s)<pos\b([^>]*)/?>"#).unwrap();
    let alphabet = NucleotideAlphabet;
    let mut motifs = Vec::new();

    for motif_caps in motif_re.captures_iter(xml) {
        let motif_xml = motif_caps.get(1).unwrap().as_str();
        let mut columns = Vec::new();

        for pos_caps in pos_re.captures_iter(motif_xml) {
            let attrs = parse_attributes(pos_caps.get(1).unwrap().as_str());
            let mut column = [0.0; 4];

            for (key, value) in attrs {
                if key.len() != 1 {
                    continue;
                }
                let code = alphabet
                    .code(key.as_bytes()[0])
                    .map_err(|_| format!("invalid letter `{key}`"))?;
                column[code as usize] = value
                    .parse()
                    .map_err(|error| format!("invalid value for `{key}`: {error}"))?;
            }
            columns.push(column);
        }

        if columns.is_empty() {
            return Err("motif has no positions".to_string());
        }
        motifs.push(ParsedMotif {
            e_value: None,
            columns,
        });
    }

    Ok(ParsedDocument { background, motifs })
}

fn transpose_columns(columns: &[[f64; 4]]) -> Vec<Vec<f64>> {
    let mut rows = vec![vec![0.0; columns.len()]; 4];
    for (j, column) in columns.iter().enumerate() {
        for i in 0..4 {
            rows[i][j] = column[i];
        }
    }
    rows
}

fn motif_to_matrix(
    config: Config,
    motif: &ParsedMotif,
    background: &[f64; 4],
) -> Result<TFMatrix, String> {
    let mut columns = motif.columns.clone();
    let alphabet_len = 4.0;

    match config.output_type {
        OutputType::Ppm => {
            if config.alpha != 0.0 {
                for column in &mut columns {
                    for i in 0..4 {
                        column[i] =
                            (column[i] + config.alpha) / (1.0 + alphabet_len * config.alpha);
                    }
                }
            }
        }
        OutputType::Pwm => {
            for column in &mut columns {
                for i in 0..4 {
                    let normalized = if config.alpha < 0.0 {
                        column[i] - background[i] * config.alpha
                    } else {
                        (column[i] + config.alpha) / (1.0 + alphabet_len * config.alpha)
                    };
                    if normalized <= 0.0 || background[i] <= 0.0 {
                        return Err(
                            "cannot compute log-odds from non-positive probabilities".to_string()
                        );
                    }
                    column[i] = (normalized / background[i]).log2();
                }
            }
        }
    }

    Ok(TFMatrix::new(transpose_columns(&columns)))
}

fn output_extension(format: OutputFormat) -> &'static str {
    match format {
        OutputFormat::Table => "table",
        OutputFormat::Jaspar => "jaspar",
    }
}

fn write_matrix(path: &str, matrix: &TFMatrix, format: OutputFormat) -> io::Result<()> {
    let mut writer = common::open_writer(Some(path))
        .map_err(|error| io::Error::new(io::ErrorKind::Other, error.to_string()))?;
    match format {
        OutputFormat::Table => matrix.write_matrix(&mut writer)?,
        OutputFormat::Jaspar => matrix.write_jaspar(&mut writer)?,
    }
    writer.flush()
}

fn run(config: Config, input: &str, basename: &str) -> Result<(), String> {
    let xml =
        fs::read_to_string(input).map_err(|error| format!("reading `{input}` failed: {error}"))?;
    let document = match config.input_format {
        InputFormat::Meme => parse_meme(&xml)?,
        InputFormat::Dreme => parse_dreme(&xml)?,
    };

    for (i, motif) in document.motifs.iter().enumerate() {
        if let Some(e_value) = motif.e_value {
            if e_value > config.filter_e_value {
                log_verbose(config, 1, &format!("Skipping motif {}...", i + 1));
                continue;
            }
        }

        log_verbose(config, 1, &format!("Parsing motif {}...", i + 1));
        let matrix = motif_to_matrix(config, motif, &document.background)?;
        let path = format!(
            "{basename}-{i:04}.{}",
            output_extension(config.output_format)
        );
        write_matrix(&path, &matrix, config.output_format)
            .map_err(|error| format!("writing `{path}` failed: {error}"))?;
    }

    Ok(())
}

fn main() {
    let matches = Command::new("meme-extract")
        .about("Extract motif matrices from MEME or DREME XML output")
        .arg(
            Arg::new("pseudo-probability")
                .long("pseudo-probability")
                .default_value("-0.00001"),
        )
        .arg(
            Arg::new("filter-e-value")
                .long("filter-e-value")
                .default_value("0.05"),
        )
        .arg(
            Arg::new("input-format")
                .long("input-format")
                .default_value("meme")
                .value_parser(["meme", "dreme"]),
        )
        .arg(
            Arg::new("output-format")
                .long("output-format")
                .default_value("table")
                .value_parser(["table", "jaspar"]),
        )
        .arg(
            Arg::new("output-type")
                .long("output-type")
                .default_value("pwm")
                .value_parser(["pwm", "ppm"]),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::Count),
        )
        .arg(Arg::new("input").required(true).index(1))
        .arg(Arg::new("basename").required(true).index(2))
        .get_matches();

    let alpha: f64 = matches
        .get_one::<String>("pseudo-probability")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid pseudo probability: {error}");
            process::exit(1);
        });
    let filter_e_value: f64 = matches
        .get_one::<String>("filter-e-value")
        .unwrap()
        .parse()
        .unwrap_or_else(|error| {
            eprintln!("invalid filter e-value: {error}");
            process::exit(1);
        });
    if filter_e_value < 0.0 {
        eprintln!("filter e-value must be non-negative");
        process::exit(1);
    }

    let config = Config {
        alpha,
        filter_e_value,
        input_format: parse_input_format(matches.get_one::<String>("input-format").unwrap()),
        output_format: parse_output_format(matches.get_one::<String>("output-format").unwrap()),
        output_type: parse_output_type(matches.get_one::<String>("output-type").unwrap()),
        verbose: matches.get_count("verbose"),
    };

    let input = matches.get_one::<String>("input").unwrap();
    let basename = matches.get_one::<String>("basename").unwrap();

    run(config, input, basename).unwrap_or_else(|error| {
        eprintln!("{error}");
        process::exit(1);
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    const MEME_XML: &str = r#"
<MEME>
  <training_set>
    <alphabet name="DNA" />
  </training_set>
  <model>
    <background_frequencies>
      <alphabet_array>
        <value letter_id="A">0.25</value>
        <value letter_id="C">0.25</value>
        <value letter_id="G">0.25</value>
        <value letter_id="T">0.25</value>
      </alphabet_array>
    </background_frequencies>
  </model>
  <motifs>
    <motif e_value="0.01">
      <probabilities>
        <alphabet_matrix>
          <alphabet_array>
            <value letter_id="A">0.7</value>
            <value letter_id="C">0.1</value>
            <value letter_id="G">0.1</value>
            <value letter_id="T">0.1</value>
          </alphabet_array>
          <alphabet_array>
            <value letter_id="A">0.1</value>
            <value letter_id="C">0.2</value>
            <value letter_id="G">0.3</value>
            <value letter_id="T">0.4</value>
          </alphabet_array>
        </alphabet_matrix>
      </probabilities>
    </motif>
    <motif e_value="0.50">
      <probabilities>
        <alphabet_matrix>
          <alphabet_array>
            <value letter_id="A">0.25</value>
            <value letter_id="C">0.25</value>
            <value letter_id="G">0.25</value>
            <value letter_id="T">0.25</value>
          </alphabet_array>
        </alphabet_matrix>
      </probabilities>
    </motif>
  </motifs>
</MEME>
"#;

    const DREME_XML: &str = r#"
<dreme>
  <model>
    <alphabet name="DNA" />
    <background A="0.30" C="0.20" G="0.20" T="0.30" />
  </model>
  <motifs>
    <motif>
      <pos A="0.8" C="0.1" G="0.05" T="0.05" />
      <pos A="0.1" C="0.2" G="0.3" T="0.4" />
    </motif>
  </motifs>
</dreme>
"#;

    #[test]
    fn parse_meme_document_and_filter() {
        let document = parse_meme(MEME_XML).unwrap();
        assert_eq!(document.background, [0.25, 0.25, 0.25, 0.25]);
        assert_eq!(document.motifs.len(), 2);
        assert_eq!(document.motifs[0].e_value, Some(0.01));

        let config = Config {
            alpha: 0.0,
            filter_e_value: 0.05,
            input_format: InputFormat::Meme,
            output_format: OutputFormat::Table,
            output_type: OutputType::Ppm,
            verbose: 0,
        };
        let matrix = motif_to_matrix(config, &document.motifs[0], &document.background).unwrap();
        assert_eq!(matrix.values[0], vec![0.7, 0.1]);
        assert_eq!(matrix.values[3], vec![0.1, 0.4]);
    }

    #[test]
    fn parse_dreme_document_and_build_pwm() {
        let document = parse_dreme(DREME_XML).unwrap();
        assert_eq!(document.background, [0.30, 0.20, 0.20, 0.30]);
        assert_eq!(document.motifs.len(), 1);

        let config = Config {
            alpha: 0.0,
            filter_e_value: 0.05,
            input_format: InputFormat::Dreme,
            output_format: OutputFormat::Table,
            output_type: OutputType::Pwm,
            verbose: 0,
        };
        let matrix = motif_to_matrix(config, &document.motifs[0], &document.background).unwrap();
        assert!(matrix.values[0][0] > 1.0);
        assert!(matrix.values[3][1] > 0.0);
    }
}
