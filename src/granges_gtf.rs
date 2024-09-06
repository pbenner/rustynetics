
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
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::collections::HashMap;

use flate2::read::GzDecoder;
use regex::Regex;

use crate::granges::GRanges;
use crate::meta::MetaData;
use crate::range::Range;

/* -------------------------------------------------------------------------- */
 
impl GRanges {

    fn read_gtf<R: BufRead>(
        &mut self,
        reader: R,
        opt_names: Vec<String>,
        opt_types: Vec<String>,
        defaults: Vec<Option<String>>,
    ) -> io::Result<()> {
        let mut type_map = HashMap::new();
        let mut gtf_opt = HashMap::new();
        let mut gtf_def = HashMap::new();

        for (i, name) in opt_names.iter().enumerate() {
            type_map.insert(name.clone(), opt_types[i].clone());
            if let Some(def) = defaults.get(i) {
                gtf_def.insert(name.clone(), def.clone());
            }
        }

        for line in reader.lines() {
            let line = line?;
            let fields: Vec<String> = Self::parse_gtf_line(&line);
            if fields.len() < 8 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "file must have at least eight columns"));
            }
            self.seqnames   .push(fields[0].clone());
            self.source     .push(fields[1].clone());
            self.feature    .push(fields[2].clone());
            self.ranges     .push(Range::new(fields[3].parse::<usize>()?, fields[4].parse::<usize>()?));
            self.score      .push(fields[5].parse::<f64>().unwrap_or(0.0));
            self.strand     .push(fields[6].chars().next().unwrap_or('*'));
            self.frame      .push(fields[7].parse::<i32>().unwrap_or(-1));

            let optional_fields = &fields[8..];
            self.parse_optional_fields(optional_fields, &type_map, &mut gtf_opt, &gtf_def, self.seqnames.len())?;
        }

        // Add meta data to GRanges
        self.meta.add("source", self.source.clone());
        self.meta.add("feature", self.feature.clone());
        self.meta.add("score", self.score.iter().map(|s| s.to_string()).collect());
        self.meta.add("frame", self.frame.iter().map(|f| f.to_string()).collect());

        for (name, values) in gtf_opt {
            self.meta.add(&name, MetaData::StringArray(values));
        }

        Ok(())
    }

    fn parse_gtf_line(line: &str) -> Vec<String> {
        let re = Regex::new(r#""[^"]*"|[^ \t;]+"#).unwrap();
        re.find_iter(line)
            .map(|mat| mat.as_str().to_string())
            .collect()
    }

    fn parse_optional_fields(
        &self,
        fields: &[String],
        type_map: &HashMap<String, String>,
        gtf_opt: &mut HashMap<String, Vec<String>>,
        gtf_def: &HashMap<String, Option<String>>,
        length: usize,
    ) -> io::Result<()> {
        for i in (0..fields.len()).step_by(2) {
            let name = &fields[i];
            let value_str = &fields[i + 1];
            if let Some(type_str) = type_map.get(name) {
                let entry = gtf_opt.entry(name.clone()).or_insert_with(Vec::new);
                entry.push(value_str.clone());
            }
        }

        for (name, values) in gtf_opt.iter_mut() {
            if values.len() < length {
                if let Some(def_val) = gtf_def.get(name) {
                    if let Some(def) = def_val {
                        values.push(def.clone());
                    } else {
                        return Err(io::Error::new(io::ErrorKind::InvalidData, format!("optional field `{}` is missing at line `{}` with no default", name, length + 1)));
                    }
                }
            }
        }

        Ok(())
    }

    fn import_gtf<P: AsRef<Path>>(
        &mut self,
        path: P,
        opt_names: Vec<String>,
        opt_types: Vec<String>,
        opt_def: Vec<Option<String>>,
    ) -> io::Result<()> {
        let file = File::open(path)?;
        let reader: Box<dyn BufRead> = if path.as_ref().to_str().unwrap().ends_with(".gz") {
            let decoder = GzDecoder::new(file);
            Box::new(BufReader::new(decoder))
        } else {
            Box::new(BufReader::new(file))
        };
        self.read_gtf(reader, opt_names, opt_types, opt_def)
    }

    fn write_gtf<W: Write>(&self, writer: W) -> io::Result<()> {
        let mut w = io::BufWriter::new(writer);

        for i in 0..self.seqnames.len() {
            writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                self.seqnames[i],
                self.source.get(i).unwrap_or(&".".to_string()),
                self.feature.get(i).unwrap_or(&".".to_string()),
                self.start[i],
                self.end[i],
                self.score.get(i).unwrap_or(&0.0).to_string(),
                self.strand[i],
                self.frame.get(i).unwrap_or(&-1).to_string()
            )?;

            for (name, values) in self.meta.iter() {
                if name == "source" || name == "feature" || name == "score" || name == "frame" {
                    continue;
                }
                write!(w, "\t{} ", name)?;
                write!(w, "\"{}\"", values[i])?;
                write!(w, ";")?;
            }
            writeln!(w)?;
        }

        Ok(())
    }

    fn export_gtf<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let file = File::create(path)?;
        self.write_gtf(file)
    }

}
