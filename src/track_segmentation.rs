// Copyright (C) 2024 Philipp Benner
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the “Software”), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufWriter, Write};

use flate2::write::GzEncoder;
use flate2::Compression;

use crate::granges::GRanges;
use crate::meta::MetaData;
use crate::track::Track;
use crate::track_generic::{GenericMutableTrack, GenericTrack};
use crate::utility::is_gzip;

/* -------------------------------------------------------------------------- */

impl GenericTrack<'_> {
    pub fn export_segmentation(
        &self,
        bed_filename: &str,
        bed_name: &str,
        bed_description: &str,
        compress: bool,
        state_names: Option<&[String]>,
        rgb_map: Option<&HashMap<String, String>>,
        scores: &[&dyn Track],
    ) -> Result<(), Box<dyn Error>> {
        if self.track.get_bin_size() == 0 {
            return Err(Box::new(io::Error::new(
                io::ErrorKind::InvalidInput,
                "track has invalid bin size",
            )));
        }

        let mut granges = self.granges("state")?;
        let states = granges
            .meta
            .get_column_float("state")
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "state column is missing"))?
            .clone();

        let mut state_ids = Vec::with_capacity(states.len());
        let mut max_state = 0usize;
        for (i, state) in states.iter().copied().enumerate() {
            if state < 0.0 || state.floor() != state {
                return Err(Box::new(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "invalid state `{}` at `{}:{}-{}`",
                        state, granges.seqnames[i], granges.ranges[i].from, granges.ranges[i].to
                    ),
                )));
            }
            let state = state as usize;
            max_state = max_state.max(state);
            state_ids.push(state);
        }

        let state_names = resolve_state_names(state_names, max_state);
        let rgb_map = resolve_rgb_map(rgb_map, &state_names);

        let mut name = Vec::with_capacity(granges.num_rows());
        let mut score = Vec::with_capacity(granges.num_rows());
        let mut thick_start = Vec::with_capacity(granges.num_rows());
        let mut thick_end = Vec::with_capacity(granges.num_rows());
        let mut item_rgb = Vec::with_capacity(granges.num_rows());

        for i in 0..granges.num_rows() {
            let state = state_ids[i];
            if state >= state_names.len() {
                return Err(Box::new(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "insufficient number of state names",
                )));
            }
            let state_name = &state_names[state];
            let color = rgb_map.get(state_name).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("RGB map is missing a color for state `{state_name}`"),
                )
            })?;

            name.push(state_name.clone());
            score.push(0i64);
            thick_start.push(granges.ranges[i].from as i64);
            thick_end.push(granges.ranges[i].to as i64);
            item_rgb.push(color.clone());

            if !scores.is_empty() {
                if state >= scores.len() {
                    return Err(Box::new(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "insufficient number of score tracks",
                    )));
                }
                let values = scores[state].get_slice(&granges.row(i))?;
                let max_score = values
                    .into_iter()
                    .map(|value| (value * 100.0) as i64)
                    .max()
                    .unwrap_or(0);
                score[i] = max_score;
            }
        }

        set_meta_column(&mut granges, "name", MetaData::StringArray(name))?;
        set_meta_column(&mut granges, "score", MetaData::IntArray(score))?;
        set_meta_column(&mut granges, "thickStart", MetaData::IntArray(thick_start))?;
        set_meta_column(&mut granges, "thickEnd", MetaData::IntArray(thick_end))?;
        set_meta_column(&mut granges, "itemRgb", MetaData::StringArray(item_rgb))?;

        write_segmentation_file(&granges, bed_filename, bed_name, bed_description, compress)
    }
}

/* -------------------------------------------------------------------------- */

impl GenericMutableTrack<'_> {
    pub fn import_segmentation(
        &mut self,
        bed_filename: &str,
        def_state: Option<&str>,
    ) -> Result<Vec<String>, Box<dyn Error>> {
        let bin_size = self.track.get_bin_size();
        if bin_size == 0 {
            return Err(Box::new(io::Error::new(
                io::ErrorKind::InvalidInput,
                "track has invalid bin size",
            )));
        }

        let mut granges = GRanges::default();
        granges.import_bed9(bed_filename, is_gzip(bed_filename))?;
        let states = granges.meta.get_column_str("name").ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid segmentation bed file: name column is missing",
            )
        })?;

        let mut state_map = HashMap::new();
        if let Some(def_state) = def_state.filter(|state| !state.is_empty()) {
            state_map.insert(def_state.to_string(), 0usize);
            for seqname in self.track.get_seq_names() {
                let mut sequence = self.track.get_sequence_mut(&seqname)?;
                for i in 0..sequence.n_bins() {
                    sequence.set_bin(i, 0.0);
                }
            }
        }

        for i in 0..granges.num_rows() {
            let state_idx = if let Some(idx) = state_map.get(&states[i]).copied() {
                idx
            } else {
                let idx = state_map.len();
                state_map.insert(states[i].clone(), idx);
                idx
            };

            let mut sequence = self.track.get_sequence_mut(&granges.seqnames[i])?;
            let from = granges.ranges[i].from;
            let to = granges.ranges[i].to;

            for position in (from..to).step_by(bin_size) {
                if position / bin_size >= sequence.n_bins() {
                    break;
                }
                sequence.set(position, state_idx as f64);
            }
        }

        let mut state_names = vec![String::new(); state_map.len()];
        for (name, idx) in state_map {
            state_names[idx] = name;
        }
        Ok(state_names)
    }
}

/* -------------------------------------------------------------------------- */

fn set_meta_column(
    granges: &mut GRanges,
    name: &str,
    data: MetaData,
) -> Result<(), Box<dyn Error>> {
    granges.meta.delete_meta(name);
    granges.meta.add(name, data)
}

fn write_segmentation_file(
    granges: &GRanges,
    bed_filename: &str,
    bed_name: &str,
    bed_description: &str,
    compress: bool,
) -> Result<(), Box<dyn Error>> {
    let file = File::create(bed_filename)?;
    let writer: Box<dyn Write> = if compress {
        Box::new(GzEncoder::new(file, Compression::default()))
    } else {
        Box::new(file)
    };
    let mut writer = BufWriter::new(writer);
    writeln!(
        writer,
        "track name=\"{}\" description=\"{}\" visibility=1 itemRgb=\"On\"",
        bed_name, bed_description
    )?;
    granges.write_bed9(&mut writer)?;
    writer.flush()?;
    Ok(())
}

fn resolve_state_names(state_names: Option<&[String]>, max_state: usize) -> Vec<String> {
    if let Some(state_names) = state_names.filter(|names| !names.is_empty()) {
        state_names.to_vec()
    } else {
        (0..=max_state).map(|i| format!("s{i}")).collect()
    }
}

fn resolve_rgb_map(
    rgb_map: Option<&HashMap<String, String>>,
    state_names: &[String],
) -> HashMap<String, String> {
    if let Some(rgb_map) = rgb_map.filter(|map| !map.is_empty()) {
        rgb_map.clone()
    } else {
        let mut unique = state_names.to_vec();
        unique.sort();
        unique.dedup();
        let chart = get_n_colors(unique.len());
        unique.into_iter().zip(chart).collect()
    }
}

fn get_n_colors(n: usize) -> Vec<String> {
    if n == 0 {
        return Vec::new();
    }
    let k = (n as f64).cbrt() as usize + 1;
    let d = 255 / k;
    (0..n)
        .map(|i| {
            let r = ((i / k / k) % k) * d;
            let g = ((i / k) % k) * d;
            let b = (i % k) * d;
            format!("{r},{g},{b}")
        })
        .collect()
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use std::collections::HashMap;
    use std::env;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::genome::Genome;
    use crate::track::Track;
    use crate::track_generic::{GenericMutableTrack, GenericTrack};
    use crate::track_simple::SimpleTrack;

    fn temp_path(name: &str) -> String {
        let suffix = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        env::temp_dir()
            .join(format!("rustynetics_{suffix}_{name}"))
            .to_string_lossy()
            .into_owned()
    }

    #[test]
    fn export_segmentation_writes_bed9_track() {
        let genome = Genome::new(vec!["chr1".into()], vec![6]);
        let track = SimpleTrack::new(
            "states".into(),
            vec![vec![0.0, 0.0, 1.0, 1.0, 2.0, 2.0]],
            genome,
            1,
        )
        .unwrap();
        let path = temp_path("segmentation.bed");

        GenericTrack::wrap(&track)
            .export_segmentation(&path, "seg", "desc", false, None, None, &[])
            .unwrap();

        let content = fs::read_to_string(&path).unwrap();
        fs::remove_file(&path).ok();

        let lines: Vec<_> = content.lines().collect();
        assert_eq!(
            lines[0],
            "track name=\"seg\" description=\"desc\" visibility=1 itemRgb=\"On\""
        );
        let row1: Vec<_> = lines[1].split('\t').collect();
        let row2: Vec<_> = lines[2].split('\t').collect();
        let row3: Vec<_> = lines[3].split('\t').collect();
        assert_eq!(
            (&row1[3], &row1[4], &row1[6], &row1[7], &row1[8]),
            (&"s0", &"0", &"0", &"2", &"0,0,0")
        );
        assert_eq!(
            (&row2[3], &row2[4], &row2[6], &row2[7], &row2[8]),
            (&"s1", &"0", &"2", &"4", &"0,0,127")
        );
        assert_eq!(
            (&row3[3], &row3[4], &row3[6], &row3[7], &row3[8]),
            (&"s2", &"0", &"4", &"6", &"0,127,0")
        );
    }

    #[test]
    fn import_segmentation_assigns_state_ids_and_default_state() {
        let path = temp_path("import_segmentation.bed");
        fs::write(
            &path,
            concat!(
                "track name=\"seg\" description=\"desc\" visibility=1 itemRgb=\"On\"\n",
                "chr1\t1\t3\tactive\t0\t.\t1\t3\t255,0,0\n",
                "chr1\t4\t5\trepressed\t0\t.\t4\t5\t0,0,255\n"
            ),
        )
        .unwrap();

        let genome = Genome::new(vec!["chr1".into()], vec![6]);
        let mut track = SimpleTrack::alloc("states".into(), genome, f64::NAN, 1);
        let state_names =
            GenericMutableTrack::wrap(&mut track).import_segmentation(&path, Some("background"));

        fs::remove_file(&path).ok();

        let state_names = state_names.unwrap();
        assert_eq!(
            state_names,
            vec![
                "background".to_string(),
                "active".to_string(),
                "repressed".to_string()
            ]
        );
        assert_eq!(
            track.get_sequence("chr1").unwrap().clone_as_vec(),
            vec![0.0, 1.0, 1.0, 0.0, 2.0, 0.0]
        );
    }

    #[test]
    fn export_segmentation_uses_custom_state_names_and_rgb_map() {
        let genome = Genome::new(vec!["chr1".into()], vec![4]);
        let track =
            SimpleTrack::new("states".into(), vec![vec![0.0, 1.0, 1.0, 0.0]], genome, 1).unwrap();
        let path = temp_path("custom_segmentation.bed");
        let state_names = vec!["background".into(), "enhancer".into()];
        let rgb_map = HashMap::from([
            ("background".to_string(), "1,2,3".to_string()),
            ("enhancer".to_string(), "4,5,6".to_string()),
        ]);

        GenericTrack::wrap(&track)
            .export_segmentation(
                &path,
                "seg",
                "desc",
                false,
                Some(&state_names),
                Some(&rgb_map),
                &[],
            )
            .unwrap();

        let content = fs::read_to_string(&path).unwrap();
        fs::remove_file(&path).ok();
        let lines: Vec<_> = content.lines().skip(1).collect();
        let row1: Vec<_> = lines[0].split('\t').collect();
        let row2: Vec<_> = lines[1].split('\t').collect();
        let row3: Vec<_> = lines[2].split('\t').collect();
        assert_eq!((&row1[3], &row1[8]), (&"background", &"1,2,3"));
        assert_eq!((&row2[3], &row2[8]), (&"enhancer", &"4,5,6"));
        assert_eq!((&row3[3], &row3[8]), (&"background", &"1,2,3"));
    }
}
