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

use std::any::Any;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::granges::GRanges;
use crate::granges_table_reader::GRangesTableReader;
use crate::granges_table_writer::GRangesTableWriter;
use crate::meta_table_reader::MetaTableReader;

/* -------------------------------------------------------------------------- */

fn read_table_all_schema(content: &str) -> io::Result<(Vec<String>, Vec<String>)> {
    let header = content
        .lines()
        .find(|line| !line.trim().is_empty())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "table input is empty"))?;
    let mut names = Vec::new();
    let mut types = Vec::new();

    for field in header.split_whitespace() {
        match field {
            "seqnames" | "from" | "to" | "start" | "end" | "strand" => {}
            _ => {
                names.push(field.to_string());
                types.push("Vec<String>".to_string());
            }
        }
    }

    Ok((names, types))
}

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub struct OptionPrintScientific(pub bool);

#[derive(Debug)]
pub struct OptionPrintStrand(pub bool);

/* -------------------------------------------------------------------------- */

impl GRanges {
    /// Reads a table from a buffered reader and populates the `GRanges` structure.
    ///
    /// This method reads from the provided `BufRead` instance, interpreting the first line
    /// as a header and subsequent lines as data. It utilizes `MetaTableReader` to read metadata
    /// and `GRangesTableReader` for genomic ranges.
    ///
    /// # Arguments
    /// - `buf_reader`: A mutable buffered reader to read from.
    /// - `names`: An array of expected column names.
    /// - `types`: An array of expected column types.
    ///
    /// # Returns
    /// An `io::Result<()>`, which will be `Ok(())` if the operation succeeds, or an error if any issues arise during reading.
    pub fn bufread_table<R: BufRead>(
        &mut self,
        mut buf_reader: R,
        names: &[&str],
        types: &[&str],
    ) -> io::Result<()> {
        let mut mreader = MetaTableReader::new(names, types);
        let mut greader = GRangesTableReader::new();

        let mut line = String::new();

        // Read first line as header
        buf_reader.read_line(&mut line)?;
        greader.read_header(&line)?;
        mreader.read_header(&line)?;

        line.clear();

        let mut line_counter = 0;

        while buf_reader.read_line(&mut line)? > 0 {
            if line.is_empty() {
                continue;
            }
            greader.read_line(&line, line_counter)?;
            mreader.read_line(&line, line_counter)?;

            line_counter += 1;

            line.clear();
        }
        greader.push(self);
        mreader.push(&mut self.meta);

        Ok(())
    }

    /// Writes the `GRanges` table to the specified writer.
    ///
    /// This method formats and writes the contents of the `GRanges` structure to the provided writer,
    /// optionally printing scientific notation or strand information based on the given arguments.
    ///
    /// # Arguments
    /// - `writer`: A mutable reference to a writer where the table will be output.
    /// - `args`: An array of dynamic arguments that may include options for scientific notation and strand.
    ///
    /// # Returns
    /// An `io::Result<()>`, which will be `Ok(())` if the operation succeeds, or an error if writing fails.
    pub fn write_table<W: Write>(&self, writer: &mut W, args: &[&dyn Any]) -> io::Result<()> {
        let mut use_scientific = false;
        let mut use_strand = false;

        for arg in args {
            if let Some(option) = arg.downcast_ref::<OptionPrintScientific>() {
                use_scientific = option.0;
            }
            if let Some(option) = arg.downcast_ref::<OptionPrintStrand>() {
                use_strand = option.0;
            }
        }

        let mut gwriter = GRangesTableWriter::new(self, use_scientific, use_strand);

        let meta_str = self.meta.print_table(args);
        let mut meta_reader = BufReader::new(meta_str.as_bytes());

        gwriter.determine_widths()?;
        gwriter.write_header(writer, &mut meta_reader)?;

        for i in 0..self.num_rows() {
            gwriter.write_row(writer, &mut meta_reader, i)?;
        }
        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

impl GRanges {
    /// Reads a table from a reader and populates the `GRanges` structure.
    ///
    /// This method wraps the buffered reading functionality, creating a `BufReader` and
    /// invoking `bufread_table` to read from it.
    ///
    /// # Arguments
    /// - `reader`: A reader to read from.
    /// - `names`: An array of expected column names.
    /// - `types`: An array of expected column types.
    ///
    /// # Returns
    /// An `io::Result<()>`, which will be `Ok(())` if the operation succeeds, or an error if reading fails.
    pub fn read_table<R: Read>(
        &mut self,
        reader: R,
        names: &[&str],
        types: &[&str],
    ) -> io::Result<()> {
        let buf_reader = BufReader::new(reader);

        self.bufread_table(buf_reader, names, types)
    }
}

/* -------------------------------------------------------------------------- */

impl GRanges {
    /// Reads a table and imports all metadata columns as string matrices.
    ///
    /// This mirrors the behavior of Go's `ReadTableAll`, which preserves all
    /// non-core metadata columns as textual data instead of inferring numeric
    /// types from their content.
    pub fn read_table_all<R: Read>(&mut self, mut reader: R) -> io::Result<()> {
        let mut content = String::new();
        reader.read_to_string(&mut content)?;

        let (names, types) = read_table_all_schema(&content)?;
        let name_refs = names.iter().map(String::as_str).collect::<Vec<_>>();
        let type_refs = types.iter().map(String::as_str).collect::<Vec<_>>();

        self.read_table(io::Cursor::new(content), &name_refs, &type_refs)
    }

    /// Imports a table from a file and preserves all metadata columns as
    /// string matrices.
    pub fn import_table_all(&mut self, filename: &str, compress: bool) -> io::Result<()> {
        let file = File::open(filename)?;
        let reader: Box<dyn Read> = if compress {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };

        self.read_table_all(reader)
    }
}

/* -------------------------------------------------------------------------- */

impl GRanges {
    /// Imports a table from a file and populates the `GRanges` structure.
    ///
    /// This method opens the specified file, optionally decompressing it if required,
    /// and reads the data into the `GRanges` instance using `bufread_table`.
    ///
    /// # Arguments
    /// - `filename`: The path to the file to import.
    /// - `names`: An array of expected column names.
    /// - `types`: An array of expected column types.
    /// - `compress`: A boolean indicating whether the file is compressed.
    ///
    /// # Returns
    /// An `io::Result<()>`, which will be `Ok(())` if the operation succeeds, or an error if reading fails.
    pub fn import_table(
        &mut self,
        filename: &str,
        names: &[&str],
        types: &[&str],
        compress: bool,
    ) -> io::Result<()> {
        let file = File::open(filename)?;
        let mut reader: Box<dyn BufRead> = if compress {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        self.bufread_table(&mut reader, names, types)?;

        Ok(())
    }

    /// Exports the `GRanges` table to a file.
    ///
    /// This method creates a file at the specified path, optionally compressing the output,
    /// and writes the `GRanges` data to it using `write_table`.
    ///
    /// # Arguments
    /// - `filename`: The path to the file to export to.
    /// - `compress`: A boolean indicating whether to compress the output file.
    /// - `args`: An array of dynamic arguments that may include options for writing.
    ///
    /// # Returns
    /// An `io::Result<()>`, which will be `Ok(())` if the operation succeeds, or an error if writing fails.
    pub fn export_table(
        &self,
        filename: &str,
        compress: bool,
        args: &[&dyn Any],
    ) -> io::Result<()> {
        let file = File::create(filename)?;
        let mut writer: Box<dyn Write> = if compress {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };

        self.write_table(&mut writer, args)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meta::MetaData;

    #[test]
    fn read_table_all_preserves_textual_columns() {
        let table = "\
seqnames from to strand score labels
chr1 1 10 + 7 a,b
chr2 5 12 - 42 nil
";
        let mut granges = GRanges::default();

        granges.read_table_all(io::Cursor::new(table)).unwrap();

        assert_eq!(
            granges.meta.get_column("score"),
            Some(&MetaData::StringMatrix(vec![
                vec!["7".to_string()],
                vec!["42".to_string()],
            ]))
        );
        assert_eq!(
            granges.meta.get_column("labels"),
            Some(&MetaData::StringMatrix(vec![
                vec!["a".to_string(), "b".to_string()],
                Vec::<String>::new(),
            ]))
        );
    }
}
