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
use std::io::{self, Cursor, Write};

use flate2::write::GzEncoder;
use flate2::Compression;

use crate::genes::Genes;
use crate::meta::Meta;

/* -------------------------------------------------------------------------- */

impl Genes {
    fn table_meta(&self) -> Meta {
        let mut meta = self.granges.meta.clone();
        meta.delete_meta("names");
        meta.delete_meta("cds");
        meta
    }

    fn table_meta_lines(&self, args: &[&dyn Any]) -> Vec<String> {
        let meta = self.table_meta();
        if meta.num_cols() == 0 {
            Vec::new()
        } else {
            meta.print_table(args).lines().map(String::from).collect()
        }
    }

    /// Writes a gene table with `names`, transcript coordinates, CDS coordinates, and extra metadata.
    pub fn write_table<W: Write>(
        &self,
        writer: &mut W,
        header: bool,
        args: &[&dyn Any],
    ) -> io::Result<()> {
        let names = self.names();
        let cds = self.cds();
        let meta_lines = self.table_meta_lines(args);

        let mut widths = [
            "names".len(),
            "seqnames".len(),
            "strand".len(),
            "txStart".len(),
            "txEnd".len(),
            "cdsStart".len(),
            "cdsEnd".len(),
        ];

        for i in 0..self.granges.num_rows() {
            widths[0] = widths[0].max(names[i].len());
            widths[1] = widths[1].max(self.granges.seqnames[i].len());
            widths[2] = widths[2].max(1);
            widths[3] = widths[3].max(self.granges.ranges[i].from.to_string().len());
            widths[4] = widths[4].max(self.granges.ranges[i].to.to_string().len());
            widths[5] = widths[5].max(cds[i].from.to_string().len());
            widths[6] = widths[6].max(cds[i].to.to_string().len());
        }

        if header {
            write!(
                writer,
                "{:>w0$} {:>w1$} {:>w2$} {:>w3$} {:>w4$} {:>w5$} {:>w6$}",
                "names",
                "seqnames",
                "strand",
                "txStart",
                "txEnd",
                "cdsStart",
                "cdsEnd",
                w0 = widths[0],
                w1 = widths[1],
                w2 = widths[2],
                w3 = widths[3],
                w4 = widths[4],
                w5 = widths[5],
                w6 = widths[6],
            )?;
            if !meta_lines.is_empty() {
                write!(writer, " {}", meta_lines[0])?;
            }
            writeln!(writer)?;
        }

        for i in 0..self.granges.num_rows() {
            if i != 0 {
                writeln!(writer)?;
            }
            write!(
                writer,
                "{:>w0$} {:>w1$} {:>w2$} {:>w3$} {:>w4$} {:>w5$} {:>w6$}",
                names[i],
                self.granges.seqnames[i],
                self.granges.strand[i],
                self.granges.ranges[i].from,
                self.granges.ranges[i].to,
                cds[i].from,
                cds[i].to,
                w0 = widths[0],
                w1 = widths[1],
                w2 = widths[2],
                w3 = widths[3],
                w4 = widths[4],
                w5 = widths[5],
                w6 = widths[6],
            )?;
            if !meta_lines.is_empty() {
                write!(writer, " {}", meta_lines[i + 1])?;
            }
        }

        Ok(())
    }

    pub fn print_table(&self, header: bool, args: &[&dyn Any]) -> String {
        let mut buffer = Vec::new();
        {
            let mut writer = Cursor::new(&mut buffer);
            self.write_table(&mut writer, header, args).unwrap();
        }
        String::from_utf8(buffer).unwrap()
    }

    pub fn export_table(
        &self,
        filename: &str,
        header: bool,
        compress: bool,
        args: &[&dyn Any],
    ) -> io::Result<()> {
        let file = File::create(filename)?;
        let mut writer: Box<dyn Write> = if compress {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };

        self.write_table(&mut writer, header, args)
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use crate::genes::Genes;
    use crate::meta::MetaData;

    #[test]
    fn genes_write_table_includes_core_columns_and_extra_meta() {
        let mut genes = Genes::new(
            vec!["gene1".into(), "gene2".into()],
            vec!["chr1".into(), "chr2".into()],
            vec![100, 200],
            vec![150, 260],
            vec![110, 210],
            vec![140, 250],
            vec!['+', '-'],
        );
        genes
            .granges
            .meta
            .add("expr", MetaData::FloatArray(vec![1.5, 2.5]))
            .unwrap();

        let table = genes.print_table(true, &[]);
        let lines: Vec<&str> = table.lines().collect();

        assert_eq!(
            lines[0].split_whitespace().collect::<Vec<_>>(),
            vec!["names", "seqnames", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "expr"]
        );
        assert_eq!(
            lines[1].split_whitespace().collect::<Vec<_>>(),
            vec!["gene1", "chr1", "+", "100", "150", "110", "140", "1.5"]
        );
        assert_eq!(
            lines[2].split_whitespace().collect::<Vec<_>>(),
            vec!["gene2", "chr2", "-", "200", "260", "210", "250", "2.5"]
        );
    }
}
