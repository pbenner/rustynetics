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
use std::io::{Read, Write, Seek, SeekFrom};
use std::error::Error;
use std::path::Path;
use std::result::Result;
use std::io;
use std::thread;
use std::sync::mpsc::{channel, Sender};
use std::ops::Range;

use async_stream::stream;
use futures::executor::block_on_stream;
use futures::executor::BlockingStream;
use futures_core::stream::Stream;
use futures_util::pin_mut;

use reqwest::blocking::{Client, Response};
use byteorder::{ByteOrder, ReadBytesExt, WriteBytesExt, BigEndian, LittleEndian};

use crate::genome::Genome;
use crate::bbi::RVertex;
use crate::bbi::BbiQueryType;
use crate::bbi::BbiFile;

/* -------------------------------------------------------------------------- */

const BIGWIG_MAGIC : u32 = 0x888FFC26;

/* -------------------------------------------------------------------------- */

// Struct for BigWigParameters
#[derive(Debug)]
struct BigWigParameters {
    block_size      : usize,
    items_per_slot  : usize,
    reduction_levels: Option<Vec<i32>>,
}

/* -------------------------------------------------------------------------- */

impl BigWigParameters {
    // Default constructor for BigWigParameters
    pub fn default() -> Self {
        BigWigParameters {
            block_size      : 256,
            items_per_slot  : 1024,
            reduction_levels: None,
        }
    }
}

/* -------------------------------------------------------------------------- */

// Wrapper for a file or HTTP stream that supports Read + Seek
enum BigWigStream {
    File(File),
    Http(HttpSeekableReader),
}

struct BigWigFile {
    stream: BigWigStream,
}

impl BigWigFile {
    fn new(stream: BigWigStream) -> Self {
        BigWigFile { stream }
    }

    fn open_file(filename: &str) -> Result<BigWigFile, Box<dyn Error>> {
        let path = Path::new(filename);

        if path.exists() && path.is_file() {
            let file = File::open(path)?;
            Ok(BigWigFile::new(BigWigStream::File(file)))
        } else {
            Err(Box::new(io::Error::new(io::ErrorKind::NotFound, "File not found")))
        }
    }

    fn open_http(url: &str) -> Result<BigWigFile, Box<dyn Error>> {
        let client = Client::new();
        let head_resp = client.head(url).send()?;
        
        if head_resp.status().is_success() {
            if let Some(content_length) = head_resp.content_length() {
                let http_reader = HttpSeekableReader::new(client, url.to_string(), content_length);
                return Ok(BigWigFile::new(BigWigStream::Http(http_reader)));
            }
        }

        Err(Box::new(io::Error::new(io::ErrorKind::InvalidInput, "HTTP request failed")))
    }

    fn open_bigwig_file(filename: &str) -> Result<BigWigFile, Box<dyn Error>> {
        if filename.starts_with("http://") || filename.starts_with("https://") {
            BigWigFile::open_http(filename)
        } else {
            BigWigFile::open_file(filename)
        }
    }

}

/* -------------------------------------------------------------------------- */

// HTTP reader that supports seeking using Range requests
struct HttpSeekableReader {
    client: Client,
    url: String,
    content_length: u64,
    current_pos: u64,
}

impl HttpSeekableReader {
    fn new(client: Client, url: String, content_length: u64) -> Self {
        HttpSeekableReader {
            client,
            url,
            content_length,
            current_pos: 0,
        }
    }

    fn get_range(&self, range: Range<u64>) -> Result<Response, reqwest::Error> {
        let range_header = format!("bytes={}-{}", range.start, range.end - 1);
        self.client
            .get(&self.url)
            .header("Range", range_header)
            .send()
    }
}

impl Read for HttpSeekableReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let range_end = (self.current_pos + buf.len() as u64).min(self.content_length);
        let response = self
            .get_range(self.current_pos..range_end)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

        let bytes = response.bytes().map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        let bytes_read = bytes.len().min(buf.len());
        buf[..bytes_read].copy_from_slice(&bytes[..bytes_read]);
        self.current_pos += bytes_read as u64;
        Ok(bytes_read)
    }
}

impl Seek for HttpSeekableReader {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        let new_pos = match pos {
            SeekFrom::Start(p) => p,
            SeekFrom::End(p) => {
                if p >= 0 {
                    self.content_length
                } else {
                    (self.content_length as i64 + p) as u64
                }
            }
            SeekFrom::Current(p) => {
                if p >= 0 {
                    self.current_pos + p as u64
                } else {
                    (self.current_pos as i64 + p) as u64
                }
            }
        };

        if new_pos <= self.content_length {
            self.current_pos = new_pos;
            Ok(self.current_pos)
        } else {
            Err(io::Error::new(io::ErrorKind::InvalidInput, "Seek position out of bounds"))
        }
    }
}

/* -------------------------------------------------------------------------- */

#[derive(PartialEq)]
enum BigWigOrder {
    LE,
    BE,
}

/* -------------------------------------------------------------------------- */

struct BigWigReader<R: Read + Seek> {
    reader: R,
    bwf   : BbiFile,
    genome: Genome,
    order : BigWigOrder
}

/* -------------------------------------------------------------------------- */

impl<R: Read + Seek> BigWigReader<R> {

    fn open(reader: &mut R) -> io::Result<(BbiFile, BigWigOrder)> {

        let mut bwf = BbiFile::new();

        // Try to open with little endian
        let r1 = bwf.open::<LittleEndian, R>(reader, BIGWIG_MAGIC);

        if let Err(err) = r1 {
            if err.kind() != io::ErrorKind::InvalidData || err.to_string() != "Invalid magic number" {
                return Err(err)
            }
        } else {
            return Ok((bwf, BigWigOrder::LE))
        }

        let r2 = bwf.open::<BigEndian, R>(reader, BIGWIG_MAGIC);

        if let Err(err) = r2 {
            return Err(err)
        } else {
            return Ok((bwf, BigWigOrder::BE))
        }
    }

    fn new(mut reader: R) -> Result<Self, Box<dyn Error>> {

        let (mut bwf, mut order) = BigWigReader::<R>::open(&mut reader)?;
    
        let mut seqnames = vec![String::new(); bwf.chrom_data.keys.len()];
        let mut lengths  = vec![0; bwf.chrom_data.keys.len()];

        for i in 0..bwf.chrom_data.keys.len() {
            if bwf.chrom_data.values[i].len() != 8 {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "invalid chromosome list",
                )));
            }
            let idx = match order {
                BigWigOrder::LE => (&bwf.chrom_data.values[i][0..4]).read_u32::<LittleEndian>()? as usize,
                BigWigOrder::BE => (&bwf.chrom_data.values[i][0..4]).read_u32::<BigEndian>()? as usize,
            };
            if idx >= bwf.chrom_data.keys.len() {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "invalid chromosome index",
                )));
            }
            seqnames[idx] = String::from_utf8_lossy(&bwf.chrom_data.keys[i]).trim_end_matches('\x00').to_string();
            lengths [idx] = match order {
                BigWigOrder::LE => (&bwf.chrom_data.values[i][4..8]).read_u32::<LittleEndian>()? as usize,
                BigWigOrder::BE => (&bwf.chrom_data.values[i][4..8]).read_u32::<BigEndian>()? as usize
            };
        }

        let genome = Genome::new(seqnames, lengths);

        Ok(BigWigReader {
            reader,
            bwf,
            genome,
            order,
        })
    }
    /*
    fn read_blocks(&mut self) -> std::sync::mpsc::Receiver<BigWigReaderType> {
        let (tx, rx) = channel();
        let bwf = self.bwf.clone();
        let mut reader = self.reader.clone();

        thread::spawn(move || {
            Self::fill_channel(tx, &mut reader, &bwf.index.root).unwrap();
        });

        rx
    }

    fn fill_channel(tx: Sender<BigWigReaderType>, reader: &mut R, vertex: &RVertex) -> Result<(), Box<dyn Error>> {
        if vertex.is_leaf != 0 {
            for i in 0..vertex.n_children as usize {
                match vertex.read_block(reader, i) {
                    Ok(block) => tx.send(BigWigReaderType { block: Some(block), error: None }).unwrap(),
                    Err(err)  => tx.send(BigWigReaderType { block: None, error: Some(Box::new(err)) }).unwrap(),
                }
            }
        } else {
            for i in 0..vertex.n_children as usize {
                Self::fill_channel(tx.clone(), reader, &vertex.children[i])?;
            }
        }
        Ok(())
    }
    */
    fn query<'a>(
        &'a mut self,
        seq_regex: &'a str,
        from     : usize,
        to       : usize,
        bin_size : usize
    ) -> impl Stream<Item = io::Result<BbiQueryType>> + 'a {

        stream! {

            let re = regex::Regex::new(&format!("^{}$", seq_regex)).unwrap();

            for seqname in &self.genome.seqnames {
                if !re.is_match(seqname) {
                    continue;
                }
                if let Some(idx) = self.genome.get_idx(seqname) {

                    let mut iterator = match self.order {
                        BigWigOrder::LE => self.bwf.query_iterator::<LittleEndian, R>(&mut self.reader, idx as u32, from as u32, to as u32, bin_size as u32),
                        BigWigOrder::BE => self.bwf.query_iterator::<BigEndian   , R>(&mut self.reader, idx as u32, from as u32, to as u32, bin_size as u32),
                    };
                    while let Some(item) = iterator.next() {

                        yield item;

                    }
                }
            }
        }
    }

    pub fn query_iterator<'a>(
        &'a mut self,
        seq_regex: &'a str,
        from     : usize,
        to       : usize,
        bin_size : usize,
    ) -> BlockingStream<impl Stream<Item = io::Result<BbiQueryType>> + 'a> {

        let s = Box::pin(self.query(seq_regex, from, to, bin_size));
        //pin_mut!(s);

        block_on_stream(s)

    }
}
