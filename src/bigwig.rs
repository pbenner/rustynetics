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
use std::io::{BufReader, BufRead, Read, Seek, SeekFrom};
use std::error::Error;
use std::path::Path;
use std::result::Result;
use std::io;

use url::Url;
use std::ops::Range;
use reqwest::blocking::{Client, Response};

/* -------------------------------------------------------------------------- */

// Constants
const BIGWIG_MAGIC: u32 = 0x888FFC26;

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
