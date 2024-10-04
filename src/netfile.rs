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

use std::io;
use std::error::Error;
use std::path::Path;
use std::result::Result;
use std::io::{Read, Seek, SeekFrom};
use std::fs::File;
use std::ops::Range;

use reqwest::blocking::{Client, Response};

/* -------------------------------------------------------------------------- */

// Wrapper for a file or HTTP stream that supports Read + Seek
#[derive(Debug)]
enum NetFileStream {
    File(File),
    Http(HttpSeekableReader),
}

#[derive(Debug)]
pub struct NetFile {
    stream: NetFileStream,
}

impl NetFile {
    fn new(stream: NetFileStream) -> Self {
        NetFile { stream }
    }

    fn open_file(filename: &str) -> Result<NetFile, Box<dyn Error>> {
        let path = Path::new(filename);

        if path.exists() && path.is_file() {
            let file = File::open(path)?;
            Ok(NetFile::new(NetFileStream::File(file)))
        } else {
            Err(Box::new(io::Error::new(io::ErrorKind::NotFound, "File not found")))
        }
    }

    fn open_http(url: &str) -> Result<NetFile, Box<dyn Error>> {
        let client    = Client::new();
        let head_resp = client.head(url).send()?;

        if !head_resp.status().is_success() {
            return Err(Box::new(io::Error::new(io::ErrorKind::InvalidInput, "HTTP request failed")));
        }

        let content_length = head_resp
            .headers()
            .get("Content-Length")
            .ok_or_else(|| Box::new(io::Error::new(io::ErrorKind::InvalidData, "Missing Content-Length header")))?
            .to_str()
            .map_err(|_| Box::new(io::Error::new(io::ErrorKind::InvalidData, "Invalid Content-Length header")))?
            .parse::<u64>()
            .map_err(|_| Box::new(io::Error::new(io::ErrorKind::InvalidData, "Invalid Content-Length header")))?;

        let http_reader = HttpSeekableReader::new(client, url.to_string(), content_length);

        Ok(NetFile::new(NetFileStream::Http(http_reader)))
    }

    pub fn open(filename: &str) -> Result<NetFile, Box<dyn Error>> {
        if filename.starts_with("http://") || filename.starts_with("https://") {
            NetFile::open_http(filename)
        } else {
            NetFile::open_file(filename)
        }
    }

}

impl Read for NetFile {

    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match &mut self.stream {
            NetFileStream::File(file) => file.read(buf),
            NetFileStream::Http(file) => file.read(buf)
        }
    }

}

impl Seek for NetFile {

    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        match &mut self.stream {
            NetFileStream::File(file) => file.seek(pos),
            NetFileStream::Http(file) => file.seek(pos)
        }
    }

}

/* -------------------------------------------------------------------------- */

// HTTP reader that supports seeking using Range requests
#[derive(Debug)]
struct HttpSeekableReader {
    client        : Client,
    url           : String,
    current_pos   : u64,
    content_length: u64,
}

impl HttpSeekableReader {

    fn new(client: Client, url: String, content_length: u64) -> Self {
        HttpSeekableReader {
            client,
            url,
            current_pos: 0,
            content_length,
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
        let range_end = self.current_pos + buf.len() as u64;
        let response  = self
            .get_range(self.current_pos..range_end)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        let bytes      = response.bytes().map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        let bytes_read = bytes.len().min(buf.len());
        buf[..bytes_read].copy_from_slice(&bytes[..bytes_read]);
        self.current_pos += bytes_read as u64;
        Ok(bytes_read)
    }

}

impl Seek for HttpSeekableReader {

    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        let new_pos = match pos {
            SeekFrom::Start(p) => {
                p
            },
            SeekFrom::Current(p) => {
                if p >= 0 {
                    self.current_pos + p as u64
                } else {
                    self.current_pos.saturating_sub((-p) as u64)
                }
            }
            SeekFrom::End(p) => {
                if p >= 0 {
                    self.content_length + p as u64
                } else {
                    self.content_length.saturating_sub((-p) as u64)
                }
            }
        };

        // Prevent seeking beyond EOF
        if new_pos > self.content_length {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "Seek position beyond file size"));
        }

        self.current_pos = new_pos;

        Ok(new_pos)
    }

}
