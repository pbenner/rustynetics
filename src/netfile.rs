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
const BUFFER_SIZE: usize = 8192; // Example buffer size (8 KB)

#[derive(Debug)]
struct HttpSeekableReader {
    client        : Client,
    url           : String,
    current_pos   : u64,
    content_length: u64,
    buffer        : Vec<u8>,   // Buffer to hold fetched data
    buffer_start  : u64,       // Start position of the buffer in the file
    buffer_end    : u64,       // End position of the buffer in the file
}

impl HttpSeekableReader {
    fn new(client: Client, url: String, content_length: u64) -> Self {
        HttpSeekableReader {
            client,
            url,
            current_pos : 0,
            content_length,
            buffer      : Vec::new(),
            buffer_start: 0,
            buffer_end  : 0,
        }
    }

    fn get_range(&self, range: Range<u64>) -> Result<Response, reqwest::Error> {
        let range_header = format!("bytes={}-{}", range.start, range.end - 1);
        self.client
            .get(&self.url)
            .header("Range", range_header)
            .send()
    }

    fn fill_buffer(&mut self) -> io::Result<()> {
        let range_end     = (self.current_pos + BUFFER_SIZE as u64).min(self.content_length);
        let response      = self
            .get_range(self.current_pos..range_end)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

        let bytes         = response.bytes().map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        self.buffer       = bytes.to_vec(); // Store the fetched data in the buffer
        self.buffer_start = self.current_pos;
        self.buffer_end   = self.buffer_start + self.buffer.len() as u64;
        Ok(())
    }
}

impl Read for HttpSeekableReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        // Refill the buffer if it's empty or current_pos is outside the buffer range
        if self.current_pos < self.buffer_start || self.current_pos >= self.buffer_end {
            self.fill_buffer()?;
        }

        // Calculate how much we can read from the buffer
        let buffer_offset   = (self.current_pos - self.buffer_start) as usize;
        let available_bytes = (self.buffer_end - self.current_pos) as usize;
        let bytes_to_read   = buf.len().min(available_bytes);

        buf[..bytes_to_read].copy_from_slice(&self.buffer[buffer_offset..buffer_offset + bytes_to_read]);
        self.current_pos += bytes_to_read as u64;

        Ok(bytes_to_read)
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
