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
use core::pin::Pin;

use futures::{Stream, StreamExt};
use async_stream::stream;

use crate::read;
use crate::infologger::Logger;

/* -------------------------------------------------------------------------- */

type ReadStreamType<'a> = Pin<Box<dyn Stream<Item = io::Result<read::Read>> + 'a>>;

/* -------------------------------------------------------------------------- */

pub enum ReadStream {}

/* -------------------------------------------------------------------------- */

impl ReadStream {

    /// Converts a stream of paired reads into a stream of single-end reads.
    ///
    /// If the `switch` parameter is set to `true`, the function will modify the input stream to
    /// indicate that all reads are single-end. If `switch` is `false`, the input stream is returned
    /// unchanged.
    ///
    /// # Parameters
    ///
    /// - `stream_in`: The input stream of reads.
    /// - `logger`: An optional logger for logging messages.
    /// - `switch`: A boolean flag to indicate whether to perform the conversion.
    ///
    /// # Returns
    ///
    /// A new stream where all reads are marked as single-end if `switch` is `true`.
    pub fn paired_as_single_end<'a>(
        mut stream_in: ReadStreamType<'a>,
        _logger      : Option<&'a Logger>,
        switch       : bool,
    ) -> ReadStreamType<'a> {

        // If PairedAsSingleEnd is false, return the input stream directly
        if !switch {
            return stream_in;
        }

        let output_stream = Box::pin(stream! {
            while let Some(item) = stream_in.next().await {

                match item {
                    Ok(mut r) => {
                        r.paired_end = false;
                        yield Ok(r)
                    },
                    Err(e) => yield Err(e),
                };
            }
        });

        Box::pin(output_stream)
    }

/* -------------------------------------------------------------------------- */

    /// Filters the input stream to include only paired-end reads.
    ///
    /// If the `switch` parameter is set to `true`, the function will filter out unpaired reads
    /// from the input stream. If `switch` is `false`, the input stream is returned unchanged.
    ///
    /// # Parameters
    ///
    /// - `stream_in`: The input stream of reads.
    /// - `logger`: An optional logger for logging the number of filtered reads.
    /// - `switch`: A boolean flag to indicate whether to filter paired-end reads.
    ///
    /// # Returns
    ///
    /// A new stream that only includes paired-end reads if `switch` is `true`.
    pub fn filter_paired_end<'a>(
        mut stream_in: ReadStreamType<'a>,
        logger       : Option<&'a Logger>,
        switch       : bool,
    ) -> ReadStreamType<'a> {

        if !switch {
            return stream_in;
        }

        let output_stream = async_stream::stream! {
            let mut n = 0;
            let mut m = 0;

            while let Some(item) = stream_in.next().await {
                match item {
                    Ok(r) => {
                        if r.paired_end {
                            yield Ok(r);
                            m += 1;
                        }
                        n += 1;
                    },
                    Err(e) => yield Err(e),
                }
            }

            if let Some(log) = logger {
                log!(log, "Filtered out {} unpaired reads ({:.2}%)", n - m, 100.0 * (n - m) as f64 / n as f64);
            }
        };

        Box::pin(output_stream)
    }

/* -------------------------------------------------------------------------- */

    /// Filters the input stream to include only single-end reads.
    ///
    /// If the `switch` parameter is set to `true`, the function will filter out paired reads
    /// from the input stream. If `switch` is `false`, the input stream is returned unchanged.
    ///
    /// # Parameters
    ///
    /// - `stream_in`: The input stream of reads.
    /// - `logger`: An optional logger for logging the number of filtered reads.
    /// - `switch`: A boolean flag to indicate whether to filter single-end reads.
    ///
    /// # Returns
    ///
    /// A new stream that only includes single-end reads if `switch` is `true`.
    pub fn filter_single_end<'a>(
        mut stream_in: ReadStreamType<'a>,
        logger       : Option<&'a Logger>,
        switch       : bool,
    ) -> ReadStreamType<'a> {

        if !switch {
            return stream_in;
        }

        let output_stream = async_stream::stream! {
            let mut n = 0;
            let mut m = 0;

            while let Some(item) = stream_in.next().await {
                match item {
                    Ok(r) => {
                        if !r.paired_end {
                            yield Ok(r);
                            m += 1;
                        }
                        n += 1;
                    },
                    Err(e) => yield Err(e),
                }
            }

            if let Some(log) = logger {
                log!(log, "Filtered out {} paired reads ({:.2}%)", n - m, 100.0 * (n - m) as f64 / n as f64);
            }
        };

        Box::pin(output_stream)
    }

/* -------------------------------------------------------------------------- */

    /// Filters the input stream to exclude duplicate reads.
    ///
    /// If the `switch` parameter is set to `true`, the function will filter out duplicates
    /// from the input stream. If `switch` is `false`, the input stream is returned unchanged.
    ///
    /// # Parameters
    ///
    /// - `stream_in`: The input stream of reads.
    /// - `logger`: An optional logger for logging the number of filtered duplicates.
    /// - `switch`: A boolean flag to indicate whether to filter duplicates.
    ///
    /// # Returns
    ///
    /// A new stream that only includes non-duplicate reads if `switch` is `true`.
    pub fn filter_duplicates<'a>(
        mut stream_in: ReadStreamType<'a>,
        logger       : Option<&'a Logger>,
        switch       : bool,
    ) -> ReadStreamType<'a> {

        if !switch {
            return stream_in;
        }

        let output_stream = async_stream::stream! {
            let mut n = 0;
            let mut m = 0;

            while let Some(item) = stream_in.next().await {
                match item {
                    Ok(r) => {
                        if !r.duplicate {
                            yield Ok(r);
                            m += 1;
                        }
                        n += 1;
                    },
                    Err(e) => yield Err(e),
                }
            }

            if let Some(log) = logger {
                log!(log, "Filtered out {} duplicates ({:.2}%)", n - m, 100.0 * (n - m) as f64 / n as f64);
            }
        };

        Box::pin(output_stream)
    }

/* -------------------------------------------------------------------------- */

    /// Filters the input stream based on strand information.
    ///
    /// If the `strand` parameter is '*', the function returns the input stream unchanged.
    /// Otherwise, it filters out reads that do not match the specified strand.
    ///
    /// # Parameters
    ///
    /// - `stream_in`: The input stream of reads.
    /// - `logger`: An optional logger for logging the number of filtered reads.
    /// - `strand`: The strand character ('+' or '-') to filter reads by.
    ///
    /// # Returns
    ///
    /// A new stream that only includes reads matching the specified strand.
    pub fn filter_strand<'a>(
        mut stream_in: ReadStreamType<'a>,
        logger       : Option<&'a Logger>,
        strand       : char,
    ) -> ReadStreamType<'a> {

        if strand == '*' {
            return stream_in;
        }

        let output_stream = async_stream::stream! {
            let mut n = 0;
            let mut m = 0;

            while let Some(item) = stream_in.next().await {
                match item {
                    Ok(r) => {
                        if r.strand == strand {
                            yield Ok(r);
                            m += 1;
                        }
                        n += 1;
                    },
                    Err(e) => yield Err(e),
                }
            }

            if let Some(log) = logger {
                log!(log, "Filtered out {} reads not on strand {} ({:.2}%)", n - m, strand, 100.0 * (n - m) as f64 / n as f64);
            }
        };

        Box::pin(output_stream)
    }

/* -------------------------------------------------------------------------- */

    /// Filters the input stream based on mapping quality.
    ///
    /// The function filters out reads with a mapping quality lower than the specified `mapq` value.
    /// If `mapq` is less than or equal to zero, the input stream is returned unchanged.
    ///
    /// # Parameters
    ///
    /// - `stream_in`: The input stream of reads.
    /// - `logger`: An optional logger for logging the number of filtered reads.
    /// - `mapq`: The minimum mapping quality to include reads.
    ///
    /// # Returns
    ///
    /// A new stream that only includes reads with mapping quality above the specified threshold.
    pub fn filter_mapq<'a>(
        mut stream_in: ReadStreamType<'a>,
        logger       : Option<&'a Logger>,
        mapq         : i64,
    ) -> ReadStreamType<'a> {

        if mapq <= 0 {
            return stream_in;
        }

        let output_stream = async_stream::stream! {
            let mut n = 0;
            let mut m = 0;

            while let Some(item) = stream_in.next().await {
                match item {
                    Ok(r) => {
                        if r.mapq >= mapq {
                            yield Ok(r);
                            m += 1;
                        }
                        n += 1;
                    },
                    Err(e) => yield Err(e),
                }
            }

            if let Some(log) = logger {
                log!(log, "Filtered out {} reads with mapping quality lower than {} ({:.2}%)", n - m, mapq, 100.0 * (n - m) as f64 / n as f64);
            }
        };

        Box::pin(output_stream)
    }

/* -------------------------------------------------------------------------- */

    /// Filters the input stream based on read length.
    ///
    /// The function filters out reads whose lengths do not fall within the specified range.
    /// If both bounds are zero, the input stream is returned unchanged.
    ///
    /// # Parameters
    ///
    /// - `stream_in`: The input stream of reads.
    /// - `logger`: An optional logger for logging the number of filtered reads.
    /// - `read_lengths`: A slice containing the minimum and maximum read lengths.
    ///
    /// # Returns
    ///
    /// A new stream that only includes reads within the specified length range.
    pub fn filter_read_length<'a>(
        mut stream_in: ReadStreamType<'a>,
        logger       : Option<&'a Logger>,
        read_lengths : &'a [usize; 2],
    ) -> ReadStreamType<'a> {

        if read_lengths[0] == 0 && read_lengths[1] == 0 {
            return stream_in;
        }

        let output_stream = async_stream::stream! {
            let mut n = 0;
            let mut m = 0;

            while let Some(item) = stream_in.next().await {
                match item {
                    Ok(r) => {
                        let len = r.range.to - r.range.from;
                        if len >= read_lengths[0] &&
                        (len <= read_lengths[1] || read_lengths[1] == 0) {
                            yield Ok(r);
                            m += 1;
                        }
                        n += 1;
                    },
                    Err(e) => yield Err(e),
                }
            }

            if let Some(log) = logger {
                log!(log, "Filtered out {} reads with non-admissible length ({:.2}%)", n - m, 100.0 * (n - m) as f64 / n as f64);
            }
        };

        Box::pin(output_stream)
    }

/* -------------------------------------------------------------------------- */

    /// Shifts the reads based on their strand.
    ///
    /// This function modifies the read positions based on the specified shift values for forward and
    /// reverse strands. If both shift values are zero, the input stream is returned unchanged.
    ///
    /// # Parameters
    ///
    /// - `stream_in`: The input stream of reads.
    /// - `logger`: An optional logger for logging the shift operation.
    /// - `shift`: An array containing the shift values for forward and reverse strands.
    ///
    /// # Returns
    ///
    /// A new stream with shifted reads based on the specified strand shifts.
    pub fn shift_reads<'a>(
        mut stream_in: ReadStreamType<'a>,
        logger       : Option<&'a Logger>,
        shift        : &'a [usize; 2]
    ) -> ReadStreamType<'a> {

        if shift[0] == 0 && shift[1] == 0 {
            return stream_in;
        }

        let output_stream = async_stream::stream! {
            while let Some(item) = stream_in.next().await {
                match item {
                    Ok(mut r) => {
                        if r.strand == '+' {
                            r.range.from += shift[0];
                            r.range.to   += shift[0];
                        } else if r.strand == '-' {
                            r.range.from += shift[1];
                            r.range.to   += shift[1];
                        }

                        r.range.to  -= r.range.from;
                        r.range.from = 0;

                        yield Ok(r);
                    },
                    Err(e) => yield Err(e),
                }
            }

            if let Some(log) = logger {
                log!(log, "Shifted reads (forward strand: {}, reverse strand: {})", shift[0], shift[1]);
            }
        };

        Box::pin(output_stream)
    }

}
