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

use std::cell::RefCell;
use std::io::{self, IsTerminal, Read, Write};
use std::rc::Rc;
use std::time::{Duration, Instant};

/* -------------------------------------------------------------------------- */

const PROGRESS_BAR_WIDTH: usize = 24;
const PROGRESS_DRAW_INTERVAL: Duration = Duration::from_millis(100);

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
struct ProgressState {
    label: String,
    total_bytes: u64,
    bytes_read: u64,
    records: usize,
    start: Instant,
    last_draw: Instant,
}

/* -------------------------------------------------------------------------- */

impl ProgressState {
    fn new(label: String, total_bytes: u64) -> Self {
        let now = Instant::now();

        Self {
            label,
            total_bytes,
            bytes_read: 0,
            records: 0,
            start: now,
            last_draw: now.checked_sub(PROGRESS_DRAW_INTERVAL).unwrap_or(now),
        }
    }

    fn inc_bytes(&mut self, bytes: usize) {
        self.bytes_read = self.bytes_read.saturating_add(bytes as u64);
        self.draw(false);
    }

    fn set_records(&mut self, records: usize) {
        self.records = records;
        self.draw(false);
    }

    fn finish(&mut self) -> io::Result<()> {
        self.draw_line(true)?;
        self.last_draw = Instant::now();
        Ok(())
    }

    fn clear(&mut self) -> io::Result<()> {
        let mut stderr = io::stderr().lock();
        write!(stderr, "\r\x1b[2K")?;
        stderr.flush()
    }

    fn draw(&mut self, force: bool) {
        if force || self.last_draw.elapsed() >= PROGRESS_DRAW_INTERVAL {
            let _ = self.draw_line(false);
            self.last_draw = Instant::now();
        }
    }

    fn draw_line(&mut self, finished: bool) -> io::Result<()> {
        let current = if finished {
            self.total_bytes.max(self.bytes_read)
        } else {
            self.bytes_read.min(self.total_bytes)
        };
        let total = self.total_bytes.max(1);
        let progress = (current as f64 / total as f64).clamp(0.0, 1.0);
        let filled = (progress * PROGRESS_BAR_WIDTH as f64).round() as usize;
        let bar = format!(
            "{}{}",
            "#".repeat(filled.min(PROGRESS_BAR_WIDTH)),
            "-".repeat(PROGRESS_BAR_WIDTH.saturating_sub(filled.min(PROGRESS_BAR_WIDTH)))
        );
        let percent = progress * 100.0;
        let elapsed = self.start.elapsed();
        let rate = if elapsed.as_secs_f64() > 0.0 {
            current as f64 / elapsed.as_secs_f64()
        } else {
            0.0
        };
        let eta = if rate > 0.0 && current < self.total_bytes {
            Duration::from_secs_f64((self.total_bytes - current) as f64 / rate)
        } else {
            Duration::default()
        };

        let mut stderr = io::stderr().lock();
        write!(
            stderr,
            "\r\x1b[2K{:<5} [{}] {:>6.2}% {}/{} {}/s eta {} reads {}",
            self.label,
            bar,
            percent,
            format_bytes(current),
            format_bytes(self.total_bytes),
            format_bytes(rate.round() as u64),
            format_duration(eta),
            self.records
        )?;
        if finished {
            writeln!(stderr)?;
        } else {
            stderr.flush()?;
        }
        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

/// A clonable handle for updating a progress bar from readers or worker code.
#[derive(Clone, Debug)]
pub struct ProgressHandle(Rc<RefCell<ProgressState>>);

/* -------------------------------------------------------------------------- */

impl ProgressHandle {
    fn new(label: String, total_bytes: u64) -> Self {
        Self(Rc::new(RefCell::new(ProgressState::new(label, total_bytes))))
    }

    /// Increments the tracked byte count.
    pub fn inc_bytes(&self, bytes: usize) {
        self.0.borrow_mut().inc_bytes(bytes);
    }

    /// Sets the current record count shown next to the progress bar.
    pub fn set_records(&self, records: usize) {
        self.0.borrow_mut().set_records(records);
    }

    fn finish(&self) -> io::Result<()> {
        self.0.borrow_mut().finish()
    }

    fn clear(&self) -> io::Result<()> {
        self.0.borrow_mut().clear()
    }
}

/* -------------------------------------------------------------------------- */

/// A scoped progress bar that draws to stderr and clears itself on early exit.
pub struct ProgressScope {
    handle: Option<ProgressHandle>,
    finished: bool,
}

/* -------------------------------------------------------------------------- */

impl ProgressScope {
    /// Creates a new progress bar that is only visible when stderr is a terminal.
    pub fn new<S: Into<String>>(label: S, total_bytes: u64) -> Self {
        Self::with_visibility(label, total_bytes, io::stderr().is_terminal())
    }

    /// Creates a new progress bar and explicitly controls whether it is visible.
    pub fn with_visibility<S: Into<String>>(label: S, total_bytes: u64, visible: bool) -> Self {
        let handle = if visible {
            Some(ProgressHandle::new(label.into(), total_bytes))
        } else {
            None
        };

        Self {
            handle,
            finished: false,
        }
    }

    /// Returns a clonable handle for use with `CountingReader` or manual updates.
    pub fn handle(&self) -> Option<ProgressHandle> {
        self.handle.clone()
    }

    /// Updates the record counter shown next to the progress bar.
    pub fn set_records(&self, records: usize) {
        if let Some(handle) = &self.handle {
            handle.set_records(records);
        }
    }

    /// Marks the progress bar as complete and prints a trailing newline.
    pub fn finish(&mut self, records: usize) {
        self.set_records(records);
        if let Some(handle) = &self.handle {
            let _ = handle.finish();
        }
        self.finished = true;
    }
}

/* -------------------------------------------------------------------------- */

impl Drop for ProgressScope {
    fn drop(&mut self) {
        if !self.finished {
            if let Some(handle) = &self.handle {
                let _ = handle.clear();
            }
        }
    }
}

/* -------------------------------------------------------------------------- */

/// A reader wrapper that forwards byte counts to a progress bar.
pub struct CountingReader<R> {
    inner: R,
    progress: Option<ProgressHandle>,
}

/* -------------------------------------------------------------------------- */

impl<R> CountingReader<R> {
    /// Wraps a reader and reports all successful reads to the provided progress handle.
    pub fn new(inner: R, progress: Option<ProgressHandle>) -> Self {
        Self { inner, progress }
    }
}

/* -------------------------------------------------------------------------- */

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let bytes_read = self.inner.read(buf)?;

        if bytes_read > 0 {
            if let Some(progress) = &self.progress {
                progress.inc_bytes(bytes_read);
            }
        }

        Ok(bytes_read)
    }
}

/* -------------------------------------------------------------------------- */

/// Formats a byte count using binary units.
pub fn format_bytes(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];

    let mut value = bytes as f64;
    let mut unit = UNITS[0];

    for next_unit in UNITS.iter().skip(1) {
        if value < 1024.0 {
            break;
        }
        value /= 1024.0;
        unit = next_unit;
    }

    if unit == "B" {
        format!("{} {}", bytes, unit)
    } else {
        format!("{value:.1} {unit}")
    }
}

/* -------------------------------------------------------------------------- */

/// Formats a duration as `MM:SS` or `HH:MM:SS`.
pub fn format_duration(duration: Duration) -> String {
    let seconds = duration.as_secs();
    let hours = seconds / 3600;
    let minutes = (seconds % 3600) / 60;
    let seconds = seconds % 60;

    if hours > 0 {
        format!("{hours:02}:{minutes:02}:{seconds:02}")
    } else {
        format!("{minutes:02}:{seconds:02}")
    }
}

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use super::{format_bytes, format_duration};

    #[test]
    fn format_duration_short() {
        assert_eq!(format_duration(std::time::Duration::from_secs(65)), "01:05");
    }

    #[test]
    fn format_bytes_binary_units() {
        assert_eq!(format_bytes(1536), "1.5 KiB");
    }
}
