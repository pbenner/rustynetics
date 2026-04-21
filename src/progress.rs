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
use std::env;
use std::io::{self, IsTerminal, Read, Write};
use std::rc::Rc;
use std::time::{Duration, Instant};

/* -------------------------------------------------------------------------- */

const PROGRESS_BAR_WIDTH: usize = 24;
const PROGRESS_BAR_WIDTH_COMPACT: usize = 16;
const PROGRESS_BAR_WIDTH_MINIMAL: usize = 12;
const PROGRESS_DRAW_INTERVAL: Duration = Duration::from_millis(100);
const DEFAULT_TERMINAL_WIDTH: usize = 80;

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
        let rendered = render_progress_line(
            &self.label,
            progress,
            current,
            self.total_bytes,
            rate.round() as u64,
            eta,
            self.records,
            terminal_width_hint(),
        );

        let mut stderr = io::stderr().lock();
        write!(stderr, "\r\x1b[2K{}", rendered)?;
        if finished {
            writeln!(stderr)?;
        } else {
            stderr.flush()?;
        }
        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

fn terminal_width_hint() -> usize {
    env::var("COLUMNS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .filter(|width| *width >= 40)
        .unwrap_or(DEFAULT_TERMINAL_WIDTH)
        .saturating_sub(1)
}

/* -------------------------------------------------------------------------- */

fn render_progress_bar(progress: f64, width: usize) -> String {
    let filled = (progress * width as f64).round() as usize;

    format!(
        "{}{}",
        "#".repeat(filled.min(width)),
        "-".repeat(width.saturating_sub(filled.min(width)))
    )
}

/* -------------------------------------------------------------------------- */

fn render_progress_line(
    label: &str,
    progress: f64,
    current: u64,
    total: u64,
    rate: u64,
    eta: Duration,
    records: usize,
    max_width: usize,
) -> String {
    let percent = format!("{:>6.2}%", progress * 100.0);
    let bytes = format!("{}/{}", format_bytes(current), format_bytes(total));
    let rate = format!("{}/s", format_bytes(rate));
    let eta = format!("eta {}", format_duration(eta));
    let reads = format!("reads {}", format_count(records));
    let candidates = [
        compose_progress_line(
            label,
            &percent,
            &render_progress_bar(progress, PROGRESS_BAR_WIDTH),
            &[&bytes, &rate, &eta, &reads],
        ),
        compose_progress_line(
            label,
            &percent,
            &render_progress_bar(progress, PROGRESS_BAR_WIDTH_COMPACT),
            &[&bytes, &reads, &rate],
        ),
        compose_progress_line(
            label,
            &percent,
            &render_progress_bar(progress, PROGRESS_BAR_WIDTH_MINIMAL),
            &[&bytes, &reads],
        ),
        format!("{} {} {}", label, percent, reads),
    ];

    for candidate in candidates {
        if candidate.chars().count() <= max_width {
            return candidate;
        }
    }

    truncate_with_ellipsis(&format!("{} {}", label, percent), max_width)
}

/* -------------------------------------------------------------------------- */

fn compose_progress_line(label: &str, percent: &str, bar: &str, fields: &[&str]) -> String {
    let mut line = format!("{label} [{bar}] {percent}");

    for field in fields {
        line.push(' ');
        line.push_str(field);
    }

    line
}

/* -------------------------------------------------------------------------- */

fn truncate_with_ellipsis(text: &str, max_width: usize) -> String {
    if text.chars().count() <= max_width {
        return text.to_string();
    }

    if max_width <= 3 {
        return ".".repeat(max_width);
    }

    let mut truncated = String::with_capacity(max_width);

    for ch in text.chars().take(max_width - 3) {
        truncated.push(ch);
    }
    truncated.push_str("...");

    truncated
}

/* -------------------------------------------------------------------------- */

fn format_count(value: usize) -> String {
    let digits = value.to_string();
    let mut formatted = String::with_capacity(digits.len() + digits.len() / 3);

    for (index, ch) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index) % 3 == 0 {
            formatted.push(',');
        }
        formatted.push(ch);
    }

    formatted
}

/* -------------------------------------------------------------------------- */

/// A clonable handle for updating a progress bar from readers or worker code.
#[derive(Clone, Debug)]
pub struct ProgressHandle(Rc<RefCell<ProgressState>>);

/* -------------------------------------------------------------------------- */

impl ProgressHandle {
    fn new(label: String, total_bytes: u64) -> Self {
        Self(Rc::new(RefCell::new(ProgressState::new(
            label,
            total_bytes,
        ))))
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

    use std::time::Duration;

    use super::{format_bytes, format_duration, render_progress_line, truncate_with_ellipsis};

    #[test]
    fn format_duration_short() {
        assert_eq!(format_duration(Duration::from_secs(65)), "01:05");
    }

    #[test]
    fn format_bytes_binary_units() {
        assert_eq!(format_bytes(1536), "1.5 KiB");
    }

    #[test]
    fn progress_line_respects_requested_width() {
        let line = render_progress_line(
            "FASTQ 12/12",
            1.0,
            2_147_483_648,
            2_147_483_648,
            123_456_789,
            Duration::from_secs(0),
            123_456_789,
            79,
        );

        assert!(line.chars().count() <= 79, "line exceeds width: {line}");
    }

    #[test]
    fn progress_line_falls_back_for_narrow_widths() {
        let line = render_progress_line(
            "FASTQ 12/12",
            0.5,
            1_073_741_824,
            2_147_483_648,
            67_108_864,
            Duration::from_secs(16),
            12_345_678,
            32,
        );

        assert!(line.chars().count() <= 32, "line exceeds width: {line}");
    }

    #[test]
    fn ellipsis_truncation_preserves_width() {
        let text = truncate_with_ellipsis("abcdefghijklmnopqrstuvwxyz", 10);

        assert_eq!(text, "abcdefg...");
        assert_eq!(text.chars().count(), 10);
    }
}
