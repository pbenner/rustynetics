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

use std::fmt;
use std::io::{self, Write};
use std::cell::RefCell;
use std::fmt::Error;

/* -------------------------------------------------------------------------- */

/// A simple logging structure that wraps a writer object.
///
/// The `Logger` struct provides a flexible mechanism for logging output, allowing the
/// target output to be directed to any type that implements `Write`. It includes methods
/// for creating loggers directed to standard output, standard error, or a "null" output
/// that discards all logs.
///
/// # Examples
///
/// ```
/// use rustynetics::infologger::Logger;
/// use rustynetics::log;
///
/// let logger = Logger::new_stdout();
/// log!(logger, "This is a log message to standard output.");
/// ```
pub struct Logger(pub Box<RefCell<dyn Write>>);

/* -------------------------------------------------------------------------- */

impl Logger {

    /// Creates a `Logger` instance that discards all output.
    ///
    /// This is useful for situations where logging is required by function
    /// signatures or interfaces but where the output is not needed.
    ///
    /// # Examples
    ///
    /// ```
    /// use rustynetics::infologger::Logger;
    /// use rustynetics::log;
    ///
    /// let logger = Logger::new_null();
    /// log!(logger, "This message will not be printed.");
    /// ```
    pub fn new_null() -> Logger {
        Logger(Box::new(RefCell::new(io::sink())))
    }

    /// Creates a `Logger` instance that writes to standard output (`stdout`).
    ///
    /// This is useful for debugging or when log messages should be visible in the console.
    ///
    /// # Examples
    ///
    /// ```
    /// use rustynetics::infologger::Logger;
    /// use rustynetics::log;
    ///
    /// let logger = Logger::new_stdout();
    /// log!(logger, "This message goes to standard output.");
    /// ```
    pub fn new_stdout() -> Logger {
        Logger(Box::new(RefCell::new(io::stdout())))
    }

    /// Creates a `Logger` instance that writes to standard error (`stderr`).
    ///
    /// This is useful for logging error or warning messages separately from standard output.
    ///
    /// # Examples
    ///
    /// ```
    /// use rustynetics::infologger::Logger;
    /// use rustynetics::log;
    ///
    /// let logger = Logger::new_stderr();
    /// log!(logger, "This message goes to standard error.");
    /// ```
    pub fn new_stderr() -> Logger {
        Logger(Box::new(RefCell::new(io::stderr())))
    }

}

/* -------------------------------------------------------------------------- */

impl fmt::Write for Logger {
    /// Writes a string slice to the logger’s internal writer.
    ///
    /// This method allows the `Logger` to implement `fmt::Write` so that it can be used
    /// in formatting macros like `write!`. It converts the string slice into bytes and
    /// writes them to the inner writer.
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the inner writer fails.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::fmt::Write;
    /// use rustynetics::infologger::Logger;
    ///
    /// let mut logger = Logger::new_stdout();
    /// write!(&mut logger, "Formatted log message: {}", 42).unwrap();
    /// ```
    fn write_str(&mut self, s: &str) -> Result<(), Error> {
        self.0.borrow_mut().write(s.as_bytes()).map_err(|_io_err| Error)?;
        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

/// A macro for logging messages to a specified `Logger` instance.
///
/// This macro provides a convenient way to log formatted messages to the provided
/// `Logger`. If logging fails, it outputs the error to standard error.
///
/// # Parameters
///
/// - `$logger`: The `Logger` instance to use for logging.
/// - `$fmt`: A format string.
/// - `$args`: Optional arguments for the format string.
///
/// # Examples
///
/// ```
/// use rustynetics::infologger::Logger;
/// use rustynetics::log;
///
/// let logger = Logger::new_stdout();
/// log!(logger, "Log message with argument: {}", 42);
/// ```
#[macro_export]
macro_rules! log {
    ($logger:expr, $fmt:expr $(, $args:expr)*) => {{
        if let Err(err) = write!($logger.0.borrow_mut(), $fmt $(, $args)*) {
            eprintln!("Logging error: {}", err); // Print error to standard error if logging fails
        }
    }};
}
