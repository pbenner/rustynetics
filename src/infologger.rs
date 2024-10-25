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

pub struct Logger(pub Box<RefCell<dyn Write>>);

/* -------------------------------------------------------------------------- */

impl Logger {

    pub fn new_null() -> Logger {
        Logger(Box::new(RefCell::new(io::sink())))
    }

    pub fn new_stdout() -> Logger {
        Logger(Box::new(RefCell::new(io::stdout())))
    }

    pub fn new_stderr() -> Logger {
        Logger(Box::new(RefCell::new(io::stderr())))
    }

}

/* -------------------------------------------------------------------------- */

impl fmt::Write for Logger {
    fn write_str(&mut self, s: &str) -> Result<(), Error> {
        self.0.borrow_mut().write(s.as_bytes()).map_err(|_io_err| Error)?;
        Ok(())
    }
}

/* -------------------------------------------------------------------------- */

#[macro_export]
macro_rules! log {
    ($logger:expr, $fmt:expr $(, $args:expr)*) => {{
        if let Err(err) = write!($logger.0.borrow_mut(), $fmt $(, $args)*) {
            eprintln!("Logging error: {}", err); // Print error to standard error if logging fails
        }
    }};
}
