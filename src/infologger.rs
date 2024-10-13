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

use std::fmt;
use std::io::{self, Write};
use std::cell::RefCell;
use std::fmt::Error;

/* -------------------------------------------------------------------------- */

pub struct Logger(Box<RefCell<dyn Write>>);

/* -------------------------------------------------------------------------- */

impl Logger {

    pub fn new_void() -> Logger {
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
