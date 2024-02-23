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
use std::io;

/* -------------------------------------------------------------------------- */

#[derive(Debug)]
pub enum Error {
    Generic(String),
    IO(io::Error)
}

/* -------------------------------------------------------------------------- */

impl From<String> for Error {
    fn from(str : String) -> Self {
        Error::Generic(str)
    }
}

/* -------------------------------------------------------------------------- */

impl From<io::Error> for Error {
    fn from(e : io::Error) -> Self {
        Error::IO(e)
    }
}

/* -------------------------------------------------------------------------- */

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::Generic(v) => f.pad(&format!("{}", v)),
            Error::IO(v)      => f.pad(&format!("{}", v))
        }
    }
}
