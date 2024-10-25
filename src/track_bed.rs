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

use std::io::{Seek, Write};
use std::error::Error;

use crate::track_generic::GenericTrack;

/* -------------------------------------------------------------------------- */

impl GenericTrack<'_> {

    pub fn write_bed<W: Write + Seek>(&self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        let r = self.granges("score")?;

        r.write_bed6(writer)?;
        Ok(())
    }


    pub fn export_bed<W: Write + Seek>(&self, filename: &str, compress: bool) -> Result<(), Box<dyn Error>> {
        let r = self.granges("score")?;

        r.export_bed6(filename, compress)?;
        Ok(())
    }
}
