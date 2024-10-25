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

/* -------------------------------------------------------------------------- */

pub trait Alphabet {
    fn bases(&self, i: u8) -> Result<Vec<u8>, String>;
    fn matching(&self, i: u8) -> Result<Vec<u8>, String>;
    fn code(&self, i: u8) -> Result<u8, String>;
    fn decode(&self, i: u8) -> Result<u8, String>;
    fn is_ambiguous(&self, i: u8) -> Result<bool, String>;
    fn is_wildcard(&self, i: u8) -> Result<bool, String>;
    fn length(&self) -> usize;
    fn length_unambiguous(&self) -> usize;
    fn string(&self) -> String;
}

pub trait ComplementableAlphabet: Alphabet {
    fn complement(&self, i: u8) -> Result<u8, String>;
    fn complement_coded(&self, i: u8) -> Result<u8, String>;
}

/* -------------------------------------------------------------------------- */

pub struct NucleotideAlphabet;

/* -------------------------------------------------------------------------- */

impl Alphabet for NucleotideAlphabet {
    fn bases(&self, i: u8) -> Result<Vec<u8>, String> {
        match i {
            b'A' | b'a' => Ok(vec![b'a']),
            b'C' | b'c' => Ok(vec![b'c']),
            b'G' | b'g' => Ok(vec![b'g']),
            b'T' | b't' => Ok(vec![b't']),
            _ => Err(format!("Bases(): `{}` is not part of the alphabet", i as char)),
        }
    }

    fn matching(&self, i: u8) -> Result<Vec<u8>, String> {
        match i {
            b'A' | b'a' => Ok(vec![b'a']),
            b'C' | b'c' => Ok(vec![b'c']),
            b'G' | b'g' => Ok(vec![b'g']),
            b'T' | b't' => Ok(vec![b't']),
            _ => Err(format!(
                "Matching(): `{}` is not a non-ambiguous letter of the alphabet",
                i as char
            )),
        }
    }

    fn code(&self, i: u8) -> Result<u8, String> {
        match i {
            b'A' | b'a' => Ok(0),
            b'C' | b'c' => Ok(1),
            b'G' | b'g' => Ok(2),
            b'T' | b't' => Ok(3),
            _ => Err(format!("Code(): `{}` is not part of the alphabet", i as char)),
        }
    }

    fn decode(&self, i: u8) -> Result<u8, String> {
        match i {
            0 => Ok(b'a'),
            1 => Ok(b'c'),
            2 => Ok(b'g'),
            3 => Ok(b't'),
            _ => Err(format!("Decode(): `{}` is not a code of the alphabet", i)),
        }
    }

    fn is_ambiguous(&self, i: u8) -> Result<bool, String> {
        match i {
            b'A' | b'a' => Ok(false),
            b'C' | b'c' => Ok(false),
            b'G' | b'g' => Ok(false),
            b'T' | b't' => Ok(false),
            _ => Err(format!("IsAmbiguous(): `{}` is not part of the alphabet", i as char)),
        }
    }

    fn is_wildcard(&self, i: u8) -> Result<bool, String> {
        match i {
            b'A' | b'a' => Ok(false),
            b'C' | b'c' => Ok(false),
            b'G' | b'g' => Ok(false),
            b'T' | b't' => Ok(false),
            _ => Err(format!("IsWildcard(): `{}` is not part of the alphabet", i as char)),
        }
    }

    fn length(&self) -> usize {
        4
    }

    fn length_unambiguous(&self) -> usize {
        4
    }

    fn string(&self) -> String {
        String::from("nucleotide alphabet")
    }
}

/* -------------------------------------------------------------------------- */

impl ComplementableAlphabet for NucleotideAlphabet {
    fn complement_coded(&self, i: u8) -> Result<u8, String> {
        match i {
            0 => Ok(3),
            1 => Ok(2),
            2 => Ok(1),
            3 => Ok(0),
            _ => Err(format!(
                "ComplementCoded(): `{}` is not a code of the alphabet",
                i
            )),
        }
    }

    fn complement(&self, i: u8) -> Result<u8, String> {
        match i {
            b'A' | b'a' => Ok(b't'),
            b'C' | b'c' => Ok(b'g'),
            b'G' | b'g' => Ok(b'c'),
            b'T' | b't' => Ok(b'a'),
            _ => Err(format!("Complement(): `{}` is not part of the alphabet", i as char)),
        }
    }
}

/* -------------------------------------------------------------------------- */

pub struct GappedNucleotideAlphabet;

/* -------------------------------------------------------------------------- */

impl Alphabet for GappedNucleotideAlphabet {
    fn bases(&self, i: u8) -> Result<Vec<u8>, String> {
        match i {
            b'A' | b'a' => Ok(vec![b'a']),
            b'C' | b'c' => Ok(vec![b'c']),
            b'G' | b'g' => Ok(vec![b'g']),
            b'T' | b't' => Ok(vec![b't']),
            b'N' | b'n' => Ok(vec![b'a', b'c', b'g', b't']),
            _ => Err(format!("Bases(): `{}` is not part of the alphabet", i as char)),
        }
    }

    fn matching(&self, i: u8) -> Result<Vec<u8>, String> {
        match i {
            b'A' | b'a' => Ok(vec![b'a', b'n']),
            b'C' | b'c' => Ok(vec![b'c', b'n']),
            b'G' | b'g' => Ok(vec![b'g', b'n']),
            b'T' | b't' => Ok(vec![b't', b'n']),
            _ => Err(format!(
                "Matching(): `{}` is not a non-ambiguous letter of the alphabet",
                i as char
            )),
        }
    }

    fn code(&self, i: u8) -> Result<u8, String> {
        match i {
            b'A' | b'a' => Ok(0),
            b'C' | b'c' => Ok(1),
            b'G' | b'g' => Ok(2),
            b'T' | b't' => Ok(3),
            b'N' | b'n' => Ok(4),
            _ => Err(format!("Code(): `{}` is not part of the alphabet", i as char)),
        }
    }

    fn decode(&self, i: u8) -> Result<u8, String> {
        match i {
            0 => Ok(b'a'),
            1 => Ok(b'c'),
            2 => Ok(b'g'),
            3 => Ok(b't'),
            4 => Ok(b'n'),
            _ => Err(format!("Decode(): `{}` is not a code of the alphabet", i)),
        }
    }

    fn is_ambiguous(&self, i: u8) -> Result<bool, String> {
        match i {
            b'A' | b'a' => Ok(false),
            b'C' | b'c' => Ok(false),
            b'G' | b'g' => Ok(false),
            b'T' | b't' => Ok(false),
            b'N' | b'n' => Ok(true),
            _ => Err(format!("is_ambiguous(): `{}` is not part of the alphabet", i as char)),
        }
    }

    fn is_wildcard(&self, i: u8) -> Result<bool, String> {
        match i {
            b'A' | b'a' => Ok(false),
            b'C' | b'c' => Ok(false),
            b'G' | b'g' => Ok(false),
            b'T' | b't' => Ok(false),
            b'N' | b'n' => Ok(true),
            _ => Err(format!("is_wildcard(): `{}` is not part of the alphabet", i as char)),
        }
    }

    fn length(&self) -> usize {
        5
    }

    fn length_unambiguous(&self) -> usize {
        4
    }

    fn string(&self) -> String {
        String::from("gapped nucleotide alphabet")
    }
}

/* -------------------------------------------------------------------------- */

impl ComplementableAlphabet for GappedNucleotideAlphabet {
    fn complement_coded(&self, i: u8) -> Result<u8, String> {
        match i {
            0 => Ok(3),
            1 => Ok(2),
            2 => Ok(1),
            3 => Ok(0),
            4 => Ok(4),
            _ => Err(format!(
                "ComplementCoded(): `{}` is not a code of the alphabet",
                i
            )),
        }
    }

    fn complement(&self, i: u8) -> Result<u8, String> {
        match i {
            b'A' | b'a' => Ok(b't'),
            b'C' | b'c' => Ok(b'g'),
            b'G' | b'g' => Ok(b'c'),
            b'T' | b't' => Ok(b'a'),
            b'N' | b'n' => Ok(b'n'),
            _ => Err(format!("Complement(): `{}` is not part of the alphabet", i as char)),
        }
    }
}

/* -------------------------------------------------------------------------- */

pub struct AmbiguousNucleotideAlphabet;

/* -------------------------------------------------------------------------- */

impl Alphabet for AmbiguousNucleotideAlphabet {
    fn bases(&self, i: u8) -> Result<Vec<u8>, String> {
        match i {
            b'A' | b'a' => Ok(vec![b'a']),
            b'C' | b'c' => Ok(vec![b'c']),
            b'G' | b'g' => Ok(vec![b'g']),
            b'T' | b't' => Ok(vec![b't']),
            b'W' | b'w' => Ok(vec![b'a', b't']),
            b'S' | b's' => Ok(vec![b'c', b'g']),
            b'M' | b'm' => Ok(vec![b'a', b'c']),
            b'K' | b'k' => Ok(vec![b'g', b't']),
            b'R' | b'r' => Ok(vec![b'a', b'g']),
            b'Y' | b'y' => Ok(vec![b'c', b't']),
            b'B' | b'b' => Ok(vec![b'c', b'g', b't']),
            b'D' | b'd' => Ok(vec![b'a', b'g', b't']),
            b'H' | b'h' => Ok(vec![b'a', b'c', b't']),
            b'V' | b'v' => Ok(vec![b'a', b'c', b'g']),
            b'N' | b'n' => Ok(vec![b'a', b'c', b'g', b't']),
            _ => Err(format!("Bases(): `{}` is not part of the alphabet", i as char)),
        }
    }

    fn matching(&self, i: u8) -> Result<Vec<u8>, String> {
        match i {
            b'A' | b'a' => Ok(vec![
                b'a', b'w', b'm', b'r', b'd', b'h', b'v', b'n',
            ]),
            b'C' | b'c' => Ok(vec![
                b'c', b's', b'm', b'y', b'b', b'h', b'v', b'n',
            ]),
            b'G' | b'g' => Ok(vec![
                b'g', b's', b'k', b'r', b'b', b'd', b'v', b'n',
            ]),
            b'T' | b't' => Ok(vec![
                b't', b'w', b'k', b'y', b'b', b'd', b'h', b'n',
            ]),
            _ => Err(format!(
                "Matching(): `{}` is not a non-ambiguous letter of the alphabet",
                i as char
            )),
        }
    }

    fn code(&self, i: u8) -> Result<u8, String> {
        match i {
            b'A' | b'a' => Ok(0),
            b'C' | b'c' => Ok(1),
            b'G' | b'g' => Ok(2),
            b'T' | b't' => Ok(3),
            b'W' | b'w' => Ok(4),
            b'S' | b's' => Ok(5),
            b'M' | b'm' => Ok(6),
            b'K' | b'k' => Ok(7),
            b'R' | b'r' => Ok(8),
            b'Y' | b'y' => Ok(9),
            b'B' | b'b' => Ok(10),
            b'D' | b'd' => Ok(11),
            b'H' | b'h' => Ok(12),
            b'V' | b'v' => Ok(13),
            b'N' | b'n' => Ok(14),
            _ => Err(format!("Code(): `{}` is not part of the alphabet", i as char)),
        }
    }

    fn decode(&self, i: u8) -> Result<u8, String> {
        match i {
            0 => Ok(b'a'),
            1 => Ok(b'c'),
            2 => Ok(b'g'),
            3 => Ok(b't'),
            4 => Ok(b'w'),
            5 => Ok(b's'),
            6 => Ok(b'm'),
            7 => Ok(b'k'),
            8 => Ok(b'r'),
            9 => Ok(b'y'),
            10 => Ok(b'b'),
            11 => Ok(b'd'),
            12 => Ok(b'h'),
            13 => Ok(b'v'),
            14 => Ok(b'n'),
            _ => Err(format!("Decode(): `{}` is not a code of the alphabet", i)),
        }
    }

    fn is_ambiguous(&self, i: u8) -> Result<bool, String> {
        match i {
            b'A' | b'a' => Ok(false),
            b'C' | b'c' => Ok(false),
            b'G' | b'g' => Ok(false),
            b'T' | b't' => Ok(false),
            b'W' | b'w' => Ok(true),
            b'S' | b's' => Ok(true),
            b'M' | b'm' => Ok(true),
            b'K' | b'k' => Ok(true),
            b'R' | b'r' => Ok(true),
            b'Y' | b'y' => Ok(true),
            b'B' | b'b' => Ok(true),
            b'D' | b'd' => Ok(true),
            b'H' | b'h' => Ok(true),
            b'V' | b'v' => Ok(true),
            b'N' | b'n' => Ok(true),
            _ => Err(format!("is_ambiguous(): `{}` is not part of the alphabet", i as char)),
        }
    }

    fn is_wildcard(&self, i: u8) -> Result<bool, String> {
        match i {
            b'A' | b'a' => Ok(false),
            b'C' | b'c' => Ok(false),
            b'G' | b'g' => Ok(false),
            b'T' | b't' => Ok(false),
            b'W' | b'w' => Ok(false),
            b'S' | b's' => Ok(false),
            b'M' | b'm' => Ok(false),
            b'K' | b'k' => Ok(false),
            b'R' | b'r' => Ok(false),
            b'Y' | b'y' => Ok(false),
            b'B' | b'b' => Ok(false),
            b'D' | b'd' => Ok(false),
            b'H' | b'h' => Ok(false),
            b'V' | b'v' => Ok(false),
            b'N' | b'n' => Ok(true),
            _ => Err(format!("is_wildcard(): `{}` is not part of the alphabet", i as char)),
        }
    }

    fn length(&self) -> usize {
        15
    }

    fn length_unambiguous(&self) -> usize {
        4
    }

    fn string(&self) -> String {
        String::from("ambiguous nucleotide alphabet")
    }
}

/* -------------------------------------------------------------------------- */

impl ComplementableAlphabet for AmbiguousNucleotideAlphabet {
    fn complement_coded(&self, i: u8) -> Result<u8, String> {
        match i {
            0 => Ok(3),
            1 => Ok(2),
            2 => Ok(1),
            3 => Ok(0),
            4 => Ok(4),
            5 => Ok(5),
            6 => Ok(7),
            7 => Ok(6),
            8 => Ok(9),
            9 => Ok(8),
            10 => Ok(13),
            11 => Ok(12),
            12 => Ok(11),
            13 => Ok(10),
            14 => Ok(14),
            _ => Err(format!(
                "ComplementCoded(): `{}` is not a code of the alphabet",
                i
            )),
        }
    }

    fn complement(&self, i: u8) -> Result<u8, String> {
        match i {
            b'A' | b'a' => Ok(b't'),
            b'C' | b'c' => Ok(b'g'),
            b'G' | b'g' => Ok(b'c'),
            b'T' | b't' => Ok(b'a'),
            b'W' | b'w' => Ok(b'w'),
            b'S' | b's' => Ok(b's'),
            b'M' | b'm' => Ok(b'k'),
            b'K' | b'k' => Ok(b'm'),
            b'R' | b'r' => Ok(b'y'),
            b'Y' | b'y' => Ok(b'r'),
            b'B' | b'b' => Ok(b'v'),
            b'D' | b'd' => Ok(b'h'),
            b'H' | b'h' => Ok(b'd'),
            b'V' | b'v' => Ok(b'b'),
            b'N' | b'n' => Ok(b'n'),
            _ => Err(format!("Complement(): `{}` is not part of the alphabet", i as char)),
        }
    }
}
