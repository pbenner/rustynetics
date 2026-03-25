use std::ffi::c_char;
use std::slice;
use std::str::{self, Utf8Error};

/// Default symbol name for Rust-oriented `bigwig-map` plugins.
pub const BIGWIG_MAP_RUST_SYMBOL: &str = "rustynetics_bigwig_map";

/// Compatibility symbol used by older shared-library mappers.
pub const BIGWIG_MAP_LEGACY_SYMBOL: &str = "bigwig_map";

/// Legacy `bigwig-map` ABI:
/// `(seqname_ptr, position, values_ptr, values_len) -> f64`.
pub type BigWigMapLegacyFn = unsafe extern "C" fn(*const c_char, usize, *const f64, usize) -> f64;

/// Rust-oriented `bigwig-map` ABI:
/// `(&BigWigMapInput) -> f64`.
pub type BigWigMapRustFn = unsafe extern "C" fn(*const BigWigMapInput) -> f64;

/// FFI payload passed to Rust-oriented `bigwig-map` plugins.
///
/// A plugin can depend on `rustynetics` and export:
///
/// ```ignore
/// use rustynetics::bigwig_map_plugin::{BigWigMapInput, BIGWIG_MAP_RUST_SYMBOL};
///
/// #[no_mangle]
/// pub unsafe extern "C" fn rustynetics_bigwig_map(input: *const BigWigMapInput) -> f64 {
///     let input = &*input;
///     let values = input.values();
///     values.iter().copied().sum()
/// }
/// ```
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct BigWigMapInput {
    pub seqname_ptr: *const c_char,
    pub seqname_len: usize,
    pub position: usize,
    pub values_ptr: *const f64,
    pub values_len: usize,
}

impl BigWigMapInput {
    /// Returns the chromosome name as raw bytes.
    ///
    /// # Safety
    /// `seqname_ptr` must point to a valid buffer of length `seqname_len`.
    pub unsafe fn seqname_bytes(&self) -> &[u8] {
        if self.seqname_len == 0 {
            &[]
        } else {
            slice::from_raw_parts(self.seqname_ptr as *const u8, self.seqname_len)
        }
    }

    /// Returns the chromosome name as UTF-8 text.
    ///
    /// # Safety
    /// `seqname_ptr` must point to a valid UTF-8 buffer of length `seqname_len`.
    pub unsafe fn seqname(&self) -> Result<&str, Utf8Error> {
        str::from_utf8(self.seqname_bytes())
    }

    /// Returns the input values for the current genomic position.
    ///
    /// # Safety
    /// `values_ptr` must point to a valid buffer of length `values_len`.
    pub unsafe fn values(&self) -> &[f64] {
        if self.values_len == 0 {
            &[]
        } else {
            slice::from_raw_parts(self.values_ptr, self.values_len)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::CString;

    #[test]
    fn input_exposes_seqname_and_values() {
        let seqname = CString::new("chr1").unwrap();
        let values = [1.0, 2.5, 3.0];
        let input = BigWigMapInput {
            seqname_ptr: seqname.as_ptr(),
            seqname_len: 4,
            position: 10,
            values_ptr: values.as_ptr(),
            values_len: values.len(),
        };

        let parsed_seqname = unsafe { input.seqname().unwrap() };
        let parsed_values = unsafe { input.values() };

        assert_eq!(parsed_seqname, "chr1");
        assert_eq!(parsed_values, &values);
    }
}
