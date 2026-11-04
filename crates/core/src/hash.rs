//! Content hashes (BLAKE3) for provenance: foil geometry, meshes and
//! polar families.  Floating-point inputs are quantized before hashing so a
//! hash does not change with the last bits of a recomputed coordinate.

use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;

/// A 32-byte BLAKE3 digest.  Serialized as a 64-character hex string.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct Hash32(pub [u8; 32]);

impl Hash32 {
    pub const ZERO: Hash32 = Hash32([0; 32]);
    pub fn is_zero(&self) -> bool {
        self.0 == [0; 32]
    }
    pub fn hex(&self) -> String {
        self.0.iter().map(|b| format!("{b:02x}")).collect()
    }
    /// The first 12 hex digits (for messages).
    pub fn short(&self) -> String {
        self.hex()[..12].to_string()
    }
    pub fn from_hex(s: &str) -> Option<Hash32> {
        let s = s.trim();
        if s.len() != 64 || !s.is_ascii() {
            return None;
        }
        let mut out = [0u8; 32];
        for (i, o) in out.iter_mut().enumerate() {
            *o = u8::from_str_radix(&s[2 * i..2 * i + 2], 16).ok()?;
        }
        Some(Hash32(out))
    }
}

impl fmt::Debug for Hash32 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Hash32({})", self.short())
    }
}

impl fmt::Display for Hash32 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.hex())
    }
}

impl Serialize for Hash32 {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(&self.hex())
    }
}

impl<'de> Deserialize<'de> for Hash32 {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let s = String::deserialize(d)?;
        Hash32::from_hex(&s).ok_or_else(|| serde::de::Error::custom(format!("not a 64-digit hex hash: '{s}'")))
    }
}

/// Incremental hasher with typed, quantized inputs.
pub struct Hasher(blake3::Hasher);

impl Default for Hasher {
    fn default() -> Self {
        Self::new()
    }
}

/// Quantum for floating-point inputs (values are rounded to 1e-9 first).
const QUANTUM: f64 = 1e9;

impl Hasher {
    pub fn new() -> Hasher {
        Hasher(blake3::Hasher::new())
    }
    /// A domain tag, so equal data of different kinds hashes differently.
    pub fn tag(&mut self, t: &str) -> &mut Self {
        self.str(t)
    }
    pub fn bytes(&mut self, b: &[u8]) -> &mut Self {
        self.0.update(&(b.len() as u64).to_le_bytes());
        self.0.update(b);
        self
    }
    pub fn str(&mut self, s: &str) -> &mut Self {
        self.bytes(s.as_bytes())
    }
    pub fn u64(&mut self, v: u64) -> &mut Self {
        self.0.update(&v.to_le_bytes());
        self
    }
    /// A float rounded to 1e-9 (NaN and infinities hash by class).
    pub fn f64(&mut self, v: f64) -> &mut Self {
        let q = if v.is_finite() {
            let r = (v * QUANTUM).round();
            // -0 and +0 are the same value
            if r == 0.0 {
                0i64
            } else {
                r as i64
            }
        } else if v.is_nan() {
            i64::MIN
        } else if v > 0.0 {
            i64::MAX
        } else {
            i64::MIN + 1
        };
        self.0.update(&q.to_le_bytes());
        self
    }
    pub fn f64s(&mut self, v: &[f64]) -> &mut Self {
        self.u64(v.len() as u64);
        for &x in v {
            self.f64(x);
        }
        self
    }
    pub fn finish(&self) -> Hash32 {
        Hash32(*self.0.finalize().as_bytes())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hex_roundtrip_and_quantization() {
        let a = Hasher::new().f64s(&[0.1, 0.2, -0.0]).finish();
        let b = Hasher::new().f64s(&[0.1 + 1e-14, 0.2, 0.0]).finish();
        assert_eq!(a, b);
        assert_ne!(a, Hasher::new().f64s(&[0.1, 0.2000001, 0.0]).finish());
        assert_eq!(Hash32::from_hex(&a.hex()), Some(a));
    }
}
