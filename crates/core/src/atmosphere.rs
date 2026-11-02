//! International Standard Atmosphere (troposphere + lower stratosphere).

use serde::{Deserialize, Serialize};

/// Fluid properties.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Fluid {
    /// density, kg/m^3
    pub rho: f64,
    /// dynamic viscosity, Pa.s
    pub mu: f64,
    /// speed of sound, m/s
    pub a: f64,
}

impl Default for Fluid {
    fn default() -> Self {
        let a = isa(0.0, 0.0);
        Fluid { rho: a.rho, mu: a.mu, a: a.a }
    }
}

impl Fluid {
    /// Kinematic viscosity, m^2/s.
    pub fn nu(&self) -> f64 {
        self.mu / self.rho
    }
    pub fn reynolds(&self, v: f64, chord: f64) -> f64 {
        v * chord / self.nu()
    }
    pub fn mach(&self, v: f64) -> f64 {
        v / self.a
    }
    /// Fresh water, 15 C (hydrofoils).
    pub fn water() -> Self {
        Fluid { rho: 999.1, mu: 1.138e-3, a: 1481.0 }
    }
    /// Standard air at altitude `h` metres.
    pub fn air_at(h: f64) -> Self {
        let a = isa(h, 0.0);
        Fluid { rho: a.rho, mu: a.mu, a: a.a }
    }
}

/// Atmosphere state at an altitude.
#[derive(Clone, Copy, Debug)]
pub struct Air {
    pub rho: f64,
    pub mu: f64,
    pub a: f64,
    /// temperature, K
    pub t: f64,
    /// pressure, Pa
    pub p: f64,
}

/// ISA properties at geometric altitude `h` metres (valid -1..20 km),
/// with a temperature offset `dt` (K) from standard.
pub fn isa(h: f64, dt: f64) -> Air {
    const G: f64 = 9.80665;
    const R: f64 = 287.052_87;
    const T0: f64 = 288.15;
    const P0: f64 = 101_325.0;
    const L: f64 = 0.0065;
    let h = h.clamp(-1000.0, 20_000.0);
    let (ts, p) = if h <= 11_000.0 {
        let t = T0 - L * h;
        (t, P0 * (t / T0).powf(G / (L * R)))
    } else {
        let t11 = T0 - L * 11_000.0;
        let p11 = P0 * (t11 / T0).powf(G / (L * R));
        (t11, p11 * (-G * (h - 11_000.0) / (R * t11)).exp())
    };
    let t = ts + dt;
    let rho = p / (R * t);
    let mu = 1.458e-6 * t.powf(1.5) / (t + 110.4); // Sutherland
    let a = (1.4 * R * t).sqrt();
    Air { rho, mu, a, t, p }
}

#[cfg(test)]
mod tests {
    #[test]
    fn sea_level() {
        let a = super::isa(0.0, 0.0);
        assert!((a.rho - 1.225).abs() < 1e-3);
        assert!((a.mu - 1.789e-5).abs() < 1e-7);
        assert!((a.a - 340.29).abs() < 0.1);
    }
}
