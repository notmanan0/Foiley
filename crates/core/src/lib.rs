//! Shared building blocks for the Foiley solvers:
//!
//! * [`Vec3`] small-vector algebra used by the 3D solvers;
//! * [`Progress`], a cheap, thread-safe progress/cancellation handle that
//!   every long-running computation takes.  The GUI and the CLI poll it from
//!   their own thread, so a running analysis never blocks the user interface;
//! * dense linear algebra ([`linalg`]) built on `faer` (blocked, parallel LU
//!   and nonsymmetric eigenvalues);
//! * interpolation helpers and a standard atmosphere.

pub mod atmosphere;
pub mod hash;
pub mod interp;
pub mod linalg;
pub mod progress;
pub mod vec3;

pub use progress::{Cancelled, LogLevel, Progress};
pub use vec3::Vec3;

pub const DEG: f64 = std::f64::consts::PI / 180.0;

/// Point-distribution laws used for paneling (chordwise and spanwise).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
pub enum Spacing {
    /// equal spacing
    Uniform,
    /// clustered at both ends (best default for wings and chords)
    #[default]
    Cosine,
    /// clustered at the end (e.g. wing tip)
    Sine,
    /// clustered at the start (e.g. wing root)
    MinusSine,
    /// smooth clustering at both ends, gentler than cosine
    Tanh,
}

impl Spacing {
    pub const ALL: [Spacing; 5] = [Spacing::Uniform, Spacing::Cosine, Spacing::Sine, Spacing::MinusSine, Spacing::Tanh];
    pub fn label(self) -> &'static str {
        match self {
            Spacing::Uniform => "Uniform",
            Spacing::Cosine => "Cosine (both ends)",
            Spacing::Sine => "Sine (outer end)",
            Spacing::MinusSine => "-Sine (inner end)",
            Spacing::Tanh => "Tanh (mild, both ends)",
        }
    }
}
