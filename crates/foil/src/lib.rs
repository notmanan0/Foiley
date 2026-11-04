//! 2D airfoil layer of Foiley.
//!
//! * [`Foil`]: airfoil geometry, file I/O, NACA generation, geometric
//!   properties and all XFoil GDES design operations (flap, TE gap, LE
//!   radius, thickness/camber scaling, ...), exposed as plain functions.
//! * [`analysis`]: operating points and polars (types 1-4) solved by the
//!   XFoil engine, with progress reporting, cancellation, panic isolation
//!   and an automatic recovery strategy for unconverged points.
//! * [`batch`]: multithreaded polar meshes (Re x Mach x Ncrit grids) used to
//!   feed the 3D viscous analyses.

pub mod analysis;
pub mod batch;
pub mod foil;
pub mod inverse;
pub mod optim;
pub mod polar;

pub use analysis::{analyze, AnalysisSettings, PointDetail, PolarType, Sweep};
pub use foil::{Foil, FoilProps};
pub use polar::{FoilPoint, FoilPolar, Recovery};

/// Error type for the 2D layer.
#[derive(Debug, Clone, PartialEq)]
pub enum FoilError {
    Io(String),
    Parse(String),
    Geometry(String),
    Cancelled,
}

impl std::fmt::Display for FoilError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FoilError::Io(s) => write!(f, "I/O error: {s}"),
            FoilError::Parse(s) => write!(f, "could not read airfoil: {s}"),
            FoilError::Geometry(s) => write!(f, "invalid airfoil geometry: {s}"),
            FoilError::Cancelled => write!(f, "cancelled by user"),
        }
    }
}
impl std::error::Error for FoilError {}

impl From<foiley_core::Cancelled> for FoilError {
    fn from(_: foiley_core::Cancelled) -> Self {
        FoilError::Cancelled
    }
}
