//! 3D layer of Foiley: planes made of wings and fuselages, meshing, and the
//! potential-flow solvers.
//!
//! * [`geom`]: plane / wing / fuselage definitions, reference dimensions,
//!   inertia.
//! * [`mesh`]: panel meshes (thin mid-surface for VLM, closed surfaces for
//!   the panel method), strips, flaps.
//! * [`llt`]: nonlinear lifting-line theory using the 2D polars.
//! * [`vlm`]: vortex-lattice (horseshoe VLM1 and ring VLM2).
//! * [`panel`]: constant source/doublet panel method (Dirichlet).
//! * [`analysis`]: polar definitions (types 1-8) and the drivers.
//! * [`stability`]: stability derivatives, state matrices, modes, time response.
//! * [`xml`]: import of xflr5 / flow5 plane XML files.

pub mod analysis;
pub mod cad;
pub mod export;
pub mod foildb;
pub mod geom;
pub mod llt;
pub mod mesh;
pub mod panel;
pub mod stability;
pub mod vlm;
pub mod xml;

pub use analysis::{run_polar, Method, PlaneOpp, PlanePolar, PlanePolarDef, PolarKind};
pub use foildb::FoilDb;
pub use geom::{Fuse, NurbsBody, NurbsFrame, Plane, PointMass, Wing, WingKind, WingSection};

/// Error type for the 3D layer.
#[derive(Debug, Clone, PartialEq)]
pub enum PlaneError {
    Geometry(String),
    MissingFoil(String),
    Solver(String),
    Input(String),
    Cancelled,
}

impl std::fmt::Display for PlaneError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PlaneError::Geometry(s) => write!(f, "geometry error: {s}"),
            PlaneError::MissingFoil(s) => write!(f, "foil '{s}' is used by the plane but is not loaded"),
            PlaneError::Solver(s) => write!(f, "solver error: {s}"),
            PlaneError::Input(s) => write!(f, "invalid analysis settings: {s}"),
            PlaneError::Cancelled => write!(f, "cancelled by user"),
        }
    }
}
impl std::error::Error for PlaneError {}

impl From<foiley_core::Cancelled> for PlaneError {
    fn from(_: foiley_core::Cancelled) -> Self {
        PlaneError::Cancelled
    }
}
impl From<foiley_core::linalg::LinalgError> for PlaneError {
    fn from(e: foiley_core::linalg::LinalgError) -> Self {
        PlaneError::Solver(e.to_string())
    }
}
