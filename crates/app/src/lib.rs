//! Application layer shared by the Foiley CLI and GUI:
//!
//! * [`project`]: the project file (foils, 2D polars, planes, 3D polars) and
//!   imports / exports;
//! * [`jobs`]: the background job manager.  Every analysis runs on a
//!   dedicated compute pool (one core is left free for the user interface),
//!   is isolated against panics, can be cancelled at any time, streams
//!   partial results, and always ends in a terminal state - so callers never
//!   wait on a solver and nothing can hang.

pub mod auto;
pub mod jobs;
pub mod legacy_io;
pub mod project;

pub use jobs::{CadFormat, CadImport, CadImportResult, InverseTarget, Job, JobId, JobKind, JobManager, JobOutput, JobStatus};
pub use project::Project;
