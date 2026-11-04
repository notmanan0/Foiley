//! Fixtures shared by the per-finding test files (docs/REMEDIATION_PLAN.md).
#![allow(dead_code)]
use foiley_core::{Progress, Vec3};
use foiley_foil::{Foil, FoilPoint, FoilPolar};
use foiley_plane::analysis::*;
use foiley_plane::geom::*;
use foiley_plane::FoilDb;

/// Root of the repository's tracked test data.
pub fn data_dir() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data")
}

pub fn naca0012() -> Foil {
    Foil::naca("0012").unwrap()
}

pub fn db_with(f: &Foil) -> FoilDb {
    let mut db = FoilDb::new();
    db.add_foil(f);
    db
}

/// The 8 m x 1 m NACA 0012 rectangle of `examples/benchmark_review.rs`.
pub fn review_rect(f: &Foil) -> Plane {
    let mut root = WingSection::new(0.0, 1.0, 0.0, &f.name);
    root.nx = 3;
    root.ny = 8;
    let tip = WingSection::new(4.0, 1.0, 0.0, &f.name);
    Plane {
        name: "review rectangle".into(),
        description: String::new(),
        wings: vec![Wing::new("wing", WingKind::Main, vec![root, tip])],
        fuses: vec![],
        point_masses: vec![],
    }
}

pub fn review_inertia() -> Inertia {
    Inertia { mass: 10.0, cog: Vec3::new(0.25, 0.0, 0.0), ixx: 10.0, iyy: 10.0, izz: 10.0, ..Default::default() }
}

/// Inviscid VLM2, relaxed wake, single point at alpha = 8 deg.
pub fn review_def() -> PlanePolarDef {
    PlanePolarDef {
        viscous: false,
        free_wake: 4,
        wake_length: 10.0,
        inertia: Some(review_inertia()),
        alpha: 8.0,
        range: Range { from: 8.0, to: 8.0, step: 1.0 },
        ..Default::default()
    }
}

pub fn run_polar1(p: &Plane, db: &FoilDb, d: &PlanePolarDef) -> PlanePolar {
    run_polar(p, db, d, &Progress::none(), &|_| {}).unwrap()
}

pub fn run1(p: &Plane, db: &FoilDb, d: &PlanePolarDef) -> PlaneOpp {
    run_polar1(p, db, d).points.remove(0)
}

/// Synthetic linear polar: cl = 0.1 alpha, constant cd.
pub fn synth_polar(foil: &str, name: &str, re: f64, ncrit: f64, cd: f64) -> FoilPolar {
    FoilPolar {
        foil: foil.into(),
        name: name.into(),
        re,
        ncrit,
        polar_type: 1,
        viscous: true,
        points: (-20..=20).map(|a| FoilPoint { alpha: a as f64, cl: 0.1 * a as f64, cd, re, ..Default::default() }).collect(),
        ..Default::default()
    }
}
