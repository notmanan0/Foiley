//! Writes the project fixtures in `tests/data/projects/` (F016): every
//! native example plane, with 2D and 3D polars, saved by the current code.
//!
//!     cargo run --release -p foiley-app --example project_fixtures -- <dir> [suffix]
//!
//! The `_v1` files in the directory were written by this generator at
//! commit 98515c9 (before format 2) with the suffix `_v1`; the others by
//! the current code with no suffix.  `INVENTORY.txt` lists them all.

use foiley_core::Progress;
use foiley_foil::{analyze, AnalysisSettings, Foil, Sweep};
use foiley_plane::analysis::{PolarKind, Range};
use foiley_plane::{run_polar, Method, Plane, PlanePolarDef};
use foiley_app::Project;
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// A 2D polar of every foil at two Reynolds numbers (cheap sweep).
fn add_foil_polars(p: &mut Project) {
    let foils: Vec<Foil> = p.foils.clone();
    for f in &foils {
        for re in [2e5, 6e5] {
            let st = AnalysisSettings { re, ..Default::default() };
            let (polar, _) = analyze(f, &st, &Sweep::Alpha { from: -4.0, to: 10.0, step: 2.0 }, &Progress::none()).expect("2D polar");
            p.add_foil_polar(polar);
        }
    }
}

/// One viscous VLM2 T1 polar on every plane.
fn add_plane_polars(p: &mut Project) {
    let db = p.foil_db();
    let planes: Vec<Plane> = p.planes.clone();
    for pl in &planes {
        let def = PlanePolarDef { kind: PolarKind::FixedSpeed, method: Method::Vlm2, velocity: 15.0, range: Range { from: 0.0, to: 6.0, step: 3.0 }, keep_cp: false, ..Default::default() };
        let polar = run_polar(pl, &db, &def, &Progress::none(), &|_| {}).expect("3D polar");
        p.add_plane_polar(polar);
    }
}

fn starter() -> Project {
    let mut p = Project::starter();
    p.name = "starter".into();
    add_foil_polars(&mut p);
    add_plane_polars(&mut p);
    p
}

fn glider_xml() -> Project {
    let mut p = Project::new("glider.xml");
    p.import(&root().join("docs/examples/glider.xml")).expect("glider.xml");
    add_foil_polars(&mut p);
    add_plane_polars(&mut p);
    p
}

/// The 8 m x 1 m NACA 0012 rectangle of `benchmark_review.rs`, geometry only.
fn review_rect() -> Project {
    let mut p = Project::new("review rectangle");
    let f = Foil::naca("0012").unwrap();
    let mut root = foiley_plane::WingSection::new(0.0, 1.0, 0.0, &f.name);
    root.nx = 3;
    root.ny = 8;
    let tip = foiley_plane::WingSection::new(4.0, 1.0, 0.0, &f.name);
    p.planes.push(Plane {
        name: "review rectangle".into(),
        description: String::new(),
        wings: vec![foiley_plane::Wing::new("wing", foiley_plane::WingKind::Main, vec![root, tip])],
        fuses: vec![],
        point_masses: vec![],
    });
    p.foils.push(f);
    p
}

/// The default glider with a hinged elevator foil and deflected 2D polar
/// families at -10 and +10 deg (format 2 only: needs FoilPolar::flap_deg).
fn flapped() -> Project {
    let mut p = Project::starter();
    p.name = "flapped elevator".into();
    let mut ef = Foil::naca("0009").unwrap();
    ef.name = "elevator 0009".into();
    ef.hinge = Some((0.7, 0.5));
    for w in p.planes[0].wings.iter_mut().filter(|w| w.kind == foiley_plane::WingKind::Elevator) {
        for s in &mut w.sections {
            s.foil_left = ef.name.clone();
            s.foil_right = ef.name.clone();
        }
    }
    p.foils.push(ef.clone());
    add_foil_polars(&mut p);
    for d in [-10.0, 10.0] {
        let g = ef.with_flap(0.7, 0.5, d).unwrap();
        for re in [2e5, 6e5] {
            let st = AnalysisSettings { re, ..Default::default() };
            let (mut polar, _) = analyze(&g, &st, &Sweep::Alpha { from: -4.0, to: 8.0, step: 2.0 }, &Progress::none()).expect("deflected 2D polar");
            polar.name = format!("{} flap {d:+}", polar.name);
            polar.foil = ef.name.clone();
            polar.flap_deg = d;
            polar.foil_hash = Some(ef.geom_hash());
            p.add_foil_polar(polar);
        }
    }
    add_plane_polars(&mut p);
    p
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let dir = PathBuf::from(args.get(1).map(|s| s.as_str()).unwrap_or("tests/data/projects"));
    let suffix = args.get(2).cloned().unwrap_or_default();
    std::fs::create_dir_all(&dir).unwrap();
    let mut written = vec![];
    type Case = (&'static str, fn() -> Project);
    let mut cases: Vec<Case> = vec![("starter", starter), ("glider_xml", glider_xml)];
    if suffix.is_empty() {
        cases.push(("review_rect", review_rect));
        cases.push(("flapped", flapped));
    }
    for (name, make) in cases {
        let path = dir.join(format!("{name}{suffix}.foiley"));
        make().save(&path).unwrap();
        written.push(path.display().to_string());
    }
    for w in written {
        println!("wrote {w}");
    }
}
