//! F009: input validation and finite output.
mod common;
use common::*;
use foiley_core::Progress;
use foiley_plane::analysis::*;
use foiley_plane::geom::Plane;

#[test]
fn f009_t8_zero_speed_rejected() {
    let d = PlanePolarDef { kind: PolarKind::List, free_wake: 0, list: vec![(5.0, 0.0, 0.0)], ..review_def() };
    let e = d.validate().unwrap_err();
    assert!(e.contains("T8 point 1"), "{e}");
    // and run_polar refuses it instead of producing a NaN row
    let f = naca0012();
    assert!(run_polar(&review_rect(&f), &db_with(&f), &d, &Progress::none(), &|_| {}).is_err());
    let nan = PlanePolarDef { list: vec![(5.0, 0.0, 10.0), (f64::NAN, 0.0, 10.0)], ..d.clone() };
    assert!(nan.validate().unwrap_err().contains("T8 point 2"));
    let ok = PlanePolarDef { list: vec![(5.0, 0.0, 10.0)], ..d };
    ok.validate().unwrap();
}

#[test]
fn f009_t6_range_crossing_zero_rejected() {
    let speed = |min, max| PlanePolarDef {
        kind: PolarKind::Control,
        free_wake: 0,
        range: Range { from: 0.0, to: 1.0, step: 0.5 },
        ranges: vec![ControlRange { target: ControlTarget::Velocity, min, max }],
        ..review_def()
    };
    assert!(speed(-5.0, 5.0).validate().is_err());
    assert!(speed(0.0, 10.0).validate().is_err());
    assert!(speed(10.0, 0.0).validate().is_err());
    speed(5.0, 15.0).validate().unwrap();
    // T4 range through zero
    let t4 = PlanePolarDef { kind: PolarKind::FixedAlpha, free_wake: 0, range: Range { from: -2.0, to: 10.0, step: 2.0 }, ..review_def() };
    assert!(t4.validate().is_err());
}

#[test]
fn f009_nonfinite_scalars_rejected() {
    let base = PlanePolarDef { free_wake: 0, ..review_def() };
    let cases: Vec<(&str, PlanePolarDef)> = vec![
        ("speed", PlanePolarDef { velocity: f64::NAN, ..base.clone() }),
        ("alpha", PlanePolarDef { alpha: f64::INFINITY, ..base.clone() }),
        ("range", PlanePolarDef { range: Range { from: 0.0, to: f64::NAN, step: 1.0 }, ..base.clone() }),
        ("density", PlanePolarDef { fluid: foiley_core::atmosphere::Fluid { rho: f64::NAN, ..base.fluid }, ..base.clone() }),
        ("reference", PlanePolarDef { reference: RefDims::Custom { area: f64::NAN, span: 1.0, chord: 1.0 }, ..base.clone() }),
        ("mass", PlanePolarDef { inertia: Some(foiley_plane::geom::Inertia { mass: f64::NAN, ..review_inertia() }), ..base.clone() }),
        ("extra drag", PlanePolarDef { extra_drag: vec![ExtraDrag { name: "x".into(), area: 1.0, cd: f64::NAN }], ..base.clone() }),
        ("ground", PlanePolarDef { ground_height: Some(f64::NAN), ..base.clone() }),
    ];
    for (what, d) in cases {
        assert!(d.validate().is_err(), "nonfinite {what} accepted");
    }
    base.validate().unwrap();
}

#[test]
fn f009_unsupported_settings_rejected() {
    // free wake: VLM only (every VLM driver, T6 and T7 included, since F007)
    let fw = PlanePolarDef { method: Method::Panel, ..review_def() };
    assert!(fw.validate().unwrap_err().contains("free-wake"));
    let t6 = PlanePolarDef { kind: PolarKind::Control, range: Range { from: 0.0, to: 0.0, step: 1.0 }, ..review_def() };
    t6.validate().unwrap();
    PlanePolarDef { kind: PolarKind::Stability, ..review_def() }.validate().unwrap();
    let llt5 = PlanePolarDef { method: Method::Llt, kind: PolarKind::Sideslip, free_wake: 0, ..review_def() };
    assert!(llt5.validate().is_err());
    // clear_unsupported makes any definition valid for its method
    for m in Method::ALL {
        for k in PolarKind::ALL {
            if !m.supports(k, Setting::PolarType) {
                continue;
            }
            let mut d = PlanePolarDef { method: m, kind: k, beta: 2.0, ground_height: Some(1.0), fuse_drag: true, free_wake: 3, list: vec![(2.0, 1.0, 10.0)], ..review_def() };
            if k == PolarKind::FixedAlpha {
                d.range = Range { from: 5.0, to: 10.0, step: 5.0 };
            }
            d.clear_unsupported();
            assert!(d.unsupported().is_empty(), "{m:?} {k:?}: {:?}", d.unsupported());
            d.validate().unwrap_or_else(|e| panic!("{m:?} {k:?}: {e}"));
        }
    }
}

/// Every point any driver returns has finite coefficients; anything else is
/// in `failed` with a reason.
#[test]
fn g06_no_nonfinite_points() {
    let f = naca0012();
    let f9 = foiley_foil::Foil::naca("0009").unwrap();
    let mut db = db_with(&f);
    db.add_foil(&f9);
    for name in [&f.name, &f9.name] {
        for re in [1e4, 1e7] {
            db.add_polar(&synth_polar(name, &format!("Re{re}"), re, 9.0, 0.01)); // for LLT
        }
    }
    let glider = Plane::default_plane(&f.name, &f9.name);
    let rect = review_rect(&f);
    let one = |a: f64| Range { from: a, to: a, step: 1.0 };
    let defs = vec![
        PlanePolarDef { viscous: false, range: Range { from: -4.0, to: 8.0, step: 4.0 }, ..Default::default() },
        PlanePolarDef { viscous: false, kind: PolarKind::FixedLift, inertia: Some(review_inertia()), range: Range { from: -4.0, to: 8.0, step: 4.0 }, ..Default::default() },
        PlanePolarDef { viscous: false, kind: PolarKind::Glide, inertia: Some(review_inertia()), range: Range { from: -4.0, to: 8.0, step: 4.0 }, ..Default::default() },
        PlanePolarDef { viscous: false, kind: PolarKind::FixedAlpha, alpha: 4.0, range: Range { from: 5.0, to: 25.0, step: 10.0 }, ..Default::default() },
        PlanePolarDef { viscous: false, kind: PolarKind::Sideslip, alpha: 4.0, range: Range { from: -5.0, to: 5.0, step: 5.0 }, ..Default::default() },
        PlanePolarDef { viscous: false, kind: PolarKind::List, list: vec![(2.0, 0.0, 10.0), (4.0, 3.0, 20.0)], ..Default::default() },
        PlanePolarDef { viscous: false, kind: PolarKind::Control, alpha: 3.0, range: Range { from: 0.0, to: 1.0, step: 0.5 }, ranges: vec![ControlRange { target: ControlTarget::Velocity, min: 5.0, max: 30.0 }], ..Default::default() },
        PlanePolarDef { viscous: false, kind: PolarKind::Stability, range: one(0.0), ..Default::default() },
        PlanePolarDef { viscous: false, method: Method::Llt, range: Range { from: -2.0, to: 6.0, step: 4.0 }, ..Default::default() },
    ];
    for d in &defs {
        for p in [&glider, &rect] {
            let d = if d.kind.needs_mass() && p.name == rect.name { PlanePolarDef { inertia: Some(review_inertia()), ..d.clone() } } else { d.clone() };
            d.validate().unwrap();
            let pol = match run_polar(p, &db, &d, &Progress::none(), &|_| {}) {
                Ok(pol) => pol,
                Err(e) => panic!("{} {:?} {:?}: {e}", p.name, d.method, d.kind),
            };
            for o in &pol.points {
                assert!(nonfinite_reason(o).is_none(), "{} {:?} {:?}: {:?}", p.name, d.method, d.kind, nonfinite_reason(o));
            }
            assert!(pol.failed.iter().all(|f| !f.is_empty()));
        }
    }
}

/// The GUI and the CLI decide what is allowed through `Method::supports`
/// (via `validate` / `uses`), not their own lists.
#[test]
fn f009_cli_and_gui_share_supports_table() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let gui = std::fs::read_to_string(root.join("crates/gui/src/plane_ui.rs")).unwrap();
    let cli = std::fs::read_to_string(root.join("crates/cli/src/main.rs")).unwrap();
    assert!(gui.contains(".supports("), "GUI does not consult Method::supports");
    assert!(gui.contains("clear_unsupported()"), "GUI does not clear unsupported settings");
    assert!(gui.contains("def.validate()"), "GUI does not validate before running");
    assert!(cli.contains("def.validate()?"), "CLI does not validate the definition");
    // and neither keeps a private method/kind rule
    for (name, src) in [("GUI", &gui), ("CLI", &cli)] {
        assert!(!src.contains("supports polar types 1-4"), "{name} has its own LLT rule");
    }
    // the table itself: validate agrees with supports for every pair
    for m in Method::ALL {
        for k in PolarKind::ALL {
            let mut d = PlanePolarDef { method: m, kind: k, free_wake: 0, inertia: Some(review_inertia()), list: vec![(2.0, 0.0, 10.0)], ..Default::default() };
            if k == PolarKind::FixedAlpha {
                d.range = Range { from: 5.0, to: 10.0, step: 5.0 };
            }
            assert_eq!(d.validate().is_ok(), m.supports(k, Setting::PolarType), "{m:?} {k:?}: {:?}", d.validate());
        }
    }
}
