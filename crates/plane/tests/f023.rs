//! F023: LLT rejects what it would ignore, and has no lateral coefficients.
mod common;
use common::*;
use foiley_plane::analysis::*;

fn llt_db() -> (foiley_plane::geom::Plane, foiley_plane::FoilDb) {
    let f = naca0012();
    let mut db = db_with(&f);
    for re in [1e4, 1e8] {
        db.add_polar(&synth_polar(&f.name, &format!("Re{re}"), re, 9.0, 0.01));
    }
    (review_rect(&f), db)
}

fn llt() -> PlanePolarDef {
    PlanePolarDef { method: Method::Llt, free_wake: 0, ..review_def() }
}

#[test]
fn f023_llt_rejects_unsupported() {
    llt().validate().unwrap();
    let cases: Vec<(&str, PlanePolarDef)> = vec![
        ("beta", PlanePolarDef { beta: 3.0, ..llt() }),
        ("ground", PlanePolarDef { ground_height: Some(1.0), ..llt() }),
        ("flaps", PlanePolarDef { flaps: vec![((0, 0), 5.0)], ..llt() }),
        ("fuse_drag", PlanePolarDef { fuse_drag: true, ..llt() }),
        ("derivatives", PlanePolarDef { derivatives: true, ..llt() }),
        ("free wake", PlanePolarDef { free_wake: 2, ..llt() }),
        ("T5", PlanePolarDef { kind: PolarKind::Sideslip, ..llt() }),
        ("T6", PlanePolarDef { kind: PolarKind::Control, ..llt() }),
        ("T7", PlanePolarDef { kind: PolarKind::Stability, ..llt() }),
        ("T8", PlanePolarDef { kind: PolarKind::List, list: vec![(2.0, 0.0, 10.0)], ..llt() }),
    ];
    let (p, db) = llt_db();
    for (what, d) in cases {
        let e = d.validate().expect_err(what);
        assert!(e.starts_with("LLT does not support"), "{what}: {e}");
        assert!(run_polar(&p, &db, &d, &foiley_core::Progress::none(), &|_| {}).is_err(), "{what} ran");
        // the same setting is fine for VLM (bar the polar types LLT lacks)
        let vlm = PlanePolarDef { method: Method::Vlm2, ..d };
        assert!(vlm.unsupported().is_empty(), "{what}: {:?}", vlm.unsupported());
    }
    // every LLT row of the table is false except the polar types 1-4
    for s in Setting::ALL {
        if s != Setting::PolarType {
            assert!(!Method::Llt.supports(PolarKind::FixedSpeed, s), "{s:?}");
        }
    }
}

#[test]
fn f023_llt_lateral_absent() {
    let (p, db) = llt_db();
    let o = run1(&p, &db, &llt());
    assert_eq!((o.cy, o.cr, o.cn), (None, None, None), "LLT reports lateral coefficients");
    // CSV leaves them empty instead of writing zeros
    let pol = run_polar1(&p, &db, &llt());
    let row: Vec<String> = pol.to_csv().lines().nth(1).unwrap().split(',').map(String::from).collect();
    let header: Vec<String> = pol.to_csv().lines().next().unwrap().split(',').map(String::from).collect();
    for col in ["CY", "Cl_roll", "Cn"] {
        let k = header.iter().position(|h| h == col).unwrap();
        assert_eq!(row[k], "", "{col}");
    }
    // a VLM point has them
    let v = run1(&p, &db, &PlanePolarDef { method: Method::Vlm2, ..llt() });
    assert!(v.cy.is_some() && v.cr.is_some() && v.cn.is_some());
}
