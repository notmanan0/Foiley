//! F008: LLT glide balance.
mod common;
use common::*;
use foiley_plane::analysis::*;

#[test]
fn f008_llt_glide_resultant_equals_weight() {
    let f = naca0012();
    let p = review_rect(&f);
    let mut db = db_with(&f);
    for re in [1e4, 1e8] {
        db.add_polar(&synth_polar(&f.name, &format!("Re{re}"), re, 9.0, 0.3));
    }
    let llt = run1(&p, &db, &PlanePolarDef { method: Method::Llt, kind: PolarKind::Glide, free_wake: 0, viscous: true, ..review_def() });
    let w = review_inertia().mass * G;
    let r = llt.lift.hypot(llt.drag) / w;
    assert!((r - 1.0).abs() < 1e-4, "resultant/weight = {r}");
    let res = llt.balance_residual.expect("LLT T3 point carries its balance residual");
    assert!(res < GLIDE_TOL, "residual {res}");
    // the drag-dependent Re is consistent: Re-dependent drag still balances
    let mut db2 = db_with(&f);
    for re in [1e5, 1e6] {
        db2.add_polar(&synth_polar(&f.name, &format!("Re{re}"), re, 9.0, 0.02 + 2e3 / re));
    }
    let o = run1(&p, &db2, &PlanePolarDef { method: Method::Llt, kind: PolarKind::Glide, free_wake: 0, viscous: true, ..review_def() });
    assert!((o.lift.hypot(o.drag) / w - 1.0).abs() < 1e-4 && o.balance_residual.unwrap() < GLIDE_TOL);
    // T2 keeps the lift-only balance
    let t2 = run1(&p, &db, &PlanePolarDef { method: Method::Llt, kind: PolarKind::FixedLift, free_wake: 0, viscous: true, ..review_def() });
    assert!((t2.lift / w - 1.0).abs() < 1e-4, "T2 lift/weight = {}", t2.lift / w);
    assert!(t2.balance_residual.is_none());
}
