//! F008: LLT glide balance.
mod common;
use common::*;
use foiley_plane::analysis::*;

#[test]
#[ignore = "F008: LLT T3 balances CL only"]
fn f008_llt_glide_resultant_equals_weight() {
    let f = naca0012();
    let p = review_rect(&f);
    let mut db = db_with(&f);
    for re in [1e4, 1e8] {
        db.add_polar(&synth_polar(&f.name, &format!("Re{re}"), re, 9.0, 0.3));
    }
    let llt = run1(&p, &db, &PlanePolarDef { method: Method::Llt, kind: PolarKind::Glide, free_wake: 0, viscous: true, ..review_def() });
    let r = llt.lift.hypot(llt.drag) / (review_inertia().mass * G);
    assert!((r - 1.0).abs() < 1e-4, "resultant/weight = {r}");
}
