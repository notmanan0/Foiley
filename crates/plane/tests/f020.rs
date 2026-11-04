//! F020: the T3 glide speed balance reports convergence and its residual.
mod common;
use common::*;
use foiley_plane::analysis::*;

fn rigid_t3() -> PlanePolarDef {
    PlanePolarDef { kind: PolarKind::Glide, free_wake: 0, ..review_def() }
}

#[test]
fn f020_t3_converged_reports_residual() {
    let f = naca0012();
    let (p, mut db) = (review_rect(&f), db_with(&f));
    for re in [1e5, 1e6] {
        db.add_polar(&synth_polar(&f.name, &format!("Re{re}"), re, 9.0, 0.01 + 1e3 / re));
    }
    let pol = run_polar1(&p, &db, &PlanePolarDef { viscous: true, ..rigid_t3() });
    assert!(pol.failed.is_empty(), "{:?}", pol.failed);
    let o = &pol.points[0];
    let r = o.balance_residual.expect("T3 point carries its residual");
    assert!(r < GLIDE_TOL, "residual {r}");
    let w = review_inertia().mass * G;
    assert!((o.lift.hypot(o.drag) / w - 1.0).abs() < GLIDE_TOL);
}

#[test]
fn f020_t3_unconverged_fails() {
    // Drag that steps up sharply across the balance Reynolds number: the
    // fixed-point speed iteration flips between the two sides and never settles.
    let f = naca0012();
    let (p, mut db) = (review_rect(&f), db_with(&f));
    let t2 = run1(&p, &db, &PlanePolarDef { kind: PolarKind::FixedLift, free_wake: 0, ..review_def() });
    let re0 = t2.velocity * 1.0 / PlanePolarDef::default().fluid.nu();
    db.add_polar(&synth_polar(&f.name, "lo", 0.97 * re0, 9.0, 0.0));
    db.add_polar(&synth_polar(&f.name, "hi", 1.03 * re0, 9.0, 2.0));
    let pol = run_polar1(&p, &db, &PlanePolarDef { viscous: true, ..rigid_t3() });
    assert!(pol.points.is_empty(), "unconverged glide reported as a point: {:?}", pol.points.first().map(|o| (o.velocity, o.balance_residual)));
    assert_eq!(pol.failed.len(), 1);
    assert!(pol.failed[0].contains("did not converge"), "{}", pol.failed[0]);
}
