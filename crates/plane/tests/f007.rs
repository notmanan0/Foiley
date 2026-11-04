//! F007 / F019: one wake-aware evaluator for every driver.
mod common;
use common::*;
use foiley_plane::analysis::*;

#[test]
#[ignore = "F007: T2 ignores the relaxed wake"]
fn f007_t2_lift_equals_weight() {
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    let t2 = run1(&p, &db, &PlanePolarDef { kind: PolarKind::FixedLift, ..review_def() });
    let r = t2.lift / (review_inertia().mass * G);
    assert!((r - 1.0).abs() < 1e-4, "lift/weight = {r}");
}

#[test]
#[ignore = "F007: T3 ignores the relaxed wake"]
fn f007_t3_relaxed_differs_from_rigid() {
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    let t3 = run1(&p, &db, &PlanePolarDef { kind: PolarKind::Glide, ..review_def() });
    let rigid = run1(&p, &db, &PlanePolarDef { kind: PolarKind::Glide, free_wake: 0, ..review_def() });
    assert!((t3.cl - rigid.cl).abs() > 1e-6, "delta_CL = {}", t3.cl - rigid.cl);
    assert!(!t3.wake.is_empty(), "relaxed T3 keeps no wake lines");
}

#[test]
#[ignore = "F019: T6 ignores the relaxed wake that T1 uses (F009 rejects free wake in T6 until F007 routes it)"]
fn f019_t1_t6_equivalent() {
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    let t1 = run1(&p, &db, &review_def());
    let t6 = run1(&p, &db, &PlanePolarDef { kind: PolarKind::Control, range: Range { from: 0.0, to: 0.0, step: 1.0 }, ..review_def() });
    assert!((t1.cl - t6.cl).abs() < 1e-9, "CL T1 {} vs T6 {}", t1.cl, t6.cl);
    assert_eq!(t1.wake.len(), t6.wake.len());
}
