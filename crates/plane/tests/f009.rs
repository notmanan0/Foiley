//! F009: input validation and finite output.
mod common;
use common::*;
use foiley_plane::analysis::*;

#[test]
#[ignore = "F009: a T8 point at zero speed passes validation"]
fn f009_t8_zero_speed_rejected() {
    let d = PlanePolarDef { kind: PolarKind::List, free_wake: 0, list: vec![(5.0, 0.0, 0.0)], ..review_def() };
    assert!(d.validate().is_err());
}
