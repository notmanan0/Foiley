//! F003: flap-deflected polars.
mod common;
use common::*;
use foiley_foil::FoilPolar;
use foiley_plane::analysis::*;
use foiley_plane::FoilDb;

/// NACA 0012 with a 75 % chord hinge, used on the whole span.
fn flapped() -> (foiley_plane::geom::Plane, FoilDb) {
    let mut f = naca0012();
    f.hinge = Some((0.75, 0.5));
    (review_rect(&f), db_with(&f))
}

fn polar_at(foil: &str, flap: f64, cd: f64) -> Vec<FoilPolar> {
    [1e5, 1e7].iter().map(|&re| FoilPolar { flap_deg: flap, ..synth_polar(foil, &format!("Re{re} flap{flap}"), re, 9.0, cd) }).collect()
}

fn def(flaps: &[((usize, usize), f64)]) -> PlanePolarDef {
    PlanePolarDef { viscous: true, free_wake: 0, flaps: flaps.to_vec(), ..review_def() }
}

#[test]
fn f003_deflected_strip_uses_deflected_polar() {
    let (p, mut db) = flapped();
    let name = p.foil_names()[0].clone();
    for (flap, cd) in [(0.0, 0.01), (10.0, 0.05), (20.0, 0.09)] {
        for q in polar_at(&name, flap, cd) {
            db.add_polar(&q);
        }
    }
    // deflected families are not ambiguous with the undeflected one
    assert_eq!(db.families_of(&name).len(), 3);
    // right flap (index 1) down 10 deg: right strips read the 10 deg family
    let o = run1(&p, &db, &def(&[((0, 1), 10.0)]));
    assert!(!o.strips.is_empty());
    for s in &o.strips {
        let want = if s.span_pos > 0.0 { 0.05 } else { 0.01 };
        assert!((s.cdv - want).abs() < 1e-12, "strip at y={:.2}: cdv {} (want {want})", s.span_pos, s.cdv);
        assert!(!s.flap_clamped);
    }
    assert!(!o.warnings.iter().any(|w| w.contains("flap angle")), "{:?}", o.warnings);
    // between two deflected families: linear in flap angle
    let o = run1(&p, &db, &def(&[((0, 0), 15.0), ((0, 1), 15.0)]));
    for s in &o.strips {
        assert!((s.cdv - 0.07).abs() < 1e-12, "15 deg strip: cdv {}", s.cdv);
    }
    // the control polar's flap range goes through the same lookup
    let t6 = PlanePolarDef { kind: PolarKind::Control, range: Range { from: 1.0, to: 1.0, step: 1.0 }, ranges: vec![ControlRange { target: ControlTarget::Flap { wing: 0, index: 1 }, min: 0.0, max: 20.0 }], ..def(&[]) };
    let o = run1(&p, &db, &t6);
    let right: Vec<f64> = o.strips.iter().filter(|s| s.span_pos > 0.0).map(|s| s.cdv).collect();
    assert!(!right.is_empty() && right.iter().all(|c| (c - 0.09).abs() < 1e-12), "T6 flap 20: {right:?}");
}

#[test]
fn f003_strict_fails_without_family() {
    let (p, mut db) = flapped();
    let name = p.foil_names()[0].clone();
    for q in polar_at(&name, 0.0, 0.01) {
        db.add_polar(&q);
    }
    // default: warn, mark the strips, use the undeflected family
    let o = run1(&p, &db, &def(&[((0, 1), 10.0)]));
    let n = o.strips.iter().filter(|s| s.flap_clamped).count();
    assert!(n > 0 && o.strips.iter().all(|s| s.flap_clamped == (s.span_pos > 0.0)), "only the deflected strips are flagged");
    let w = o.warnings.iter().find(|w| w.contains("flap angle")).expect("warning about the missing deflected family");
    assert!(w.starts_with(&n.to_string()), "{w}");
    // strict: the point fails, naming the cause
    let strict = PlanePolarDef { strict_envelope: true, ..def(&[((0, 1), 10.0)]) };
    let pol = run_polar1(&p, &db, &strict);
    assert!(pol.points.is_empty(), "strict mode kept the point");
    assert!(pol.failed.iter().any(|f| f.contains("flap angle")), "{:?}", pol.failed);
    // beyond the deflected families is flagged too
    for q in polar_at(&name, 5.0, 0.03) {
        db.add_polar(&q);
    }
    let o = run1(&p, &db, &def(&[((0, 1), 10.0)]));
    assert!(o.strips.iter().filter(|s| s.span_pos > 0.0).all(|s| s.flap_clamped && (s.cdv - 0.03).abs() < 1e-12));
    // and an undeflected run is unaffected
    let o = run1(&p, &db, &PlanePolarDef { strict_envelope: true, ..def(&[]) });
    assert!(o.strips.iter().all(|s| !s.flap_clamped && (s.cdv - 0.01).abs() < 1e-12));
}
