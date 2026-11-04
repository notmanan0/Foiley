//! F012: viscous interpolation reports which envelope limit it hit.
mod common;
use common::*;
use foiley_plane::analysis::*;
use foiley_plane::foildb::Envelope;

fn db_res(res: &[f64]) -> (foiley_plane::geom::Plane, foiley_plane::FoilDb) {
    let f = naca0012();
    let mut db = db_with(&f);
    for &re in res {
        db.add_polar(&synth_polar(&f.name, &format!("Re{re}"), re, 9.0, 0.01));
    }
    (review_rect(&f), db)
}

fn viscous() -> PlanePolarDef {
    PlanePolarDef { viscous: true, free_wake: 0, ..review_def() }
}

#[test]
fn f012_envelope_kinds_reported() {
    let (_, db) = db_res(&[1e5, 1e6]);
    let name = &naca0012().name;
    assert_eq!(db.at_re_cl(name, 5e5, 0.5).unwrap().1, Envelope::Inside);
    assert_eq!(db.at_re_cl(name, 1e4, 0.5).unwrap().1, Envelope::ClampedRe);
    assert_eq!(db.at_re_cl(name, 5e7, 0.5).unwrap().1, Envelope::ClampedRe);
    assert_eq!(db.at_re_cl(name, 5e5, 5.0).unwrap().1, Envelope::ClampedCl);
    assert_eq!(db.at_re_cl(name, 1e4, 5.0).unwrap().1, Envelope::ClampedBoth);
    assert_eq!(db.at_re_alpha(name, 5e5, 40.0).unwrap().1, Envelope::ClampedCl);
    assert!(db.at_re_cl("no such foil", 5e5, 0.5).is_none());

    // in a 3D point: flight Re ~ 7e5 against polars far above it
    let (p, hi) = db_res(&[1e8, 2e8]);
    let o = run1(&p, &hi, &viscous());
    assert!(o.strips.iter().all(|s| Envelope::from_u8(s.envelope) == Envelope::ClampedRe && s.out_of_polar));
    let w = o.warnings.iter().find(|w| w.contains("clamped")).expect("clamp warning");
    assert!(w.contains(&format!("{} in Re", o.strips.len())), "{w}");
    // inside: no warning, every strip Inside
    let (p, ok) = db_res(&[1e5, 1e7]);
    let o = run1(&p, &ok, &viscous());
    assert!(o.strips.iter().all(|s| s.envelope == Envelope::Inside as u8 && !s.out_of_polar));
    assert!(!o.warnings.iter().any(|w| w.contains("clamped")), "{:?}", o.warnings);
    // no polar at all
    let (p, none) = db_res(&[]);
    let o = run1(&p, &none, &viscous());
    assert!(o.strips.iter().all(|s| Envelope::from_u8(s.envelope) == Envelope::Missing));
}

#[test]
fn f012_strict_envelope_fails_point() {
    let strict = PlanePolarDef { strict_envelope: true, ..viscous() };
    let (p, hi) = db_res(&[1e8, 2e8]);
    let pol = run_polar1(&p, &hi, &strict);
    assert!(pol.points.is_empty());
    assert!(pol.failed[0].contains("strict envelope") && pol.failed[0].contains("clamped in Re"), "{:?}", pol.failed);
    // missing polars also fail in strict mode
    let (p, none) = db_res(&[]);
    assert!(run_polar1(&p, &none, &strict).failed[0].contains("no polar"));
    // inside the envelope strict mode changes nothing
    let (p, ok) = db_res(&[1e5, 1e7]);
    let a = run1(&p, &ok, &strict);
    let b = run1(&p, &ok, &viscous());
    assert_eq!(a.cd, b.cd);
    // and LLT honours it too
    let llt = PlanePolarDef { method: Method::Llt, ..strict };
    assert!(run_polar1(&p, &hi, &llt).points.is_empty());
}
