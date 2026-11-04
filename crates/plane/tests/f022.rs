//! F022: pure-Rust STEP B-rep import (trimming, seams, B-splines,
//! watertight tessellation, mesh checks).  Fixtures: tests/data/cad/*.stp,
//! written by tests/data/cad/make_fixtures.py with closed-form areas and
//! volumes.

use foiley_core::Vec3;
use foiley_plane::cad::{Domain, Model, StepOptions};
use std::f64::consts::{PI, TAU};

const CAD_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/cad");

fn text(name: &str) -> String {
    std::fs::read_to_string(std::path::Path::new(CAD_DIR).join(name)).unwrap_or_else(|e| panic!("{name}: {e}"))
}

fn model(name: &str) -> Model {
    Model::from_step(&text(name), &StepOptions::default()).unwrap_or_else(|e| panic!("{name}: {e}"))
}

/// Every loop point maps back onto its 3D vertex; returns the domain.
fn checked_domain(m: &Model, fi: usize) -> Domain {
    let d = m.face_domain(fi).unwrap_or_else(|e| panic!("face {fi}: {e}"));
    let s = &m.faces[fi].surface;
    for l in &d.loops {
        assert_eq!(l.uv.len(), l.ids.len());
        assert!(l.uv.len() >= 3, "degenerate loop");
        for (uv, &id) in l.uv.iter().zip(&l.ids) {
            let r = (s.eval(uv[0], uv[1]) - m.vertex(&d, id)).norm();
            assert!(r < 1e-9, "{}: residual {r:e} at {uv:?}", s.kind());
        }
    }
    d
}

fn face_of(m: &Model, kind: &str) -> usize {
    m.faces.iter().position(|f| f.surface.kind() == kind).unwrap_or_else(|| panic!("no {kind} face"))
}

fn span(d: &Domain) -> [f64; 2] {
    let (lo, hi) = d.loops[0].bbox();
    [hi[0] - lo[0], hi[1] - lo[1]]
}

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * (1.0 + b.abs())
}

/// Edge inversion without PCURVEs at surface seams: the domain of a closed
/// cylinder, a surface of revolution, a cone (apex), a sphere (poles), a
/// torus (two seams), an extrusion, a patch crossing the seam and a band
/// without a seam edge are all simple polygons of the right size, every
/// point maps back to its edge sample, and the seam's two uses lie one
/// period apart.
#[test]
fn f022_seam_inversion_cylinder_and_revolution() {
    // cylinder along x with a tilted seam: [0, 2pi] x [0, 2]
    for (name, pcurves) in [("cyl_seam.stp", false), ("cyl_pcurve.stp", true)] {
        let m = model(name);
        assert_eq!(m.faces.len(), 3, "{name}");
        let fi = face_of(&m, "CYLINDRICAL_SURFACE");
        let d = checked_domain(&m, fi);
        assert_eq!(d.loops.len(), 1);
        let sp = span(&d);
        assert!(close(sp[0], TAU, 1e-12) && close(sp[1], 2.0, 1e-12), "{name}: span {sp:?}");
        assert!(close(d.area(), TAU * 2.0, 1e-12), "{name}: area {}", d.area());
        assert_eq!(d.used_pcurves > 0, pcurves, "{name}: pcurves used {}", d.used_pcurves);
        // the seam vertices appear twice, a period apart
        let seam_edge = m.edges.iter().find(|e| e.curve.kind() == "LINE").unwrap();
        for &(_, vid) in &seam_edge.samples {
            let us: Vec<f64> = d.loops[0].ids.iter().zip(&d.loops[0].uv).filter(|(i, _)| **i == vid).map(|(_, uv)| uv[0]).collect();
            assert!(us.len() >= 2, "seam vertex used {} times", us.len());
            let (lo, hi) = us.iter().fold((f64::MAX, f64::MIN), |a, &u| (a.0.min(u), a.1.max(u)));
            assert!(close(hi - lo, TAU, 1e-12), "seam uses {us:?}");
        }
        // caps are discs of radius 0.5
        for (k, f) in m.faces.iter().enumerate() {
            if f.surface.kind() == "PLANE" {
                let d = checked_domain(&m, k);
                let a = d.area().abs();
                assert!(a < PI * 0.25 && a > PI * 0.25 * 0.99, "cap area {a}");
            }
        }
    }
    // surface of revolution (frustum): u in [0, 2pi], v = generatrix parameter in [0, 1]
    let m = model("frustum_rev.stp");
    let d = checked_domain(&m, face_of(&m, "SURFACE_OF_REVOLUTION"));
    let sp = span(&d);
    assert!(close(sp[0], TAU, 1e-12) && close(sp[1], 1.0, 1e-9), "revolution span {sp:?}");
    assert!(close(d.area(), TAU, 1e-9), "revolution area {}", d.area());
    // patch from -60 to +60 degrees crosses u = 0 without wrapping
    let m = model("cyl_partial.stp");
    let d = checked_domain(&m, 0);
    let sp = span(&d);
    assert!(close(sp[0], TAU / 3.0, 1e-12) && close(sp[1], 1.0, 1e-12), "partial span {sp:?}");
    assert!(close(d.area(), TAU / 3.0, 1e-12));
    // band bounded by two circles and no seam edge: an artificial seam joins them
    let m = model("cyl_annulus.stp");
    let d = checked_domain(&m, 0);
    assert_eq!((d.loops.len(), d.artificial_seams), (1, 1));
    assert!(close(d.area(), TAU, 1e-12), "band area {}", d.area());
    // cone with its apex: rectangle u in [0, 2pi], v from the apex (-2/3) to the base (1)
    let m = model("cone_apex.stp");
    let d = checked_domain(&m, face_of(&m, "CONICAL_SURFACE"));
    assert!(close(d.area(), TAU * 5.0 / 3.0, 1e-9), "cone area {}", d.area());
    // sphere bounded by one meridian seam between the poles
    let m = model("sphere.stp");
    let d = checked_domain(&m, 0);
    assert!(close(d.area(), TAU * PI, 1e-12), "sphere area {}", d.area());
    // torus: seams in u and v through one vertex
    let m = model("torus.stp");
    let d = checked_domain(&m, 0);
    assert!(close(d.area(), TAU * TAU, 1e-12), "torus area {}", d.area());
    // linear extrusion of an ellipse: periodic in the curve parameter
    let m = model("ellipse_prism.stp");
    let d = checked_domain(&m, face_of(&m, "SURFACE_OF_LINEAR_EXTRUSION"));
    assert!(close(d.area(), TAU, 1e-12), "extrusion area {}", d.area());
}

/// Analytic faces with holes: a box with a round hole (millimetre file),
/// every analytic surface type, and the strict / lenient handling of
/// unsupported entities.
#[test]
fn f022_analytic_faces_with_holes_parse() {
    let fine = StepOptions { chord_tol: 1e-7, ..Default::default() };
    let m = Model::from_step(&text("box_hole.stp"), &fine).unwrap();
    assert!(close(m.length_scale, 1e-3, 1e-15), "unit {}", m.length_scale);
    assert_eq!(m.faces.len(), 7);
    assert_eq!(m.n_shells, 1);
    assert!(m.unsupported.is_empty(), "{:?}", m.unsupported);
    let with_holes: Vec<usize> = (0..7).filter(|&i| m.faces[i].bounds.len() == 2).collect();
    assert_eq!(with_holes.len(), 2, "top and bottom faces carry the hole");
    let (r, h) = (0.2, 0.5);
    let mut planar = 0.0;
    for fi in 0..7 {
        let d = checked_domain(&m, fi);
        let f = &m.faces[fi];
        match f.surface.kind() {
            // planes are isometric: domain area = face area (outer minus hole)
            "PLANE" => {
                let a = d.area();
                assert!(a > 0.0, "outward planar faces are counter-clockwise");
                if f.bounds.len() == 2 {
                    assert!(close(a, 2.0 - PI * r * r, 1e-6), "face with hole: {a}");
                    let hole = d.loops.iter().map(|l| l.signed_area()).fold(f64::MAX, f64::min);
                    assert!(hole < 0.0 && close(-hole, PI * r * r, 1e-6), "hole loop {hole}");
                }
                planar += a;
            }
            // the hole wall faces the axis (same_sense false): clockwise, |area| * r = 2 pi r h
            "CYLINDRICAL_SURFACE" => {
                assert!(!f.same_sense);
                assert!(d.area() < 0.0);
                assert!(close(-d.area() * r, TAU * r * h, 1e-12), "wall {}", d.area());
            }
            k => panic!("unexpected {k}"),
        }
    }
    assert!(close(planar, 2.0 * (2.0 - PI * r * r) + 2.0 * 2.0 * h + 2.0 * 1.0 * h, 1e-6), "planar area {planar}");

    // every analytic surface and curve type is recognised
    let mut surfaces = std::collections::BTreeSet::new();
    let mut curves = std::collections::BTreeSet::new();
    for name in ["box_hole.stp", "cyl_seam.stp", "frustum_rev.stp", "cone_apex.stp", "sphere.stp", "torus.stp", "ellipse_prism.stp"] {
        let m = model(name);
        assert!(m.unsupported.is_empty(), "{name}: {:?}", m.unsupported);
        surfaces.extend(m.surface_kinds().into_keys());
        curves.extend(m.curve_kinds().into_keys());
        for fi in 0..m.faces.len() {
            checked_domain(&m, fi);
        }
    }
    for k in ["PLANE", "CYLINDRICAL_SURFACE", "CONICAL_SURFACE", "SPHERICAL_SURFACE", "TOROIDAL_SURFACE", "SURFACE_OF_REVOLUTION", "SURFACE_OF_LINEAR_EXTRUSION"] {
        assert!(surfaces.contains(k), "{k} missing from {surfaces:?}");
    }
    for k in ["LINE", "CIRCLE", "ELLIPSE"] {
        assert!(curves.contains(k), "{k} missing from {curves:?}");
    }

    // an unsupported surface is listed in the notes, and strict mode fails on it
    let t = text("cyl_seam.stp").replace("CYLINDRICAL_SURFACE('',#4,0.5)", "OFFSET_SURFACE('',#4,0.5,.F.)");
    let m = Model::from_step(&t, &StepOptions::default()).unwrap();
    assert_eq!(m.faces.len(), 2);
    assert_eq!(m.unsupported.get("OFFSET_SURFACE"), Some(&1));
    assert!(m.notes.iter().any(|n| n.contains("OFFSET_SURFACE")), "{:?}", m.notes);
    let e = Model::from_step(&t, &StepOptions { strict: true, ..Default::default() }).err().expect("strict fails");
    assert!(e.contains("OFFSET_SURFACE"), "{e}");
    // the old reader's surfaces still come out (legacy API)
    let (tris, notes) = foiley_plane::cad::read_step(&text("box_hole.stp"), 3).unwrap();
    assert!(!tris.is_empty(), "{notes:?}");
    let _ = Vec3::ZERO;
}

fn bspline_of(m: &Model, fi: usize) -> (&foiley_plane::cad::NurbsSurface, [bool; 2]) {
    match &m.faces[fi].surface {
        foiley_plane::cad::Surface::BSpline { s, closed } => (s, *closed),
        s => panic!("not a B-spline: {}", s.kind()),
    }
}

/// B-spline surface evaluation against closed forms: the bicubic patch
/// S = (u, v, u^2 v) (simple STEP entity) and the rational loft of exact
/// circles (complex rational entity), with first derivatives.
#[test]
fn f022_bspline_surface_eval_matches_reference() {
    let m = model("bspline_patch.stp");
    let (s, closed) = bspline_of(&m, 0);
    assert_eq!((s.pu, s.pv), (3, 3));
    assert_eq!(closed, [false, false]);
    for i in 0..=8 {
        for j in 0..=8 {
            let (u, v) = (i as f64 / 8.0, j as f64 / 8.0);
            let (p, su, sv) = s.eval_derivs(u, v);
            let want = Vec3::new(u, v, u * u * v);
            assert!((p - want).norm() < 1e-14, "S({u},{v}) = {p:?}");
            assert!((s.eval(u, v) - want).norm() < 1e-14);
            assert!((su - Vec3::new(1.0, 0.0, 2.0 * u * v)).norm() < 1e-12, "Su({u},{v}) = {su:?}");
            assert!((sv - Vec3::new(0.0, 1.0, u * u)).norm() < 1e-12, "Sv({u},{v}) = {sv:?}");
        }
    }
    // rational loft: degree 1 along x through (x, r) stations, exact circles around
    let m = model("loft_fuselage.stp");
    let (s, closed) = bspline_of(&m, 0);
    assert_eq!(closed, [false, true], "closed around, open along");
    assert_eq!(m.faces[0].surface.poles(0), vec![0.0, 1.0], "nose and tail are poles");
    let st = [(0.0, 0.0), (0.3, 0.1), (1.0, 0.1), (1.4, 0.0)];
    let radius = |x: f64| {
        let k = st.windows(2).position(|w| x <= w[1].0 + 1e-15).unwrap_or(2);
        let (a, b) = (st[k], st[k + 1]);
        a.1 + (b.1 - a.1) * (x - a.0) / (b.0 - a.0)
    };
    for i in 0..=12 {
        for j in 0..=16 {
            let (u, v) = (i as f64 / 12.0, j as f64 / 16.0);
            let (p, su, sv) = s.eval_derivs(u, v);
            // x is linear in u on each third
            let k = ((u * 3.0).floor() as usize).min(2);
            let x = st[k].0 + (st[k + 1].0 - st[k].0) * (u * 3.0 - k as f64);
            assert!((p.x - x).abs() < 1e-13, "x({u}) = {} vs {x}", p.x);
            assert!((p.y.hypot(p.z) - radius(x)).abs() < 1e-13, "radius at u {u} v {v}");
            // derivatives against central differences (inside a span)
            let h = 1e-6;
            if (u * 3.0 - (u * 3.0).round()).abs() > 1e-3 {
                let fd = (s.eval(u + h, v) - s.eval(u - h, v)) / (2.0 * h);
                assert!((su - fd).norm() < 1e-6 * (1.0 + fd.norm()), "Su {su:?} vs {fd:?}");
            }
            if (v * 4.0 - (v * 4.0).round()).abs() > 1e-3 {
                let fd = (s.eval(u, v + h) - s.eval(u, v - h)) / (2.0 * h);
                assert!((sv - fd).norm() < 1e-6 * (1.0 + fd.norm()), "Sv {sv:?} vs {fd:?}");
            }
        }
    }
    // rational B-spline curve: the hole of the plate is an exact circle
    let m = model("bspline_plate_hole.stp");
    let e = m.edges.iter().find(|e| e.curve.kind() == "B_SPLINE_CURVE").expect("NURBS circle edge");
    for k in 0..=64 {
        let (t0, t1) = (e.samples[0].0, e.samples.last().unwrap().0);
        let t = t0 + (t1 - t0) * k as f64 / 64.0;
        let p = e.curve.eval(t);
        let r = (p - Vec3::new(0.5, 0.5, 0.0)).norm();
        assert!((r - 0.2).abs() < 1e-14, "circle radius {r}");
        assert!(e.curve.deriv(t).dot(p - Vec3::new(0.5, 0.5, 0.0)).abs() < 1e-12, "tangent");
    }
}

/// Edge inversion onto B-spline surfaces without PCURVEs: a curve across
/// the patch (diagonal S(t, t)), a NURBS-circle hole in a B-spline plate,
/// and the closed loft whose only edge is its seam between two poles.
#[test]
fn f022_bspline_edge_inversion_residual_1e9() {
    // direct inversion of surface points from the coarse-grid seed
    let m = model("bspline_patch.stp");
    let s = &m.faces[0].surface;
    for i in 0..=10 {
        for j in 0..=10 {
            let (u, v) = (0.03 + 0.094 * i as f64, 0.01 + 0.098 * j as f64);
            let p = s.eval(u, v);
            let uv = s.invert(p, None);
            assert!((s.eval(uv[0], uv[1]) - p).norm() < 1e-12, "residual at {u},{v}");
            assert!((uv[0] - u).abs() < 1e-9 && (uv[1] - v).abs() < 1e-9, "{uv:?} vs {u},{v}");
        }
    }
    // the two triangles of the patch, split along the diagonal curve
    let diag = m.edges.iter().find(|e| e.curve.kind() == "B_SPLINE_CURVE" && e.samples.len() > 2 && (m.vertices[e.v[1] as usize] - Vec3::new(1.0, 1.0, 1.0)).norm() < 1e-12 && (m.vertices[e.v[0] as usize]).norm() < 1e-12).expect("diagonal edge");
    let on_diag: std::collections::HashSet<u32> = diag.samples.iter().map(|s| s.1).collect();
    for fi in 0..2 {
        let d = checked_domain(&m, fi);
        assert!(close(d.area().abs(), 0.5, 1e-12), "triangle {fi}: area {}", d.area());
        assert!(d.area() > 0.0);
        for (uv, id) in d.loops[0].uv.iter().zip(&d.loops[0].ids) {
            if on_diag.contains(id) {
                assert!((uv[0] - uv[1]).abs() < 1e-9, "diagonal point at {uv:?}");
            }
        }
    }
    // plate with a NURBS-circle hole: area 1 - pi r^2
    let fine = StepOptions { chord_tol: 1e-7, ..Default::default() };
    let m = Model::from_step(&text("bspline_plate_hole.stp"), &fine).unwrap();
    let d = checked_domain(&m, 0);
    assert_eq!(d.loops.len(), 2);
    assert!(close(d.area(), 1.0 - PI * 0.04, 1e-6), "plate area {}", d.area());
    // closed loft: the seam's two uses lie one period (1) apart in v, the poles are duplicated
    let m = model("loft_fuselage.stp");
    let d = checked_domain(&m, 0);
    assert!(close(d.area(), 1.0, 1e-12), "loft domain area {}", d.area());
    let sp = span(&d);
    assert!(close(sp[0], 1.0, 1e-12) && close(sp[1], 1.0, 1e-12), "loft span {sp:?}");
}
