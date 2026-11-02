//! Physics validation of the 3D solvers against theory and consistency
//! checks.  Run with `cargo test --release -p foiley-plane`.

use foiley_core::{Progress, Spacing, Vec3};
use foiley_foil::{analyze, AnalysisSettings, Foil, Sweep};
use foiley_plane::analysis::*;
use foiley_plane::geom::*;
use foiley_plane::FoilDb;
use std::sync::OnceLock;

fn rect(foil: &str, ar: f64, nx: usize, ny: usize) -> Plane {
    let mut a = WingSection::new(0.0, 1.0, 0.0, foil);
    a.nx = nx;
    a.ny = ny;
    a.ydist = Spacing::Uniform;
    let b = WingSection::new(ar / 2.0, 1.0, 0.0, foil);
    Plane { name: "rect".into(), description: String::new(), wings: vec![Wing::new("w", WingKind::Main, vec![a, b])], fuses: vec![], point_masses: vec![] }
}

/// Foils with polar meshes (computed once with the 2D engine).
fn db() -> &'static FoilDb {
    static DB: OnceLock<FoilDb> = OnceLock::new();
    DB.get_or_init(|| {
        let mut db = FoilDb::new();
        let mut thin = Foil::naca("0012").unwrap().scaled_tc(0.01, 1.0).unwrap();
        thin.name = "thin".into();
        db.add_foil(&thin);
        for code in ["2412", "0009", "0012"] {
            let f = Foil::naca(code).unwrap();
            db.add_foil(&f);
            for re in [6e4, 1.2e5, 2.5e5, 5e5, 1e6] {
                let st = AnalysisSettings { re, ..Default::default() };
                let (p, _) = analyze(&f, &st, &Sweep::Alpha { from: -8.0, to: 14.0, step: 0.5 }, &Progress::none()).unwrap();
                db.add_polar(&p);
            }
        }
        db
    })
}

fn run(p: &Plane, def: &PlanePolarDef) -> PlanePolar {
    run_polar(p, db(), def, &Progress::none(), &|_| {}).unwrap()
}

fn inviscid(method: Method, alpha: f64) -> PlanePolarDef {
    PlanePolarDef { method, viscous: false, reference: RefDims::Planform, range: Range { from: alpha, to: alpha, step: 1.0 }, ..Default::default() }
}

#[test]
fn vlm_rectangular_wing_matches_lifting_surface_theory() {
    let p = rect("thin", 8.0, 8, 40);
    for m in [Method::Vlm1, Method::Vlm2] {
        let r = run(&p, &inviscid(m, 5.0));
        let o = &r.points[0];
        let cla = o.cl / 5f64.to_radians().sin();
        // lifting-surface value for AR 8 rectangular ~4.6/rad (Helmbold 4.9 is an upper bound)
        assert!((4.5..4.75).contains(&cla), "{m:?} CLa = {cla}");
        // rectangular wing span efficiency ~0.97-0.99
        assert!((0.95..1.0).contains(&o.efficiency), "{m:?} e = {}", o.efficiency);
        assert!(o.cy.abs() < 1e-9 && o.cn.abs() < 1e-9, "symmetric flight must have no side force / yaw");
    }
}

#[test]
fn panel_method_thickness_effect_and_consistency() {
    let pv = run(&rect("thin", 8.0, 8, 20), &inviscid(Method::Vlm2, 5.0));
    let pp = run(&rect("NACA 0012", 8.0, 16, 20), &inviscid(Method::Panel, 5.0));
    let (cv, cp) = (pv.points[0].cl, pp.points[0].cl);
    // a 12% thick wing lifts a few percent more than a flat plate in potential flow
    assert!(cp > cv && cp < 1.10 * cv, "panel {cp} vs VLM {cv}");
    let r0 = run(&rect("NACA 0012", 8.0, 12, 12), &inviscid(Method::Panel, 0.0));
    assert!(r0.points[0].cl.abs() < 1e-6, "symmetric section at 0 deg must have no lift");
}

#[test]
fn ground_effect_reduces_induced_drag() {
    let p = rect("thin", 8.0, 6, 20);
    let free = run(&p, &inviscid(Method::Vlm2, 4.0)).points[0].clone();
    let mut d = inviscid(Method::Vlm2, 4.0);
    d.ground_height = Some(0.3);
    let ge = run(&p, &d).points[0].clone();
    assert!(ge.cdi < 0.8 * free.cdi, "CDi free {} vs ground {}", free.cdi, ge.cdi);
    assert!(ge.cl > free.cl, "ground effect should increase lift");
}

#[test]
fn dihedral_gives_stable_roll_due_to_sideslip() {
    let mut p = rect("thin", 8.0, 4, 12);
    p.wings[0].sections[0].dihedral = 6.0;
    let mut d = inviscid(Method::Vlm2, 4.0);
    d.kind = PolarKind::Sideslip;
    d.alpha = 4.0;
    d.range = Range { from: -5.0, to: 5.0, step: 5.0 };
    let r = run(&p, &d);
    let (m, z, pl) = (&r.points[0], &r.points[1], &r.points[2]);
    assert!(z.cr.abs() < 1e-9 && z.cy.abs() < 1e-9);
    // antisymmetric response
    assert!((m.cr + pl.cr).abs() < 1e-6 * (1.0 + m.cr.abs()));
    assert!(m.cr != 0.0 && m.cy != 0.0);
}

#[test]
fn default_plane_viscous_t1_t2_t3() {
    let p = Plane::default_plane("NACA 2412", "NACA 0009");
    let mut d = PlanePolarDef { range: Range { from: -2.0, to: 8.0, step: 2.0 }, ..Default::default() };
    let r = run(&p, &d);
    assert_eq!(r.points.len(), 6, "failed: {:?}", r.failed);
    for o in &r.points {
        assert!(o.cdv > 0.004 && o.cdv < 0.05, "CDv {}", o.cdv);
        assert!(o.cd > o.cdi);
    }
    // Cm decreasing with alpha (statically stable default plane)
    assert!(r.points.last().unwrap().cm < r.points[0].cm);
    let ld = r.points.iter().map(|o| o.ld()).fold(0.0, f64::max);
    assert!((10.0..40.0).contains(&ld), "L/D max {ld}");
    // T2: lift = weight
    d.kind = PolarKind::FixedLift;
    let r2 = run(&p, &d);
    let w = r2.points[0].mass * G;
    for o in &r2.points {
        if o.cl > 0.05 {
            assert!((o.lift - w).abs() / w < 1e-6, "lift {} weight {}", o.lift, w);
        }
    }
    // T3: aerodynamic force = weight
    d.kind = PolarKind::Glide;
    let r3 = run(&p, &d);
    for o in &r3.points {
        let ca = (o.cl * o.cl + o.cd * o.cd).sqrt();
        let fa = 0.5 * 1.225 * o.velocity * o.velocity * r3.ref_dims.0 * ca;
        assert!((fa - w).abs() / w < 0.01, "aero force {fa} vs weight {w}");
        assert!(o.vz() > 0.0 && o.vz() < 2.0);
    }
}

#[test]
fn stability_polar_trims_and_finds_modes() {
    let p = Plane::default_plane("NACA 2412", "NACA 0009");
    let d = PlanePolarDef { kind: PolarKind::Stability, range: Range { from: 0.0, to: 0.0, step: 1.0 }, ..Default::default() };
    let r = run(&p, &d);
    assert_eq!(r.points.len(), 1, "{:?}", r.failed);
    let o = &r.points[0];
    assert!(o.cm.abs() < 0.01, "trimmed Cm {} (viscous moment only)", o.cm);
    let s = o.stability.as_ref().expect("stability result");
    assert!(s.static_margin > 0.02 && s.static_margin < 0.6, "static margin {}", s.static_margin);
    assert!(s.nondim.cma < 0.0 && s.nondim.cmq < 0.0 && s.nondim.clp < 0.0 && s.nondim.cnb > 0.0, "{:?}", s.nondim);
    let names: Vec<&str> = s.modes.iter().map(|m| m.name.as_str()).collect();
    for want in ["Short period", "Phugoid", "Dutch roll", "Roll damping", "Spiral"] {
        assert!(names.iter().any(|n| n.starts_with(want)), "missing {want}: {names:?}");
    }
    let sp = s.modes.iter().find(|m| m.name == "Short period").unwrap();
    assert!(sp.stable && sp.damping > 0.3, "short period {:?}", sp);
    let ph = s.modes.iter().find(|m| m.name == "Phugoid").unwrap();
    assert!(ph.freq_hz < sp.freq_hz && ph.freq_hz > 0.02);
    let roll = s.modes.iter().find(|m| m.name == "Roll damping").unwrap();
    assert!(roll.stable);
    // phugoid frequency vs the quasi-steady (short-period residualised)
    // approximation w^2 = -g Zu / (m u0 + Zq - Zw Mq / Mw)
    let d = &s.dim;
    let m = o.mass;
    let w2 = -G * d.zu / (m * s.u0 + d.zq - d.zw * d.mq / d.mw);
    let w_ph = ph.eigenvalue.abs();
    assert!((w_ph / w2.sqrt() - 1.0).abs() < 0.15, "phugoid {w_ph} rad/s vs approximation {}", w2.sqrt());
    assert!((d.zu + 2.0 * m * G / s.u0).abs() < 0.02 * d.zu.abs(), "Zu = -2mg/u0");
    let resp = foiley_plane::stability::time_response(&s.a_long, [1.0, 0.0, 0.0, 0.0], 5.0, 0.01);
    assert!(resp.len() > 400);
}

#[test]
fn llt_close_to_vlm_in_linear_range() {
    // chord 1 m at 15 m/s: Re ~ 1e6, where the NACA 0012 polar is linear
    // (at Re ~ 1e5 laminar bubbles make it strongly nonlinear in XFoil)
    let p = rect("NACA 0012", 8.0, 8, 20);
    let mut d = PlanePolarDef { method: Method::Llt, velocity: 15.0, range: Range { from: 4.0, to: 4.0, step: 1.0 }, reference: RefDims::Planform, ..Default::default() };
    let l = run(&p, &d);
    assert_eq!(l.points.len(), 1, "{:?}", l.failed);
    d.method = Method::Vlm2;
    let v = run(&p, &d);
    let (cl_l, cl_v) = (l.points[0].cl, v.points[0].cl);
    // lifting line vs lifting surface: within ~10% in the linear range
    assert!((cl_l / cl_v - 1.0).abs() < 0.12, "LLT {cl_l} VLM {cl_v}");
    assert!((0.9..1.0).contains(&l.points[0].efficiency), "LLT e = {}", l.points[0].efficiency);
    assert!(l.points[0].cdv > 0.0 && l.points[0].cdi > 0.0);
}

#[test]
fn control_polar_flap_changes_moment() {
    let mut p = Plane::default_plane("NACA 2412", "NACA 0009");
    // make the elevator foil a control surface
    let mut ff = Foil::naca("0009").unwrap();
    ff.name = "elev".into();
    ff.hinge = Some((0.7, 0.5));
    let mut dbx = db().clone();
    dbx.add_foil(&ff);
    for s in &mut p.wings[1].sections {
        s.foil_left = "elev".into();
        s.foil_right = "elev".into();
    }
    let d = PlanePolarDef {
        kind: PolarKind::Control,
        viscous: false,
        alpha: 2.0,
        range: Range { from: 0.0, to: 1.0, step: 1.0 },
        ranges: vec![ControlRange { target: ControlTarget::Flap { wing: 1, index: 0 }, min: 0.0, max: 10.0 }, ControlRange { target: ControlTarget::Flap { wing: 1, index: 1 }, min: 0.0, max: 10.0 }],
        ..Default::default()
    };
    let r = run_polar(&p, &dbx, &d, &Progress::none(), &|_| {}).unwrap();
    assert_eq!(r.points.len(), 2);
    // elevator TE down -> more tail lift -> nose-down moment
    assert!(r.points[1].cm < r.points[0].cm - 0.02, "{} -> {}", r.points[0].cm, r.points[1].cm);
    let _ = Vec3::ZERO;
}

#[test]
fn panel_method_with_fuselage_runs() {
    let mut p = Plane::default_plane("NACA 2412", "NACA 0009");
    p.fuses.push(Fuse::default_pod());
    for w in &mut p.wings {
        for s in &mut w.sections {
            s.ny = s.ny.min(8);
            s.nx = s.nx.min(8);
        }
    }
    let d = PlanePolarDef { method: Method::Panel, range: Range { from: 0.0, to: 4.0, step: 4.0 }, ..Default::default() };
    let r = run(&p, &d);
    assert_eq!(r.points.len(), 2, "{:?}", r.failed);
    assert!(r.points[1].cl > r.points[0].cl + 0.2);
    assert!(r.points.iter().all(|o| o.cp.iter().all(|c| c.is_finite())));
}

const FLOW5_EX: &str = "../../refsrc/flow5/API_examples";

#[test]
fn triangle_panels_match_quad_panels() {
    let p = rect("NACA 0012", 8.0, 16, 20);
    let q = run(&p, &inviscid(Method::Panel, 5.0)).points[0].cl;
    let t = run(&p, &inviscid(Method::TriPanel, 5.0)).points[0].cl;
    assert!((t / q - 1.0).abs() < 0.04, "tri {t} vs quad {q}");
    assert!(run(&p, &inviscid(Method::TriPanel, 0.0)).points[0].cl.abs() < 1e-4);
}

#[test]
fn nurbs_fuselage_and_cad_imports() {
    // NURBS body converted from the elliptic pod meshes and analyses
    let mut p = Plane::default_plane("NACA 2412", "NACA 0009");
    let mut f = Fuse::default_pod();
    f.nurbs = f.to_nurbs();
    let a_ell = Fuse::default_pod().wetted_area();
    let a_nrb = f.wetted_area();
    assert!((a_nrb / a_ell - 1.0).abs() < 0.15, "NURBS area {a_nrb} vs elliptic {a_ell}");
    p.fuses.push(f);
    for w in &mut p.wings {
        for s in &mut w.sections {
            s.nx = s.nx.min(8);
            s.ny = s.ny.min(8);
        }
    }
    let r = run(&p, &PlanePolarDef { method: Method::Panel, viscous: false, range: Range { from: 2.0, to: 2.0, step: 1.0 }, ..Default::default() });
    assert_eq!(r.points.len(), 1);
    // CAD fuselages (flow5 example files): STL and STEP of the same body
    let stl = std::path::Path::new(FLOW5_EX).join("PlaneRun5/stl/pw5_fuse.stl");
    let stp = std::path::Path::new(FLOW5_EX).join("PlaneRun6/step/pw5_2024-06-27T20_19_49.010Z.stp");
    if stl.exists() && stp.exists() {
        let t1 = foiley_plane::cad::read_stl(&std::fs::read(&stl).unwrap()).unwrap();
        let (t2, notes) = foiley_plane::cad::read_step(&std::fs::read_to_string(&stp).unwrap(), 3).unwrap();
        let (a, b) = (foiley_plane::cad::bbox(&t1), foiley_plane::cad::bbox(&t2));
        let (l1, l2) = ((a.1 - a.0).x, (b.1 - b.0).x);
        println!("STL {} tris, length {l1}; STEP {} tris, length {l2}; {:?}", t1.len(), t2.len(), notes);
        assert!(t2.len() > 500);
        let fz = Fuse::from_triangles("stl", t1, 1.0);
        let vol: f64 = fz.triangles.iter().map(|t| t[0].dot(t[1].cross(t[2]))).sum::<f64>() / 6.0;
        assert!(vol > 0.0, "outward orientation");
        assert!(fz.section_at(fz.position.x + 0.5 * fz.length()).is_some());
    }
}

#[test]
fn free_wake_close_to_rigid_and_rolls_down() {
    let p = rect("thin", 8.0, 4, 16);
    let rigid = run(&p, &inviscid(Method::Vlm2, 6.0)).points[0].clone();
    let mut d = inviscid(Method::Vlm2, 6.0);
    d.free_wake = 4;
    let free = run(&p, &d).points[0].clone();
    assert!((free.cl / rigid.cl - 1.0).abs() < 0.05, "free {} rigid {}", free.cl, rigid.cl);
    assert!(!free.wake.is_empty());
    // wake trails downward relative to the freestream direction behind a lifting wing
    let mid = &free.wake[free.wake.len() / 2];
    let (a, b) = (mid[0], mid[mid.len() / 3]);
    let dz = (b[2] - a[2]) as f64 - ((b[0] - a[0]) as f64) * 6f64.to_radians().tan();
    assert!(dz < 0.0, "downwash should deflect the wake down: {dz}");
}

#[test]
fn control_derivatives_and_forced_response() {
    let mut p = Plane::default_plane("NACA 2412", "NACA 0009");
    let mut ff = Foil::naca("0009").unwrap();
    ff.name = "elev".into();
    ff.hinge = Some((0.7, 0.5));
    let mut dbx = db().clone();
    dbx.add_foil(&ff);
    for s in &mut p.wings[1].sections {
        s.foil_left = "elev".into();
        s.foil_right = "elev".into();
    }
    let d = PlanePolarDef {
        kind: PolarKind::Stability,
        viscous: false,
        range: Range { from: 0.0, to: 0.0, step: 1.0 },
        gains: vec![(ControlTarget::Flap { wing: 1, index: 0 }, 5.0), (ControlTarget::Flap { wing: 1, index: 1 }, 5.0)],
        ..Default::default()
    };
    let r = run_polar(&p, &dbx, &d, &Progress::none(), &|_| {}).unwrap();
    let s = r.points[0].stability.as_ref().unwrap();
    let c = s.control.expect("control derivatives");
    assert!(c[4] < 0.0, "elevator TE down must pitch nose down: M_delta = {}", c[4]);
    let resp = foiley_plane::stability::forced_response(&s.a_long, &s.b_long, foiley_plane::stability::ControlInput::Step, 0.2, 3.0, 0.005);
    let theta_end = resp.last().unwrap().1[3];
    assert!(theta_end < 0.0, "step down-elevator pitches nose down: theta = {theta_end}");
}

#[test]
fn tri_panels_with_fuselage_consistent_with_quads() {
    let mut p = Plane::default_plane("NACA 2412", "NACA 0009");
    p.wings.truncate(2);
    p.fuses.push(Fuse::default_pod());
    for w in &mut p.wings {
        for s in &mut w.sections {
            s.nx = s.nx.min(8);
            s.ny = s.ny.min(8);
        }
    }
    let mut d = PlanePolarDef { method: Method::Panel, viscous: false, range: Range { from: 0.0, to: 4.0, step: 2.0 }, ..Default::default() };
    let q = run(&p, &d);
    d.method = Method::TriPanel;
    let t = run(&p, &d);
    for (a, b) in q.points.iter().zip(&t.points) {
        println!("alpha {}: quad CL {:.4} CDi {:.5} Cn {:.1e} | tri CL {:.4} CDi {:.5} Cn {:.1e}", a.alpha, a.cl, a.cdi, a.cn, b.cl, b.cdi, b.cn);
        assert!((b.cl - a.cl).abs() < 0.03 + 0.05 * a.cl.abs());
        assert!(b.cn.abs() < 1e-3 && a.cn.abs() < 1e-3, "symmetric plane must have no yaw");
        assert!(b.cdi >= 0.0 && a.cdi >= 0.0);
    }
}
