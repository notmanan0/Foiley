//! F005: stability derivatives from the full (viscous + extra drag) loads.
mod common;
use common::*;
use foiley_core::{Progress, Vec3};
use foiley_foil::FoilPolar;
use foiley_plane::analysis::*;
use foiley_plane::evaluator::Evaluator;
use foiley_plane::geom::Plane;
use foiley_plane::stability::StabilityModel;
use foiley_plane::FoilDb;
use std::collections::HashMap;

/// The default glider (NACA 0012 wing, 0009 tail) and a database with both foils.
fn glider() -> (Plane, FoilDb) {
    let f = naca0012();
    let f9 = foiley_foil::Foil::naca("0009").unwrap();
    let mut db = db_with(&f);
    db.add_foil(&f9);
    (Plane::default_plane(&f.name, &f9.name), db)
}

/// Synthetic polars covering Re 1e4..1e7 with a cd that falls with Re and a
/// section Cm `cm`.
fn add_polars(db: &mut FoilDb, foil: &str, cm: f64) {
    for re in [1e4, 1e5, 1e6, 1e7] {
        let mut p: FoilPolar = synth_polar(foil, &format!("Re{re}"), re, 9.0, 0.005 + 50.0 / re);
        for pt in &mut p.points {
            pt.cm = cm;
        }
        db.add_polar(&p);
    }
}

fn sd() -> PlanePolarDef {
    PlanePolarDef { viscous: false, derivatives: true, range: Range { from: 4.0, to: 4.0, step: 1.0 }, ..Default::default() }
}

fn drag() -> Vec<ExtraDrag> {
    vec![ExtraDrag { name: "audit drag".into(), area: 1.0, cd: 0.1 }]
}

#[test]
fn f005_fullload_xu_moves_with_drag() {
    let (glider, db) = glider();
    let sd = PlanePolarDef { stability_model: StabilityModel::FullLoad, ..sd() };
    let clean = run1(&glider, &db, &sd);
    let draggy = run1(&glider, &db, &PlanePolarDef { extra_drag: drag(), ..sd });
    assert!(draggy.cd > clean.cd);
    let (a, b) = (clean.stability.unwrap().dim.xu, draggy.stability.unwrap().dim.xu);
    assert!((a - b).abs() > 1e-6, "Xu {a} vs {b}");
    // a constant-CD drag D = q S CD gives dX/du = -2 D / u0
    let u0 = draggy.velocity;
    let d = 0.5 * PlanePolarDef::default().fluid.rho * u0 * u0 * 0.1;
    assert!(((b - a) / (-2.0 * d / u0) - 1.0).abs() < 1e-6, "dXu {} vs {}", b - a, -2.0 * d / u0);
}

#[test]
fn f005_inviscid_xu_unchanged() {
    let (glider, db) = glider();
    let inv = run1(&glider, &db, &sd());
    let inv_drag = run1(&glider, &db, &PlanePolarDef { extra_drag: drag(), ..sd() });
    let (a, b) = (inv.stability.unwrap(), inv_drag.stability.unwrap());
    // the default model ignores the extra drag, bit for bit
    assert_eq!(a.dim, b.dim);
    // and without any non-potential load FullLoad reproduces it
    let full = run1(&glider, &db, &PlanePolarDef { stability_model: StabilityModel::FullLoad, ..sd() }).stability.unwrap();
    for (x, y, n) in [(a.dim.xu, full.dim.xu, "Xu"), (a.dim.zw, full.dim.zw, "Zw"), (a.dim.mq, full.dim.mq, "Mq"), (a.dim.lp, full.dim.lp, "Lp"), (a.dim.nr, full.dim.nr, "Nr")] {
        assert!((x - y).abs() <= 1e-9 * x.abs().max(1e-9), "{n}: inviscid {x} vs FullLoad {y}");
    }
}

#[test]
fn f005_fullload_trim_includes_viscous_cm() {
    let (glider, mut db) = glider();
    add_polars(&mut db, "NACA 0012", -0.01);
    add_polars(&mut db, "NACA 0009", 0.0);
    let t7 = PlanePolarDef { kind: PolarKind::Stability, viscous: true, range: Range { from: 0.0, to: 0.0, step: 1.0 }, ..Default::default() };
    let inv = run1(&glider, &db, &t7);
    let full = run1(&glider, &db, &PlanePolarDef { stability_model: StabilityModel::FullLoad, ..t7.clone() });
    // the FullLoad trim zeroes the reduced Cm, section Cm and drag moments included
    assert!(full.cm.abs() < 1e-6, "FullLoad trimmed Cm = {}", full.cm);
    assert!(full.strips.iter().any(|s| s.cm_c4 != 0.0));
    // the nose-down section Cm moves the (statically stable) trim to a lower angle of attack
    assert!(full.alpha < inv.alpha - 0.05, "trim alpha FullLoad {} vs inviscid {}", full.alpha, inv.alpha);
    // the inviscid trim leaves the full Cm unbalanced
    let (p, d) = (&glider, PlanePolarDef { kind: PolarKind::FixedSpeed, stability_model: StabilityModel::FullLoad, velocity: inv.velocity, range: Range { from: inv.alpha, to: inv.alpha, step: 1.0 }, ..t7 });
    let at_inv = run1(p, &db, &d);
    assert!(at_inv.cm.abs() > 1e-3, "full Cm at the inviscid trim = {}", at_inv.cm);
    let (si, sf) = (inv.stability.unwrap(), full.stability.unwrap());
    assert!((si.dim.xu - sf.dim.xu).abs() > 1e-6, "Xu {} vs {}", si.dim.xu, sf.dim.xu);
}

#[test]
fn f005_reduce_uses_local_velocity() {
    // a roll rate speeds up the strips of one wing and slows the other:
    // their Reynolds numbers follow V - omega x r
    let (glider, mut db) = glider();
    add_polars(&mut db, "NACA 0012", 0.0);
    add_polars(&mut db, "NACA 0009", 0.0);
    let def = PlanePolarDef { viscous: true, ..Default::default() };
    let ctx = Ctx::new(&glider, &db, &def).unwrap();
    let mesh = foiley_plane::mesh::build(&glider, &db, &mesh_options(&def, &HashMap::new(), &HashMap::new())).unwrap();
    let cog = ctx.inertia.cog;
    let ev = Evaluator::new(make_system(&def, &mesh, cog, ctx.c_ref, &Progress::none()).unwrap(), &def, ctx.c_ref);
    let (alpha, v) = (4.0, 10.0);
    let wd = wind_dir(alpha, 0.0);
    // yaw rate r about z: local air velocity V - r z x (c4 - cog) changes along y
    let omega = Vec3::new(0.0, 0.0, 0.5);
    let (raw, _) = ev.loads(wd, v, omega, &Progress::none()).unwrap();
    let plain = ctx.reduce(&mesh, &raw, cog, alpha, 0.0, v, ctx.inertia.mass, None);
    let local = ctx.reduce(&mesh, &raw, cog, alpha, 0.0, v, ctx.inertia.mass, Some(omega));
    let nu = def.fluid.nu();
    for (k, s) in mesh.strips.iter().enumerate() {
        let vl = (wd * v - omega.cross(s.c4 - cog)).norm();
        assert!((local.strips[k].re / (vl * s.chord / nu) - 1.0).abs() < 1e-12, "strip {k}");
        assert!((plain.strips[k].re / (v * s.chord / nu) - 1.0).abs() < 1e-12);
    }
    assert!(local.cdv != plain.cdv);
    // yawing right speeds up the left wing: its strips drag more
    let side = |o: &PlaneOpp, left: bool| -> f64 { mesh.strips.iter().zip(&o.strips).filter(|(s, _)| s.wing == 0 && s.left == left).map(|(_, r)| r.re).sum() };
    assert!((side(&local, true) - side(&local, false)).abs() > 1e-3 * side(&local, true));
}
