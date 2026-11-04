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

// ---------------------------------------------------------------- F005b / F006

/// Projected derivatives of the evaluator's inviscid loads by central
/// differences (step `h`), in the stability axes of (alpha, beta):
/// rows X Y Z L M N, columns u v w p q r.
fn fd_derivs(p: &Plane, db: &FoilDb, def: &PlanePolarDef, alpha: f64, beta: f64, u0: f64, h: f64) -> [[f64; 6]; 6] {
    let ctx = Ctx::new(p, db, def).unwrap();
    let cog = ctx.inertia.cog;
    let mesh = foiley_plane::mesh::build(p, db, &mesh_options(def, &HashMap::new(), &HashMap::new())).unwrap();
    let ev = Evaluator::new(make_system(def, &mesh, cog, ctx.c_ref, &Progress::none()).unwrap(), def, ctx.c_ref);
    let rho = def.fluid.rho;
    let is = -wind_dir(alpha, beta);
    let ks = -wind_normal(alpha);
    let js = ks.cross(is);
    let e = [is, js, ks];
    let v0 = -is * u0;
    let fm = |vinf: Vec3, om: Vec3| {
        let (r, _) = ev.loads_at(vinf, om, &Progress::none()).unwrap();
        let (f, m) = ((r.force + r.induced_force) * rho, r.moment * rho);
        [f.dot(is), f.dot(js), f.dot(ks), m.dot(is), m.dot(js), m.dot(ks)]
    };
    let hw = h / ctx.c_ref;
    let mut out = [[0.0; 6]; 6];
    for k in 0..3 {
        // translation: air velocity v0 - e dh when the aircraft speeds up along e
        let (a, b) = (fm(v0 - e[k] * h, Vec3::ZERO), fm(v0 + e[k] * h, Vec3::ZERO));
        let (c, d) = (fm(v0, e[k] * hw), fm(v0, -e[k] * hw));
        for r in 0..6 {
            out[r][k] = (a[r] - b[r]) / (2.0 * h);
            out[r][3 + k] = (c[r] - d[r]) / (2.0 * hw);
        }
    }
    out
}

/// The derivatives of a result in the same layout as `fd_derivs`.
fn dims(s: &foiley_plane::stability::StabilityResult) -> [[f64; 6]; 6] {
    let d = &s.dim;
    [
        [d.xu, d.xv, d.xw, d.xp, d.xq, d.xr],
        [d.yu, d.yv, d.yw, d.yp, d.yq, d.yr],
        [d.zu, d.zv, d.zw, d.zp, d.zq, d.zr],
        [d.lu, d.lv, d.lw, d.lp, d.lq, d.lr],
        [d.mu, d.mv, d.mw, d.mp, d.mq, d.mr],
        [d.nu, d.nv, d.nw, d.np, d.nq, d.nr],
    ]
}

fn assert_fd_match(s: &foiley_plane::stability::StabilityResult, fd: &[[f64; 6]; 6], what: &str) {
    let dm = dims(s);
    let scale = dm.iter().flatten().fold(0.0f64, |a, v| a.max(v.abs()));
    for r in 0..6 {
        for c in 0..6 {
            let (a, b) = (dm[r][c], fd[r][c]);
            assert!((a - b).abs() <= 1e-6 * a.abs().max(1e-6 * scale), "{what}: derivative [{}][{}] {a} vs finite difference {b}", "XYZLMN".as_bytes()[r] as char, ["u", "v", "w", "p", "q", "r"][c]);
        }
    }
}

fn t1_derivs(alpha: f64, beta: f64) -> PlanePolarDef {
    PlanePolarDef { derivatives: true, free_wake: 0, alpha, beta, range: Range { from: alpha, to: alpha, step: 1.0 }, ..review_def() }
}

#[test]
fn f006_symmetric_beta0_block_diagonal_matches_4x4() {
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    let s = run1(&p, &db, &t1_derivs(4.0, 0.0)).stability.unwrap();
    // off-diagonal blocks vanish against their rows' diagonal blocks
    for r in 0..8 {
        let blk = if r < 4 { 0 } else { 4 };
        let diag = (blk..blk + 4).map(|c| s.a[r][c].abs()).fold(0.0, f64::max);
        for c in (0..8).filter(|c| (c / 4) * 4 != blk) {
            assert!(s.a[r][c].abs() <= 1e-12 * diag, "a[{r}][{c}] = {} vs row scale {diag}", s.a[r][c]);
        }
    }
    // the blocks are the previous 4x4 assembly from the same derivatives
    let (d, m, u0, g) = (s.dim, review_inertia().mass, s.u0, G);
    let [ixx, iyy, izz, izx] = s.inertia_s;
    let a_long = [[d.xu / m, d.xw / m, 0.0, -g], [d.zu / m, d.zw / m, (d.zq + m * u0) / m, 0.0], [d.mu / iyy, d.mw / iyy, d.mq / iyy, 0.0], [0.0, 0.0, 1.0, 0.0]];
    let den = ixx * izz - izx * izx;
    let (ixp, izp, izxp) = (den / izz, den / ixx, izx / den);
    let a_lat = [
        [d.yv / m, d.yp / m, d.yr / m - u0, g],
        [d.lv / ixp + izxp * d.nv, d.lp / ixp + izxp * d.np, d.lr / ixp + izxp * d.nr, 0.0],
        [izxp * d.lv + d.nv / izp, izxp * d.lp + d.np / izp, izxp * d.lr + d.nr / izp, 0.0],
        [0.0, 1.0, 0.0, 0.0],
    ];
    for (old, new, off) in [(&a_long, &s.a_long, 0), (&a_lat, &s.a_lat, 4)] {
        for r in 0..4 {
            let sc = old[r].iter().fold(0.0f64, |a, v| a.max(v.abs()));
            for c in 0..4 {
                assert!((old[r][c] - new[r][c]).abs() <= 1e-12 * sc, "block {off} [{r}][{c}]: {} vs 4x4 {}", new[r][c], old[r][c]);
                assert_eq!(new[r][c], s.a[off + r][off + c]);
            }
        }
    }
    // and the modes are those of the two 4x4 matrices
    let eig = |a: &[[f64; 4]; 4]| foiley_core::linalg::eigen(&a.iter().map(|r| r.to_vec()).collect::<Vec<_>>()).unwrap().0;
    for (blk, long) in [(&a_long, true), (&a_lat, false)] {
        for l in eig(blk).into_iter().filter(|l| l.im >= 0.0) {
            let hit = s.modes.iter().filter(|m| m.longitudinal == long).any(|m| (m.eigenvalue.re - l.re).abs() + (m.eigenvalue.im - l.im).abs() < 1e-9 * l.abs().max(1.0));
            // the lateral classifier names roll and spiral, not every real root
            let named = long || l.im > 0.0;
            assert!(hit || !named, "eigenvalue {l:?} (longitudinal {long}) missing from {:?}", s.modes.iter().map(|m| (&m.name, m.eigenvalue)).collect::<Vec<_>>());
        }
    }
    for m in &s.modes {
        let f = m.longitudinal_fraction;
        assert!(if m.longitudinal { f > 0.999 } else { f < 0.001 }, "{}: fraction {f}", m.name);
    }
}

#[test]
fn f006_asymmetric_beta0_coupling_matches_fd() {
    // CoG off the plane of symmetry, with a product of inertia
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    let inertia = foiley_plane::geom::Inertia { cog: Vec3::new(0.25, 0.4, 0.0), ixy: 1.5, ..review_inertia() };
    let def = PlanePolarDef { inertia: Some(inertia), ..t1_derivs(4.0, 0.0) };
    let o = run1(&p, &db, &def);
    let s = o.stability.unwrap();
    assert_fd_match(&s, &fd_derivs(&p, &db, &def, 4.0, 0.0, o.velocity, 1e-4 * o.velocity), "asymmetric CoG");
    // real coupling terms: rolling moment with speed and with w
    let d = s.dim;
    assert!(d.lu.abs() > 1e-3 * d.mw.abs() && d.lw.abs() > 1e-3 * d.mw.abs(), "Lu {} Lw {} vs Mw {}", d.lu, d.lw, d.mw);
    // and they reach the coupled matrix through the inverse inertia tensor
    let (is, ks) = (-wind_dir(4.0, 0.0), -wind_normal(4.0));
    let js = ks.cross(is);
    let jm = [[inertia.ixx, -inertia.ixy, -inertia.ixz], [-inertia.ixy, inertia.iyy, -inertia.iyz], [-inertia.ixz, -inertia.iyz, inertia.izz]];
    let e = [is, js, ks];
    let jb = |a: Vec3, b: Vec3| {
        let (a, b) = ([a.x, a.y, a.z], [b.x, b.y, b.z]);
        (0..3).map(|i| (0..3).map(|k| a[i] * jm[i][k] * b[k]).sum::<f64>()).sum::<f64>()
    };
    let jj: [[f64; 3]; 3] = std::array::from_fn(|r| std::array::from_fn(|c| jb(e[r], e[c])));
    let det3 = |m: [[f64; 3]; 3]| m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    // J x = (L, M, N), by Cramer
    let solve = |rhs: [f64; 3]| -> [f64; 3] {
        std::array::from_fn(|k| {
            let mut m = jj;
            for r in 0..3 {
                m[r][k] = rhs[r];
            }
            det3(m) / det3(jj)
        })
    };
    for (col, (l, m_, n)) in [(0usize, (d.lu, d.mu, d.nu)), (1, (d.lw, d.mw, d.nw)), (4, (d.lv, d.mv, d.nv))] {
        let w = solve([l, m_, n]);
        for (row, k) in [(5usize, 0usize), (2, 1), (6, 2)] {
            assert!((s.a[row][col] - w[k]).abs() <= 1e-10 * w[k].abs().max(1e-12), "a[{row}][{col}] = {} vs {}", s.a[row][col], w[k]);
        }
    }
    let off = (0..4).flat_map(|r| (4..8).map(move |c| (r, c))).map(|(r, c)| s.a[r][c].abs().max(s.a[c][r].abs())).fold(0.0, f64::max);
    assert!(off > 1e-3, "coupling dropped: max off-diagonal {off}");
}

#[test]
fn f006_nonzero_beta_derivs_match_fd() {
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    let def = t1_derivs(4.0, 5.0);
    let o = run1(&p, &db, &def);
    let s = o.stability.unwrap();
    assert_eq!(s.beta, 5.0);
    assert_fd_match(&s, &fd_derivs(&p, &db, &def, 4.0, 5.0, o.velocity, 1e-4 * o.velocity), "beta = 5");
    // the sideslip couples the sets: these are not the beta = 0 derivatives
    let s0 = run1(&p, &db, &t1_derivs(4.0, 0.0)).stability.unwrap();
    assert!((s.dim.xv - s0.dim.xv).abs() > 1e-3 * s.dim.xu.abs(), "Xv beta 5 {} vs beta 0 {}", s.dim.xv, s0.dim.xv);
    assert!((s.dim.lu - s0.dim.lu).abs() > 1e-6, "Lu beta 5 {} vs beta 0 {}", s.dim.lu, s0.dim.lu);
}

#[test]
fn g34_phugoid_lanchester() {
    // Lanchester's phugoid: wn = sqrt(2) g / u0, damping 1 / (sqrt(2) L/D).
    // The approximation needs a large relative density and a short period
    // well clear of the phugoid, so the default glider is ballasted 100x
    // (inertias too), its CoG moved 4 cm forward (static margin ~0.21), and
    // flown with constant-CD sections; FullLoad puts that drag in Xu.
    let (glider, mut db) = glider();
    for foil in ["NACA 0012", "NACA 0009"] {
        for re in [1e4, 1e5, 1e6, 1e7, 1e8] {
            db.add_polar(&synth_polar(foil, &format!("Re{re}"), re, 9.0, 0.01));
        }
    }
    let b = glider.inertia(&db);
    let k = 100.0;
    let inertia = foiley_plane::geom::Inertia { cog: Vec3::new(b.cog.x - 0.04, b.cog.y, b.cog.z), mass: b.mass * k, ixx: b.ixx * k, iyy: b.iyy * k, izz: b.izz * k, ixz: b.ixz * k, ..b };
    let t7 = PlanePolarDef { kind: PolarKind::Stability, viscous: true, stability_model: StabilityModel::FullLoad, inertia: Some(inertia), range: Range { from: 0.0, to: 0.0, step: 1.0 }, ..Default::default() };
    let o = run1(&glider, &db, &t7);
    let s = o.stability.as_ref().unwrap();
    let ph = s.modes.iter().find(|m| m.name == "Phugoid").unwrap_or_else(|| panic!("no phugoid: {:?}", s.modes.iter().map(|m| &m.name).collect::<Vec<_>>()));
    let (zeta, wn) = (1.0 / (2f64.sqrt() * o.ld()), 2f64.sqrt() * G / s.u0);
    assert!(o.ld() > 10.0 && o.ld() < 30.0, "L/D {}", o.ld());
    assert!((ph.damping / zeta - 1.0).abs() < 0.05, "phugoid damping {} vs Lanchester {zeta}", ph.damping);
    assert!((ph.eigenvalue.abs() / wn - 1.0).abs() < 0.05, "phugoid wn {} vs Lanchester {wn}", ph.eigenvalue.abs());
}
