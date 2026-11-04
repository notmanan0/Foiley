//! F013: numerical diagnostics (condition estimate, residuals, balances).
mod common;
use common::*;
use foiley_core::linalg::{row_residuals, DenseMatrix, Diagnostics, Lu};
use foiley_core::{Progress, Vec3};
use foiley_plane::analysis::*;
use foiley_plane::mesh::{build, MeshOptions};
use foiley_plane::vlm::{Lattice, VlmSystem};

/// ||M||_1 (max column sum).
fn norm1(m: &DenseMatrix) -> f64 {
    (0..m.cols).map(|j| m.col(j).iter().map(|v| v.abs()).sum::<f64>()).fold(0.0, f64::max)
}

/// The VLM influence matrix, assembled independently of the solver.
fn vlm_matrix(lat: &Lattice) -> DenseMatrix {
    let mesh = lat.mesh;
    let n = mesh.panels.len();
    let mut a = DenseMatrix::zeros(n, n);
    for i in 0..n {
        let p = &mesh.panels[i];
        for j in 0..n {
            a.set(i, j, lat.element_velocity(j, p.ctrl, false).dot(p.normal));
        }
    }
    a
}

fn rect200() -> (foiley_plane::geom::Plane, foiley_plane::FoilDb) {
    let f = naca0012();
    let mut p = review_rect(&f);
    p.wings[0].sections[0].nx = 5;
    p.wings[0].sections[0].ny = 20;
    (p, db_with(&f))
}

#[test]
fn f013_condition_estimate_within_10x() {
    let (p, db) = rect200();
    let mesh = build(&p, &db, &MeshOptions::default()).unwrap();
    assert_eq!(mesh.panels.len(), 200);
    let sys = VlmSystem::build_with(&mesh, true, Vec3::ZERO, None, &Progress::none(), Diagnostics::Cheap).unwrap();
    let est = sys.report.cond1_estimate.unwrap();
    // exact kappa_1 from an explicit inverse
    let a = vlm_matrix(&Lattice::new(&mesh, true, None));
    let n = a.rows;
    let na = norm1(&a);
    let mut inv = DenseMatrix::zeros(n, n);
    for i in 0..n {
        inv.set(i, i, 1.0);
    }
    Lu::factor(a).unwrap().solve_many(&mut inv);
    let exact = na * norm1(&inv);
    assert!(est <= exact * (1.0 + 1e-8), "estimate {est} above exact {exact}");
    assert!(est >= exact / 10.0, "estimate {est} not within 10x of exact {exact}");
}

#[test]
fn f013_condition_estimate_hard_matrix() {
    // Hilbert-like matrix: kappa grows fast; the estimate must track it
    let n = 10;
    let a = DenseMatrix::from_rows(&(0..n).map(|i| (0..n).map(|j| 1.0 / (i + j + 1) as f64).collect()).collect::<Vec<_>>());
    let na = norm1(&a);
    let mut inv = DenseMatrix::zeros(n, n);
    (0..n).for_each(|i| inv.set(i, i, 1.0));
    let lu = Lu::factor(a).unwrap();
    lu.solve_many(&mut inv);
    let exact = na * norm1(&inv);
    let est = lu.cond1_estimate();
    assert!(exact > 1e12 && est <= exact * 1.001 && est >= exact / 10.0, "est {est} exact {exact}");
}

#[test]
fn f013_full_residual_uses_f64_operator() {
    // A rounded (f32) copy of A gives an exactly zero residual for the x that
    // solves the rounded system, while the true f64 operator does not.
    let n = 40;
    let a64 = DenseMatrix::from_rows(&(0..n).map(|i| (0..n).map(|j| 1.0 / (i + j + 1) as f64 + if i == j { 1.0 } else { 0.0 }).collect()).collect::<Vec<_>>());
    let mut a32 = a64.clone();
    a32.data.iter_mut().for_each(|v| *v = *v as f32 as f64);
    let x = DenseMatrix { rows: n, cols: 1, data: (0..n).map(|i| ((i * 7) % 17) as f64 - 8.0).collect() };
    let b = DenseMatrix { rows: n, cols: 1, data: a32.mul_vec(&x.data) };
    let all: Vec<usize> = (0..n).collect();
    let r32 = row_residuals(&a32, &all, &x, &b);
    let r64 = row_residuals(&a64, &all, &x, &b);
    assert!(r32.r[0].iter().all(|v| *v == 0.0), "rounded operator should hide the error");
    assert!(r64.r[0].iter().any(|v| v.abs() > 1e-9), "f64 operator must expose it");

    // Full mode on a real system: every row checked, and the residual equals
    // the one from an independently assembled f64 matrix
    let (p, db) = rect200();
    let mesh = build(&p, &db, &MeshOptions::default()).unwrap();
    let sys = VlmSystem::build_with(&mesh, true, Vec3::new(1.0, 0.0, 0.0), None, &Progress::none(), Diagnostics::Full).unwrap();
    let rep = &sys.report;
    assert_eq!(rep.residuals.rows.len(), mesh.panels.len());
    let a = vlm_matrix(&Lattice::new(&mesh, true, None));
    let m = mesh.panels.len();
    let x = DenseMatrix { rows: m, cols: 6, data: sys.unit.iter().flatten().copied().collect() };
    let bm = DenseMatrix { rows: m, cols: 6, data: (0..6).flat_map(|k| rep.residuals.b[k].clone()).collect() };
    let rows: Vec<usize> = (0..m).collect();
    let indep = row_residuals(&a, &rows, &x, &bm);
    for k in 0..6 {
        for (u, v) in indep.r[k].iter().zip(&rep.residuals.r[k]) {
            assert_eq!(u, v, "residual not from the exact f64 operator");
        }
    }
}

#[test]
fn f013_refinement_never_raises_residual() {
    let (p, db) = rect200();
    let mesh = build(&p, &db, &MeshOptions::default()).unwrap();
    let cheap = VlmSystem::build_with(&mesh, true, Vec3::ZERO, None, &Progress::none(), Diagnostics::Cheap).unwrap();
    let full = VlmSystem::build_with(&mesh, true, Vec3::ZERO, None, &Progress::none(), Diagnostics::Full).unwrap();
    let off = VlmSystem::build(&mesh, true, Vec3::ZERO, None, &Progress::none()).unwrap();
    // cheap diagnostics leave the solution bit-identical
    assert_eq!(cheap.unit, off.unit);
    let a = vlm_matrix(&Lattice::new(&mesh, true, None));
    let m = mesh.panels.len();
    let rows: Vec<usize> = (0..m).collect();
    let b = DenseMatrix { rows: m, cols: 6, data: (0..6).flat_map(|k| full.report.residuals.b[k].clone()).collect() };
    let res = |u: &[Vec<f64>; 6]| row_residuals(&a, &rows, &DenseMatrix { rows: m, cols: 6, data: u.iter().flatten().copied().collect() }, &b);
    let (r0, r1) = (res(&off.unit), res(&full.unit));
    let inf = |v: &[f64]| v.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    for k in 0..6 {
        assert!(inf(&r1.r[k]) <= inf(&r0.r[k]), "refinement raised residual {k}");
    }
}

#[test]
fn g05_backward_error_reported() {
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    for method in [Method::Vlm1, Method::Vlm2, Method::Panel, Method::TriPanel, Method::TriLinear] {
        let d = PlanePolarDef { method, viscous: false, alpha: 4.0, range: Range { from: 4.0, to: 4.0, step: 1.0 }, ..Default::default() };
        let off = run1(&p, &db, &d);
        assert!(off.diagnostics.is_none(), "diagnostics must be off by default");
        for level in [Diagnostics::Cheap, Diagnostics::Full] {
            let o = run1(&p, &db, &PlanePolarDef { diagnostics: level, ..d.clone() });
            let g = o.diagnostics.as_ref().unwrap_or_else(|| panic!("{method:?} {level:?}: no diagnostics"));
            let be = g.backward_error.unwrap();
            assert!(be < 1e-12, "{method:?} {level:?}: backward error {be}");
            assert!(g.cond1_estimate.unwrap() >= 1.0);
            assert!(g.diag_ratio > 0.0 && g.diag_ratio <= 1.0);
            if level == Diagnostics::Full {
                assert_eq!(g.rows_checked, g.unknowns);
            } else {
                assert!(g.rows_checked >= 64.min(g.unknowns));
            }
            assert!(g.cd_near_minus_trefftz.is_finite() && g.cl_pressure_minus_wake.is_finite());
            if level == Diagnostics::Cheap {
                if method == Method::TriLinear {
                    // its RHS is a parallel reduction: last-bit run-to-run noise (F015)
                    assert!((o.cl / off.cl - 1.0).abs() < 1e-12, "{method:?}: cheap diagnostics changed the answer");
                } else {
                    assert_eq!(o.cl, off.cl, "{method:?}: cheap diagnostics changed the answer");
                }
            }
        }
    }
    // a relaxed-wake point is not checked against the rigid operator
    let o = run1(&p, &db, &PlanePolarDef { diagnostics: Diagnostics::Cheap, ..review_def() });
    assert!(o.diagnostics.unwrap().backward_error.is_none());
}

/// The near-field vs wake balances expose a coarse panel mesh and shrink
/// as it is refined; on a VLM the bound and shed Kutta-Joukowski lifts are
/// the same circulation, and the freestream near field carries no induced drag.
#[test]
fn f013_balances_track_mesh_quality() {
    let f = naca0012();
    let db = db_with(&f);
    let at = |method, nx| {
        let mut p = review_rect(&f);
        p.wings[0].sections[0].nx = nx;
        let d = PlanePolarDef { method, viscous: false, alpha: 4.0, range: Range { from: 4.0, to: 4.0, step: 1.0 }, diagnostics: Diagnostics::Cheap, ..Default::default() };
        let o = run1(&p, &db, &d);
        (o.cl, o.cdi, o.diagnostics.unwrap())
    };
    let (_, _, coarse) = at(Method::Panel, 3);
    let (cl, _, fine) = at(Method::Panel, 16);
    assert!(coarse.cl_pressure_minus_wake.abs() > 1.0, "coarse panel mesh should show a large lift imbalance");
    assert!(fine.cl_pressure_minus_wake.abs() < 0.05 * cl, "fine mesh imbalance {} of CL {cl}", fine.cl_pressure_minus_wake);
    assert!(fine.cl_pressure_minus_wake.abs() < coarse.cl_pressure_minus_wake.abs() / 10.0);
    let (_, cdi, vlm) = at(Method::Vlm2, 8);
    assert!(vlm.cl_pressure_minus_wake.abs() < 1e-9, "{}", vlm.cl_pressure_minus_wake);
    assert!((vlm.cd_near_minus_trefftz + cdi).abs() < 1e-9, "{} vs -{cdi}", vlm.cd_near_minus_trefftz);
}
