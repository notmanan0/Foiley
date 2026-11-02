//! Regression tests: Foiley polars must reproduce the reference polars
//! produced by double-precision Fortran XFoil 6.97 (tests/data/*.pol) to
//! every printed digit, for both the sequential (pipelined) sweep and the
//! plain point-by-point path.  The parallel sweep must agree to within the
//! Newton convergence tolerance.

use foiley::fmt::f_fmt;
use foiley::sweep::{run_parallel, run_sequential, sequence_specs, solve_point};
use foiley::xfoil::XFoil;

fn setup(naca: i32, re: f64) -> XFoil {
    foiley::set_quiet(true);
    let mut xf = XFoil::new();
    xf.naca(naca);
    xf.lvisc = true;
    xf.reinf1 = re;
    let mut r = 0.0;
    xf.mrcl(1.0, &mut r);
    xf.itmax = 100;
    xf
}

/// XFoil polar data line: 1X,F7.3,F9.4,F10.5,F10.5,F9.4,F9.4,F9.4
fn polar_text(pts: &[foiley::sweep::SeqPoint]) -> Vec<String> {
    pts.iter()
        .filter(|p| p.converged)
        .map(|p| {
            let q = &p.point;
            format!(
                " {}{}{}{}{}{}{}",
                f_fmt(q.alpha, 7, 3),
                f_fmt(q.cl, 9, 4),
                f_fmt(q.cd, 10, 5),
                f_fmt(q.cdp, 10, 5),
                f_fmt(q.cm, 9, 4),
                f_fmt(q.xtr_top, 9, 4),
                f_fmt(q.xtr_bot, 9, 4)
            )
        })
        .collect()
}

fn reference_lines(file: &str) -> Vec<String> {
    let txt = std::fs::read_to_string(format!("{}/tests/data/{}", env!("CARGO_MANIFEST_DIR"), file)).unwrap();
    let mut out = Vec::new();
    let mut data = false;
    for l in txt.lines() {
        if l.trim_start().starts_with("---") {
            data = true;
            continue;
        }
        if data && !l.trim().is_empty() {
            out.push(l.trim_end().to_string());
        }
    }
    out
}

fn check(naca: i32, re: f64, a: f64, b: f64, d: f64, file: &str) {
    let refl = reference_lines(file);
    let specs = sequence_specs(true, a, b, d);

    // pipelined sequential sweep
    let mut xf = setup(naca, re);
    let pts = run_sequential(&mut xf, true, &specs, 105, |_, _| {});
    assert_eq!(polar_text(&pts), refl, "sequential sweep differs from Fortran reference ({file})");

    // plain point-by-point loop
    let mut xf = setup(naca, re);
    let pts: Vec<_> = specs.iter().map(|&s| solve_point(&mut xf, true, s, 105)).collect();
    assert_eq!(polar_text(&pts), refl, "point loop differs from Fortran reference ({file})");

    // parallel sweep: within convergence tolerance
    let mut xf = setup(naca, re);
    let par = run_parallel(&mut xf, true, &specs, 105, 4, |_, _| {});
    let r: Vec<Vec<f64>> = refl.iter().map(|l| l.split_whitespace().map(|t| t.parse().unwrap()).collect()).collect();
    let p: Vec<_> = par.iter().filter(|p| p.converged).collect();
    assert_eq!(p.len(), r.len());
    for (u, v) in p.iter().zip(r.iter()) {
        assert!((u.point.alpha - v[0]).abs() < 1e-9);
        assert!((u.point.cl - v[1]).abs() < 2e-4, "CL {} vs {}", u.point.cl, v[1]);
        assert!((u.point.cd - v[2]).abs() < 2e-5, "CD {} vs {}", u.point.cd, v[2]);
        assert!((u.point.cm - v[4]).abs() < 2e-4, "CM {} vs {}", u.point.cm, v[4]);
    }
}

#[test]
fn naca0012_re1e6() {
    check(12, 1e6, 0.0, 10.0, 1.0, "ref_naca0012_re1000000.pol");
}

#[test]
fn naca2412_re3e5() {
    check(2412, 3e5, -4.0, 8.0, 1.0, "ref_naca2412_re300000.pol");
}

#[test]
fn naca4415_re3e6() {
    check(4415, 3e6, 0.0, 12.0, 2.0, "ref_naca4415_re3000000.pol");
}
