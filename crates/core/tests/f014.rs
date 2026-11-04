//! F014: the blocked LU is cancellable within milliseconds and as accurate
//! as the single faer call.
use foiley_core::linalg::{DenseMatrix, Lu};
use foiley_core::Progress;
use std::time::{Duration, Instant};

/// Dense nonsymmetric test matrix with a moderately dominant diagonal
/// (pivoting still happens off the diagonal for small n).
fn matrix(n: usize, dominance: f64) -> DenseMatrix {
    let mut m = DenseMatrix::zeros(n, n);
    m.fill_columns_par(|j, c| {
        let mut s = (j as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        for (i, v) in c.iter_mut().enumerate() {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            *v = (s >> 11) as f64 / (1u64 << 53) as f64 - 0.5 + if i == j { dominance } else { 0.0 };
        }
    });
    m
}

fn inf(v: &[f64]) -> f64 {
    v.iter().fold(0.0, |m, x| m.max(x.abs()))
}

/// ||b - A x||_inf / (||A||_inf ||x||_inf)
fn backward_error(a: &DenseMatrix, x: &[f64], b: &[f64]) -> f64 {
    let ax = a.mul_vec(x);
    let r: Vec<f64> = b.iter().zip(&ax).map(|(b, y)| b - y).collect();
    let a_inf = (0..a.rows).map(|i| (0..a.cols).map(|j| a.at(i, j).abs()).sum::<f64>()).fold(0.0, f64::max);
    inf(&r) / (a_inf * inf(x))
}

#[test]
fn f014_blocked_solves_like_faer() {
    // small sizes and odd panel widths exercise every edge of the blocking
    for (n, nb) in [(1, 4), (7, 3), (64, 64), (300, 17), (513, 128)] {
        let a = matrix(n, 0.0);
        let b: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin() + 1.0).collect();
        let mut x0 = b.clone();
        Lu::factor(a.clone()).unwrap().solve(&mut x0);
        let mut x1 = b.clone();
        Lu::factor_blocked(a.clone(), &Progress::none(), nb, 1e6).unwrap().solve(&mut x1);
        let (e0, e1) = (backward_error(&a, &x0, &b), backward_error(&a, &x1, &b));
        assert!(e1 <= 2.0 * e0.max(1e-17), "n {n} nb {nb}: blocked {e1:e} vs faer {e0:e}");
        // transpose solve uses the same factors
        let mut y = b.clone();
        let lu = Lu::factor_blocked(a.clone(), &Progress::none(), nb, 1e6).unwrap();
        lu.solve_transpose(&mut y);
        let mut at = DenseMatrix::zeros(n, n);
        (0..n).for_each(|i| (0..n).for_each(|j| at.set(i, j, a.at(j, i))));
        assert!(backward_error(&at, &y, &b) < 1e-13, "n {n}: transpose");
    }
}

#[test]
fn f014_backward_error_matches_faer() {
    let n = 2600;
    let a = matrix(n, 0.0);
    let b: Vec<f64> = (0..n).map(|i| (i as f64 * 0.11).cos()).collect();
    let mono = Lu::factor_blocked(a.clone(), &Progress::none(), n, f64::INFINITY).unwrap(); // one panel = one faer call
    let blocked = Lu::factor_with(a.clone(), &Progress::none()).unwrap();
    let (mut x0, mut x1) = (b.clone(), b.clone());
    mono.solve(&mut x0);
    blocked.solve(&mut x1);
    let (e0, e1) = (backward_error(&a, &x0, &b), backward_error(&a, &x1, &b));
    assert!(e1 <= 2.0 * e0, "blocked {e1:e} vs monolithic {e0:e}");
    assert!(e1 < 1e-13, "{e1:e}");
}

/// Cancel 1 s into a 14,100-unknown factorization (the size of the
/// TriLinear fine wing) ten times; the 95th-percentile time from cancel()
/// to the call returning must be at most 100 ms.
#[test]
fn f014_cancel_latency() {
    let n = 14_100;
    let base = matrix(n, 0.1 * n as f64);
    let mut lat = Vec::new();
    for _ in 0..10 {
        let a = base.clone();
        let prog = Progress::new();
        let p2 = prog.clone();
        let h = std::thread::spawn(move || matches!(Lu::factor_with(a, &p2), Err(foiley_core::linalg::LinalgError::Cancelled)));
        std::thread::sleep(Duration::from_secs(1));
        let t = Instant::now();
        prog.cancel();
        assert!(h.join().unwrap(), "factorization finished or failed before the cancel");
        lat.push(t.elapsed().as_secs_f64() * 1e3);
    }
    lat.sort_by(f64::total_cmp);
    let p95 = lat[(0.95 * (lat.len() - 1) as f64).ceil() as usize];
    println!("cancel latency ms: {lat:.1?} (p95 {p95:.1})");
    assert!(p95 <= 100.0, "p95 cancel latency {p95:.1} ms: {lat:?}");
}
