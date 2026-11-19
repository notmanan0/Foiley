//! Cubic splines used by the NACA 6-series thickness forms: the
//! Forsythe-Malcolm-Moler end conditions and piecewise-cubic Hermite
//! evaluation of NACA456's `SplineProcedures` module (public domain).

/// Spline slopes `yp` at the knots `x` (strictly increasing) through `y`,
/// with FMM end conditions: the third derivative at each end matches the
/// cubic through the four end points.
pub(super) fn fmm_slopes(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    assert_eq!(n, y.len());
    if n < 2 {
        return vec![0.0; n];
    }
    let dx: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
    let delta: Vec<f64> = (0..n - 1).map(|i| (y[i + 1] - y[i]) / dx[i]).collect();
    if n == 2 {
        return vec![delta[0]; 2];
    }
    let dd: Vec<f64> = (0..n - 2).map(|i| delta[i + 1] - delta[i]).collect();
    if n == 3 {
        let d2 = dd[0] / (x[2] - x[0]);
        let d1 = delta[0] - d2 * dx[0];
        // p'(x) = f[x0,x1] + f[x0,x1,x2] (2x - x0 - x1)
        return vec![d1, d1 + 2.0 * d2 * dx[0], d1 + 2.0 * d2 * (x[2] - x[0])];
    }
    let mut alpha = vec![0.0; n];
    let mut beta = vec![0.0; n];
    let mut sigma = vec![0.0; n];
    alpha[0] = -dx[0];
    for i in 1..n - 1 {
        alpha[i] = 2.0 * (dx[i - 1] + dx[i]);
    }
    for i in 1..n - 1 {
        alpha[i] -= dx[i - 1] * dx[i - 1] / alpha[i - 1];
    }
    alpha[n - 1] = -dx[n - 2] - dx[n - 2] * dx[n - 2] / alpha[n - 2];

    beta[0] = dd[1] / (x[3] - x[1]) - dd[0] / (x[2] - x[0]);
    beta[0] = beta[0] * dx[0] * dx[0] / (x[3] - x[0]);
    beta[1..(n - 1)].copy_from_slice(&dd[..(n - 2)]);
    beta[n - 1] = dd[n - 3] / (x[n - 1] - x[n - 3]) - dd[n - 4] / (x[n - 2] - x[n - 4]);
    beta[n - 1] = -beta[n - 1] * dx[n - 2] * dx[n - 2] / (x[n - 1] - x[n - 4]);
    for i in 1..n {
        beta[i] -= dx[i - 1] * beta[i - 1] / alpha[i - 1];
    }
    sigma[n - 1] = beta[n - 1] / alpha[n - 1];
    for i in (0..n - 1).rev() {
        sigma[i] = (beta[i] - dx[i] * sigma[i + 1]) / alpha[i];
    }
    let mut yp = vec![0.0; n];
    for i in 0..n - 1 {
        yp[i] = delta[i] - dx[i] * (sigma[i] + sigma[i] + sigma[i + 1]);
    }
    yp[n - 1] = yp[n - 2] + dx[n - 2] * 3.0 * (sigma[n - 1] + sigma[n - 2]);
    yp
}

/// Value, first and second derivative of the Hermite cubic on [a, b].
fn hermite(a: f64, fa: f64, fpa: f64, b: f64, fb: f64, fpb: f64, u: f64) -> (f64, f64, f64) {
    let h = b - a;
    let (r1, r2, r3, r4) = (fa, fb, fpa * h, fpb * h);
    // coefficients of t^3, t^2, t, 1 (NACA456 MAGIC matrix)
    let c1 = 2.0 * r1 - 2.0 * r2 + r3 + r4;
    let c2 = -3.0 * r1 + 3.0 * r2 - 2.0 * r3 - r4;
    let c3 = r3;
    let c4 = r1;
    let hi = 1.0 / h;
    let t = (u - a) * hi;
    let f = c4 + t * (c3 + t * (c2 + t * c1));
    let fp = hi * (c3 + t * (2.0 * c2 + t * 3.0 * c1));
    let fpp = hi * hi * (2.0 * c2 + t * 6.0 * c1);
    (f, fp, fpp)
}

/// Interval k with x[k] <= u < x[k+1], clamped to the table.
fn interval(x: &[f64], u: f64) -> usize {
    let n = x.len();
    let k = x.partition_point(|&v| v <= u);
    k.saturating_sub(1).min(n - 2)
}

/// Evaluate the spline (x, y, yp) at u: value, slope and curvature.
pub(super) fn eval(x: &[f64], y: &[f64], yp: &[f64], u: f64) -> (f64, f64, f64) {
    let k = interval(x, u);
    hermite(x[k], y[k], yp[k], x[k + 1], y[k + 1], yp[k + 1], u)
}

/// Parameter u where the spline (s, f, fp) takes the value `fbar`
/// (first crossing), or `None` if it never does.
pub(super) fn solve(s: &[f64], f: &[f64], fp: &[f64], fbar: f64) -> Option<f64> {
    let n = s.len();
    for k in 0..n {
        if (f[k] - fbar).abs() < 1e-14 {
            return Some(s[k]);
        }
    }
    let k = (1..n).find(|&k| (f[k - 1] - fbar) * (f[k] - fbar) < 0.0)?;
    let (a, b) = (s[k - 1], s[k]);
    // safeguarded Newton on the bracketing cubic
    let (mut lo, mut hi) = (a, b);
    let mut glo = f[k - 1] - fbar;
    let mut u = a + (b - a) * (fbar - f[k - 1]) / (f[k] - f[k - 1]);
    for _ in 0..100 {
        let (v, d, _) = hermite(a, f[k - 1], fp[k - 1], b, f[k], fp[k], u);
        let v = v - fbar;
        if v == 0.0 {
            return Some(u);
        }
        if (v < 0.0) == (glo < 0.0) {
            lo = u;
            glo = v;
        } else {
            hi = u;
        }
        let mut next = if d != 0.0 { u - v / d } else { f64::NAN };
        if !(next > lo && next < hi) {
            next = 0.5 * (lo + hi);
        }
        if (next - u).abs() <= 1e-15 * (1.0 + u.abs()) || (hi - lo) <= 1e-15 * (1.0 + u.abs()) {
            return Some(next);
        }
        u = next;
    }
    Some(u)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fmm_spline_reproduces_a_cubic() {
        let x: Vec<f64> = (0..9).map(|i| (i as f64 * 0.37).powf(1.3)).collect();
        let f = |t: f64| 1.0 - 2.0 * t + 0.5 * t * t - 0.25 * t * t * t;
        let y: Vec<f64> = x.iter().map(|&t| f(t)).collect();
        let yp = fmm_slopes(&x, &y);
        for &u in &[0.05, 0.4, 1.1, 2.2] {
            let (v, d, _) = eval(&x, &y, &yp, u);
            assert!((v - f(u)).abs() < 1e-12, "{u}: {v} vs {}", f(u));
            assert!((d - (-2.0 + u - 0.75 * u * u)).abs() < 1e-10);
        }
        let u = solve(&x, &y, &yp, f(1.3)).unwrap();
        assert!((u - 1.3).abs() < 1e-12, "{u}");
        // three knots: the quadratic through them
        let (x, y) = ([0.0, 1.0, 2.0], [0.0, 1.0, 4.0]);
        let yp = fmm_slopes(&x, &y);
        assert_eq!(yp, vec![0.0, 2.0, 4.0]);
        assert!((eval(&x, &y, &yp, 0.5).0 - 0.25).abs() < 1e-15);
    }
}
