//! Interpolation helpers.

/// Linear interpolation in a table sorted by `xs` (clamped at the ends).
pub fn lerp_table(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let (i0, i1, t) = bracket(xs, x);
    if xs.is_empty() {
        return 0.0;
    }
    ys[i0] + t * (ys[i1] - ys[i0])
}

/// Bracket `x` in sorted `xs`: returns (i0, i1, t) with x = (1-t) xs[i0] + t xs[i1]
/// (clamped, so t in [0,1]).
pub fn bracket(xs: &[f64], x: f64) -> (usize, usize, f64) {
    let n = xs.len();
    if n <= 1 || x <= xs[0] {
        return (0, 0, 0.0);
    }
    if x >= xs[n - 1] {
        return (n - 1, n - 1, 0.0);
    }
    // binary search for first xs[i] >= x
    let i = xs.partition_point(|v| *v < x).clamp(1, n - 1);
    let (x0, x1) = (xs[i - 1], xs[i]);
    (i - 1, i, if x1 != x0 { (x - x0) / (x1 - x0) } else { 0.0 })
}

/// Natural cubic spline through (xs, ys) with `xs` strictly increasing.
#[derive(Clone, Debug, Default)]
pub struct CubicSpline {
    xs: Vec<f64>,
    ys: Vec<f64>,
    m: Vec<f64>,
}

impl CubicSpline {
    pub fn new(xs: &[f64], ys: &[f64]) -> Self {
        let n = xs.len();
        let mut m = vec![0.0; n];
        if n > 2 {
            let mut a = vec![0.0; n];
            let mut b = vec![1.0; n];
            let mut c = vec![0.0; n];
            let mut d = vec![0.0; n];
            for i in 1..n - 1 {
                let h0 = xs[i] - xs[i - 1];
                let h1 = xs[i + 1] - xs[i];
                a[i] = h0;
                b[i] = 2.0 * (h0 + h1);
                c[i] = h1;
                d[i] = 6.0 * ((ys[i + 1] - ys[i]) / h1 - (ys[i] - ys[i - 1]) / h0);
            }
            for i in 1..n {
                let w = a[i] / b[i - 1];
                b[i] -= w * c[i - 1];
                d[i] -= w * d[i - 1];
            }
            m[n - 1] = d[n - 1] / b[n - 1];
            for i in (0..n - 1).rev() {
                m[i] = (d[i] - c[i] * m[i + 1]) / b[i];
            }
        }
        CubicSpline { xs: xs.to_vec(), ys: ys.to_vec(), m }
    }
    pub fn eval(&self, x: f64) -> f64 {
        let n = self.xs.len();
        if n == 0 {
            return 0.0;
        }
        if n == 1 {
            return self.ys[0];
        }
        let i1 = self.xs.partition_point(|v| *v < x).clamp(1, n - 1);
        let i0 = i1 - 1;
        let h = self.xs[i1] - self.xs[i0];
        let a = (self.xs[i1] - x) / h;
        let b = (x - self.xs[i0]) / h;
        a * self.ys[i0] + b * self.ys[i1] + ((a * a * a - a) * self.m[i0] + (b * b * b - b) * self.m[i1]) * h * h / 6.0
    }
}

/// `n+1` points distributed over [0,1] by a named law.
pub fn distribution(kind: crate::Spacing, n: usize) -> Vec<f64> {
    use std::f64::consts::PI;
    let n = n.max(1);
    (0..=n)
        .map(|i| {
            let t = i as f64 / n as f64;
            match kind {
                crate::Spacing::Uniform => t,
                crate::Spacing::Cosine => 0.5 * (1.0 - (t * PI).cos()),
                crate::Spacing::Sine => (t * PI / 2.0).sin(),
                crate::Spacing::MinusSine => 1.0 - ((1.0 - t) * PI / 2.0).sin(),
                crate::Spacing::Tanh => {
                    let k = 2.0;
                    0.5 * (1.0 + (k * (2.0 * t - 1.0)).tanh() / k.tanh())
                }
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn table() {
        let xs = [0.0, 1.0, 2.0];
        let ys = [0.0, 10.0, 0.0];
        assert_eq!(lerp_table(&xs, &ys, 0.5), 5.0);
        assert_eq!(lerp_table(&xs, &ys, 5.0), 0.0);
        assert_eq!(lerp_table(&xs, &ys, -1.0), 0.0);
        let s = CubicSpline::new(&[0.0, 1.0, 2.0, 3.0], &[0.0, 1.0, 4.0, 9.0]);
        assert!((s.eval(1.5) - 2.25).abs() < 0.1);
        let d = distribution(crate::Spacing::Cosine, 4);
        assert!((d[2] - 0.5).abs() < 1e-12 && d[4] == 1.0);
    }
}
