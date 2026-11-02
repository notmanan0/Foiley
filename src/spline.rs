//! Cubic spline routines (port of XFoil spline.f).
//!
//! All arrays are 1-based: element 0 is unused, valid data in `1..=n`.

/// Tri-diagonal solve. a: diagonal, b: sub-diagonal, c: super-diagonal, d: rhs (overwritten).
/// Slices are 1-based with length >= kk+1.
pub fn trisol(a: &mut [f64], b: &[f64], c: &mut [f64], d: &mut [f64], kk: usize) {
    for k in 2..=kk {
        let km = k - 1;
        c[km] /= a[km];
        d[km] /= a[km];
        a[k] -= b[k] * c[km];
        d[k] -= b[k] * d[km];
    }
    d[kk] /= a[kk];
    for k in (1..kk).rev() {
        d[k] -= c[k] * d[k + 1];
    }
}

/// Spline with zero-2nd-derivative end conditions.
pub fn spline(x: &[f64], xs: &mut [f64], s: &[f64], n: usize) {
    splind(x, xs, s, n, 999.0, 999.0);
}

/// Spline with specified end conditions (999 = zero 2nd deriv, -999 = zero 3rd deriv).
pub fn splind(x: &[f64], xs: &mut [f64], s: &[f64], n: usize, xs1: f64, xs2: f64) {
    let mut a = vec![0.0; n + 1];
    let mut b = vec![0.0; n + 1];
    let mut c = vec![0.0; n + 1];
    for i in 2..n {
        let dsm = s[i] - s[i - 1];
        let dsp = s[i + 1] - s[i];
        b[i] = dsp;
        a[i] = 2.0 * (dsm + dsp);
        c[i] = dsm;
        xs[i] = 3.0 * ((x[i + 1] - x[i]) * dsm / dsp + (x[i] - x[i - 1]) * dsp / dsm);
    }
    if xs1 == 999.0 {
        a[1] = 2.0;
        c[1] = 1.0;
        xs[1] = 3.0 * (x[2] - x[1]) / (s[2] - s[1]);
    } else if xs1 == -999.0 {
        a[1] = 1.0;
        c[1] = 1.0;
        xs[1] = 2.0 * (x[2] - x[1]) / (s[2] - s[1]);
    } else {
        a[1] = 1.0;
        c[1] = 0.0;
        xs[1] = xs1;
    }
    if xs2 == 999.0 {
        b[n] = 1.0;
        a[n] = 2.0;
        xs[n] = 3.0 * (x[n] - x[n - 1]) / (s[n] - s[n - 1]);
    } else if xs2 == -999.0 {
        b[n] = 1.0;
        a[n] = 1.0;
        xs[n] = 2.0 * (x[n] - x[n - 1]) / (s[n] - s[n - 1]);
    } else {
        a[n] = 1.0;
        b[n] = 0.0;
        xs[n] = xs2;
    }
    if n == 2 && xs1 == -999.0 && xs2 == -999.0 {
        b[n] = 1.0;
        a[n] = 2.0;
        xs[n] = 3.0 * (x[n] - x[n - 1]) / (s[n] - s[n - 1]);
    }
    trisol(&mut a, &b, &mut c, xs, n);
}

/// Simple averaging spline (non-oscillatory).
pub fn splina(x: &[f64], xs: &mut [f64], s: &[f64], n: usize) {
    let mut lend = true;
    let mut xs1 = 0.0;
    let mut xs2 = 0.0;
    for i in 1..n {
        let ds = s[i + 1] - s[i];
        if ds == 0.0 {
            xs[i] = xs1;
            lend = true;
        } else {
            let dx = x[i + 1] - x[i];
            xs2 = dx / ds;
            if lend {
                xs[i] = xs2;
                lend = false;
            } else {
                xs[i] = 0.5 * (xs1 + xs2);
            }
        }
        xs1 = xs2;
    }
    xs[n] = xs1;
}

#[inline]
fn find_interval(ss: f64, s: &[f64], n: usize) -> usize {
    let mut ilow = 1usize;
    let mut i = n;
    while i - ilow > 1 {
        let imid = (i + ilow) / 2;
        if ss < s[imid] {
            i = imid;
        } else {
            ilow = imid;
        }
    }
    i
}

pub fn seval(ss: f64, x: &[f64], xs: &[f64], s: &[f64], n: usize) -> f64 {
    let i = find_interval(ss, s, n);
    let ds = s[i] - s[i - 1];
    let t = (ss - s[i - 1]) / ds;
    let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    let cx2 = ds * xs[i] - x[i] + x[i - 1];
    t * x[i] + (1.0 - t) * x[i - 1] + (t - t * t) * ((1.0 - t) * cx1 - t * cx2)
}

pub fn deval(ss: f64, x: &[f64], xs: &[f64], s: &[f64], n: usize) -> f64 {
    let i = find_interval(ss, s, n);
    let ds = s[i] - s[i - 1];
    let t = (ss - s[i - 1]) / ds;
    let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    let cx2 = ds * xs[i] - x[i] + x[i - 1];
    let d = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    d / ds
}

pub fn d2val(ss: f64, x: &[f64], xs: &[f64], s: &[f64], n: usize) -> f64 {
    let i = find_interval(ss, s, n);
    let ds = s[i] - s[i - 1];
    let t = (ss - s[i - 1]) / ds;
    let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    let cx2 = ds * xs[i] - x[i] + x[i - 1];
    ((6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2) / (ds * ds)
}

pub fn curv(ss: f64, x: &[f64], xs: &[f64], y: &[f64], ys: &[f64], s: &[f64], n: usize) -> f64 {
    let i = find_interval(ss, s, n);
    let ds = s[i] - s[i - 1];
    let t = (ss - s[i - 1]) / ds;
    let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    let cx2 = ds * xs[i] - x[i] + x[i - 1];
    let xd = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    let xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
    let cy1 = ds * ys[i - 1] - y[i] + y[i - 1];
    let cy2 = ds * ys[i] - y[i] + y[i - 1];
    let yd = y[i] - y[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2;
    let ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2;
    let sd = (xd * xd + yd * yd).sqrt().max(0.001 * ds);
    (xd * ydd - yd * xdd) / (sd * sd * sd)
}

pub fn curvs(ss: f64, x: &[f64], xs: &[f64], y: &[f64], ys: &[f64], s: &[f64], n: usize) -> f64 {
    let i = find_interval(ss, s, n);
    let ds = s[i] - s[i - 1];
    let t = (ss - s[i - 1]) / ds;
    let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    let cx2 = ds * xs[i] - x[i] + x[i - 1];
    let xd = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    let xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
    let xddd = 6.0 * cx1 + 6.0 * cx2;
    let cy1 = ds * ys[i - 1] - y[i] + y[i - 1];
    let cy2 = ds * ys[i] - y[i] + y[i - 1];
    let yd = y[i] - y[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2;
    let ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2;
    let yddd = 6.0 * cy1 + 6.0 * cy2;
    let sd = (xd * xd + yd * yd).sqrt().max(0.001 * ds);
    let bot = sd * sd * sd;
    let dbotdt = 3.0 * sd * (xd * xdd + yd * ydd);
    let top = xd * ydd - yd * xdd;
    let dtopdt = xd * yddd - yd * xddd;
    (dtopdt * bot - dbotdt * top) / (bot * bot)
}

/// Inverse spline: find s such that x(s) = xi, starting from guess `si`.
pub fn sinvrt(si: &mut f64, xi: f64, x: &[f64], xs: &[f64], s: &[f64], n: usize) {
    let sisav = *si;
    for _ in 0..10 {
        let res = seval(*si, x, xs, s, n) - xi;
        let resp = deval(*si, x, xs, s, n);
        let ds = -res / resp;
        *si += ds;
        if (ds / (s[n] - s[1])).abs() < 1.0e-5 {
            return;
        }
    }
    crate::outln!("SINVRT: spline inversion failed. Input value returned.");
    *si = sisav;
}

/// Arc length array.
pub fn scalc(x: &[f64], y: &[f64], s: &mut [f64], n: usize) {
    s[1] = 0.0;
    for i in 2..=n {
        s[i] = s[i - 1] + ((x[i] - x[i - 1]).powi(2) + (y[i] - y[i - 1]).powi(2)).sqrt();
    }
}

/// Segmented spline: derivative discontinuities where successive s values are identical.
pub fn segspl(x: &[f64], xs: &mut [f64], s: &[f64], n: usize) {
    segspld(x, xs, s, n, -999.0, -999.0);
}

pub fn segspld(x: &[f64], xs: &mut [f64], s: &[f64], n: usize, xs1: f64, xs2: f64) {
    if s[1] == s[2] {
        panic!("SEGSPL:  First input point duplicated");
    }
    if s[n] == s[n - 1] {
        panic!("SEGSPL:  Last  input point duplicated");
    }
    let mut iseg0 = 1usize;
    for iseg in 2..=n.saturating_sub(2) {
        if s[iseg] == s[iseg + 1] {
            let nseg = iseg - iseg0 + 1;
            let o = iseg0 - 1;
            splind(&x[o..], &mut xs[o..], &s[o..], nseg, xs1, xs2);
            iseg0 = iseg + 1;
        }
    }
    let nseg = n - iseg0 + 1;
    let o = iseg0 - 1;
    splind(&x[o..], &mut xs[o..], &s[o..], nseg, xs1, xs2);
}

/// Geometrically stretched array (SETEXP). s is 1-based local, s[1] = 0.
pub fn setexp(s: &mut [f64], ds1: f64, smax: f64, nn: usize) {
    let sigma = smax / ds1;
    let nex = nn - 1;
    let rnex = nex as f64;
    let rni = 1.0 / rnex;
    let aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0;
    let bbb = rnex * (rnex - 1.0) / 2.0;
    let ccc = rnex - sigma;
    let disc = (bbb * bbb - 4.0 * aaa * ccc).max(0.0);
    let mut ratio = if nex <= 1 {
        panic!("SETEXP: Cannot fill array.  N too small.");
    } else if nex == 2 {
        -ccc / bbb + 1.0
    } else {
        (-bbb + disc.sqrt()) / (2.0 * aaa) + 1.0
    };
    if ratio != 1.0 {
        let mut ok = false;
        for _ in 0..100 {
            let sigman = (ratio.powi(nex as i32) - 1.0) / (ratio - 1.0);
            let res = sigman.powf(rni) - sigma.powf(rni);
            let dresdr = rni * sigman.powf(rni) * (rnex * ratio.powi(nex as i32 - 1) - sigman) / (ratio.powi(nex as i32) - 1.0);
            let dratio = -res / dresdr;
            ratio += dratio;
            if dratio.abs() < 1.0e-5 {
                ok = true;
                break;
            }
        }
        if !ok {
            crate::outln!(" SETEXP: Convergence failed.  Continuing anyway ...");
        }
    }
    s[1] = 0.0;
    let mut ds = ds1;
    for n in 2..=nn {
        s[n] = s[n - 1] + ds;
        ds *= ratio;
    }
}
