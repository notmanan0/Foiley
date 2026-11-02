//! Geometry helpers used by GDES (ports of routines from xgeom.f, xtcam.f,
//! sort.f, spline.f and xfoil.f that are not already in geom.rs/spline.rs).
//!
//! All arrays are 1-based like the Fortran (index 0 unused).

use crate::geom::{atanc, lefind, sopps};
use crate::spline::{curv, d2val, deval, scalc, segspl, seval};

// ---------------------------------------------------------------------------
// Fortran-style formatting helpers
// ---------------------------------------------------------------------------

/// Fortran Iw (stars on overflow).
pub fn i_fmt(v: i64, w: usize) -> String {
    let s = format!("{}", v);
    if s.len() > w {
        "*".repeat(w)
    } else {
        format!("{:>w$}", s, w = w)
    }
}

/// gfortran list-directed INTEGER item (leading separator included).
pub fn lint(v: i64) -> String {
    format!("{:>12}", v)
}

/// gfortran list-directed REAL*8 item (leading separator included, 26 chars).
pub fn lreal(v: f64) -> String {
    if !v.is_finite() {
        return format!("{:>26}", if v.is_nan() { "NaN" } else if v > 0.0 { "Infinity" } else { "-Infinity" });
    }
    if v == 0.0 {
        let s = if v.is_sign_negative() { "-0.0000000000000000" } else { "0.0000000000000000" };
        return format!("{:>21}     ", s);
    }
    let a = v.abs();
    if (0.1..1.0e17).contains(&a) {
        // G editing with 17 significant digits: F form with 17-k decimals
        let mut k = a.log10().floor() as i32 + 1;
        let dec = (17 - k).max(0) as usize;
        let r = format!("{:.*}", dec, a);
        if r.parse::<f64>().unwrap_or(a) >= 10f64.powi(k) {
            k += 1;
        }
        if k <= 17 {
            let dec = (17 - k).max(0) as usize;
            let s = format!("{:.*}", dec, v);
            return format!("{:>21}     ", s);
        }
    }
    let e = format!("{:.16e}", v); // d.dddde-14
    let (m, ex) = e.split_once('e').unwrap();
    let exv: i32 = ex.parse().unwrap_or(0);
    let s = format!("{}E{}{:03}", m, if exv < 0 { '-' } else { '+' }, exv.abs());
    format!("{:>26}", s)
}

// ---------------------------------------------------------------------------
// xfoil.f
// ---------------------------------------------------------------------------

/// ROTATE: rotate coordinates about origin by angle ALFA (rad).
pub fn rotate(x: &mut [f64], y: &mut [f64], n: usize, alfa: f64) {
    let sa = alfa.sin();
    let ca = alfa.cos();
    let xoff = 0.0;
    let yoff = 0.0;
    for i in 1..=n {
        let xt = x[i];
        let yt = y[i];
        x[i] = ca * xt + sa * yt + xoff;
        y[i] = ca * yt - sa * xt + yoff;
    }
}

// ---------------------------------------------------------------------------
// xgeom.f
// ---------------------------------------------------------------------------

/// XLFIND: leftmost (minimum x) point location.
pub fn xlfind(x: &[f64], xp: &[f64], _y: &[f64], _yp: &[f64], s: &[f64], n: usize) -> f64 {
    let dslen = s[n] - s[1];
    let dseps = (s[n] - s[1]) * 1.0e-5;
    let mut i = 3;
    while i + 2 <= n {
        let dx = x[i + 1] - x[i];
        if dx > 0.0 {
            break;
        }
        i += 1;
    }
    let mut sle = s[i];
    if s[i] == s[i - 1] {
        return sle;
    }
    for _ in 0..50 {
        let dxds = deval(sle, x, xp, s, n);
        let dxdd = d2val(sle, x, xp, s, n);
        let res = dxds;
        let ress = dxdd;
        let mut dsle = -res / ress;
        dsle = dsle.max(-0.01 * dslen.abs());
        dsle = dsle.min(0.01 * dslen.abs());
        sle += dsle;
        if dsle.abs() < dseps {
            return sle;
        }
    }
    crate::outln!(" XLFIND:  Left point not found.  Continuing...");
    s[i]
}

/// YSYM: make airfoil symmetric about chord line using side ISIDE.
/// Returns the new number of points.
#[allow(clippy::too_many_arguments)]
pub fn ysym(x: &mut [f64], xp: &mut [f64], y: &mut [f64], yp: &mut [f64], s: &mut [f64], nx: usize, n: usize, iside: i32) -> usize {
    let sle = lefind(x, xp, y, yp, s, n);
    let xle = seval(sle, x, xp, s, n);
    let yle = seval(sle, y, yp, s, n);
    let xte = 0.5 * (x[1] + x[n]);
    let yte = 0.5 * (y[1] + y[n]);
    let chsq = (xte - xle).powi(2) + (yte - yle).powi(2);
    let dxc = (xte - xle) / chsq.sqrt();
    let dyc = (yte - yle) / chsq.sqrt();

    let mut i = 2;
    while i <= n {
        let ds = s[i] - s[i - 1];
        if s[i] - sle >= -0.01 * ds {
            break;
        }
        i += 1;
    }
    let ile = i - 1;

    let ds = s[ile + 1] - s[ile];
    let (ile1, ile2) = if sle - s[ile - 1] < 0.1 * ds {
        (ile - 1, ile + 1)
    } else if s[ile + 1] - sle < 0.1 * ds {
        (ile, ile + 2)
    } else {
        (ile, ile + 1)
    };

    let (ig1, ig2, igdir): (i64, i64, i64) = if iside == 1 { (1, ile1 as i64, 1) } else { (n as i64, ile2 as i64, -1) };

    let nnew = (2 * ((ig2 - ig1).abs() + 1) + 1) as usize;
    if nnew > nx {
        panic!("YSYM:  Array overflow on passed arrays.");
    }
    let mut xnew = vec![0.0; nnew + 2];
    let mut ynew = vec![0.0; nnew + 2];

    let mut i = ig1;
    loop {
        let iu = i as usize;
        let xbar = (x[iu] - xle) * dxc + (y[iu] - yle) * dyc;
        let ybar = (y[iu] - yle) * dxc - (x[iu] - xle) * dyc;
        let i1 = (1 + (i - ig1) * igdir) as usize;
        let i2 = (nnew as i64 - (i - ig1) * igdir) as usize;
        xnew[i1] = xle + xbar * dxc - ybar * dyc;
        xnew[i2] = xle + xbar * dxc + ybar * dyc;
        ynew[i1] = yle + ybar * dxc + xbar * dyc;
        ynew[i2] = yle - ybar * dxc + xbar * dyc;
        if i == ig2 {
            break;
        }
        i += igdir;
    }
    xnew[nnew / 2 + 1] = xle;
    ynew[nnew / 2 + 1] = yle;

    let n = nnew;
    for ig in 1..=n {
        let i = if igdir == 1 { ig } else { n - ig + 1 };
        x[i] = xnew[ig];
        y[i] = ynew[ig];
    }
    scalc(x, y, s, n);
    segspl(x, xp, s, n);
    segspl(y, yp, s, n);
    n
}

/// LERSCL: scale LE radius by RFAC with blending distance DOC.
#[allow(clippy::too_many_arguments)]
pub fn lerscl(x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64], n: usize, doc: f64, rfac: f64, xnew: &mut [f64], ynew: &mut [f64]) {
    let sle = lefind(x, xp, y, yp, s, n);
    let xle = seval(sle, x, xp, s, n);
    let yle = seval(sle, y, yp, s, n);
    let xte = 0.5 * (x[1] + x[n]);
    let yte = 0.5 * (y[1] + y[n]);
    let chord = ((xte - xle).powi(2) + (yte - yle).powi(2)).sqrt();
    let dxc = (xte - xle) / chord;
    let dyc = (yte - yle) / chord;
    let srfac = rfac.abs().sqrt();
    for i in 1..=n {
        let xbar = (x[i] - xle) * dxc + (y[i] - yle) * dyc;
        let ybar = (y[i] - yle) * dxc - (x[i] - xle) * dyc;
        let sopp = sopps(s[i], x, xp, y, yp, s, n, sle);
        let xopp = seval(sopp, x, xp, s, n);
        let yopp = seval(sopp, y, yp, s, n);
        let ybarop = (yopp - yle) * dxc - (xopp - xle) * dyc;
        let xoc = xbar / chord;
        let arg = (xoc / doc).min(15.0);
        let tfac = 1.0 - (1.0 - srfac) * (-arg).exp();
        let ybarct = 0.5 * (ybar + ybarop) + tfac * 0.5 * (ybar - ybarop);
        xnew[i] = xle + xbar * dxc - ybarct * dyc;
        ynew[i] = yle + ybarct * dxc + xbar * dyc;
    }
}

/// SSS: arc-length points S1,S2 at flap surface break locations.
#[allow(clippy::too_many_arguments)]
pub fn sss(ss: f64, del: f64, xbf: f64, ybf: f64, x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64], n: usize, iside: i32) -> (f64, f64) {
    let eps = 1.0e-5;
    let stot = (s[n] - s[1]).abs();
    let sind = (0.5 * del.abs()).sin();
    let ssgn = if iside == 1 { -1.0 } else { 1.0 };

    let rsq = (seval(ss, x, xp, s, n) - xbf).powi(2) + (seval(ss, y, yp, s, n) - ybf).powi(2);
    let mut s1 = ss - (sind * rsq.sqrt() + eps * stot) * ssgn;
    let mut s2 = ss + (sind * rsq.sqrt() + eps * stot) * ssgn;

    let mut conv = false;
    for _ in 0..10 {
        let x1 = seval(s1, x, xp, s, n);
        let x1p = deval(s1, x, xp, s, n);
        let y1 = seval(s1, y, yp, s, n);
        let y1p = deval(s1, y, yp, s, n);
        let x2 = seval(s2, x, xp, s, n);
        let x2p = deval(s2, x, xp, s, n);
        let y2 = seval(s2, y, yp, s, n);
        let y2p = deval(s2, y, yp, s, n);

        let r1sq = (x1 - xbf).powi(2) + (y1 - ybf).powi(2);
        let r2sq = (x2 - xbf).powi(2) + (y2 - ybf).powi(2);
        let r1 = r1sq.sqrt();
        let r2 = r2sq.sqrt();
        let rrsq = (x1 - x2).powi(2) + (y1 - y2).powi(2);
        let rr = rrsq.sqrt();

        if r1 <= eps * stot || r2 <= eps * stot {
            return (ss, ss);
        }
        let r1_s1 = (x1p * (x1 - xbf) + y1p * (y1 - ybf)) / r1;
        let r2_s2 = (x2p * (x2 - xbf) + y2p * (y2 - ybf)) / r2;

        let (rs1, a11, a12, rs2, a21, a22);
        if sind > 0.01 {
            if rr == 0.0 {
                return (s1, s2);
            }
            let rr_s1 = (x1p * (x1 - x2) + y1p * (y1 - y2)) / rr;
            let rr_s2 = -(x2p * (x1 - x2) + y2p * (y1 - y2)) / rr;
            rs1 = ((xbf - x1) * (x2 - x1) + (ybf - y1) * (y2 - y1)) / rr - sind * r1;
            a11 = ((xbf - x1) * (-x1p) + (ybf - y1) * (-y1p)) / rr + ((-x1p) * (x2 - x1) + (-y1p) * (y2 - y1)) / rr
                - ((xbf - x1) * (x2 - x1) + (ybf - y1) * (y2 - y1)) * rr_s1 / rrsq
                - sind * r1_s1;
            a12 = ((xbf - x1) * (x2p) + (ybf - y1) * (y2p)) / rr - ((xbf - x1) * (x2 - x1) + (ybf - y1) * (y2 - y1)) * rr_s2 / rrsq;
            rs2 = r1 - r2;
            a21 = r1_s1;
            a22 = -r2_s2;
        } else {
            rs1 = (r1 + r2) * sind + (s1 - s2) * ssgn;
            a11 = r1_s1 * sind + ssgn;
            a12 = r2_s2 * sind - ssgn;
            let x1pp = d2val(s1, x, xp, s, n);
            let y1pp = d2val(s1, y, yp, s, n);
            let x2pp = d2val(s2, x, xp, s, n);
            let y2pp = d2val(s2, y, yp, s, n);
            let xtot = x1 + x2 - 2.0 * xbf;
            let ytot = y1 + y2 - 2.0 * ybf;
            rs2 = xtot * (x1p + x2p) + ytot * (y1p + y2p);
            a21 = x1p * (x1p + x2p) + y1p * (y1p + y2p) + xtot * x1pp + ytot * y1pp;
            a22 = x2p * (x1p + x2p) + y2p * (y1p + y2p) + xtot * x2pp + ytot * y2pp;
        }
        let det = a11 * a22 - a12 * a21;
        let mut ds1 = -(rs1 * a22 - a12 * rs2) / det;
        let mut ds2 = -(a11 * rs2 - rs1 * a21) / det;
        ds1 = ds1.min(0.01 * stot);
        ds1 = ds1.max(-0.01 * stot);
        ds2 = ds2.min(0.01 * stot);
        ds2 = ds2.max(-0.01 * stot);
        s1 += ds1;
        s2 += ds2;
        if ds1.abs() + ds2.abs() < eps * stot {
            conv = true;
            break;
        }
    }
    if !conv {
        crate::outln!(" SSS: failed to converge subtending angle points");
        s1 = ss;
        s2 = ss;
    }
    if del == 0.0 {
        s1 = 0.5 * (s1 + s2);
        s2 = s1;
    }
    (s1, s2)
}

/// CLIS: list curvatures at nodes.
pub fn clis(x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64], n: usize) {
    let pi = 4.0 * 1.0f64.atan();
    let mut cmax = 0.0f64;
    let mut imax = 1;
    crate::outln!("\n  i        x         y         s       theta        curv");
    let mut arad = 0.0;
    for i in 1..=n {
        if i == 1 {
            arad = (-yp[i]).atan2(-xp[i]);
        } else {
            arad = atanc(-yp[i], -xp[i], arad);
        }
        let adeg = arad * 180.0 / pi;
        let cv = curv(s[i], x, xp, y, yp, s, n);
        crate::outln!(
            " {}{:10.5}{:10.5}{:10.5}{}{}",
            i_fmt(i as i64, 3),
            x[i],
            y[i],
            s[i],
            crate::fmt::f_fmt(adeg, 11, 3),
            crate::fmt::f_fmt(cv, 12, 3)
        );
        if cv.abs() > cmax.abs() {
            cmax = cv;
            imax = i;
        }
    }
    crate::outln!(
        "\n Maximum curvature ={}   at  i,x,y,s  = {}{}{}{}",
        crate::fmt::f_fmt(cmax, 14, 3),
        i_fmt(imax as i64, 3),
        crate::fmt::f_fmt(x[imax], 9, 4),
        crate::fmt::f_fmt(y[imax], 9, 4),
        crate::fmt::f_fmt(s[imax], 9, 4)
    );
}

/// CANG with printing option. IPRINT=2 all angles, 1 max only, 0 none.
pub fn cang_p(x: &[f64], y: &[f64], n: usize, iprint: i32) -> (usize, f64) {
    use crate::fmt::f_fmt;
    let mut amax: f64 = 0.0;
    let mut imax = 1;
    if iprint == 2 {
        crate::outln!("\n  i       x        y      angle");
    }
    for i in 2..n {
        let (mut dx1, mut dy1) = (x[i] - x[i - 1], y[i] - y[i - 1]);
        let (mut dx2, mut dy2) = (x[i] - x[i + 1], y[i] - y[i + 1]);
        if dx1 == 0.0 && dy1 == 0.0 {
            let im2 = i.saturating_sub(2);
            dx1 = x[i] - x[im2];
            dy1 = y[i] - y[im2];
        }
        if dx2 == 0.0 && dy2 == 0.0 {
            let ip2 = (i + 2).min(x.len() - 1);
            dx2 = x[i] - x[ip2];
            dy2 = y[i] - y[ip2];
        }
        let crossp = (dx2 * dy1 - dy2 * dx1) / ((dx1 * dx1 + dy1 * dy1) * (dx2 * dx2 + dy2 * dy2)).sqrt();
        let angl = crossp.asin() * (180.0 / 3.1415926);
        if iprint == 2 {
            crate::outln!(" {}{}{}{}", i_fmt(i as i64, 3), f_fmt(x[i], 9, 4), f_fmt(y[i], 9, 4), f_fmt(angl, 9, 3));
        }
        if angl.abs() > amax.abs() {
            amax = angl;
            imax = i;
        }
    }
    if iprint >= 1 {
        crate::outln!(
            "\n Maximum panel corner angle ={}   at  i,x,y  = {}{}{}",
            f_fmt(amax, 7, 3),
            i_fmt(imax as i64, 3),
            f_fmt(x[imax], 9, 4),
            f_fmt(y[imax], 9, 4)
        );
    }
    (imax, amax)
}

/// IJSECT results.
#[derive(Default, Debug, Clone, Copy)]
pub struct IjSect {
    pub area: f64,
    pub slen: f64,
    pub xmin: f64,
    pub xmax: f64,
    pub ymin: f64,
    pub ymax: f64,
    pub xc: f64,
    pub yc: f64,
    pub xct: f64,
    pub yct: f64,
    pub aixx: f64,
    pub aiyy: f64,
    pub aixxt: f64,
    pub aiyyt: f64,
    pub aj: f64,
    pub ajt: f64,
}

/// IJSECT: section bending/torsion properties.
pub fn ijsect(n: usize, x: &[f64], y: &[f64]) -> IjSect {
    let mut r = IjSect { xmin: x[1], xmax: x[1], ymin: y[1], ymax: y[1], ..Default::default() };

    let mut dx = x[1] - x[n];
    let mut dy = y[1] - y[n];
    let mut ds = (dx * dx + dy * dy).sqrt();
    let mut xavg = 0.5 * (x[1] + x[n]);
    let mut yavg = 0.5 * (y[1] + y[n]);

    let mut x_dy = dy * xavg;
    let mut xx_dy = dy * xavg.powi(2);
    let mut xxx_dy = dy * xavg.powi(3);
    let mut x_ds = ds * xavg;
    let mut xx_ds = ds * xavg.powi(2);

    let mut y_dx = dx * yavg;
    let mut yy_dx = dx * yavg.powi(2);
    let mut yyy_dx = dx * yavg.powi(3);
    let mut y_ds = ds * yavg;
    let mut yy_ds = ds * yavg.powi(2);

    let mut c_ds = ds;

    for i in 2..=n {
        dx = x[i] - x[i - 1];
        dy = y[i] - y[i - 1];
        ds = (dx * dx + dy * dy).sqrt();
        xavg = 0.5 * (x[i] + x[i - 1]);
        yavg = 0.5 * (y[i] + y[i - 1]);

        x_dy += dy * xavg;
        xx_dy += dy * xavg.powi(2);
        xxx_dy += dy * xavg.powi(3);
        x_ds += ds * xavg;
        xx_ds += ds * xavg.powi(2);

        y_dx += dx * yavg;
        yy_dx += dx * yavg.powi(2);
        yyy_dx += dx * yavg.powi(3);
        y_ds += ds * yavg;
        yy_ds += ds * yavg.powi(2);

        c_ds += ds;

        r.xmin = r.xmin.min(x[i]);
        r.xmax = r.xmax.max(x[i]);
        r.ymin = r.ymin.min(y[i]);
        r.ymax = r.ymax.max(y[i]);
    }

    r.area = -y_dx;
    r.slen = c_ds;
    if r.area == 0.0 {
        return r;
    }

    r.xc = xx_dy / (2.0 * x_dy);
    r.xct = x_ds / c_ds;
    r.aiyy = xxx_dy / 3.0 - xx_dy * r.xc + x_dy * r.xc.powi(2);
    r.aiyyt = xx_ds - x_ds * r.xct * 2.0 + c_ds * r.xct.powi(2);

    r.yc = yy_dx / (2.0 * y_dx);
    r.yct = y_ds / c_ds;
    r.aixx = -yyy_dx / 3.0 + yy_dx * r.yc - y_dx * r.yc.powi(2);
    r.aixxt = yy_ds - y_ds * r.yct * 2.0 + c_ds * r.yct.powi(2);

    // find IMID (as in the Fortran, IMID = N if no x increase is found)
    let mut i = 1;
    while i < n {
        if x[i + 1] >= x[i] {
            break;
        }
        i += 1;
    }
    let imid = i;

    let mut aj = 0.0;
    for i in 2..=imid {
        let xavg = 0.5 * (x[i] + x[i - 1]);
        let yavg = 0.5 * (y[i] + y[i - 1]);
        let dx = x[i - 1] - x[i];
        let mut yopp = 0.0;
        if xavg > x[n] {
            yopp = y[n];
        } else if xavg <= x[imid] {
            yopp = y[imid];
        } else {
            let mut j = n;
            while j >= imid {
                if xavg > x[j - 1] && xavg <= x[j] {
                    let frac = (xavg - x[j - 1]) / (x[j] - x[j - 1]);
                    yopp = y[j - 1] + (y[j] - y[j - 1]) * frac;
                    break;
                }
                j -= 1;
            }
        }
        aj += (yavg - yopp).abs().powi(3) * dx / 3.0;
    }
    r.aj = aj;
    r.ajt = 4.0 * r.area.powi(2) / r.slen;
    r
}

/// BENDUMP: display structural properties.
pub fn bendump(n: usize, x: &[f64], y: &[f64]) {
    use crate::fmt::g_fmt;
    let r = ijsect(n, x, y);
    let g = |v: f64| g_fmt(v, 14, 6);
    let p = |a: &str, v: f64| crate::outln!(" {}{}", a, g(v));
    crate::outln!("");
    p("Area =", r.area);
    p("Slen =", r.slen);
    crate::outln!("");
    crate::outln!(" X-bending parameters(solid):");
    p("        Xc =", r.xc);
    p("  max X-Xc =", r.xmax - r.xc);
    p("  min X-Xc =", r.xmin - r.xc);
    p("       Iyy =", r.aiyy);
    let xbar = (r.xmax - r.xc).abs().max((r.xmin - r.xc).abs());
    p(" Iyy/(X-Xc)=", r.aiyy / xbar);
    crate::outln!("");
    crate::outln!(" Y-bending parameters(solid):");
    p("        Yc =", r.yc);
    p("  max Y-Yc =", r.ymax - r.yc);
    p("  min Y-Yc =", r.ymin - r.yc);
    p("       Ixx =", r.aixx);
    let ybar = (r.ymax - r.yc).abs().max((r.ymin - r.yc).abs());
    p(" Ixx/(Y-Yc)=", r.aixx / ybar);
    crate::outln!("");
    p("       J   =", r.aj);
    crate::outln!("");
    crate::outln!("");
    crate::outln!(" X-bending parameters(skin):");
    p("         Xc =", r.xct);
    p("   max X-Xc =", r.xmax - r.xct);
    p("   min X-Xc =", r.xmin - r.xct);
    p("      Iyy/t =", r.aiyyt);
    let xbart = (r.xmax - r.xct).abs().max((r.xmin - r.xct).abs());
    p(" Iyy/t(X-Xc)=", r.aiyyt / xbart);
    crate::outln!("");
    crate::outln!(" Y-bending parameters(skin):");
    p("         Yc =", r.yct);
    p("   max Y-Yc =", r.ymax - r.yct);
    p("   min Y-Yc =", r.ymin - r.yct);
    p("      Ixx/t =", r.aixxt);
    let ybart = (r.ymax - r.yct).abs().max((r.ymin - r.yct).abs());
    p(" Ixx/t(Y-Yc)=", r.aixxt / ybart);
    crate::outln!("");
    p("      J/t   =", r.ajt);
}

/// AREFINE: add points around corners exceeding ATOL degrees, within X1..X2.
/// Returns NNEW (0 on overflow).
#[allow(clippy::too_many_arguments)]
pub fn arefine(
    x: &[f64],
    y: &[f64],
    s: &[f64],
    xs: &[f64],
    ys: &[f64],
    n: usize,
    atol: f64,
    ndim: usize,
    xnew: &mut [f64],
    ynew: &mut [f64],
    x1: f64,
    x2: f64,
) -> usize {
    let atolr = atol * 3.14159 / 180.0;
    let overflow = || {
        crate::outln!(" SDOUBLE:  Arrays will overflow.  No action taken.");
        0
    };
    let mut k = 1;
    xnew[k] = x[1];
    ynew[k] = y[1];
    for i in 2..n {
        let dxm = x[i] - x[i - 1];
        let dym = y[i] - y[i - 1];
        let dxp = x[i + 1] - x[i];
        let dyp = y[i + 1] - y[i];
        let crsp = dxm * dyp - dym * dxp;
        let dotp = dxm * dxp + dym * dyp;
        let aseg = if crsp == 0.0 && dotp == 0.0 { 0.0 } else { crsp.atan2(dotp) };
        let lref = aseg.abs() > atolr;
        if lref {
            let smid = s[i] - 0.3333 * (s[i] - s[i - 1]);
            let xk = seval(smid, x, xs, s, n);
            let yk = seval(smid, y, ys, s, n);
            if xk >= x1 && xk <= x2 {
                k += 1;
                if k > ndim {
                    return overflow();
                }
                xnew[k] = xk;
                ynew[k] = yk;
            }
        }
        k += 1;
        if k > ndim {
            return overflow();
        }
        xnew[k] = x[i];
        ynew[k] = y[i];
        if lref {
            let smid = s[i] + 0.3333 * (s[i + 1] - s[i]);
            let xk = seval(smid, x, xs, s, n);
            let yk = seval(smid, y, ys, s, n);
            if xk >= x1 && xk <= x2 {
                k += 1;
                if k > ndim {
                    return overflow();
                }
                xnew[k] = xk;
                ynew[k] = yk;
            }
        }
    }
    k += 1;
    if k > ndim {
        return overflow();
    }
    xnew[k] = x[n];
    ynew[k] = y[n];
    k
}

/// SCHECK: remove tiny segments. Returns LCHANGE.
pub fn scheck(x: &mut [f64], y: &mut [f64], n: &mut usize, stol: f64) -> bool {
    let mut lchange = false;
    if stol > 0.3 {
        crate::outln!(" SCHECK:  Bad value for small panels (STOL > 0.3)");
        return false;
    }
    'restart: loop {
        let nn = *n;
        let mut i = 2;
        while i + 2 <= nn {
            let dxm1 = x[i] - x[i - 1];
            let dym1 = y[i] - y[i - 1];
            let dsm1 = (dxm1 * dxm1 + dym1 * dym1).sqrt();
            let dxp1 = x[i + 1] - x[i];
            let dyp1 = y[i + 1] - y[i];
            let dsp1 = (dxp1 * dxp1 + dyp1 * dyp1).sqrt();
            let dxp2 = x[i + 2] - x[i + 1];
            let dyp2 = y[i + 2] - y[i + 1];
            let dsp2 = (dxp2 * dxp2 + dyp2 * dyp2).sqrt();
            if dsp1 != 0.0 && (dsp1 < stol * dsm1 || dsp1 < stol * dsp2) {
                x[i] = 0.5 * (x[i] + x[i + 1]);
                y[i] = 0.5 * (y[i] + y[i + 1]);
                for l in i + 1..=nn {
                    x[l] = x.get(l + 1).copied().unwrap_or(0.0);
                    y[l] = y.get(l + 1).copied().unwrap_or(0.0);
                }
                *n = nn - 1;
                lchange = true;
                crate::outln!(" SCHECK segment removed at {}", lint(i as i64));
                continue 'restart;
            }
            i += 1;
        }
        break;
    }
    lchange
}

/// INSIDE: true if point (xf,yf) is inside contour x,y.
pub fn inside(x: &[f64], y: &[f64], n: usize, xf: f64, yf: f64) -> bool {
    let mut angle = 0.0;
    for i in 1..=n {
        let ip = if i == n { 1 } else { i + 1 };
        let xb1 = x[i] - xf;
        let yb1 = y[i] - yf;
        let xb2 = x[ip] - xf;
        let yb2 = y[ip] - yf;
        angle += (xb1 * yb2 - yb1 * xb2) / ((xb1 * xb1 + yb1 * yb1) * (xb2 * xb2 + yb2 * yb2)).sqrt();
    }
    angle.abs() > 1.0
}

// ---------------------------------------------------------------------------
// sort.f
// ---------------------------------------------------------------------------

/// SORTOL: sort by S, removing near-duplicate (S,W) points within TOL.
pub fn sortol(tol: f64, kk: &mut usize, s: &mut [f64], w: &mut [f64]) {
    let mut sorted = false;
    for _ in 0..1234 {
        let mut done = true;
        for n in 1..*kk {
            let np = n + 1;
            if s[np] < s[n] {
                s.swap(np, n);
                w.swap(np, n);
                done = false;
            }
        }
        if done {
            sorted = true;
            break;
        }
    }
    if !sorted {
        crate::outln!(" Sort failed");
    }
    loop {
        let kks = *kk;
        let mut done = true;
        for k in 1..=kks {
            if k >= *kk {
                continue;
            }
            let dsq = (s[k] - s[k + 1]).powi(2) + (w[k] - w[k + 1]).powi(2);
            if dsq >= tol * tol {
                continue;
            }
            *kk -= 1;
            for kt in k + 1..=*kk {
                s[kt] = s[kt + 1];
                w[kt] = w[kt + 1];
            }
            done = false;
        }
        if done {
            break;
        }
    }
}

/// SORTDUP: sort arrays by S keeping duplicates.
pub fn sortdup(kk: usize, s: &mut [f64], w: &mut [f64]) {
    for _ in 0..1234 {
        let mut done = true;
        for n in 1..kk {
            let np = n + 1;
            if s[np] >= s[n] {
                continue;
            }
            s.swap(np, n);
            w.swap(np, n);
            done = false;
        }
        if done {
            return;
        }
    }
    crate::outln!(" Sort failed");
}

/// FIXDUP: remove leading/trailing duplicates and triple duplicates.
pub fn fixdup(kk: &mut usize, s: &mut [f64], w: &mut [f64]) {
    if s[2] == s[1] {
        for n in 1..*kk {
            s[n] = s[n + 1];
            w[n] = w[n + 1];
        }
        *kk -= 1;
    }
    if s[*kk] == s[*kk - 1] {
        s[*kk - 1] = s[*kk];
        w[*kk - 1] = w[*kk];
        *kk -= 1;
    }
    'again: loop {
        let k = *kk;
        for n in 1..=k.saturating_sub(2) {
            if s[n] == s[n + 1] && s[n] == s[n + 2] {
                for i in n..k {
                    s[i] = s[i + 1];
                    w[i] = w[i + 1];
                }
                *kk -= 1;
                continue 'again;
            }
        }
        break;
    }
}

// ---------------------------------------------------------------------------
// xtcam.f
// ---------------------------------------------------------------------------

/// GETCAM: camber and thickness distributions of an airfoil.
#[allow(clippy::too_many_arguments)]
pub fn getcam(
    xcm: &mut [f64],
    ycm: &mut [f64],
    ncm: &mut usize,
    xtk: &mut [f64],
    ytk: &mut [f64],
    ntk: &mut usize,
    x: &[f64],
    xp: &[f64],
    y: &[f64],
    yp: &[f64],
    s: &[f64],
    n: usize,
) {
    let sl = xlfind(x, xp, y, yp, s, n);
    let xl = seval(sl, x, xp, s, n);
    let yl = seval(sl, y, yp, s, n);
    for i in 1..=n {
        let sopp = sopps(s[i], x, xp, y, yp, s, n, sl);
        let xopp = seval(sopp, x, xp, s, n);
        let yopp = seval(sopp, y, yp, s, n);
        xcm[i] = 0.5 * (x[i] + xopp);
        ycm[i] = 0.5 * (y[i] + yopp);
        xtk[i] = 0.5 * (x[i] + xopp);
        ytk[i] = 0.5 * (y[i] - yopp);
        ytk[i] = ytk[i].abs();
    }
    let tol = 1.0e-5 * (s[n] - s[1]);
    *ncm = n + 1;
    xcm[n + 1] = xl;
    ycm[n + 1] = yl;
    sortol(tol, ncm, xcm, ycm);
    let yof = ycm[1];
    for i in 1..=*ncm {
        ycm[i] -= yof;
    }
    *ntk = n + 1;
    xtk[n + 1] = xl;
    ytk[n + 1] = 0.0;
    sortol(tol, ntk, xtk, ytk);
}

/// GETMAX: highpoint of y(x) and its location. Returns (xmax, ymax).
pub fn getmax(x: &[f64], y: &[f64], yp: &mut [f64], n: usize) -> (f64, f64) {
    let xlen = x[n] - x[1];
    let xtol = xlen * 1.0e-5;
    segspl(y, yp, x, n);
    let mut ymax0 = y[1];
    let mut xmax0 = x[1];
    let mut ddx = 0.0;
    for i in 2..=n {
        if y[i].abs() > ymax0.abs() {
            ymax0 = y[i];
            xmax0 = 0.5 * (x[i - 1] + x[i]);
            // Fortran reads X(N+1) here when I=N (left-over array content)
            let xip = x.get(i + 1).copied().unwrap_or(0.0);
            ddx = 0.5 * (xip - x[i - 1]).abs();
        }
    }
    let mut xmax = xmax0;
    let mut ymax;
    for _ in 0..10 {
        ymax = seval(xmax, y, yp, x, n);
        let res = deval(xmax, y, yp, x, n);
        let resp = d2val(xmax, y, yp, x, n);
        if (xlen * resp).abs() < 1.0e-6 {
            return (xmax, ymax);
        }
        let mut dx = -res / resp;
        dx = (0.5 * ddx).min(dx.abs()).copysign(dx);
        xmax += dx;
        if dx.abs() < xtol {
            return (xmax, ymax);
        }
    }
    crate::outln!(" GETMAX: Newton iteration for max camber/thickness failed.");
    (xmax0, ymax0)
}

/// CPCAM: camber line y(x) from loading DCp(x).
pub fn cpcam(n: usize, x: &[f64], y: &mut [f64], dydx: &mut [f64], p: &[f64], dpdx: &[f64]) {
    let qopi = 7.9577471545948E-02;
    let p0 = p[1];
    let p1 = p[n];
    let x0 = x[1];
    let x1 = x[n];
    let ysing = |xt: f64| {
        qopi * p1 * ((xt - x1) * ((x1 - xt) / (x1 - x0)).max(1.0e-6).ln() - xt) - qopi * p0 * ((xt - x0) * ((xt - x0) / (x1 - x0)).max(1.0e-6).ln() - xt)
    };
    for i in 1..=n {
        crate::outln!(" {}{:10.4}{:10.4}{:10.4}", i_fmt(i as i64, 4), x[i], p[i], dpdx[i]);
        dydx[i] = 0.0;
        let j = 1;
        let mut yp1 = if i == j { dpdx[j] } else { (p[j] - p[i]) / (x[j] - x[i]) };
        for j in 2..=n {
            let yp2 = if i == j { dpdx[j] } else { (p[j] - p[i]) / (x[j] - x[i]) };
            dydx[i] += 0.5 * (yp1 + yp2) * (x[j] - x[j - 1]);
            yp1 = yp2;
        }
        dydx[i] *= qopi;
        if i != 1 {
            dydx[i] -= qopi * (p[i] - p0) * (x[i] - x0).ln();
        }
        if i != n {
            dydx[i] += qopi * (p[i] - p1) * (x1 - x[i]).ln();
        }
    }
    y[1] = 0.0;
    for i in 2..=n {
        y[i] = y[i - 1] + 0.5 * (dydx[i] + dydx[i - 1]) * (x[i] - x[i - 1]);
    }
    for i in 1..=n {
        y[i] += ysing(x[i]);
    }
    let y0 = y[1];
    let y1 = y[n];
    for i in 1..=n {
        y[i] = y[i] - y0 * (x1 - x[i]) / (x1 - x0) - y1 * (x[i] - x0) / (x1 - x0);
    }
}

// ---------------------------------------------------------------------------
// spline.f
// ---------------------------------------------------------------------------

/// SPLNXY: spline x(s),y(s) with true arc length parameter.
pub fn splnxy(x: &[f64], xs: &mut [f64], y: &[f64], ys: &mut [f64], s: &mut [f64], n: usize) {
    const KMAX: usize = 32;
    let kk = KMAX;
    let npass = 10;
    let mut xt = [0.0f64; KMAX + 1];
    let mut yt = [0.0f64; KMAX + 1];
    scalc(x, y, s, n);
    segspl(x, xs, s, n);
    segspl(y, ys, s, n);
    for ipass in 1..=npass {
        let mut serr: f64 = 0.0;
        let mut ds = s[2] - s[1];
        for i in 2..=n {
            let dx = x[i] - x[i - 1];
            let dy = y[i] - y[i - 1];
            let cx1 = ds * xs[i - 1] - dx;
            let cx2 = ds * xs[i] - dx;
            let cy1 = ds * ys[i - 1] - dy;
            let cy2 = ds * ys[i] - dy;
            xt[0] = 0.0;
            yt[0] = 0.0;
            for k in 1..kk {
                let t = k as f64 / kk as f64;
                xt[k] = t * dx + (t - t * t) * ((1.0 - t) * cx1 - t * cx2);
                yt[k] = t * dy + (t - t * t) * ((1.0 - t) * cy1 - t * cy2);
            }
            xt[kk] = dx;
            yt[kk] = dy;
            let mut sint1 = 0.0;
            for k in 1..=kk {
                sint1 += ((xt[k] - xt[k - 1]).powi(2) + (yt[k] - yt[k - 1]).powi(2)).sqrt();
            }
            let mut sint2 = 0.0;
            let mut k = 2;
            while k <= kk {
                sint2 += ((xt[k] - xt[k - 2]).powi(2) + (yt[k] - yt[k - 2]).powi(2)).sqrt();
                k += 2;
            }
            let sint = (4.0 * sint1 - sint2) / 3.0;
            if (sint - ds).abs() > serr.abs() {
                serr = sint - ds;
            }
            if i < n {
                ds = s[i + 1] - s[i];
            }
            s[i] = s[i - 1] + sint.sqrt();
        }
        serr /= s[n] - s[1];
        crate::outln!("{}{}", lint(ipass as i64), lreal(serr));
        segspl(x, xs, s, n);
        segspl(y, ys, s, n);
        if serr.abs() < 1.0e-7 {
            return;
        }
    }
}
