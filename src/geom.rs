//! Geometry routines: NACA generation, paneling, LE finding, geometric
//! parameters, airfoil file reading (ports of naca.f, xgeom.f, PANGEN,
//! ABCOPY, LOAD, AREAD).

use crate::panel::ncalc;
use crate::spline::{curv, d2val, deval, scalc, segspl, seval, trisol};
use crate::xfoil::{XFoil, IBX, IQX};

/// ATAN2 with branch-cut continuity relative to THOLD (ATANC).
pub fn atanc(y: f64, x: f64, thold: f64) -> f64 {
    let pi = std::f64::consts::PI;
    let tpi = 2.0 * pi;
    let thnew = y.atan2(x);
    let dthet = thnew - thold;
    let s = if dthet >= 0.0 { pi } else { -pi };
    let dtcorr = dthet - tpi * ((dthet + s) / tpi).trunc();
    thold + dtcorr
}

/// NACA 4-digit. Returns (xb, yb, nb, name) with 1-based arrays.
pub fn naca4(ides: i32, nside: usize) -> (Vec<f64>, Vec<f64>, usize, String) {
    let an = 1.5;
    let n4 = ides / 1000;
    let n3 = (ides - n4 * 1000) / 100;
    let n2 = (ides - n4 * 1000 - n3 * 100) / 10;
    let n1 = ides - n4 * 1000 - n3 * 100 - n2 * 10;
    let m = n4 as f64 / 100.0;
    let p = n3 as f64 / 10.0;
    let t = (n2 * 10 + n1) as f64 / 100.0;
    let anp = an + 1.0;
    let mut xx = vec![0.0; nside + 1];
    let mut yt = vec![0.0; nside + 1];
    let mut yc = vec![0.0; nside + 1];
    for i in 1..=nside {
        let frac = (i - 1) as f64 / (nside - 1) as f64;
        xx[i] = if i == nside { 1.0 } else { 1.0 - anp * frac * (1.0 - frac).powf(an) - (1.0 - frac).powf(anp) };
        let x = xx[i];
        yt[i] = (0.29690 * x.sqrt() - 0.12600 * x - 0.35160 * x * x + 0.28430 * x * x * x - 0.10150 * x.powi(4)) * t / 0.20;
        yc[i] = if x < p {
            m / (p * p) * (2.0 * p * x - x * x)
        } else {
            m / ((1.0 - p) * (1.0 - p)) * ((1.0 - 2.0 * p) + 2.0 * p * x - x * x)
        };
    }
    let (xb, yb, nb) = assemble(&xx, &yt, &yc, nside);
    (xb, yb, nb, format!("NACA {}{}{}{}", n4, n3, n2, n1))
}

fn assemble(xx: &[f64], yt: &[f64], yc: &[f64], nside: usize) -> (Vec<f64>, Vec<f64>, usize) {
    let mut xb = vec![0.0; 2 * nside + 1];
    let mut yb = vec![0.0; 2 * nside + 1];
    let mut ib = 0;
    for i in (1..=nside).rev() {
        ib += 1;
        xb[ib] = xx[i];
        yb[ib] = yc[i] + yt[i];
    }
    for i in 2..=nside {
        ib += 1;
        xb[ib] = xx[i];
        yb[ib] = yc[i] - yt[i];
    }
    (xb, yb, ib)
}

/// NACA 5-digit. Returns None for illegal designation.
pub fn naca5(ides: i32, nside: usize) -> Option<(Vec<f64>, Vec<f64>, usize, String)> {
    let an = 1.5;
    let n5 = ides / 10000;
    let n4 = (ides - n5 * 10000) / 1000;
    let n3 = (ides - n5 * 10000 - n4 * 1000) / 100;
    let n2 = (ides - n5 * 10000 - n4 * 1000 - n3 * 100) / 10;
    let n1 = ides - n5 * 10000 - n4 * 1000 - n3 * 100 - n2 * 10;
    let n543 = 100 * n5 + 10 * n4 + n3;
    let (m, c) = match n543 {
        210 => (0.0580, 361.4),
        220 => (0.1260, 51.64),
        230 => (0.2025, 15.957),
        240 => (0.2900, 6.643),
        250 => (0.3910, 3.230),
        _ => {
            crate::outln!(" Illegal 5-digit designation");
            crate::outln!(" First three digits must be 210, 220, ... 250");
            return None;
        }
    };
    let t = (n2 * 10 + n1) as f64 / 100.0;
    let anp = an + 1.0;
    let mut xx = vec![0.0; nside + 1];
    let mut yt = vec![0.0; nside + 1];
    let mut yc = vec![0.0; nside + 1];
    for i in 1..=nside {
        let frac = (i - 1) as f64 / (nside - 1) as f64;
        xx[i] = if i == nside { 1.0 } else { 1.0 - anp * frac * (1.0 - frac).powf(an) - (1.0 - frac).powf(anp) };
        let x = xx[i];
        yt[i] = (0.29690 * x.sqrt() - 0.12600 * x - 0.35160 * x * x + 0.28430 * x * x * x - 0.10150 * x.powi(4)) * t / 0.20;
        yc[i] = if x < m {
            (c / 6.0) * (x * x * x - 3.0 * m * x * x + m * m * (3.0 - m) * x)
        } else {
            (c / 6.0) * m * m * m * (1.0 - x)
        };
    }
    let (xb, yb, nb) = assemble(&xx, &yt, &yc, nside);
    Some((xb, yb, nb, format!("NACA {}{}{}{}{}", n5, n4, n3, n2, n1)))
}

/// Leading edge spline parameter (LEFIND).
pub fn lefind(x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64], n: usize) -> f64 {
    let dseps = (s[n] - s[1]) * 1.0e-5;
    let xte = 0.5 * (x[1] + x[n]);
    let yte = 0.5 * (y[1] + y[n]);
    let mut i = 3;
    while i <= n - 2 {
        let dxte = x[i] - xte;
        let dyte = y[i] - yte;
        let dx = x[i + 1] - x[i];
        let dy = y[i + 1] - y[i];
        if dxte * dx + dyte * dy < 0.0 {
            break;
        }
        i += 1;
    }
    let mut sle = s[i];
    if s[i] == s[i - 1] {
        return sle;
    }
    for _ in 0..50 {
        let xle = seval(sle, x, xp, s, n);
        let yle = seval(sle, y, yp, s, n);
        let dxds = deval(sle, x, xp, s, n);
        let dyds = deval(sle, y, yp, s, n);
        let dxdd = d2val(sle, x, xp, s, n);
        let dydd = d2val(sle, y, yp, s, n);
        let xchord = xle - xte;
        let ychord = yle - yte;
        let res = xchord * dxds + ychord * dyds;
        let ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd;
        let mut dsle = -res / ress;
        dsle = dsle.max(-0.02 * (xchord + ychord).abs());
        dsle = dsle.min(0.02 * (xchord + ychord).abs());
        sle += dsle;
        if dsle.abs() < dseps {
            return sle;
        }
    }
    crate::outln!(" LEFIND:  LE point not found.  Continuing...");
    s[i]
}

/// Opposite point on airfoil with same chordwise x (SOPPS).
#[allow(clippy::too_many_arguments)]
pub fn sopps(si: f64, x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64], n: usize, sle: f64) -> f64 {
    let slen = s[n] - s[1];
    let xle = seval(sle, x, xp, s, n);
    let yle = seval(sle, y, yp, s, n);
    let xte = 0.5 * (x[1] + x[n]);
    let yte = 0.5 * (y[1] + y[n]);
    let chord = ((xte - xle).powi(2) + (yte - yle).powi(2)).sqrt();
    let dxc = (xte - xle) / chord;
    let dyc = (yte - yle) / chord;
    let (inn, inopp) = if si < sle { (1, n) } else { (n, 1) };
    let sfrac = (si - sle) / (s[inn] - sle);
    let mut sopp = sle + sfrac * (s[inopp] - sle);
    if sfrac.abs() <= 1.0e-5 {
        return sle;
    }
    let xi = seval(si, x, xp, s, n);
    let yi = seval(si, y, yp, s, n);
    let xbar = (xi - xle) * dxc + (yi - yle) * dyc;
    for _ in 0..12 {
        let xopp = seval(sopp, x, xp, s, n);
        let yopp = seval(sopp, y, yp, s, n);
        let xoppd = deval(sopp, x, xp, s, n);
        let yoppd = deval(sopp, y, yp, s, n);
        let res = (xopp - xle) * dxc + (yopp - yle) * dyc - xbar;
        let resd = xoppd * dxc + yoppd * dyc;
        if res.abs() / slen < 1.0e-5 {
            return sopp;
        }
        if resd == 0.0 {
            break;
        }
        let dsopp = -res / resd;
        sopp += dsopp;
        if dsopp.abs() / slen < 1.0e-5 {
            return sopp;
        }
    }
    crate::outln!(" SOPPS: Opposite-point location failed. Continuing...");
    sle + sfrac * (s[inopp] - sle)
}

/// Scale coordinates to unit chord (NORM).
pub fn norm(x: &mut [f64], xp: &mut [f64], y: &mut [f64], yp: &mut [f64], s: &mut [f64], n: usize) {
    scalc(x, y, s, n);
    segspl(x, xp, s, n);
    segspl(y, yp, s, n);
    let sle = lefind(x, xp, y, yp, s, n);
    let xmax = 0.5 * (x[1] + x[n]);
    let xmin = seval(sle, x, xp, s, n);
    let ymin = seval(sle, y, yp, s, n);
    let fudge = 1.0 / (xmax - xmin);
    for i in 1..=n {
        x[i] = (x[i] - xmin) * fudge;
        y[i] = (y[i] - ymin) * fudge;
        s[i] *= fudge;
    }
}

/// Section properties (AECALC). Returns (area, xcen, ycen, ei11, ei22, apx1, apx2).
pub fn aecalc(n: usize, x: &[f64], y: &[f64], t: &[f64], itype: i32) -> (f64, f64, f64, f64, f64, f64, f64) {
    let pi = std::f64::consts::PI;
    let (mut sint, mut aint, mut xint, mut yint, mut xxint, mut xyint, mut yyint) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for io in 1..=n {
        let ip = if io == n { 1 } else { io + 1 };
        let dx = x[io] - x[ip];
        let dy = y[io] - y[ip];
        let xa = (x[io] + x[ip]) * 0.5;
        let ya = (y[io] + y[ip]) * 0.5;
        let ta = (t[io] + t[ip]) * 0.5;
        let ds = (dx * dx + dy * dy).sqrt();
        sint += ds;
        if itype == 1 {
            let da = ya * dx;
            aint += da;
            xint += xa * da;
            yint += ya * da / 2.0;
            xxint += xa * xa * da;
            xyint += xa * ya * da / 2.0;
            yyint += ya * ya * da / 3.0;
        } else {
            let da = ta * ds;
            aint += da;
            xint += xa * da;
            yint += ya * da;
            xxint += xa * xa * da;
            xyint += xa * ya * da;
            yyint += ya * ya * da;
        }
    }
    let area = aint;
    if aint == 0.0 {
        return (area, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0f64.atan2(0.0));
    }
    let xcen = xint / aint;
    let ycen = yint / aint;
    let eixx = yyint - ycen * ycen * aint;
    let eixy = xyint - xcen * ycen * aint;
    let eiyy = xxint - xcen * xcen * aint;
    let eisq = 0.25 * (eixx - eiyy).powi(2) + eixy * eixy;
    let sgn = if eiyy - eixx >= 0.0 { 1.0 } else { -1.0 };
    let ei11 = 0.5 * (eixx + eiyy) - sgn * eisq.sqrt();
    let ei22 = 0.5 * (eixx + eiyy) + sgn * eisq.sqrt();
    let (mut apx1, mut apx2);
    if ei11 == 0.0 || ei22 == 0.0 || eisq / (ei11 * ei22) < (0.001 * sint).powi(4) {
        apx1 = 0.0;
        apx2 = 1.0f64.atan2(0.0);
    } else {
        let c1 = eixy;
        let s1 = eixx - ei11;
        let c2 = eixy;
        let s2 = eixx - ei22;
        if s1.abs() > s2.abs() {
            apx1 = s1.atan2(c1);
            apx2 = apx1 + 0.5 * pi;
        } else {
            apx2 = s2.atan2(c2);
            apx1 = apx2 - 0.5 * pi;
        }
        if apx1 < -0.5 * pi {
            apx1 += pi;
        }
        if apx1 > 0.5 * pi {
            apx1 -= pi;
        }
        if apx2 < -0.5 * pi {
            apx2 += pi;
        }
        if apx2 > 0.5 * pi {
            apx2 -= pi;
        }
    }
    (area, xcen, ycen, ei11, ei22, apx1, apx2)
}

/// Max thickness & camber at discrete points (TCCALC). Returns (thick, xthick, cambr, xcambr).
pub fn tccalc(x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64], n: usize) -> (f64, f64, f64, f64) {
    let sle = lefind(x, xp, y, yp, s, n);
    let xle = seval(sle, x, xp, s, n);
    let yle = seval(sle, y, yp, s, n);
    let xte = 0.5 * (x[1] + x[n]);
    let yte = 0.5 * (y[1] + y[n]);
    let chord = ((xte - xle).powi(2) + (yte - yle).powi(2)).sqrt();
    let dxc = (xte - xle) / chord;
    let dyc = (yte - yle) / chord;
    let (mut thick, mut xthick, mut cambr, mut xcambr) = (0.0f64, 0.0, 0.0f64, 0.0);
    for i in 1..=n {
        let ybar = (y[i] - yle) * dxc - (x[i] - xle) * dyc;
        let sopp = sopps(s[i], x, xp, y, yp, s, n, sle);
        let xopp = seval(sopp, x, xp, s, n);
        let yopp = seval(sopp, y, yp, s, n);
        let ybarop = (yopp - yle) * dxc - (xopp - xle) * dyc;
        let yc = 0.5 * (ybar + ybarop);
        let yt = (ybar - ybarop).abs();
        if yc.abs() > cambr.abs() {
            cambr = yc;
            xcambr = xopp;
        }
        if yt.abs() > thick.abs() {
            thick = yt;
            xthick = xopp;
        }
    }
    (thick, xthick, cambr, xcambr)
}

/// Geometric parameters of airfoil shape (GEOPAR).
#[derive(Clone, Copy, Debug, Default)]
pub struct GeoPar {
    pub sle: f64,
    pub chord: f64,
    pub area: f64,
    pub radle: f64,
    pub angte: f64,
    pub ei11a: f64,
    pub ei22a: f64,
    pub apx1a: f64,
    pub apx2a: f64,
    pub ei11t: f64,
    pub ei22t: f64,
    pub apx1t: f64,
    pub apx2t: f64,
    pub thick: f64,
    pub xthick: f64,
    pub cambr: f64,
    pub xcambr: f64,
}

pub fn geopar(x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64], n: usize, verbose: bool) -> GeoPar {
    let sle = lefind(x, xp, y, yp, s, n);
    let xle = seval(sle, x, xp, s, n);
    let yle = seval(sle, y, yp, s, n);
    let xte = 0.5 * (x[1] + x[n]);
    let yte = 0.5 * (y[1] + y[n]);
    let chord = ((xte - xle).powi(2) + (yte - yle).powi(2)).sqrt();
    let curvle = curv(sle, x, xp, y, yp, s, n);
    let radle = if curvle.abs() > 0.001 * (s[n] - s[1]) { 1.0 / curvle } else { 0.0 };
    let ang1 = (-yp[1]).atan2(-xp[1]);
    let ang2 = atanc(yp[n], xp[n], ang1);
    let angte = ang2 - ang1;
    let t = vec![1.0; n + 1];
    let (area, _, _, ei11a, ei22a, apx1a, apx2a) = aecalc(n, x, y, &t, 1);
    let (_, _, _, ei11t, ei22t, apx1t, apx2t) = aecalc(n, x, y, &t, 2);
    let (thick, xthick, cambr, xcambr) = tccalc(x, xp, y, yp, s, n);
    if verbose {
        crate::outln!(" Max thickness = {:12.6}  at x = {:7.3}", thick, xthick);
        crate::outln!(" Max camber    = {:12.6}  at x = {:7.3}", cambr, xcambr);
    }
    GeoPar {
        sle,
        chord,
        area,
        radle,
        angte,
        ei11a,
        ei22a,
        apx1a,
        apx2a,
        ei11t,
        ei22t,
        apx1t,
        apx2t,
        thick,
        xthick,
        cambr,
        xcambr,
    }
}

/// Max panel corner angle (CANG). Returns (imax, amax).
pub fn cang(x: &[f64], y: &[f64], n: usize) -> (usize, f64) {
    let mut amax: f64 = 0.0;
    let mut imax = 1;
    for i in 2..n {
        let (mut dx1, mut dy1) = (x[i] - x[i - 1], y[i] - y[i - 1]);
        let (mut dx2, mut dy2) = (x[i] - x[i + 1], y[i] - y[i + 1]);
        if dx1 == 0.0 && dy1 == 0.0 && i >= 3 {
            dx1 = x[i] - x[i - 2];
            dy1 = y[i] - y[i - 2];
        }
        if dx2 == 0.0 && dy2 == 0.0 && i + 2 <= n {
            dx2 = x[i] - x[i + 2];
            dy2 = y[i] - y[i + 2];
        }
        let crossp = (dx2 * dy1 - dy2 * dx1) / ((dx1 * dx1 + dy1 * dy1) * (dx2 * dx2 + dy2 * dy2)).sqrt();
        let angl = crossp.asin() * (180.0 / 3.1415926);
        if angl.abs() > amax.abs() {
            amax = angl;
            imax = i;
        }
    }
    (imax, amax)
}

/// Parse up to `nmax` floats from a string in Fortran list-directed style
/// (blanks/commas/tabs separate). Returns None on a parse error.
pub fn getflt(line: &str, nmax: usize) -> Option<Vec<f64>> {
    let mut out = Vec::new();
    for tok in line.split(|c: char| c == ' ' || c == ',' || c == '\t').filter(|t| !t.is_empty()) {
        if out.len() >= nmax {
            break;
        }
        let t = tok.replace(['d', 'D'], "e");
        match t.parse::<f64>() {
            Ok(v) => out.push(v),
            Err(_) => {
                // Fortran list-directed read accepts integers like "1." etc. — anything else is an error
                return None;
            }
        }
    }
    Some(out)
}

/// Result of reading an airfoil file (AREAD).
pub struct AreadResult {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub n: usize,
    pub name: String,
    pub ispars: String,
    pub itype: i32,
}

/// Read airfoil coordinate file (AREAD). itype: 0 error, 1 plain, 2 labeled, 3 MSES.
pub fn aread(fname: &str, nmax: usize) -> AreadResult {
    let fail = |msg: &str| {
        crate::outln!("{}", msg);
        crate::outln!(" *** LOAD NOT COMPLETED ***");
        AreadResult { x: vec![0.0], y: vec![0.0], n: 0, name: String::new(), ispars: String::new(), itype: 0 }
    };
    let text = match std::fs::read_to_string(fname.trim()) {
        Ok(t) => t,
        Err(_) => return fail(&format!("\n File OPEN error.  Nonexistent file:  {}", fname.trim())),
    };
    let lines: Vec<&str> = text.lines().collect();
    let is_comment = |l: &str| l.starts_with('#') || l.starts_with('!');
    let mut idx = 0;
    while idx < lines.len() && is_comment(lines[idx]) {
        idx += 1;
    }
    if idx >= lines.len() {
        return fail("\n File READ error.  Unrecognizable file format");
    }
    let line1 = lines[idx];
    let mut idx2 = idx + 1;
    while idx2 < lines.len() && is_comment(lines[idx2]) {
        idx2 += 1;
    }
    if idx2 >= lines.len() {
        return fail("\n File READ error.  Unrecognizable file format");
    }
    let line2 = lines[idx2];
    let mut name = String::new();
    let mut ispars = String::new();
    let itype;
    let start;
    match getflt(line1, 2) {
        Some(a) if a.len() >= 2 => {
            crate::outln!("\n Plain airfoil file");
            itype = 1;
            start = 0;
        }
        _ => {
            name = line1.to_string();
            match getflt(line2, 4) {
                None => return fail("\n File READ error.  Unrecognizable file format"),
                Some(a) if a.len() < 2 => return fail("\n File READ error.  Unrecognizable file format"),
                Some(a) if a.len() == 2 => {
                    crate::outln!("\n Labeled airfoil file.  Name:  {}", name);
                    itype = 2;
                    start = idx + 1;
                }
                Some(_) => {
                    crate::outln!("\n MSES airfoil file.  Name:  {}", name);
                    itype = 3;
                    ispars = line2.to_string();
                    start = idx2 + 1;
                }
            }
        }
    }
    let mut x = vec![0.0];
    let mut y = vec![0.0];
    for l in &lines[start..] {
        if is_comment(l) {
            continue;
        }
        let a = match getflt(l, 2) {
            Some(a) => a,
            None => return fail("\n File READ error.  Unrecognizable file format"),
        };
        if a.len() < 2 {
            continue;
        }
        if a[0] == 999.0 && a[1] == 999.0 {
            break; // first element only
        }
        if x.len() > nmax {
            return fail(&format!("\n Buffer array size exceeded\n Maximum number of points: {:4}", nmax));
        }
        x.push(a[0]);
        y.push(a[1]);
    }
    let n = x.len() - 1;
    AreadResult { x, y, n, name, ispars, itype }
}

impl XFoil {
    fn buffer_geopar(&mut self) {
        let nb = self.nb;
        scalc(&self.xb, &self.yb, &mut self.sb, nb);
        segspl(&self.xb, &mut self.xbp, &self.sb, nb);
        segspl(&self.yb, &mut self.ybp, &self.sb, nb);
        let g = geopar(&self.xb, &self.xbp, &self.yb, &self.ybp, &self.sb, nb, true);
        self.sble = g.sle;
        self.chordb = g.chord;
        self.areab = g.area;
        self.radble = g.radle;
        self.angbte = g.angte;
        self.ei11ba = g.ei11a;
        self.ei22ba = g.ei22a;
        self.apx1ba = g.apx1a;
        self.apx2ba = g.apx2a;
        self.ei11bt = g.ei11t;
        self.ei22bt = g.ei22t;
        self.apx1bt = g.apx1t;
        self.apx2bt = g.apx2t;
        self.thickb = g.thick;
        self.cambrb = g.cambr;
    }

    /// Set buffer airfoil from external coordinates (1-based arrays).
    pub fn set_buffer(&mut self, xb: &[f64], yb: &[f64], nb: usize, name: &str) {
        self.nb = nb;
        for i in 1..=nb {
            self.xb[i] = xb[i];
            self.yb[i] = yb[i];
        }
        self.name = name.to_string();
        self.lclock = false;
        self.xbf = 0.0;
        self.ybf = 0.0;
        self.lbflap = false;
        self.buffer_geopar();
    }

    /// NACA 4/5 digit airfoil to buffer, then panel it (NACA). Returns false if illegal.
    pub fn naca(&mut self, ides: i32) -> bool {
        let nside = IQX / 3;
        let r = if ides <= 9999 {
            Some(naca4(ides, nside))
        } else if ides <= 25099 {
            naca5(ides, nside)
        } else {
            crate::outln!(" This designation not implemented.");
            return false;
        };
        let (xb, yb, nb, name) = match r {
            Some(r) => r,
            None => return false,
        };
        self.set_buffer(&xb, &yb, nb, &name);
        crate::outln!("\n Buffer airfoil set using{:4} points", nb);
        self.pangen(true);
        true
    }

    /// Read buffer airfoil from file (LOAD). Returns itype (0 = failure).
    pub fn load(&mut self, fname: &str) -> i32 {
        let r = aread(fname, IBX);
        if r.itype == 0 {
            return 0;
        }
        self.nb = r.n;
        for i in 1..=r.n {
            self.xb[i] = r.x[i];
            self.yb[i] = r.y[i];
        }
        self.name = r.name.trim_end().to_string();
        if r.itype == 3 {
            self.ispars = r.ispars;
        }
        let f = fname.trim();
        self.prefix = match f.find('.') {
            Some(k) => f[..k].to_string(),
            None => f.to_string(),
        };
        let nb = self.nb;
        let mut area = 0.0;
        for i in 1..=nb {
            let ip = if i == nb { 1 } else { i + 1 };
            area += 0.5 * (self.yb[i] + self.yb[ip]) * (self.xb[i] - self.xb[ip]);
        }
        if area >= 0.0 {
            self.lclock = false;
            crate::outln!("\n Number of input coordinate points:{:4}\n Counterclockwise ordering", nb);
        } else {
            self.lclock = true;
            crate::outln!("\n Number of input coordinate points:{:4}\n Clockwise ordering", nb);
            self.xb[1..=nb].reverse();
            self.yb[1..=nb].reverse();
        }
        if self.lnorm {
            norm(&mut self.xb, &mut self.xbp, &mut self.yb, &mut self.ybp, &mut self.sb, nb);
            crate::outln!("\n Airfoil has been normalized");
        }
        self.buffer_geopar();
        let xble = seval(self.sble, &self.xb, &self.xbp, &self.sb, nb);
        let yble = seval(self.sble, &self.yb, &self.ybp, &self.sb, nb);
        let xbte = 0.5 * (self.xb[1] + self.xb[nb]);
        let ybte = 0.5 * (self.yb[1] + self.yb[nb]);
        crate::outln!("\n  LE  x,y  ={:10.5}{:10.5}  |   Chord ={:10.5}", xble, yble, self.chordb);
        crate::outln!("  TE  x,y  ={:10.5}{:10.5}  |", xbte, ybte);
        if r.itype <= 2 {
            let ch = self.chordb;
            let rnd = |v: f64| ((20.0 * (v / ch).abs() + 0.5).trunc() / 20.0) * ch.abs() * v.signum();
            let xinl = rnd(xble - 2.0 * ch);
            let xout = rnd(xble + 3.0 * ch);
            let ybot = rnd(yble - 2.5 * ch);
            let ytop = rnd(yble + 3.5 * ch);
            self.ispars = format!(" {:8.2}{:8.2}{:8.2}{:8.2}", xinl, xout, ybot, ytop);
        }
        self.xbf = 0.0;
        self.ybf = 0.0;
        self.lbflap = false;
        r.itype
    }

    /// Reset solution flags after geometry change.
    pub fn invalidate_geometry(&mut self) {
        self.qw_valid = false;
        self.lgamu = false;
        self.lqinu = false;
        self.lwake = false;
        self.lqaij = false;
        self.ladij = false;
        self.lwdij = false;
        self.lipan = false;
        self.lvconv = false;
        self.lscini = false;
    }

    /// Copy buffer airfoil directly to current airfoil (ABCOPY).
    pub fn abcopy(&mut self, lconf: bool) {
        if self.nb <= 1 {
            crate::outln!(" ABCOPY: Buffer airfoil not available.");
            return;
        }
        if self.nb > IQX - 5 {
            crate::outln!(" Maximum number of panel nodes  : {}", IQX - 5);
            crate::outln!(" Number of buffer airfoil points: {}", self.nb);
            crate::outln!(" Current airfoil cannot be set.");
            crate::outln!(" Try executing PANE at Top Level instead.");
            return;
        }
        if self.n != self.nb {
            self.lblini = false;
        }
        self.n = self.nb;
        for i in 1..=self.n {
            self.x[i] = self.xb[i];
            self.y[i] = self.yb[i];
        }
        self.lgsame = true;
        if self.lbflap {
            self.xof = self.xbf;
            self.yof = self.ybf;
            self.lflap = true;
        }
        // strip out doubled points
        let mut i = 1;
        loop {
            i += 1;
            if self.x[i - 1] == self.x[i] && self.y[i - 1] == self.y[i] {
                for j in i..self.n {
                    self.x[j] = self.x[j + 1];
                    self.y[j] = self.y[j + 1];
                }
                self.n -= 1;
            }
            if i >= self.n {
                break;
            }
        }
        self.finish_current_airfoil();
        self.invalidate_geometry();
        if lconf {
            crate::outln!("\n Current airfoil nodes set from buffer airfoil nodes ({:4} )", self.n);
        }
    }

    /// Spline current airfoil, find LE/TE, normals, panel angles.
    pub fn finish_current_airfoil(&mut self) {
        let n = self.n;
        scalc(&self.x, &self.y, &mut self.s, n);
        segspl(&self.x, &mut self.xp, &self.s, n);
        segspl(&self.y, &mut self.yp, &self.s, n);
        ncalc(&self.x, &self.y, &self.s, n, &mut self.nx, &mut self.ny);
        self.sle = lefind(&self.x, &self.xp, &self.y, &self.yp, &self.s, n);
        self.xle = seval(self.sle, &self.x, &self.xp, &self.s, n);
        self.yle = seval(self.sle, &self.y, &self.yp, &self.s, n);
        self.xte = 0.5 * (self.x[1] + self.x[n]);
        self.yte = 0.5 * (self.y[1] + self.y[n]);
        self.chord = ((self.xte - self.xle).powi(2) + (self.yte - self.yle).powi(2)).sqrt();
        self.tecalc();
        self.apcalc();
    }

    /// Curvature-based paneling from buffer airfoil (PANGEN).
    pub fn pangen(&mut self, shopar: bool) {
        let nb = self.nb;
        if nb < 2 {
            crate::outln!(" PANGEN: Buffer airfoil not available.");
            self.n = 0;
            return;
        }
        let ipfac = 5;
        self.n = self.npan;
        let n0 = self.npan;
        scalc(&self.xb, &self.yb, &mut self.sb, nb);
        segspl(&self.xb, &mut self.xbp, &self.sb, nb);
        segspl(&self.yb, &mut self.ybp, &self.sb, nb);
        let (xb, yb, xbp, ybp, sb) = (&self.xb, &self.yb, &self.xbp, &self.ybp, &self.sb);
        let sbref = 0.5 * (sb[nb] - sb[1]);
        let mut w1 = vec![0.0; 5 * IBX + 2];
        let mut w2 = vec![0.0; 5 * IBX + 2];
        let mut w3 = vec![0.0; 5 * IBX + 2];
        let mut w4 = vec![0.0; 5 * IBX + 2];
        let mut w5 = vec![0.0; nb + 2];
        let mut w6 = vec![0.0; nb + 2];
        for i in 1..=nb {
            w5[i] = curv(sb[i], xb, xbp, yb, ybp, sb, nb).abs() * sbref;
        }
        let sble = lefind(xb, xbp, yb, ybp, sb, nb);
        let cvle = curv(sble, xb, xbp, yb, ybp, sb, nb).abs() * sbref;
        let mut ible = 0;
        for i in 1..nb {
            if sble == sb[i] && sble == sb[i + 1] {
                ible = i;
                crate::outln!("\n Sharp leading edge");
                break;
            }
        }
        let xble = seval(sble, xb, xbp, sb, nb);
        let yble = seval(sble, yb, ybp, sb, nb);
        let xbte = 0.5 * (xb[1] + xb[nb]);
        let ybte = 0.5 * (yb[1] + yb[nb]);
        let chbsq = (xbte - xble).powi(2) + (ybte - yble).powi(2);
        let nk = 3;
        let mut cvsum = 0.0;
        for k in -nk..=nk {
            let frac = k as f64 / nk as f64;
            let sbk = sble + frac * sbref / cvle.max(20.0);
            let cvk = curv(sbk, xb, xbp, yb, ybp, sb, nb).abs() * sbref;
            cvsum += cvk;
        }
        let mut cvavg = cvsum / (2 * nk + 1) as f64;
        if ible != 0 {
            cvavg = 10.0;
        }
        let cc = 6.0 * self.cvpar;
        let cvte = cvavg * self.cterat;
        w5[1] = cvte;
        w5[nb] = cvte;
        let smool = (1.0 / cvavg.max(20.0)).max(0.25 / (self.npan / 2) as f64);
        let smoosq = (smool * sbref).powi(2);
        w2[1] = 1.0;
        w3[1] = 0.0;
        for i in 2..nb {
            let dsm = sb[i] - sb[i - 1];
            let dsp = sb[i + 1] - sb[i];
            let dso = 0.5 * (sb[i + 1] - sb[i - 1]);
            if dsm == 0.0 || dsp == 0.0 {
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
            } else {
                w1[i] = smoosq * (-1.0 / dsm) / dso;
                w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
                w3[i] = smoosq * (-1.0 / dsp) / dso;
            }
        }
        w1[nb] = 0.0;
        w2[nb] = 1.0;
        for i in 2..nb {
            if sb[i] == sble || i == ible || i == ible + 1 {
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
                w5[i] = cvle;
            } else if sb[i - 1] < sble && sb[i] > sble {
                let dsm = sb[i - 1] - sb[i - 2];
                let dsp = sble - sb[i - 1];
                let dso = 0.5 * (sble - sb[i - 2]);
                w1[i - 1] = smoosq * (-1.0 / dsm) / dso;
                w2[i - 1] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
                w3[i - 1] = 0.0;
                w5[i - 1] += smoosq * cvle / (dsp * dso);
                let dsm = sb[i] - sble;
                let dsp = sb[i + 1] - sb[i];
                let dso = 0.5 * (sb[i + 1] - sble);
                w1[i] = 0.0;
                w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
                w3[i] = smoosq * (-1.0 / dsp) / dso;
                w5[i] += smoosq * cvle / (dsm * dso);
                break;
            }
        }
        for i in 2..nb {
            let xoc = ((xb[i] - xble) * (xbte - xble) + (yb[i] - yble) * (ybte - yble)) / chbsq;
            if sb[i] < sble {
                if xoc > self.xsref1 && xoc < self.xsref2 {
                    w1[i] = 0.0;
                    w2[i] = 1.0;
                    w3[i] = 0.0;
                    w5[i] = cvle * self.ctrrat;
                }
            } else if xoc > self.xpref1 && xoc < self.xpref2 {
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
                w5[i] = cvle * self.ctrrat;
            }
        }
        if ible == 0 {
            trisol(&mut w2, &w1, &mut w3, &mut w5, nb);
        } else {
            trisol(&mut w2, &w1, &mut w3, &mut w5, ible);
            let o = ible;
            trisol(&mut w2[o..], &w1[o..], &mut w3[o..], &mut w5[o..], nb - ible);
        }
        let mut cvmax: f64 = 0.0;
        for i in 1..=nb {
            cvmax = cvmax.max(w5[i].abs());
        }
        for i in 1..=nb {
            w5[i] /= cvmax;
        }
        segspl(&w5, &mut w6, sb, nb);

        let nn = ipfac * (n0 - 1) + 1;
        let rdste = 0.667;
        let rtf = (rdste - 1.0) * 2.0 + 1.0;
        let mut snew = vec![0.0; nn + 2];
        let mut nn1 = 0;
        if ible == 0 {
            let dsavg = (sb[nb] - sb[1]) / ((nn - 3) as f64 + 2.0 * rtf);
            snew[1] = sb[1];
            for i in 2..nn {
                snew[i] = sb[1] + dsavg * ((i - 2) as f64 + rtf);
            }
            snew[nn] = sb[nb];
        } else {
            let nfrac1 = (n0 * ible) / nb;
            nn1 = ipfac * (nfrac1 - 1) + 1;
            let dsavg1 = (sble - sb[1]) / ((nn1 - 2) as f64 + rtf);
            snew[1] = sb[1];
            for i in 2..=nn1 {
                snew[i] = sb[1] + dsavg1 * ((i - 2) as f64 + rtf);
            }
            let nn2 = nn - nn1 + 1;
            let dsavg2 = (sb[nb] - sble) / ((nn2 - 2) as f64 + rtf);
            for i in 2..nn2 {
                snew[i - 1 + nn1] = sble + dsavg2 * ((i - 2) as f64 + rtf);
            }
            snew[nn] = sb[nb];
        }

        let mut ok = false;
        for _iter in 0..20 {
            let mut cv1 = seval(snew[1], &w5, &w6, sb, nb);
            let mut cv2 = seval(snew[2], &w5, &w6, sb, nb);
            let mut cvs1 = deval(snew[1], &w5, &w6, sb, nb);
            let mut cvs2 = deval(snew[2], &w5, &w6, sb, nb);
            let _ = (cv1, cvs1);
            let mut cavm = (cv1 * cv1 + cv2 * cv2).sqrt();
            let (mut cavm_s1, mut cavm_s2) = if cavm == 0.0 { (0.0, 0.0) } else { (cvs1 * cv1 / cavm, cvs2 * cv2 / cavm) };
            for i in 2..nn {
                let dsm = snew[i] - snew[i - 1];
                let dsp = snew[i] - snew[i + 1];
                let cv3 = seval(snew[i + 1], &w5, &w6, sb, nb);
                let cvs3 = deval(snew[i + 1], &w5, &w6, sb, nb);
                let cavp = (cv3 * cv3 + cv2 * cv2).sqrt();
                let (cavp_s2, cavp_s3) = if cavp == 0.0 { (0.0, 0.0) } else { (cvs2 * cv2 / cavp, cvs3 * cv3 / cavp) };
                let fm = cc * cavm + 1.0;
                let fp = cc * cavp + 1.0;
                let rez = dsp * fp + dsm * fm;
                w1[i] = -fm + cc * dsm * cavm_s1;
                w2[i] = fp + fm + cc * (dsp * cavp_s2 + dsm * cavm_s2);
                w3[i] = -fp + cc * dsp * cavp_s3;
                w4[i] = -rez;
                cv1 = cv2;
                cv2 = cv3;
                cvs1 = cvs2;
                cvs2 = cvs3;
                cavm = cavp;
                cavm_s1 = cavp_s2;
                cavm_s2 = cavp_s3;
            }
            w2[1] = 1.0;
            w3[1] = 0.0;
            w4[1] = 0.0;
            w1[nn] = 0.0;
            w2[nn] = 1.0;
            w4[nn] = 0.0;
            if rtf != 1.0 {
                let i = 2;
                w4[i] = -((snew[i] - snew[i - 1]) + rtf * (snew[i] - snew[i + 1]));
                w1[i] = -1.0;
                w2[i] = 1.0 + rtf;
                w3[i] = -rtf;
                let i = nn - 1;
                w4[i] = -((snew[i] - snew[i + 1]) + rtf * (snew[i] - snew[i - 1]));
                w3[i] = -1.0;
                w2[i] = 1.0 + rtf;
                w1[i] = -rtf;
            }
            if ible != 0 {
                let i = nn1;
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
                w4[i] = sble - snew[i];
            }
            trisol(&mut w2, &w1, &mut w3, &mut w4, nn);
            let mut rlx = 1.0;
            let mut dmax: f64 = 0.0;
            for i in 1..nn {
                let ds = snew[i + 1] - snew[i];
                let dds = w4[i + 1] - w4[i];
                let dsrat = 1.0 + rlx * dds / ds;
                if dsrat > 4.0 {
                    rlx = (4.0 - 1.0) * ds / dds;
                }
                if dsrat < 0.2 {
                    rlx = (0.2 - 1.0) * ds / dds;
                }
                dmax = dmax.max(w4[i].abs());
            }
            for i in 2..nn {
                snew[i] += rlx * w4[i];
            }
            if dmax.abs() < 1.0e-3 {
                ok = true;
                break;
            }
        }
        if !ok {
            crate::outln!(" Paneling convergence failed.  Continuing anyway...");
        }
        let mut n = n0;
        for i in 1..=n {
            let ind = ipfac * (i - 1) + 1;
            self.s[i] = snew[ind];
            self.x[i] = seval(snew[ind], xb, xbp, sb, nb);
            self.y[i] = seval(snew[ind], yb, ybp, sb, nb);
        }
        // corners (double points) in buffer airfoil
        for ib in 1..nb {
            if sb[ib] == sb[ib + 1] {
                let xbcorn = xb[ib];
                let ybcorn = yb[ib];
                let sbcorn = sb[ib];
                for i in 1..=n {
                    if self.s[i] <= sbcorn {
                        continue;
                    }
                    for j in (i..=n).rev() {
                        self.x[j + 1] = self.x[j];
                        self.y[j + 1] = self.y[j];
                        self.s[j + 1] = self.s[j];
                    }
                    n += 1;
                    if n > IQX - 1 {
                        panic!("PANEL: Too many panels. Increase IQX");
                    }
                    self.x[i] = xbcorn;
                    self.y[i] = ybcorn;
                    self.s[i] = sbcorn;
                    if i >= 3 {
                        self.s[i - 1] = 0.5 * (self.s[i] + self.s[i - 2]);
                        self.x[i - 1] = seval(self.s[i - 1], xb, xbp, sb, nb);
                        self.y[i - 1] = seval(self.s[i - 1], yb, ybp, sb, nb);
                    }
                    if i + 2 <= n {
                        self.s[i + 1] = 0.5 * (self.s[i] + self.s[i + 2]);
                        self.x[i + 1] = seval(self.s[i + 1], xb, xbp, sb, nb);
                        self.y[i + 1] = seval(self.s[i + 1], yb, ybp, sb, nb);
                    }
                    break;
                }
            }
        }
        self.n = n;
        self.finish_current_airfoil();
        self.invalidate_geometry();
        self.lblini = false;
        self.lgsame = false;
        self.lqspec = false;
        if self.lbflap {
            self.xof = self.xbf;
            self.yof = self.ybf;
            self.lflap = true;
        }
        if self.sharp {
            crate::outln!("\n Sharp trailing edge");
        } else {
            let gap = ((self.x[1] - self.x[n]).powi(2) + (self.y[1] - self.y[n]).powi(2)).sqrt();
            crate::outln!("\n Blunt trailing edge.  Gap ={:9.5}", gap);
        }
        if shopar {
            crate::outln!(
                "\n Paneling parameters used...\n   Number of panel nodes      {:4}\n   Panel bunching parameter   {:6.3}\n   TE/LE panel density ratio  {:6.3}\n   Refined-area/LE panel density ratio   {:6.3}\n   Top    side refined area x/c limits {:6.3}{:6.3}\n   Bottom side refined area x/c limits {:6.3}{:6.3}",
                self.npan, self.cvpar, self.cterat, self.ctrrat, self.xsref1, self.xsref2, self.xpref1, self.xpref2
            );
        }
    }
}
