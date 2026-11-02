//! 2D polar storage, interpolation and I/O.

use serde::{Deserialize, Serialize};

/// One converged 2D operating point.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct FoilPoint {
    /// angle of attack, degrees
    pub alpha: f64,
    pub cl: f64,
    /// total drag coefficient
    pub cd: f64,
    /// pressure drag
    pub cdp: f64,
    /// friction drag
    pub cdf: f64,
    /// pitching moment about the moment reference point (default c/4)
    pub cm: f64,
    /// transition x/c, upper and lower surfaces
    pub xtr_top: f64,
    pub xtr_bot: f64,
    /// minimum pressure coefficient
    pub cpmin: f64,
    /// Reynolds and Mach number actually used (vary with CL for type 2/3)
    pub re: f64,
    pub mach: f64,
    /// Newton iterations used; 0 for inviscid
    pub iterations: u32,
    /// solved after a recovery step (re-initialised BL or sub-stepping)
    pub recovered: bool,
    /// trailing-edge flap deflection, degrees (type 6 polars)
    #[serde(default)]
    pub flap: f64,
}

impl FoilPoint {
    pub fn ld(&self) -> f64 {
        if self.cd.abs() > 1e-12 {
            self.cl / self.cd
        } else {
            0.0
        }
    }
}

/// A 2D polar: points sorted by alpha (or Re for type-4 polars).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct FoilPolar {
    pub name: String,
    pub foil: String,
    /// 1 fixed speed, 2 fixed lift, 3 rubber chord, 4 fixed alpha
    pub polar_type: u8,
    pub re: f64,
    pub mach: f64,
    pub ncrit: f64,
    pub xtr_top: f64,
    pub xtr_bot: f64,
    pub viscous: bool,
    pub points: Vec<FoilPoint>,
    /// operating points that could not be converged (alpha or CL or Re spec)
    #[serde(default)]
    pub failed: Vec<f64>,
}

impl FoilPolar {
    /// Insert or replace a point, keeping the list sorted by alpha
    /// (by Re for type 4).
    pub fn insert(&mut self, p: FoilPoint) {
        let key = |q: &FoilPoint| match self.polar_type {
            4 => q.re,
            6 => q.flap,
            _ => q.alpha,
        };
        let k = key(&p);
        let tol = if self.polar_type == 4 { 1e-6 * k.abs().max(1.0) } else { 1e-4 };
        if let Some(i) = self.points.iter().position(|q| (key(q) - k).abs() < tol) {
            self.points[i] = p;
            return;
        }
        let i = self.points.partition_point(|q| key(q) < k);
        self.points.insert(i, p);
    }

    pub fn alpha_range(&self) -> Option<(f64, f64)> {
        Some((self.points.first()?.alpha, self.points.last()?.alpha))
    }

    pub fn cl_max(&self) -> Option<FoilPoint> {
        self.points.iter().cloned().max_by(|a, b| a.cl.total_cmp(&b.cl))
    }

    pub fn ld_max(&self) -> Option<FoilPoint> {
        self.points.iter().cloned().filter(|p| p.cd > 0.0).max_by(|a, b| a.ld().total_cmp(&b.ld()))
    }

    /// Zero-lift angle (degrees), by linear interpolation.
    pub fn alpha0(&self) -> Option<f64> {
        let p = &self.points;
        for w in p.windows(2) {
            if (w[0].cl <= 0.0 && w[1].cl >= 0.0) || (w[0].cl >= 0.0 && w[1].cl <= 0.0) {
                let d = w[1].cl - w[0].cl;
                if d.abs() < 1e-12 {
                    return Some(w[0].alpha);
                }
                return Some(w[0].alpha - w[0].cl * (w[1].alpha - w[0].alpha) / d);
            }
        }
        None
    }

    /// Lift-curve slope (per degree) from a least-squares fit over the
    /// linear range (|alpha - alpha0| < 6 deg).
    pub fn cl_alpha(&self) -> Option<f64> {
        let a0 = self.alpha0().unwrap_or(0.0);
        let pts: Vec<&FoilPoint> = self.points.iter().filter(|p| (p.alpha - a0).abs() < 6.0).collect();
        if pts.len() < 2 {
            return None;
        }
        let n = pts.len() as f64;
        let sx: f64 = pts.iter().map(|p| p.alpha).sum();
        let sy: f64 = pts.iter().map(|p| p.cl).sum();
        let sxx: f64 = pts.iter().map(|p| p.alpha * p.alpha).sum();
        let sxy: f64 = pts.iter().map(|p| p.alpha * p.cl).sum();
        let d = n * sxx - sx * sx;
        if d.abs() < 1e-12 {
            return None;
        }
        Some((n * sxy - sx * sy) / d)
    }

    /// Interpolate a point at a given CL, walking out from the zero-lift
    /// region so the pre-stall branch is used.  Returns `(point, out_of_range)`.
    pub fn at_cl(&self, cl: f64) -> Option<(FoilPoint, bool)> {
        let p = &self.points;
        if p.is_empty() {
            return None;
        }
        if p.len() == 1 {
            return Some((p[0], true));
        }
        // start from the point closest to cl=0 and walk in the direction of cl
        let i0 = (0..p.len()).min_by(|&a, &b| p[a].cl.abs().total_cmp(&p[b].cl.abs())).unwrap();
        let up = cl >= p[i0].cl;
        let mut i = i0;
        loop {
            let j = if up { i + 1 } else { i.wrapping_sub(1) };
            if j >= p.len() {
                return Some((p[i], true));
            }
            let (a, b) = (&p[i], &p[j]);
            let (lo, hi) = if a.cl < b.cl { (a.cl, b.cl) } else { (b.cl, a.cl) };
            if cl >= lo && cl <= hi {
                let t = if (b.cl - a.cl).abs() > 1e-12 { (cl - a.cl) / (b.cl - a.cl) } else { 0.0 };
                return Some((lerp_point(a, b, t), false));
            }
            // past stall (CL decreasing): stop at the max
            if (up && b.cl < a.cl) || (!up && b.cl > a.cl) {
                return Some((*a, true));
            }
            i = j;
        }
    }

    /// Interpolate a point at a given alpha (clamped, flags out-of-range).
    pub fn at_alpha(&self, alpha: f64) -> Option<(FoilPoint, bool)> {
        let p = &self.points;
        let first = p.first()?;
        let last = p.last()?;
        if alpha <= first.alpha {
            return Some((*first, alpha < first.alpha - 1e-9));
        }
        if alpha >= last.alpha {
            return Some((*last, alpha > last.alpha + 1e-9));
        }
        let j = p.partition_point(|q| q.alpha < alpha);
        let (a, b) = (&p[j - 1], &p[j]);
        let t = (alpha - a.alpha) / (b.alpha - a.alpha);
        Some((lerp_point(a, b, t), false))
    }

    /// CSV text with a header line.
    pub fn to_csv(&self) -> String {
        let mut s = String::from("alpha,CL,CD,CDp,CDf,Cm,Xtr_top,Xtr_bot,Cpmin,Re,Mach,CL/CD,iterations,recovered\n");
        for p in &self.points {
            s.push_str(&format!(
                "{:.3},{:.5},{:.6},{:.6},{:.6},{:.5},{:.4},{:.4},{:.4},{:.0},{:.4},{:.3},{},{}\n",
                p.alpha, p.cl, p.cd, p.cdp, p.cdf, p.cm, p.xtr_top, p.xtr_bot, p.cpmin, p.re, p.mach, p.ld(), p.iterations, p.recovered
            ));
        }
        s
    }

    /// XFoil-style polar text (readable by XFoil/xflr5 importers).
    pub fn to_xfoil_text(&self) -> String {
        let mut s = String::new();
        s.push_str("\n       Foiley polar (XFoil 6.97 format)\n\n");
        s.push_str(&format!(" Calculated polar for: {}\n\n", self.foil));
        s.push_str(&format!(" {} {} Reynolds number fixed          Mach number fixed\n\n", self.polar_type.min(3), self.polar_type.min(3)));
        s.push_str(&format!(" xtrf = {:7.3} (top) {:9.3} (bottom)\n", self.xtr_top, self.xtr_bot));
        s.push_str(&format!(
            " Mach = {:7.3}     Re = {:9.3} e 6     Ncrit = {:7.3}\n\n",
            self.mach,
            self.re / 1e6,
            self.ncrit
        ));
        s.push_str("  alpha    CL        CD       CDp       CM     Top_Xtr  Bot_Xtr\n");
        s.push_str(" ------ -------- --------- --------- -------- -------- --------\n");
        for p in &self.points {
            s.push_str(&format!(
                " {:7.3}{:9.4}{:10.5}{:10.5}{:9.4}{:9.4}{:9.4}\n",
                p.alpha, p.cl, p.cd, p.cdp, p.cm, p.xtr_top, p.xtr_bot
            ));
        }
        s
    }

    /// Parse an XFoil/xflr5 polar text file (.txt/.pol).
    pub fn from_xfoil_text(text: &str, foil_fallback: &str) -> Option<FoilPolar> {
        let mut pol = FoilPolar { polar_type: 1, viscous: true, ncrit: 9.0, xtr_top: 1.0, xtr_bot: 1.0, foil: foil_fallback.into(), ..Default::default() };
        let mut in_table = false;
        for line in text.lines() {
            let l = line.trim();
            if let Some(r) = l.strip_prefix("Calculated polar for:") {
                pol.foil = r.trim().to_string();
            }
            if l.starts_with("Mach =") {
                let t: Vec<&str> = l.split_whitespace().collect();
                // Mach = 0.000     Re = 1.000 e 6     Ncrit = 9.000
                if let Some(m) = t.get(2).and_then(|v| v.parse().ok()) {
                    pol.mach = m;
                }
                if let Some(i) = t.iter().position(|v| *v == "Re") {
                    let mant: f64 = t.get(i + 2).and_then(|v| v.parse().ok()).unwrap_or(0.0);
                    let exp: i32 = t.get(i + 4).and_then(|v| v.parse().ok()).unwrap_or(0);
                    pol.re = mant * 10f64.powi(exp);
                }
                if let Some(i) = t.iter().position(|v| *v == "Ncrit") {
                    pol.ncrit = t.get(i + 2).and_then(|v| v.parse().ok()).unwrap_or(9.0);
                }
            }
            if l.starts_with("xtrf") {
                let nums: Vec<f64> = l.split_whitespace().filter_map(|v| v.parse().ok()).collect();
                if nums.len() >= 2 {
                    pol.xtr_top = nums[0];
                    pol.xtr_bot = nums[1];
                }
            }
            if l.starts_with("---") || l.starts_with("------") {
                in_table = true;
                continue;
            }
            if in_table {
                let v: Vec<f64> = l.split_whitespace().filter_map(|x| x.parse().ok()).collect();
                if v.len() >= 5 {
                    pol.points.push(FoilPoint {
                        alpha: v[0],
                        cl: v[1],
                        cd: v[2],
                        cdp: v[3],
                        cdf: v[2] - v[3],
                        cm: v[4],
                        xtr_top: *v.get(5).unwrap_or(&1.0),
                        xtr_bot: *v.get(6).unwrap_or(&1.0),
                        re: pol.re,
                        mach: pol.mach,
                        ..Default::default()
                    });
                }
            }
        }
        if pol.points.is_empty() {
            return None;
        }
        pol.points.sort_by(|a, b| a.alpha.total_cmp(&b.alpha));
        pol.name = format!("T1_Re{:.3}_M{:.2}_N{:.1}", pol.re / 1e6, pol.mach, pol.ncrit);
        Some(pol)
    }
}

pub fn lerp_point(a: &FoilPoint, b: &FoilPoint, t: f64) -> FoilPoint {
    let l = |x: f64, y: f64| x + t * (y - x);
    FoilPoint {
        alpha: l(a.alpha, b.alpha),
        cl: l(a.cl, b.cl),
        cd: l(a.cd, b.cd),
        cdp: l(a.cdp, b.cdp),
        cdf: l(a.cdf, b.cdf),
        cm: l(a.cm, b.cm),
        xtr_top: l(a.xtr_top, b.xtr_top),
        xtr_bot: l(a.xtr_bot, b.xtr_bot),
        cpmin: l(a.cpmin, b.cpmin),
        re: l(a.re, b.re),
        mach: l(a.mach, b.mach),
        iterations: a.iterations.max(b.iterations),
        recovered: a.recovered || b.recovered,
        flap: l(a.flap, b.flap),
    }
}
