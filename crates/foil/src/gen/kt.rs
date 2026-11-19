//! Karman-Trefftz sections: the conformal map
//!
//! ```text
//! z = n ((zeta+1)^n + (zeta-1)^n) / ((zeta+1)^n - (zeta-1)^n),   n = 2 - tau/pi
//! ```
//!
//! of a circle through zeta = 1 centred at (-x_c, y_c).  tau is the
//! trailing-edge angle (tau = 0 gives the Joukowsky section).  Because the
//! potential flow about the circle is known, every section comes with its
//! exact inviscid lift and pressure distribution ([`KarmanTrefftz::exact`]),
//! which makes it a reference case for the panel solvers.

use super::{check_contour, Recipe, NSIDE};
use crate::{Foil, FoilError};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, PartialEq)]
struct C(f64, f64);

impl C {
    fn polar(r: f64, th: f64) -> C {
        C(r * th.cos(), r * th.sin())
    }
    fn abs(self) -> f64 {
        self.0.hypot(self.1)
    }
    fn arg(self) -> f64 {
        self.1.atan2(self.0)
    }
    fn powf(self, n: f64) -> C {
        C::polar(self.abs().powf(n), self.arg() * n)
    }
    fn inv(self) -> C {
        let d = self.0 * self.0 + self.1 * self.1;
        C(self.0 / d, -self.1 / d)
    }
}
impl std::ops::Add for C {
    type Output = C;
    fn add(self, o: C) -> C {
        C(self.0 + o.0, self.1 + o.1)
    }
}
impl std::ops::Sub for C {
    type Output = C;
    fn sub(self, o: C) -> C {
        C(self.0 - o.0, self.1 - o.1)
    }
}
impl std::ops::Mul for C {
    type Output = C;
    fn mul(self, o: C) -> C {
        C(self.0 * o.0 - self.1 * o.1, self.0 * o.1 + self.1 * o.0)
    }
}
impl std::ops::Mul<f64> for C {
    type Output = C;
    fn mul(self, k: f64) -> C {
        C(self.0 * k, self.1 * k)
    }
}
impl std::ops::Div for C {
    type Output = C;
    fn div(self, o: C) -> C {
        let d = o.0 * o.0 + o.1 * o.1;
        C((self.0 * o.0 + self.1 * o.1) / d, (self.1 * o.0 - self.0 * o.1) / d)
    }
}

/// Karman-Trefftz parameters.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct KarmanTrefftz {
    /// Leftward offset of the circle centre (sets the thickness, > 0).
    pub x_center: f64,
    /// Upward offset of the circle centre (sets the camber).
    pub y_center: f64,
    /// Trailing-edge angle, degrees (0 = Joukowsky, cusped).
    pub te_angle: f64,
}

impl Default for KarmanTrefftz {
    fn default() -> Self {
        KarmanTrefftz { x_center: 0.08, y_center: 0.08, te_angle: 10.0 }
    }
}

/// Exact inviscid solution at one angle of attack.
#[derive(Clone, Debug, PartialEq)]
pub struct KtSolution {
    /// Lift coefficient (chord of the normalized section).
    pub cl: f64,
    /// Moment coefficient about the quarter chord.
    pub cm: f64,
    /// (x, y, Cp) on the normalized section, in the order of the foil's points
    /// (the trailing-edge stagnation point excluded).
    pub cp: Vec<(f64, f64, f64)>,
}

/// Normalization of the mapped section: z_norm = (z - le) rot / chord.
struct Frame {
    le: C,
    rot: C,
    chord: f64,
}

impl KarmanTrefftz {
    fn n(&self) -> f64 {
        2.0 - self.te_angle.to_radians() / PI
    }
    fn center(&self) -> C {
        C(-self.x_center, self.y_center)
    }
    fn radius(&self) -> f64 {
        (C(1.0, 0.0) - self.center()).abs()
    }
    /// Angle of the trailing edge (zeta = 1) seen from the centre.
    fn theta_te(&self) -> f64 {
        (C(1.0, 0.0) - self.center()).arg()
    }

    pub fn validate(&self) -> Result<(), FoilError> {
        if ![self.x_center, self.y_center, self.te_angle].iter().all(|v| v.is_finite()) {
            return Err(FoilError::Geometry("Karman-Trefftz: non-finite parameter".into()));
        }
        if !(0.001..=0.5).contains(&self.x_center) || self.y_center.abs() > 0.5 || !(0.0..=60.0).contains(&self.te_angle) {
            return Err(FoilError::Geometry("Karman-Trefftz: need centre x offset 0.001..0.5, |y offset| <= 0.5, TE angle 0..60 deg".into()));
        }
        Ok(())
    }

    fn map(&self, zeta: C) -> C {
        // w = (zeta-1)/(zeta+1) is never real-negative on or outside the
        // circle (that would put zeta on the slit -1..1 inside it), so the
        // principal power is continuous around the contour.
        let n = self.n();
        let r = ((zeta - C(1.0, 0.0)) / (zeta + C(1.0, 0.0))).powf(n);
        (C(1.0, 0.0) + r) / (C(1.0, 0.0) - r) * n
    }

    fn dz_dzeta(&self, zeta: C) -> C {
        let n = self.n();
        let r = ((zeta - C(1.0, 0.0)) / (zeta + C(1.0, 0.0))).powf(n);
        let omr = C(1.0, 0.0) - r;
        r * (4.0 * n * n) / (omr * omr * (zeta * zeta - C(1.0, 0.0)))
    }

    /// Circle angles of the contour points, trailing edge over the upper side
    /// and back (2 NSIDE - 1 points, bunched at both ends of the section).
    fn thetas(&self) -> Vec<f64> {
        let m = 2 * (NSIDE - 1);
        let t0 = self.theta_te();
        (0..=m).map(|i| t0 + 2.0 * PI * i as f64 / m as f64).collect()
    }

    fn frame(&self) -> Frame {
        let c = self.center();
        let r = self.radius();
        let te = self.map(C(1.0, 0.0));
        // leading edge: the contour point farthest from the trailing edge
        let dist = |th: f64| (self.map(c + C::polar(r, th)) - te).abs();
        let t0 = self.theta_te();
        let samples = 2000;
        let (mut best, mut bd) = (t0 + PI, 0.0);
        for i in 1..samples {
            let th = t0 + 2.0 * PI * i as f64 / samples as f64;
            let d = dist(th);
            if d > bd {
                bd = d;
                best = th;
            }
        }
        let (mut a, mut b) = (best - 2.0 * PI / samples as f64, best + 2.0 * PI / samples as f64);
        let g = 0.5 * (5f64.sqrt() - 1.0);
        for _ in 0..100 {
            let (u, v) = (b - g * (b - a), a + g * (b - a));
            if dist(u) > dist(v) {
                b = v;
            } else {
                a = u;
            }
        }
        let le = self.map(c + C::polar(r, 0.5 * (a + b)));
        let d = te - le;
        let chord = d.abs();
        Frame { le, rot: C(d.0 / chord, -d.1 / chord), chord }
    }

    pub fn generate(&self, name: &str) -> Result<Foil, FoilError> {
        self.validate()?;
        let c = self.center();
        let r = self.radius();
        let fr = self.frame();
        let th = self.thetas();
        let last = th.len() - 1;
        let mut x = Vec::with_capacity(th.len());
        let mut y = Vec::with_capacity(th.len());
        for (i, &t) in th.iter().enumerate() {
            // both ends are the trailing edge itself
            let z = if i == 0 || i == last { self.map(C(1.0, 0.0)) } else { self.map(c + C::polar(r, t)) };
            let p = (z - fr.le) * fr.rot * (1.0 / fr.chord);
            x.push(p.0);
            y.push(p.1);
        }
        let mut f = Foil::new(name, x, y)?;
        check_contour(&f, "Karman-Trefftz")?;
        f.note = format!("Karman-Trefftz (centre -{}, {}; TE {} deg)", self.x_center, self.y_center, self.te_angle);
        f.recipe = Some(Recipe::KarmanTrefftz(*self));
        Ok(f)
    }

    /// Complex velocity dW/dzeta about the circle, unit free stream at
    /// `a` radians (mapped plane), circulation `gamma` (clockwise positive).
    fn circle_velocity(&self, a: f64, gamma: f64, zeta: C) -> C {
        let zr = zeta - self.center();
        let r = self.radius();
        C::polar(1.0, -a) - C::polar(1.0, a) * (r * r) * (zr * zr).inv() + C(0.0, gamma / (2.0 * PI)) * zr.inv()
    }

    /// Exact inviscid lift, quarter-chord moment (nose up positive) and Cp
    /// at `alpha` degrees, relative to the chord of the normalized section.
    ///
    /// Lift is Kutta-Joukowski's; the force and moment come from the
    /// Blasius integrals, evaluated on a circle of twice the radius where
    /// the integrand is analytic, so the periodic trapezoidal rule converges
    /// to machine precision.
    pub fn exact(&self, alpha: f64) -> Result<KtSolution, FoilError> {
        self.validate()?;
        if !alpha.is_finite() {
            return Err(FoilError::Geometry("Karman-Trefftz: the angle of attack must be finite".into()));
        }
        // equivalent angles give the same flow: reduce to (-180, 180] first
        // (rem_euclid is exact; subtracting first would round 1e20 - 180)
        let alpha = match alpha.rem_euclid(360.0) {
            r if r > 180.0 => r - 360.0,
            r => r,
        };
        let c = self.center();
        let r = self.radius();
        let fr = self.frame();
        // free stream angle in the mapped (unrotated) plane; z ~ zeta far away
        let a = alpha.to_radians() - fr.rot.arg();
        let beta = -self.theta_te();
        let gamma = 4.0 * PI * r * (a + beta).sin();
        let cl = 2.0 * gamma / fr.chord;

        // Blasius: X - iY = (i/2) oint (dW/dz)^2 dz, M0 = -(1/2) Re oint z (dW/dz)^2 dz
        // (rho = U = 1, counter-clockwise moment about the origin of z)
        let m = 4096;
        let (r2, mut f_int, mut m_int) = (2.0 * r, C(0.0, 0.0), C(0.0, 0.0));
        for k in 0..m {
            let e = C::polar(1.0, 2.0 * PI * k as f64 / m as f64);
            let zeta = c + e * r2;
            let dzeta = C(0.0, 1.0) * e * (r2 * 2.0 * PI / m as f64);
            let w = self.circle_velocity(a, gamma, zeta);
            let g = w * w * self.dz_dzeta(zeta).inv() * dzeta;
            f_int = f_int + g;
            m_int = m_int + self.map(zeta) * g;
        }
        let xy = C(0.0, 0.5) * f_int;
        let (fx, fy) = (xy.0, -xy.1);
        let m0 = -0.5 * m_int.0;
        let te = self.map(C(1.0, 0.0));
        let qc = fr.le + (te - fr.le) * 0.25;
        let m_qc = m0 - (qc.0 * fy - qc.1 * fx);
        // counter-clockwise is nose down for flow from the leading edge
        let cm = -m_qc / (0.5 * fr.chord * fr.chord);

        let th = self.thetas();
        let mut cp = Vec::with_capacity(th.len());
        for &t in &th[1..th.len() - 1] {
            let zeta = c + C::polar(r, t);
            let v = (self.circle_velocity(a, gamma, zeta) / self.dz_dzeta(zeta)).abs();
            let p = (self.map(zeta) - fr.le) * fr.rot * (1.0 / fr.chord);
            cp.push((p.0, p.1, 1.0 - v * v));
        }
        Ok(KtSolution { cl, cm, cp })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_section_is_normalized_and_closed() {
        for kt in [KarmanTrefftz::default(), KarmanTrefftz { x_center: 0.1, y_center: 0.0, te_angle: 0.0 }, KarmanTrefftz { x_center: 0.05, y_center: -0.05, te_angle: 25.0 }] {
            let f = kt.generate("kt").unwrap();
            let n = f.n();
            assert_eq!((f.x[0], f.y[0]), (f.x[n - 1], f.y[n - 1]));
            assert!((f.x[0] - 1.0).abs() < 1e-12 && f.y[0].abs() < 1e-12, "{:?}", (f.x[0], f.y[0]));
            // the leading edge (farthest from the trailing edge) is at the origin
            let ile = f.le_index();
            assert!(f.x[ile] > -1e-9 && f.x[ile] < 1e-3, "{} {}", f.x[ile], f.y[ile]);
            // counter-clockwise: upper surface first
            assert!(f.y[n / 4] > f.y[3 * n / 4]);
        }
    }

    #[test]
    fn exact_lift_matches_thin_airfoil_theory_for_a_thin_section() {
        // nearly flat thin Joukowsky section: CL ~ 2 pi alpha
        let kt = KarmanTrefftz { x_center: 0.001, y_center: 0.0, te_angle: 0.0 };
        let s = kt.exact(4.0).unwrap();
        let want = 2.0 * PI * 4f64.to_radians();
        assert!((s.cl - want).abs() / want < 0.01, "{} {want}", s.cl);
    }

    /// Lift and quarter-chord moment by integrating the exact pressure over
    /// `n` contour points (an independent check of the Blasius integrals).
    fn pressure_integral(kt: &KarmanTrefftz, alpha: f64, n: usize) -> (f64, f64) {
        let c = kt.center();
        let r = kt.radius();
        let fr = kt.frame();
        let a = alpha.to_radians() - fr.rot.arg();
        let gamma = 4.0 * PI * r * (a - kt.theta_te()).sin();
        let t0 = kt.theta_te();
        let pts: Vec<(f64, f64, f64)> = (1..n)
            .map(|k| {
                let zeta = c + C::polar(r, t0 + 2.0 * PI * k as f64 / n as f64);
                let v = (kt.circle_velocity(a, gamma, zeta) / kt.dz_dzeta(zeta)).abs();
                let p = (kt.map(zeta) - fr.le) * fr.rot * (1.0 / fr.chord);
                (p.0, p.1, 1.0 - v * v)
            })
            .collect();
        let (mut fx, mut fy, mut cm) = (0.0, 0.0, 0.0);
        for w in pts.windows(2) {
            let ((x0, y0, p0), (x1, y1, p1)) = (w[0], w[1]);
            let pm = 0.5 * (p0 + p1);
            let (gx, gy) = (-pm * (y1 - y0), pm * (x1 - x0));
            fx += gx;
            fy += gy;
            cm += -(0.5 * (x0 + x1) - 0.25) * gy + 0.5 * (y0 + y1) * gx;
        }
        let al = alpha.to_radians();
        (fy * al.cos() - fx * al.sin(), cm)
    }

    #[test]
    fn blasius_moment_matches_a_converged_pressure_integral() {
        for kt in [KarmanTrefftz::default(), KarmanTrefftz { x_center: 0.001, y_center: 0.5, te_angle: 0.0 }, KarmanTrefftz { x_center: 0.05, y_center: -0.04, te_angle: 25.0 }] {
            for alpha in [-3.0, 4.0, 10.0] {
                let s = kt.exact(alpha).unwrap();
                let (cl, cm) = pressure_integral(&kt, alpha, 400_000);
                assert!((cl - s.cl).abs() < 2e-4 * s.cl.abs().max(1.0), "{kt:?} {alpha}: CL {cl} vs {}", s.cl);
                assert!((cm - s.cm).abs() < 2e-4, "{kt:?} {alpha}: Cm {cm} vs {}", s.cm);
            }
        }
        // a flat plate has no moment about its quarter chord
        let flat = KarmanTrefftz { x_center: 0.001, y_center: 0.0, te_angle: 0.0 };
        assert!(flat.exact(15.0).unwrap().cm.abs() < 1e-3);
        assert!(flat.exact(f64::NAN).is_err());
        // equivalent angles give the same answer
        let k = KarmanTrefftz::default();
        let (a, b) = (k.exact(-80.0).unwrap(), k.exact(280.0).unwrap());
        assert!((a.cl - b.cl).abs() < 1e-9 && (a.cm - b.cm).abs() < 1e-9);
        // 1e20 is exactly 280 (mod 360)
        assert!((k.exact(1e20).unwrap().cl - a.cl).abs() < 1e-9);
    }

    #[test]
    fn pressure_integrates_to_the_kutta_joukowsky_lift() {
        let kt = KarmanTrefftz::default();
        let alpha = 3.0;
        let s = kt.exact(alpha).unwrap();
        let (mut fx, mut fy) = (0.0, 0.0);
        for i in 0..s.cp.len() - 1 {
            let (x0, y0, p0) = s.cp[i];
            let (x1, y1, p1) = s.cp[i + 1];
            let pm = 0.5 * (p0 + p1);
            fx += -pm * (y1 - y0);
            fy += pm * (x1 - x0);
        }
        let a = alpha.to_radians();
        let lift = fy * a.cos() - fx * a.sin();
        assert!((lift - s.cl).abs() / s.cl < 0.02, "{lift} {}", s.cl);
    }
}
