//! Curves and surfaces of STEP B-reps (ISO 10303-42 parametrisations):
//! evaluation, first derivatives, periodicity, singular (pole) lines, and
//! point inversion (closest parameter) by Newton from a seed.

use super::{NurbsCurve, NurbsSurface};
use foiley_core::Vec3;
use std::f64::consts::{FRAC_PI_2, TAU};

/// Right-handed orthonormal frame (AXIS2_PLACEMENT_3D).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame {
    pub o: Vec3,
    pub x: Vec3,
    pub y: Vec3,
    pub z: Vec3,
}

impl Frame {
    /// From an origin, an axis and an optional reference direction (made
    /// orthogonal to the axis; an arbitrary perpendicular if missing).
    pub fn new(o: Vec3, z: Vec3, xref: Option<Vec3>) -> Frame {
        let z = z.normalized();
        let mut x = xref.map(|x| x - z * x.dot(z)).unwrap_or(Vec3::ZERO);
        if x.norm() < 1e-12 {
            let a = if z.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
            x = a - z * a.dot(z);
        }
        let x = x.normalized();
        Frame { o, x, y: z.cross(x), z }
    }
    fn local(&self, p: Vec3) -> Vec3 {
        let d = p - self.o;
        Vec3::new(d.dot(self.x), d.dot(self.y), d.dot(self.z))
    }
}

/// Rotation of `w` about the unit axis `a` by `ang` (Rodrigues).
fn rotate(w: Vec3, a: Vec3, ang: f64) -> Vec3 {
    let (s, c) = ang.sin_cos();
    w * c + a.cross(w) * s + a * (a.dot(w) * (1.0 - c))
}

/// Minimise `f` on [lo, hi]: `n` samples, then golden section around the best.
pub(crate) fn min1d(f: impl Fn(f64) -> f64, lo: f64, hi: f64, n: usize) -> f64 {
    let n = n.max(2);
    let mut best = (f64::INFINITY, lo);
    for i in 0..=n {
        let t = lo + (hi - lo) * i as f64 / n as f64;
        let v = f(t);
        if v < best.0 {
            best = (v, t);
        }
    }
    let h = (hi - lo) / n as f64;
    let (mut a, mut b) = ((best.1 - h).max(lo), (best.1 + h).min(hi));
    let g = 0.5 * (5f64.sqrt() - 1.0);
    let (mut c, mut d) = (b - g * (b - a), a + g * (b - a));
    let (mut fc, mut fd) = (f(c), f(d));
    for _ in 0..80 {
        if (b - a).abs() <= 1e-15 * (1.0 + a.abs() + b.abs()) {
            break;
        }
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - g * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + g * (b - a);
            fd = f(d);
        }
    }
    let t = 0.5 * (a + b);
    if f(t) <= best.0 {
        t
    } else {
        best.1
    }
}

/// A 3D (or, for parameter-space curves, 2D with z = 0) curve.
#[derive(Clone, Debug)]
pub enum Curve {
    /// C(t) = p + t d
    Line { p: Vec3, d: Vec3 },
    /// C(t) = o + r (cos t x + sin t y)
    Circle { f: Frame, r: f64 },
    /// C(t) = o + a cos t x + b sin t y
    Ellipse { f: Frame, a: f64, b: f64 },
    /// piecewise linear through the points, t in [0, n - 1]
    Polyline(Vec<Vec3>),
    BSpline(NurbsCurve),
}

impl Curve {
    pub fn kind(&self) -> &'static str {
        match self {
            Curve::Line { .. } => "LINE",
            Curve::Circle { .. } => "CIRCLE",
            Curve::Ellipse { .. } => "ELLIPSE",
            Curve::Polyline(_) => "POLYLINE",
            Curve::BSpline(_) => "B_SPLINE_CURVE",
        }
    }
    pub fn eval(&self, t: f64) -> Vec3 {
        match self {
            Curve::Line { p, d } => *p + *d * t,
            Curve::Circle { f, r } => {
                let (s, c) = t.sin_cos();
                f.o + (f.x * c + f.y * s) * *r
            }
            Curve::Ellipse { f, a, b } => {
                let (s, c) = t.sin_cos();
                f.o + f.x * (a * c) + f.y * (b * s)
            }
            Curve::Polyline(p) => {
                let n = p.len();
                if n == 1 {
                    return p[0];
                }
                let t = t.clamp(0.0, (n - 1) as f64);
                let i = (t.floor() as usize).min(n - 2);
                p[i].lerp(p[i + 1], t - i as f64)
            }
            Curve::BSpline(c) => c.eval(t),
        }
    }
    pub fn deriv(&self, t: f64) -> Vec3 {
        match self {
            Curve::Line { d, .. } => *d,
            Curve::Circle { f, r } => {
                let (s, c) = t.sin_cos();
                (f.y * c - f.x * s) * *r
            }
            Curve::Ellipse { f, a, b } => {
                let (s, c) = t.sin_cos();
                f.y * (b * c) - f.x * (a * s)
            }
            Curve::Polyline(p) => {
                let n = p.len();
                if n < 2 {
                    return Vec3::ZERO;
                }
                let i = (t.clamp(0.0, (n - 1) as f64).floor() as usize).min(n - 2);
                p[i + 1] - p[i]
            }
            Curve::BSpline(c) => c.eval_deriv(t).1,
        }
    }
    /// Parameter period of a closed periodic curve.
    pub fn period(&self) -> Option<f64> {
        match self {
            Curve::Circle { .. } | Curve::Ellipse { .. } => Some(TAU),
            _ => None,
        }
    }
    /// Parameter range of a bounded curve.
    pub fn range(&self) -> Option<(f64, f64)> {
        match self {
            Curve::Polyline(p) => Some((0.0, p.len().saturating_sub(1) as f64)),
            Curve::BSpline(c) => Some(c.domain()),
            _ => None,
        }
    }
    /// Natural number of pieces (for the initial split when sampling).
    pub fn pieces(&self) -> usize {
        match self {
            Curve::Polyline(p) => p.len().saturating_sub(1).max(1),
            Curve::BSpline(c) => c.spans(),
            _ => 1,
        }
    }
    /// Parameter of the closest point to `p` (global search).
    pub fn param_of(&self, p: Vec3) -> f64 {
        match self {
            Curve::Line { p: o, d } => (p - *o).dot(*d) / d.norm2().max(1e-300),
            Curve::Circle { f, .. } => {
                let l = f.local(p);
                l.y.atan2(l.x)
            }
            Curve::Ellipse { f, a, b } => {
                let l = f.local(p);
                let mut t = (l.y / b).atan2(l.x / a);
                // Newton on (C(t) - p) . C'(t) = 0
                for _ in 0..30 {
                    let (s, c) = t.sin_cos();
                    let (dx, dy) = (a * c - l.x, b * s - l.y);
                    let g = -dx * a * s + dy * b * c;
                    let h = a * a * s * s + b * b * c * c - dx * a * c - dy * b * s;
                    if h.abs() < 1e-300 {
                        break;
                    }
                    let step = g / h;
                    t -= step;
                    if step.abs() < 1e-15 {
                        break;
                    }
                }
                t
            }
            Curve::Polyline(pts) => {
                let mut best = (f64::INFINITY, 0.0);
                for i in 0..pts.len().saturating_sub(1) {
                    let (a, b) = (pts[i], pts[i + 1]);
                    let ab = b - a;
                    let s = ((p - a).dot(ab) / ab.norm2().max(1e-300)).clamp(0.0, 1.0);
                    let d = (a + ab * s - p).norm2();
                    if d < best.0 {
                        best = (d, i as f64 + s);
                    }
                }
                best.1
            }
            Curve::BSpline(c) => {
                let (a, b) = c.domain();
                let mut t = min1d(|t| (c.eval(t) - p).norm2(), a, b, 8 * c.spans());
                for _ in 0..30 {
                    let (q, d) = c.eval_deriv(t);
                    let step = (q - p).dot(d) / d.norm2().max(1e-300);
                    let nt = (t - step).clamp(a, b);
                    if (c.eval(nt) - p).norm2() > (q - p).norm2() {
                        break;
                    }
                    let done = (nt - t).abs() <= 1e-15 * (1.0 + t.abs());
                    t = nt;
                    if done {
                        break;
                    }
                }
                t
            }
        }
    }
    /// Minimise `f(t)` over the curve parameter: sampled over the period or
    /// range, or around `seeds` for unbounded curves.
    fn search(&self, f: impl Fn(f64) -> f64, seeds: &[f64], span: f64) -> f64 {
        if let Some(per) = self.period() {
            return min1d(&f, 0.0, per, 48);
        }
        if let Some((a, b)) = self.range() {
            return min1d(&f, a, b, 8 * self.pieces());
        }
        let mut best = (f64::INFINITY, 0.0);
        for &s in seeds {
            let t = min1d(&f, s - span, s + span, 64);
            let v = f(t);
            if v < best.0 {
                best = (v, t);
            }
        }
        best.1
    }
}

/// A STEP surface.
#[derive(Clone, Debug)]
pub enum Surface {
    /// S = o + u x + v y
    Plane(Frame),
    /// S = o + r (cos u x + sin u y) + v z
    Cylinder { f: Frame, r: f64 },
    /// S = o + (r + v tan a)(cos u x + sin u y) + v z
    Cone { f: Frame, r: f64, tan_a: f64 },
    /// S = o + r cos v (cos u x + sin u y) + r sin v z
    Sphere { f: Frame, r: f64 },
    /// S = o + (big + small cos v)(cos u x + sin u y) + small sin v z
    Torus { f: Frame, big: f64, small: f64 },
    /// the curve C(v) rotated by u about the axis (o, a)
    Revolution { curve: Box<Curve>, o: Vec3, a: Vec3 },
    /// S = C(u) + v d
    Extrusion { curve: Box<Curve>, d: Vec3 },
    /// `closed[k]`: the surface closes on itself in that parameter (first
    /// and last control rows coincide), which is then treated as periodic
    BSpline { s: NurbsSurface, closed: [bool; 2] },
}

impl Surface {
    pub fn kind(&self) -> &'static str {
        match self {
            Surface::Plane(_) => "PLANE",
            Surface::Cylinder { .. } => "CYLINDRICAL_SURFACE",
            Surface::Cone { .. } => "CONICAL_SURFACE",
            Surface::Sphere { .. } => "SPHERICAL_SURFACE",
            Surface::Torus { .. } => "TOROIDAL_SURFACE",
            Surface::Revolution { .. } => "SURFACE_OF_REVOLUTION",
            Surface::Extrusion { .. } => "SURFACE_OF_LINEAR_EXTRUSION",
            Surface::BSpline { .. } => "B_SPLINE_SURFACE",
        }
    }

    pub fn eval(&self, u: f64, v: f64) -> Vec3 {
        match self {
            Surface::Plane(f) => f.o + f.x * u + f.y * v,
            Surface::Cylinder { f, r } => {
                let (s, c) = u.sin_cos();
                f.o + (f.x * c + f.y * s) * *r + f.z * v
            }
            Surface::Cone { f, r, tan_a } => {
                let (s, c) = u.sin_cos();
                f.o + (f.x * c + f.y * s) * (r + v * tan_a) + f.z * v
            }
            Surface::Sphere { f, r } => {
                let (s, c) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                f.o + (f.x * c + f.y * s) * (r * cv) + f.z * (r * sv)
            }
            Surface::Torus { f, big, small } => {
                let (s, c) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                f.o + (f.x * c + f.y * s) * (big + small * cv) + f.z * (small * sv)
            }
            Surface::Revolution { curve, o, a } => *o + rotate(curve.eval(v) - *o, *a, u),
            Surface::Extrusion { curve, d } => curve.eval(u) + *d * v,
            Surface::BSpline { s, closed } => {
                let (u, v) = Self::wrap_closed(s, closed, u, v);
                s.eval(u, v)
            }
        }
    }

    /// (S, dS/du, dS/dv)
    pub fn derivs(&self, u: f64, v: f64) -> (Vec3, Vec3, Vec3) {
        match self {
            Surface::Plane(f) => (self.eval(u, v), f.x, f.y),
            Surface::Cylinder { f, r } => {
                let (s, c) = u.sin_cos();
                (self.eval(u, v), (f.y * c - f.x * s) * *r, f.z)
            }
            Surface::Cone { f, r, tan_a } => {
                let (s, c) = u.sin_cos();
                let rad = f.x * c + f.y * s;
                (self.eval(u, v), (f.y * c - f.x * s) * (r + v * tan_a), rad * *tan_a + f.z)
            }
            Surface::Sphere { f, r } => {
                let (s, c) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                let rad = f.x * c + f.y * s;
                (self.eval(u, v), (f.y * c - f.x * s) * (r * cv), rad * (-r * sv) + f.z * (r * cv))
            }
            Surface::Torus { f, big, small } => {
                let (s, c) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                let rad = f.x * c + f.y * s;
                (self.eval(u, v), (f.y * c - f.x * s) * (big + small * cv), rad * (-small * sv) + f.z * (small * cv))
            }
            Surface::Revolution { curve, o, a } => {
                let p = self.eval(u, v);
                (p, a.cross(p - *o), rotate(curve.deriv(v), *a, u))
            }
            Surface::Extrusion { curve, d } => (self.eval(u, v), curve.deriv(u), *d),
            Surface::BSpline { s, closed } => {
                let (u, v) = Self::wrap_closed(s, closed, u, v);
                s.eval_derivs(u, v)
            }
        }
    }

    /// Parameter periods (u, v).
    pub fn periods(&self) -> [Option<f64>; 2] {
        match self {
            Surface::Plane(_) => [None, None],
            Surface::BSpline { s, closed } => {
                let ((u0, u1), (v0, v1)) = s.domain();
                [closed[0].then_some(u1 - u0), closed[1].then_some(v1 - v0)]
            }
            Surface::Cylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } => [Some(TAU), None],
            Surface::Torus { .. } => [Some(TAU), Some(TAU)],
            Surface::Revolution { curve, .. } => [Some(TAU), curve.period()],
            Surface::Extrusion { curve, .. } => [curve.period(), None],
        }
    }

    /// Parameter value of the natural seam of a periodic direction.
    pub fn seam_origin(&self, dir: usize) -> f64 {
        match self {
            Surface::BSpline { s, .. } => {
                let ((u0, _), (v0, _)) = s.domain();
                if dir == 0 {
                    u0
                } else {
                    v0
                }
            }
            _ => 0.0,
        }
    }

    /// Finite parameter bounds of non-periodic directions (for clamping).
    pub fn bounds(&self) -> [Option<(f64, f64)>; 2] {
        match self {
            Surface::Sphere { .. } => [None, Some((-FRAC_PI_2, FRAC_PI_2))],
            Surface::Revolution { curve, .. } => [None, curve.range()],
            Surface::Extrusion { curve, .. } => [curve.range(), None],
            Surface::BSpline { s, closed } => {
                let (a, b) = s.domain();
                [(!closed[0]).then_some(a), (!closed[1]).then_some(b)]
            }
            _ => [None, None],
        }
    }

    /// The whole parameter domain, when it is bounded (a face without
    /// bounds covers it).
    pub fn natural_domain(&self) -> Option<[(f64, f64); 2]> {
        if let Surface::BSpline { s, .. } = self {
            let (a, b) = s.domain();
            return Some([a, b]);
        }
        let p = self.periods();
        let b = self.bounds();
        let d = |i: usize| p[i].map(|per| (0.0, per)).or(b[i]);
        Some([d(0)?, d(1)?])
    }

    /// A B-spline surface from its net, detecting closure in u and v.
    pub fn bspline(s: NurbsSurface) -> Surface {
        let nu = s.ctrl.len();
        let nv = s.ctrl.first().map(|r| r.len()).unwrap_or(0);
        let (mut lo, mut hi) = (Vec3::new(f64::MAX, f64::MAX, f64::MAX), Vec3::new(f64::MIN, f64::MIN, f64::MIN));
        for p in s.ctrl.iter().flatten() {
            lo = Vec3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Vec3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        let tol = 1e-10 * (hi - lo).norm().max(1e-300);
        let closed_u = nu > 2 && (0..nv).all(|j| (s.ctrl[0][j] - s.ctrl[nu - 1][j]).norm() <= tol && (s.w[0][j] - s.w[nu - 1][j]).abs() <= 1e-12);
        let closed_v = nv > 2 && (0..nu).all(|i| (s.ctrl[i][0] - s.ctrl[i][nv - 1]).norm() <= tol && (s.w[i][0] - s.w[i][nv - 1]).abs() <= 1e-12);
        Surface::BSpline { s, closed: [closed_u, closed_v] }
    }

    fn wrap_closed(s: &NurbsSurface, closed: &[bool; 2], u: f64, v: f64) -> (f64, f64) {
        if !closed[0] && !closed[1] {
            return (u, v);
        }
        let ((u0, u1), (v0, v1)) = s.domain();
        let w = |x: f64, a: f64, b: f64| if x < a || x > b { a + (x - a).rem_euclid(b - a) } else { x };
        (if closed[0] { w(u, u0, u1) } else { u }, if closed[1] { w(v, v0, v1) } else { v })
    }

    /// Parameter values v where the whole u line collapses to a point
    /// (cone apex, sphere poles), for closing a face around a pole.
    pub fn poles_v(&self) -> Vec<f64> {
        self.poles(1)
    }

    /// Values of parameter `dir` (0 = u, 1 = v) at which the iso-line of
    /// the other parameter collapses to a point.
    pub fn poles(&self, dir: usize) -> Vec<f64> {
        if let Surface::BSpline { s, .. } = self {
            let nu = s.ctrl.len();
            let nv = s.ctrl.first().map(|r| r.len()).unwrap_or(0);
            let ((u0, u1), (v0, v1)) = s.domain();
            let size = s.ctrl.iter().flatten().fold(0.0f64, |m, p| m.max(p.norm())).max(1e-300);
            let collapsed = |pts: Vec<Vec3>| pts.iter().all(|p| (*p - pts[0]).norm() <= 1e-10 * size);
            let mut out = Vec::new();
            if dir == 0 {
                if collapsed(s.ctrl[0].clone()) {
                    out.push(u0);
                }
                if collapsed(s.ctrl[nu - 1].clone()) {
                    out.push(u1);
                }
            } else {
                if collapsed((0..nu).map(|i| s.ctrl[i][0]).collect()) {
                    out.push(v0);
                }
                if collapsed((0..nu).map(|i| s.ctrl[i][nv - 1]).collect()) {
                    out.push(v1);
                }
            }
            return out;
        }
        if dir == 0 {
            return vec![];
        }
        match self {
            Surface::Cone { r, tan_a, .. } if tan_a.abs() > 1e-300 => vec![-r / tan_a],
            Surface::Sphere { .. } => vec![-FRAC_PI_2, FRAC_PI_2],
            Surface::Revolution { curve, o, a } => {
                let mut out = Vec::new();
                if let Some((t0, t1)) = curve.range() {
                    for t in [t0, t1] {
                        let q = curve.eval(t) - *o;
                        if (q - *a * q.dot(*a)).norm() < 1e-9 * (1.0 + q.norm()) {
                            out.push(t);
                        }
                    }
                }
                out
            }
            _ => vec![],
        }
    }

    /// Is the u (resp. v) parameter undefined at this point (|dS/du| = 0)?
    pub fn singular(&self, u: f64, v: f64) -> [bool; 2] {
        let (_, su, sv) = self.derivs(u, v);
        let (a, b) = (su.norm(), sv.norm());
        let m = a.max(b).max(1e-300);
        [a < 1e-8 * m, b < 1e-8 * m]
    }

    /// Parameters of the point of the surface closest to `p`: a closed form
    /// or a coarse-grid seed, then Newton.  `seed` (a nearby solution, e.g.
    /// the previous point of a loop) is used where the inversion is not
    /// unique (B-splines).
    pub fn invert(&self, p: Vec3, seed: Option<[f64; 2]>) -> [f64; 2] {
        let uv = match self {
            Surface::Plane(f) => {
                let l = f.local(p);
                return [l.x, l.y];
            }
            Surface::Cylinder { f, .. } => {
                let l = f.local(p);
                return [l.y.atan2(l.x), l.z];
            }
            Surface::Cone { f, .. } => {
                let l = f.local(p);
                [l.y.atan2(l.x), l.z]
            }
            Surface::Sphere { f, .. } => {
                let l = f.local(p);
                return [l.y.atan2(l.x), l.z.atan2(l.x.hypot(l.y))];
            }
            Surface::Torus { f, big, .. } => {
                let l = f.local(p);
                return [l.y.atan2(l.x), l.z.atan2(l.x.hypot(l.y) - big)];
            }
            Surface::Revolution { curve, o, a } => {
                let d = p - *o;
                let z = d.dot(*a);
                let rp = d - *a * z;
                let r = rp.norm();
                let meridian = |t: f64| {
                    let q = curve.eval(t) - *o;
                    let zc = q.dot(*a);
                    let rc = (q - *a * zc).norm();
                    (rc - r).powi(2) + (zc - z).powi(2)
                };
                let span = 2.0 * (d.norm() + (curve.eval(0.0) - *o).norm()) / curve.deriv(0.0).norm().max(1e-300);
                let seeds = [curve.param_of(p), curve.param_of(*o + *a * z), seed.map(|s| s[1]).unwrap_or(0.0)];
                let t = curve.search(meridian, &seeds, span);
                let q = curve.eval(t) - *o;
                let qr = q - *a * q.dot(*a);
                let u = if qr.norm() > 1e-12 * (1.0 + q.norm()) && r > 1e-12 * (1.0 + d.norm()) {
                    let e1 = qr.normalized();
                    let e2 = a.cross(e1);
                    rp.dot(e2).atan2(rp.dot(e1))
                } else {
                    seed.map(|s| s[0]).unwrap_or(0.0)
                };
                [u, t]
            }
            Surface::Extrusion { curve, d } => {
                let dn = d.normalized();
                let perp = |t: f64| {
                    let w = p - curve.eval(t);
                    (w - dn * w.dot(dn)).norm2()
                };
                let span = 2.0 * ((p - curve.eval(0.0)).norm()) / curve.deriv(0.0).norm().max(1e-300) + 1.0;
                let seeds = [curve.param_of(p), seed.map(|s| s[0]).unwrap_or(0.0)];
                let t = curve.search(perp, &seeds, span);
                [t, (p - curve.eval(t)).dot(*d) / d.norm2().max(1e-300)]
            }
            Surface::BSpline { s, .. } => match seed {
                Some(sd) => sd,
                None => {
                    let ((u0, u1), (v0, v1)) = s.domain();
                    let (nu, nv) = ((s.ctrl.len() * 3).clamp(8, 64), (s.ctrl.first().map(|r| r.len()).unwrap_or(2) * 3).clamp(8, 64));
                    let mut best = (f64::INFINITY, [u0, v0]);
                    for i in 0..=nu {
                        for j in 0..=nv {
                            let uv = [u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64];
                            let d = (s.eval(uv[0], uv[1]) - p).norm2();
                            if d < best.0 {
                                best = (d, uv);
                            }
                        }
                    }
                    best.1
                }
            },
        };
        self.newton(p, uv)
    }

    /// Gauss-Newton on |S(u, v) - p|^2 with backtracking, clamped to the
    /// bounds of non-periodic directions.  Handles singular points (one
    /// derivative zero) by stepping in the other parameter only.
    pub fn newton(&self, p: Vec3, mut uv: [f64; 2]) -> [f64; 2] {
        let bounds = self.bounds();
        let clamp = |mut w: [f64; 2]| {
            for k in 0..2 {
                if let Some((a, b)) = bounds[k] {
                    w[k] = w[k].clamp(a, b);
                }
            }
            w
        };
        uv = clamp(uv);
        let mut r2 = (self.eval(uv[0], uv[1]) - p).norm2();
        for _ in 0..60 {
            if r2 == 0.0 {
                break;
            }
            let (s, su, sv) = self.derivs(uv[0], uv[1]);
            let r = p - s;
            let (a11, a12, a22) = (su.norm2(), su.dot(sv), sv.norm2());
            let (b1, b2) = (su.dot(r), sv.dot(r));
            let det = a11 * a22 - a12 * a12;
            let (du, dv) = if det > 1e-12 * a11 * a22 && det > 0.0 {
                ((a22 * b1 - a12 * b2) / det, (a11 * b2 - a12 * b1) / det)
            } else if a11 >= a22 && a11 > 0.0 {
                (b1 / a11, 0.0)
            } else if a22 > 0.0 {
                (0.0, b2 / a22)
            } else {
                break;
            };
            let mut lam = 1.0;
            let mut moved = false;
            for _ in 0..30 {
                let cand = clamp([uv[0] + lam * du, uv[1] + lam * dv]);
                let rc = (self.eval(cand[0], cand[1]) - p).norm2();
                if rc < r2 || (rc == r2 && cand != uv) {
                    let small = (cand[0] - uv[0]).abs() <= 1e-15 * (1.0 + uv[0].abs()) && (cand[1] - uv[1]).abs() <= 1e-15 * (1.0 + uv[1].abs());
                    uv = cand;
                    r2 = rc;
                    moved = !small;
                    break;
                }
                lam *= 0.5;
            }
            if !moved {
                break;
            }
        }
        uv
    }
}

/// Wrap `x` into (-p/2, p/2].
pub(crate) fn wrap(x: f64, p: f64) -> f64 {
    x - p * (x / p).round()
}
