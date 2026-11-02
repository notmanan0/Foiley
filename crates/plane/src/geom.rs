//! Plane geometry: wings defined by sections, fuselages defined by frames,
//! point masses; reference dimensions and inertia.
//!
//! Axes: x aft, y to the right (starboard), z up.  Angles in degrees.

use crate::foildb::FoilDb;
use foiley_core::interp::{distribution, CubicSpline};
use foiley_core::{Spacing, Vec3, DEG};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum WingKind {
    #[default]
    Main,
    Elevator,
    Fin,
    Other,
}

impl WingKind {
    pub fn label(self) -> &'static str {
        match self {
            WingKind::Main => "Main wing",
            WingKind::Elevator => "Elevator",
            WingKind::Fin => "Fin",
            WingKind::Other => "Other wing",
        }
    }
    pub const ALL: [WingKind; 4] = [WingKind::Main, WingKind::Elevator, WingKind::Fin, WingKind::Other];
}

/// One wing section.  The panel between section i and i+1 uses the
/// dihedral, panel counts and distributions of section i.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct WingSection {
    /// span position measured along the developed span (ignoring dihedral), m
    pub y: f64,
    pub chord: f64,
    /// leading-edge x offset relative to the wing position, m
    pub offset: f64,
    /// dihedral of the panel outboard of this section, degrees
    pub dihedral: f64,
    /// twist of this section about its quarter chord, degrees (+ = nose up)
    pub twist: f64,
    pub foil_left: String,
    pub foil_right: String,
    /// chordwise panels (outboard panel)
    pub nx: usize,
    /// spanwise panels (outboard panel)
    pub ny: usize,
    pub xdist: Spacing,
    pub ydist: Spacing,
}

impl WingSection {
    pub fn new(y: f64, chord: f64, offset: f64, foil: &str) -> Self {
        WingSection { y, chord, offset, dihedral: 0.0, twist: 0.0, foil_left: foil.into(), foil_right: foil.into(), nx: 13, ny: 11, xdist: Spacing::Cosine, ydist: Spacing::Uniform }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PointMass {
    pub tag: String,
    pub mass: f64,
    pub position: Vec3,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Wing {
    pub name: String,
    pub kind: WingKind,
    pub sections: Vec<WingSection>,
    /// position of the root section's reference point in plane axes, m
    pub position: Vec3,
    /// rotation about x, degrees (90 makes a one-sided wing point up: a fin)
    pub rx: f64,
    /// tilt (incidence) about y, degrees (+ = nose up)
    pub ry: f64,
    /// build both halves (mirror about the wing's xz plane)
    pub two_sided: bool,
    /// left foils follow right foils
    pub symmetric: bool,
    /// close the root side with a panel patch (thick panel method)
    pub close_inner: bool,
    /// structural mass, kg (distributed over the wing volume)
    pub mass: f64,
    #[serde(default)]
    pub point_masses: Vec<PointMass>,
}

/// A section placed in space: quarter-chord point, chord vector (LE->TE,
/// length = chord), up vector (length = chord).
#[derive(Clone, Copy, Debug)]
pub struct SectionFrame {
    pub q4: Vec3,
    pub cx: Vec3,
    pub cu: Vec3,
    pub chord: f64,
}

impl SectionFrame {
    /// Point at chord fraction `xi` and height `eta` (fraction of chord).
    #[inline]
    pub fn point(&self, xi: f64, eta: f64) -> Vec3 {
        self.q4 + self.cx * (xi - 0.25) + self.cu * eta
    }
    pub fn le(&self) -> Vec3 {
        self.point(0.0, 0.0)
    }
    pub fn te(&self) -> Vec3 {
        self.point(1.0, 0.0)
    }
}

impl Wing {
    pub fn new(name: &str, kind: WingKind, sections: Vec<WingSection>) -> Self {
        Wing { name: name.into(), kind, sections, position: Vec3::ZERO, rx: 0.0, ry: 0.0, two_sided: true, symmetric: true, close_inner: false, mass: 0.0, point_masses: vec![] }
    }

    /// Basic validity checks with plain-language messages.
    pub fn validate(&self) -> Result<(), String> {
        let n = self.sections.len();
        if n < 2 {
            return Err(format!("wing '{}' needs at least 2 sections", self.name));
        }
        for (i, s) in self.sections.iter().enumerate() {
            if !(s.chord > 0.0) || !s.chord.is_finite() {
                return Err(format!("wing '{}' section {}: chord must be > 0", self.name, i + 1));
            }
            if s.nx == 0 || s.ny == 0 || s.nx > 200 || s.ny > 500 {
                return Err(format!("wing '{}' section {}: panel counts must be 1..200 (x) and 1..500 (y)", self.name, i + 1));
            }
            if i > 0 && s.y < self.sections[i - 1].y - 1e-12 {
                return Err(format!("wing '{}': section {} span position {} is less than the previous one", self.name, i + 1, s.y));
            }
        }
        if self.sections[n - 1].y - self.sections[0].y < 1e-6 {
            return Err(format!("wing '{}' has zero span", self.name));
        }
        Ok(())
    }

    pub fn foil(&self, i: usize, left: bool) -> &str {
        let s = &self.sections[i];
        if left && !self.symmetric {
            &s.foil_left
        } else {
            &s.foil_right
        }
    }

    /// Junction "up" directions (bisectors of adjacent panel normals).
    fn junction_ups(&self) -> Vec<Vec3> {
        let n = self.sections.len();
        let pn: Vec<Vec3> = (0..n - 1)
            .map(|i| {
                let d = self.sections[i].dihedral * DEG;
                Vec3::new(0.0, -d.sin(), d.cos())
            })
            .collect();
        let mut ups = Vec::with_capacity(n);
        for i in 0..n {
            let u = if i == 0 {
                if self.two_sided && self.sections[0].y.abs() < 1e-9 {
                    Vec3::Z
                } else {
                    pn[0]
                }
            } else if i == n - 1 {
                pn[n - 2]
            } else {
                (pn[i - 1] + pn[i]).normalized()
            };
            ups.push(u);
        }
        ups
    }

    /// Section frames of the right half in wing-local axes (before rx/ry/position).
    pub fn local_frames(&self) -> Vec<SectionFrame> {
        let n = self.sections.len();
        let ups = self.junction_ups();
        let mut le = Vec3::new(self.sections[0].offset, self.sections[0].y, 0.0);
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let s = &self.sections[i];
            if i > 0 {
                let p = &self.sections[i - 1];
                let d = p.dihedral * DEG;
                let l = s.y - p.y;
                le = Vec3::new(s.offset, le.y + l * d.cos(), le.z + l * d.sin());
            }
            let u = ups[i];
            let span_axis = u.cross(Vec3::X).normalized();
            let tw = s.twist * DEG;
            let x = Vec3::X.rot_axis(span_axis, tw);
            let up = u.rot_axis(span_axis, tw);
            let q4 = le + Vec3::X * (0.25 * s.chord);
            out.push(SectionFrame { q4, cx: x * s.chord, cu: up * s.chord, chord: s.chord });
        }
        out
    }

    /// Wing-local -> plane axes.
    #[inline]
    pub fn to_plane(&self, p: Vec3) -> Vec3 {
        self.position + p.rot_x(self.rx * DEG).rot_y(self.ry * DEG)
    }
    #[inline]
    pub fn dir_to_plane(&self, d: Vec3) -> Vec3 {
        d.rot_x(self.rx * DEG).rot_y(self.ry * DEG)
    }

    /// Section frames in plane axes, for one side (`left` mirrors y).
    pub fn frames(&self, left: bool) -> Vec<SectionFrame> {
        self.local_frames()
            .into_iter()
            .map(|f| {
                let m = |v: Vec3| if left { Vec3::new(v.x, -v.y, v.z) } else { v };
                SectionFrame { q4: self.to_plane(m(f.q4)), cx: self.dir_to_plane(m(f.cx)), cu: self.dir_to_plane(m(f.cu)), chord: f.chord }
            })
            .collect()
    }

    /// Half-wing planform quantities: (area, projected area, projected span, integral of c^2).
    fn half_sums(&self) -> (f64, f64, f64, f64) {
        let (mut s, mut sp, mut bp, mut c2) = (0.0, 0.0, 0.0, 0.0);
        for w in self.sections.windows(2) {
            let l = w[1].y - w[0].y;
            let d = w[0].dihedral * DEG;
            let a = l * 0.5 * (w[0].chord + w[1].chord);
            s += a;
            sp += a * d.cos();
            bp += l * d.cos();
            let (c1, c2v) = (w[0].chord, w[1].chord);
            c2 += l * (c1 * c1 + c1 * c2v + c2v * c2v) / 3.0;
        }
        (s, sp, bp, c2)
    }

    pub fn planform_area(&self) -> f64 {
        let k = if self.two_sided { 2.0 } else { 1.0 };
        k * self.half_sums().0
    }
    pub fn projected_area(&self) -> f64 {
        let k = if self.two_sided { 2.0 } else { 1.0 };
        k * self.half_sums().1
    }
    pub fn planform_span(&self) -> f64 {
        let tip = self.sections.last().map(|s| s.y).unwrap_or(0.0);
        if self.two_sided {
            2.0 * tip
        } else {
            tip - self.sections[0].y
        }
    }
    pub fn projected_span(&self) -> f64 {
        let (_, _, bp, _) = self.half_sums();
        let root = if self.two_sided { self.sections[0].y } else { 0.0 };
        if self.two_sided {
            2.0 * (bp + root)
        } else {
            bp
        }
    }
    /// Mean aerodynamic chord.
    pub fn mac(&self) -> f64 {
        let (s, _, _, c2) = self.half_sums();
        if s > 0.0 {
            c2 / s
        } else {
            0.0
        }
    }
    pub fn aspect_ratio(&self) -> f64 {
        let a = self.planform_area();
        if a > 0.0 {
            self.planform_span().powi(2) / a
        } else {
            0.0
        }
    }
    pub fn taper(&self) -> f64 {
        let r = self.sections.first().map(|s| s.chord).unwrap_or(0.0);
        let t = self.sections.last().map(|s| s.chord).unwrap_or(0.0);
        if r > 1e-9 {
            t / r
        } else {
            0.0
        }
    }
    /// Average quarter-chord sweep, degrees.
    pub fn sweep(&self) -> f64 {
        let (r, t) = (&self.sections[0], self.sections.last().unwrap());
        let dx = (t.offset + 0.25 * t.chord) - (r.offset + 0.25 * r.chord);
        let dy = t.y - r.y;
        dx.atan2(dy) / DEG
    }
    /// Spanwise y position (wing-local, developed) of the MAC, and its LE x.
    pub fn mac_position(&self) -> (f64, f64) {
        let (s, ..) = self.half_sums();
        let mut y_int = 0.0;
        let mut x_int = 0.0;
        for w in self.sections.windows(2) {
            // integrate c*y and c*x_le over the panel with 20-point midpoint rule
            let n = 20;
            for k in 0..n {
                let t = (k as f64 + 0.5) / n as f64;
                let y = w[0].y + t * (w[1].y - w[0].y);
                let c = w[0].chord + t * (w[1].chord - w[0].chord);
                let x = w[0].offset + t * (w[1].offset - w[0].offset);
                let dy = (w[1].y - w[0].y) / n as f64;
                y_int += c * y * dy;
                x_int += c * x * dy;
            }
        }
        if s > 0.0 {
            (y_int / s, x_int / s)
        } else {
            (0.0, 0.0)
        }
    }

    /// Insert a section halfway between section `i` and `i+1`.
    pub fn split_section(&mut self, i: usize) {
        if i + 1 >= self.sections.len() {
            return;
        }
        let (a, b) = (&self.sections[i], &self.sections[i + 1]);
        let mut s = a.clone();
        s.y = 0.5 * (a.y + b.y);
        s.chord = 0.5 * (a.chord + b.chord);
        s.offset = 0.5 * (a.offset + b.offset);
        s.twist = 0.5 * (a.twist + b.twist);
        s.ny = (a.ny / 2).max(1);
        self.sections[i].ny = (a.ny - s.ny).max(1);
        self.sections.insert(i + 1, s);
    }

    /// Default main wing (flow5 default dimensions).
    pub fn default_main(foil: &str) -> Wing {
        let mut a = WingSection::new(0.0, 0.30, 0.0, foil);
        a.nx = 13;
        a.ny = 19;
        a.ydist = Spacing::Sine;
        let mut b = WingSection::new(1.5, 0.19, 0.085, foil);
        b.nx = 13;
        b.ny = 19;
        let mut w = Wing::new("Main wing", WingKind::Main, vec![a, b]);
        w.position = Vec3::new(0.0, 0.0, 0.0);
        w.mass = 0.8;
        w
    }
    pub fn default_elevator(foil: &str) -> Wing {
        let mut a = WingSection::new(0.0, 0.20, 0.0, foil);
        a.nx = 7;
        a.ny = 7;
        let mut b = WingSection::new(0.32, 0.11, 0.05, foil);
        b.nx = 7;
        b.ny = 7;
        let mut w = Wing::new("Elevator", WingKind::Elevator, vec![a, b]);
        w.position = Vec3::new(0.95, 0.0, 0.025);
        w.ry = -1.5;
        w.mass = 0.1;
        w
    }
    pub fn default_fin(foil: &str) -> Wing {
        let mut a = WingSection::new(0.0, 0.20, 0.0, foil);
        a.nx = 7;
        a.ny = 6;
        let mut b = WingSection::new(0.29, 0.12, 0.08, foil);
        b.nx = 7;
        b.ny = 6;
        let mut w = Wing::new("Fin", WingKind::Fin, vec![a, b]);
        w.position = Vec3::new(0.95, 0.0, 0.05);
        w.rx = 90.0;
        w.two_sided = false;
        w.close_inner = true;
        w.mass = 0.05;
        w
    }
}

/// Fuselage frame: an elliptic / super-elliptic cross section.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct FuseFrame {
    /// x relative to the fuse position
    pub x: f64,
    /// centre height relative to the fuse position
    pub z: f64,
    pub half_width: f64,
    pub half_height: f64,
}

/// A frame of NURBS control points (right half, y >= 0, ordered from the
/// top of the section to the bottom), relative to the fuselage position.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NurbsFrame {
    pub x: f64,
    pub points: Vec<Vec3>,
}

/// xflr5 / flow5 style NURBS fuselage: a tensor-product B-spline surface
/// through the frames' control points, mirrored about y = 0.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NurbsBody {
    pub frames: Vec<NurbsFrame>,
    /// degree along the body (x)
    pub degree_u: usize,
    /// degree around the hoop
    pub degree_v: usize,
}

impl NurbsBody {
    pub fn surface(&self, left: bool) -> crate::cad::NurbsSurface {
        let nu = self.frames.len();
        let nv = self.frames.iter().map(|f| f.points.len()).min().unwrap_or(0);
        let pu = self.degree_u.clamp(1, nu.saturating_sub(1).max(1));
        let pv = self.degree_v.clamp(1, nv.saturating_sub(1).max(1));
        let ctrl: Vec<Vec<Vec3>> = self.frames.iter().map(|f| f.points[..nv].iter().map(|p| if left { Vec3::new(p.x, -p.y, p.z) } else { *p }).collect()).collect();
        crate::cad::NurbsSurface { pu, pv, ku: crate::cad::clamped_uniform_knots(nu, pu), kv: crate::cad::clamped_uniform_knots(nv, pv), w: vec![vec![1.0; nv]; nu], ctrl }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Fuse {
    pub name: String,
    pub position: Vec3,
    /// elliptic / super-elliptic frames (used when `nurbs` is None and `triangles` is empty)
    pub frames: Vec<FuseFrame>,
    /// cross-section exponent (2 = ellipse, larger = boxier)
    pub exponent: f64,
    /// panels along x
    pub nx: usize,
    /// panels around each half (top to bottom)
    pub nh: usize,
    /// structural mass, kg (distributed over the skin)
    pub mass: f64,
    /// NURBS definition (xflr5 / flow5 bodies)
    #[serde(default)]
    pub nurbs: Option<NurbsBody>,
    /// triangulated surface from STL / STEP, relative to `position`
    #[serde(default)]
    pub triangles: Vec<[Vec3; 3]>,
}

impl Fuse {
    pub fn default_pod() -> Fuse {
        let f = |x: f64, w: f64, h: f64, z: f64| FuseFrame { x, z, half_width: w, half_height: h };
        Fuse {
            name: "Fuselage".into(),
            position: Vec3::new(-0.25, 0.0, -0.01),
            frames: vec![f(0.0, 0.004, 0.004, 0.0), f(0.05, 0.035, 0.04, 0.0), f(0.2, 0.045, 0.055, 0.005), f(0.45, 0.03, 0.035, 0.01), f(0.9, 0.012, 0.012, 0.03), f(1.25, 0.006, 0.006, 0.035)],
            exponent: 2.0,
            nx: 30,
            nh: 8,
            mass: 0.3,
            nurbs: None,
            triangles: vec![],
        }
    }

    /// A fuselage from a triangulated surface (STL / STEP), in plane axes.
    pub fn from_triangles(name: &str, mut tris: Vec<[Vec3; 3]>, mass: f64) -> Fuse {
        crate::cad::orient(&mut tris);
        let (lo, _) = crate::cad::bbox(&tris);
        let pos = Vec3::new(lo.x, 0.0, 0.0);
        let tris = tris.into_iter().map(|t| [t[0] - pos, t[1] - pos, t[2] - pos]).collect();
        Fuse { name: name.into(), position: pos, frames: vec![], exponent: 2.0, nx: 0, nh: 0, mass, nurbs: None, triangles: tris }
    }

    pub fn kind_label(&self) -> &'static str {
        if !self.triangles.is_empty() {
            "triangulated surface (STL / STEP)"
        } else if self.nurbs.is_some() {
            "NURBS (frames of control points)"
        } else {
            "elliptic frames"
        }
    }

    /// Convert elliptic frames to an equivalent NURBS body (5 points per frame).
    pub fn to_nurbs(&self) -> Option<NurbsBody> {
        if self.frames.len() < 2 {
            return None;
        }
        let k = 0.77;
        let frames = self
            .frames
            .iter()
            .map(|f| {
                let pts = vec![
                    Vec3::new(f.x, 0.0, f.z + f.half_height),
                    Vec3::new(f.x, f.half_width * k, f.z + f.half_height * k),
                    Vec3::new(f.x, f.half_width, f.z),
                    Vec3::new(f.x, f.half_width * k, f.z - f.half_height * k),
                    Vec3::new(f.x, 0.0, f.z - f.half_height),
                ];
                NurbsFrame { x: f.x, points: pts }
            })
            .collect();
        Some(NurbsBody { frames, degree_u: 3, degree_v: 2 })
    }

    pub fn length(&self) -> f64 {
        if !self.triangles.is_empty() {
            let (lo, hi) = crate::cad::bbox(&self.triangles);
            return hi.x - lo.x;
        }
        if let Some(n) = &self.nurbs {
            return match (n.frames.first(), n.frames.last()) {
                (Some(a), Some(b)) => b.x - a.x,
                _ => 0.0,
            };
        }
        match (self.frames.first(), self.frames.last()) {
            (Some(a), Some(b)) => b.x - a.x,
            _ => 0.0,
        }
    }

    pub fn validate(&self) -> Result<(), String> {
        if !self.triangles.is_empty() {
            if self.triangles.len() < 4 {
                return Err(format!("fuselage '{}': too few triangles", self.name));
            }
            return Ok(());
        }
        if let Some(n) = &self.nurbs {
            if n.frames.len() < 2 {
                return Err(format!("fuselage '{}' needs at least 2 NURBS frames", self.name));
            }
            if n.frames.iter().any(|f| f.points.len() < 2) {
                return Err(format!("fuselage '{}': every NURBS frame needs at least 2 points", self.name));
            }
        } else {
            if self.frames.len() < 2 {
                return Err(format!("fuselage '{}' needs at least 2 frames", self.name));
            }
            for w in self.frames.windows(2) {
                if w[1].x <= w[0].x {
                    return Err(format!("fuselage '{}': frame x positions must increase", self.name));
                }
            }
        }
        if self.nx < 2 || self.nh < 2 {
            return Err(format!("fuselage '{}': at least 2 panels in each direction", self.name));
        }
        Ok(())
    }

    /// Interpolated elliptic frames with cosine spacing.
    pub fn stations(&self) -> Vec<FuseFrame> {
        let xs: Vec<f64> = self.frames.iter().map(|f| f.x).collect();
        let sp = |g: fn(&FuseFrame) -> f64| CubicSpline::new(&xs, &self.frames.iter().map(g).collect::<Vec<_>>());
        let (sz, sw, sh) = (sp(|f| f.z), sp(|f| f.half_width), sp(|f| f.half_height));
        let (x0, x1) = (xs[0], xs[xs.len() - 1]);
        let minr = 1e-3 * (x1 - x0);
        distribution(Spacing::Cosine, self.nx)
            .into_iter()
            .map(|t| {
                let x = x0 + t * (x1 - x0);
                FuseFrame { x, z: sz.eval(x), half_width: sw.eval(x).max(minr), half_height: sh.eval(x).max(minr) }
            })
            .collect()
    }

    /// Surface point of an elliptic frame: hoop angle phi in [0, pi] from top
    /// to bottom on the right side (negative y for `left`).
    pub fn point(&self, f: &FuseFrame, phi: f64, left: bool) -> Vec3 {
        let e = 2.0 / self.exponent.max(1.0);
        let (s, c) = phi.sin_cos();
        let y = f.half_width * s.abs().powf(e) * s.signum();
        let z = f.z + f.half_height * c.abs().powf(e) * c.signum();
        self.position + Vec3::new(f.x, if left { -y } else { y }, z)
    }

    /// Structured surface grid of one half, plane axes: [i along x][k top->bottom].
    /// None for triangulated fuselages.
    pub fn grid(&self, nx: usize, nh: usize, left: bool) -> Option<Vec<Vec<Vec3>>> {
        if !self.triangles.is_empty() {
            return None;
        }
        let nh = nh.max(2);
        if let Some(n) = &self.nurbs {
            let s = n.surface(left);
            let ((u0, u1), (v0, v1)) = s.domain();
            let us = distribution(Spacing::Cosine, nx.max(2));
            let g = us.iter().map(|&u| (0..=nh).map(|k| self.position + s.eval(u0 + u * (u1 - u0), v0 + (v1 - v0) * k as f64 / nh as f64)).collect()).collect();
            return Some(g);
        }
        let mut f = self.clone();
        f.nx = nx.max(2);
        let st = f.stations();
        Some(st.iter().map(|fr| (0..=nh).map(|k| self.point(fr, std::f64::consts::PI * k as f64 / nh as f64, left)).collect()).collect())
    }

    /// The whole surface as triangles in plane axes.
    pub fn world_triangles(&self, nx: usize, nh: usize) -> Vec<[Vec3; 3]> {
        if !self.triangles.is_empty() {
            return self.triangles.iter().map(|t| [t[0] + self.position, t[1] + self.position, t[2] + self.position]).collect();
        }
        let mut out = Vec::new();
        for left in [false, true] {
            if let Some(g) = self.grid(nx, nh, left) {
                for i in 0..g.len() - 1 {
                    for k in 0..g[i].len() - 1 {
                        let (a, b, c, d) = (g[i][k], g[i + 1][k], g[i + 1][k + 1], g[i][k + 1]);
                        out.push([a, b, c]);
                        out.push([a, c, d]);
                    }
                }
            }
        }
        out
    }

    /// Half width, half height and centre z at plane x (None outside the fuse).
    pub fn section_at(&self, x_plane: f64) -> Option<(f64, f64, f64)> {
        if self.triangles.is_empty() && self.nurbs.is_none() {
            let x = x_plane - self.position.x;
            let (x0, x1) = (self.frames.first()?.x, self.frames.last()?.x);
            if x < x0 || x > x1 {
                return None;
            }
            let xs: Vec<f64> = self.frames.iter().map(|f| f.x).collect();
            let w = CubicSpline::new(&xs, &self.frames.iter().map(|f| f.half_width).collect::<Vec<_>>()).eval(x);
            let h = CubicSpline::new(&xs, &self.frames.iter().map(|f| f.half_height).collect::<Vec<_>>()).eval(x);
            let z = CubicSpline::new(&xs, &self.frames.iter().map(|f| f.z).collect::<Vec<_>>()).eval(x);
            return Some((w, h, z + self.position.z));
        }
        // generic: intersect the triangulated surface with the plane x = x_plane
        let tris = self.world_triangles(40, 12);
        let mut ymax: f64 = 0.0;
        let (mut zmin, mut zmax) = (f64::MAX, f64::MIN);
        let mut hit = false;
        for t in &tris {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if (a.x - x_plane) * (b.x - x_plane) <= 0.0 && (a.x - b.x).abs() > 1e-15 {
                    let s = (x_plane - a.x) / (b.x - a.x);
                    let p = a.lerp(b, s);
                    ymax = ymax.max(p.y.abs());
                    zmin = zmin.min(p.z);
                    zmax = zmax.max(p.z);
                    hit = true;
                }
            }
        }
        if !hit {
            return None;
        }
        Some((ymax, 0.5 * (zmax - zmin), 0.5 * (zmax + zmin)))
    }

    /// Wetted area.
    pub fn wetted_area(&self) -> f64 {
        self.world_triangles(48, 24).iter().map(|t| 0.5 * (t[1] - t[0]).cross(t[2] - t[0]).norm()).sum()
    }

    /// Max equivalent diameter / length (for form factor).
    pub fn fineness_inverse(&self) -> f64 {
        let l = self.length();
        if l <= 0.0 {
            return 0.0;
        }
        let tris = self.world_triangles(40, 12);
        let (lo, _) = crate::cad::bbox(&tris);
        let mut d: f64 = 0.0;
        for k in 1..40 {
            if let Some((w, h, _)) = self.section_at(lo.x + l * k as f64 / 40.0) {
                d = d.max(2.0 * (w * h).sqrt());
            }
        }
        d / l
    }
}

/// Mass properties: mass, CoG, inertia tensor about the CoG.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Inertia {
    pub mass: f64,
    pub cog: Vec3,
    pub ixx: f64,
    pub iyy: f64,
    pub izz: f64,
    pub ixz: f64,
    pub ixy: f64,
    pub iyz: f64,
}

/// Accumulates lumped masses and returns mass, CoG and inertia about the CoG.
#[derive(Default)]
struct MassAcc {
    m: Vec<(f64, Vec3)>,
}
impl MassAcc {
    fn add(&mut self, m: f64, p: Vec3) {
        if m > 0.0 && p.is_finite() {
            self.m.push((m, p));
        }
    }
    fn finish(&self) -> Inertia {
        let mass: f64 = self.m.iter().map(|a| a.0).sum();
        if mass <= 0.0 {
            return Inertia::default();
        }
        let cog = self.m.iter().fold(Vec3::ZERO, |a, (m, p)| a + *p * *m) / mass;
        let mut i = Inertia { mass, cog, ..Default::default() };
        for (m, p) in &self.m {
            let d = *p - cog;
            i.ixx += m * (d.y * d.y + d.z * d.z);
            i.iyy += m * (d.x * d.x + d.z * d.z);
            i.izz += m * (d.x * d.x + d.y * d.y);
            i.ixz += m * d.x * d.z;
            i.ixy += m * d.x * d.y;
            i.iyz += m * d.y * d.z;
        }
        i
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Plane {
    pub name: String,
    #[serde(default)]
    pub description: String,
    pub wings: Vec<Wing>,
    #[serde(default)]
    pub fuses: Vec<Fuse>,
    #[serde(default)]
    pub point_masses: Vec<PointMass>,
}

impl Plane {
    /// A conventional glider-like default (wing, elevator, fin, pod).
    pub fn default_plane(wing_foil: &str, tail_foil: &str) -> Plane {
        Plane {
            name: "Plane".into(),
            description: String::new(),
            wings: vec![Wing::default_main(wing_foil), Wing::default_elevator(tail_foil), Wing::default_fin(tail_foil)],
            fuses: vec![],
            point_masses: vec![PointMass { tag: "ballast / battery".into(), mass: 0.6, position: Vec3::new(-0.08, 0.0, 0.0) }],
        }
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.wings.is_empty() {
            return Err("the plane has no wing".into());
        }
        for w in &self.wings {
            w.validate()?;
        }
        for f in &self.fuses {
            f.validate()?;
        }
        Ok(())
    }

    /// The reference wing (first main wing, or the first wing).
    pub fn main_wing(&self) -> Option<&Wing> {
        self.wings.iter().find(|w| w.kind == WingKind::Main).or(self.wings.first())
    }

    /// Foil names used by the wings.
    pub fn foil_names(&self) -> Vec<String> {
        let mut v: Vec<String> = Vec::new();
        for w in &self.wings {
            for s in &w.sections {
                for f in [&s.foil_left, &s.foil_right] {
                    if !v.contains(f) {
                        v.push(f.clone());
                    }
                }
            }
        }
        v
    }

    /// Tail volumes (horizontal, vertical) relative to the main wing.
    pub fn tail_volumes(&self) -> (f64, f64) {
        let Some(m) = self.main_wing() else { return (0.0, 0.0) };
        let xm = m.position.x + m.sections[0].offset + 0.25 * m.sections[0].chord;
        let (sm, cm, bm) = (m.projected_area(), m.mac(), m.projected_span());
        let mut h = 0.0;
        let mut v = 0.0;
        for w in &self.wings {
            let l = w.position.x + w.sections[0].offset + 0.25 * w.sections[0].chord - xm;
            match w.kind {
                WingKind::Elevator if sm * cm > 0.0 => h += w.projected_area() * l / (sm * cm),
                WingKind::Fin if sm * bm > 0.0 => v += w.planform_area() * l / (sm * bm),
                _ => {}
            }
        }
        (h, v)
    }

    /// Mass properties.  Wing structural mass is spread over the wing
    /// volume (by local foil thickness), fuselage mass over its skin.
    pub fn inertia(&self, db: &FoilDb) -> Inertia {
        let mut acc = MassAcc::default();
        for w in &self.wings {
            for pm in &w.point_masses {
                acc.add(pm.mass, w.position + pm.position);
            }
            if w.mass > 0.0 {
                let mut elems: Vec<(f64, Vec3)> = Vec::new();
                for left in if w.two_sided { vec![false, true] } else { vec![false] } {
                    let fr = w.frames(left);
                    for i in 0..fr.len() - 1 {
                        let (sa, sb) = (db.shape(w.foil(i, left)), db.shape(w.foil(i + 1, left)));
                        let ta = sa.as_ref().map(|s| s.area).unwrap_or(0.0);
                        let tb = sb.as_ref().map(|s| s.area).unwrap_or(0.0);
                        let ns = 10;
                        let nc = 8;
                        let span = (fr[i + 1].q4 - fr[i].q4).norm();
                        for k in 0..ns {
                            let t = (k as f64 + 0.5) / ns as f64;
                            let c = fr[i].chord + t * (fr[i + 1].chord - fr[i].chord);
                            let thick_area = (ta + t * (tb - ta)).max(0.01) * c * c;
                            for l in 0..nc {
                                let xi = (l as f64 + 0.5) / nc as f64;
                                // thickness-weighted chordwise distribution ~ sqrt(xi)(1-xi)
                                let wgt = xi.sqrt() * (1.0 - xi) + 0.05;
                                let p = fr[i].point(xi, 0.0).lerp(fr[i + 1].point(xi, 0.0), t);
                                elems.push((thick_area * span / ns as f64 * wgt, p));
                            }
                        }
                    }
                }
                let tot: f64 = elems.iter().map(|e| e.0).sum();
                if tot > 0.0 {
                    for (v, p) in elems {
                        acc.add(w.mass * v / tot, p);
                    }
                }
            }
        }
        for f in &self.fuses {
            if f.mass > 0.0 {
                let tris = f.world_triangles(30, 12);
                let tot: f64 = tris.iter().map(|t| 0.5 * (t[1] - t[0]).cross(t[2] - t[0]).norm()).sum();
                if tot > 0.0 {
                    for t in &tris {
                        let a = 0.5 * (t[1] - t[0]).cross(t[2] - t[0]).norm();
                        acc.add(f.mass * a / tot, (t[0] + t[1] + t[2]) / 3.0);
                    }
                }
            }
        }
        for pm in &self.point_masses {
            acc.add(pm.mass, pm.position);
        }
        acc.finish()
    }

    /// Reference dimensions from the main wing: (area, span, chord).
    pub fn reference(&self, projected: bool) -> (f64, f64, f64) {
        match self.main_wing() {
            Some(w) => {
                if projected {
                    (w.projected_area(), w.projected_span(), w.mac())
                } else {
                    (w.planform_area(), w.planform_span(), w.mac())
                }
            }
            None => (1.0, 1.0, 1.0),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rectangle_dims() {
        let mut w = Wing::new("w", WingKind::Main, vec![WingSection::new(0.0, 0.2, 0.0, "f"), WingSection::new(1.0, 0.2, 0.0, "f")]);
        assert!((w.planform_area() - 0.4).abs() < 1e-12);
        assert!((w.mac() - 0.2).abs() < 1e-12);
        assert!((w.aspect_ratio() - 10.0).abs() < 1e-12);
        w.sections[0].dihedral = 60.0;
        assert!((w.projected_span() - 1.0).abs() < 1e-12);
        let fr = w.frames(false);
        assert!((fr[1].le().z - 60f64.to_radians().sin()).abs() < 1e-12);
        let fl = w.frames(true);
        assert!((fl[1].le().y + 0.5).abs() < 1e-12);
    }
    #[test]
    fn twist_and_fin() {
        let mut w = Wing::new("w", WingKind::Main, vec![WingSection::new(0.0, 1.0, 0.0, "f"), WingSection::new(1.0, 1.0, 0.0, "f")]);
        w.sections[1].twist = 10.0;
        let fr = w.frames(false);
        // nose up: TE below LE
        assert!(fr[1].te().z < fr[1].le().z);
        // twisted about quarter chord
        assert!((fr[1].point(0.25, 0.0) - Vec3::new(0.25, 1.0, 0.0)).norm() < 1e-12);
        let f = Wing::default_fin("f");
        let tip = f.frames(false)[1].le();
        assert!(tip.z > 0.25 && tip.y.abs() < 1e-9, "{tip:?}");
    }
    #[test]
    fn taper_mac() {
        let w = Wing::new("w", WingKind::Main, vec![WingSection::new(0.0, 0.3, 0.0, "f"), WingSection::new(1.5, 0.15, 0.0, "f")]);
        // MAC = 2/3 cr (1+l+l^2)/(1+l)
        let l: f64 = 0.5;
        assert!((w.mac() - 2.0 / 3.0 * 0.3 * (1.0 + l + l * l) / (1.0 + l)).abs() < 1e-12);
    }
}
