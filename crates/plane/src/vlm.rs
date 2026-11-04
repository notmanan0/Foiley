//! Vortex-lattice methods on thin (camber-surface) meshes.
//!
//! * VLM1: one horseshoe vortex per panel (bound at 1/4, legs to the far
//!   field along +x).
//! * VLM2: one vortex ring per panel (front edge at the panel's 1/4 chord,
//!   rear edge at the next panel's 1/4 chord); the trailing ring sheds
//!   semi-infinite legs.
//!
//! The influence matrix is factored once; the six unit solutions (unit
//! freestream along x, y, z and unit rotations p, q, r about the CoG) are
//! superposed for any operating point.

use crate::mesh::{Mesh, PanelKind};
use crate::PlaneError;
use foiley_core::linalg::{diagnose, DenseMatrix, Diagnostics, FnRows, Lu, SolveReport};
use foiley_core::{Progress, Vec3};
use rayon::prelude::*;
use std::f64::consts::PI;

pub const FOUR_PI: f64 = 4.0 * PI;

/// Velocity induced at `p` by a straight vortex segment a->b of unit
/// circulation (Biot-Savart), with a small viscous core `rc`.
#[inline]
pub fn segment(p: Vec3, a: Vec3, b: Vec3, rc: f64) -> Vec3 {
    let r0 = b - a;
    let r1 = p - a;
    let r2 = p - b;
    let c = r1.cross(r2);
    let c2 = c.norm2();
    let n1 = r1.norm();
    let n2 = r2.norm();
    let l0 = r0.norm2();
    if n1 < 1e-10 || n2 < 1e-10 || c2 < 1e-24 * l0.max(1e-30) {
        return Vec3::ZERO;
    }
    // distance from the line, for the core
    let h2 = c2 / l0.max(1e-30);
    let damp = if rc > 0.0 { h2 / (h2 + rc * rc) } else { 1.0 };
    let k = r0.dot(r1 / n1 - r2 / n2) / (FOUR_PI * c2) * damp;
    c * k
}

/// Velocity of a semi-infinite vortex starting at `a` going along unit `d`.
#[inline]
pub fn semi_infinite(p: Vec3, a: Vec3, d: Vec3, rc: f64) -> Vec3 {
    let r = p - a;
    let c = d.cross(r);
    let c2 = c.norm2();
    if c2 < 1e-24 {
        return Vec3::ZERO;
    }
    let damp = if rc > 0.0 { c2 / (c2 + rc * rc) } else { 1.0 };
    let cosb = r.dot(d) / r.norm().max(1e-30);
    // V = G/(4 pi h) (1 + cos b) e,  e = d x r / |d x r|
    c * ((1.0 + cosb) / (FOUR_PI * c2) * damp)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Ground {
    None,
    /// ground plane at z = -height (mirror image with reversed circulation)
    Ground { height_mm: i64 },
}

/// The vortex elements of a thin mesh.
pub struct Lattice<'a> {
    pub mesh: &'a Mesh,
    pub rings: bool,
    /// downstream neighbour in the strip (None = trailing panel)
    next: Vec<Option<usize>>,
    pub rc: f64,
    /// ground plane height below the origin (None = free flight)
    pub ground: Option<f64>,
}

impl<'a> Lattice<'a> {
    pub fn new(mesh: &'a Mesh, rings: bool, ground: Option<f64>) -> Self {
        let mut next = vec![None; mesh.panels.len()];
        for s in &mesh.strips {
            for w in s.panels.windows(2) {
                next[w[0]] = Some(w[1]);
            }
        }
        let size = mesh.panels.iter().map(|p| p.size).fold(0.0, f64::max);
        Lattice { mesh, rings, next, rc: 1e-6 * size.max(1e-3), ground }
    }

    /// Velocity at `p` induced by element `j` with unit circulation.
    /// `wake_only` keeps only the trailing legs (Trefftz-type evaluations).
    pub fn element_velocity(&self, j: usize, p: Vec3, wake_only: bool) -> Vec3 {
        let v = self.element_velocity_free(j, p, wake_only);
        match self.ground {
            None => v,
            Some(h) => {
                let pm = Vec3::new(p.x, p.y, -p.z - 2.0 * h);
                let vm = self.element_velocity_free(j, pm, wake_only);
                v + Vec3::new(vm.x, vm.y, -vm.z)
            }
        }
    }

    fn element_velocity_free(&self, j: usize, p: Vec3, wake_only: bool) -> Vec3 {
        let pj = &self.mesh.panels[j];
        let (va, vb) = (pj.va, pj.vb);
        let d = Vec3::X;
        let rc = self.rc;
        if !self.rings {
            let legs = semi_infinite(p, vb, d, rc) - semi_infinite(p, va, d, rc);
            return if wake_only { legs } else { legs + segment(p, va, vb, rc) };
        }
        match self.next[j] {
            Some(k) => {
                if wake_only {
                    return Vec3::ZERO;
                }
                let (na, nb) = (self.mesh.panels[k].va, self.mesh.panels[k].vb);
                segment(p, va, vb, rc) + segment(p, vb, nb, rc) + segment(p, nb, na, rc) + segment(p, na, va, rc)
            }
            None => {
                // close a quarter panel behind the TE, then shed legs
                let [la, ta, tb, lb] = pj.nodes;
                let ca = ta + (ta - la) * 0.25;
                let cb = tb + (tb - lb) * 0.25;
                let legs = semi_infinite(p, cb, d, rc) - semi_infinite(p, ca, d, rc);
                if wake_only {
                    return legs + segment(p, vb, cb, rc) + segment(p, ca, va, rc);
                }
                segment(p, va, vb, rc) + segment(p, vb, cb, rc) + segment(p, ca, va, rc) + legs
            }
        }
    }

    /// Net bound circulation on the 1/4-chord segment of each panel.
    pub fn bound_circulation(&self, gamma: &[f64]) -> Vec<f64> {
        if !self.rings {
            return gamma.to_vec();
        }
        let mut b = gamma.to_vec();
        for (j, nx) in self.next.iter().enumerate() {
            if let Some(k) = nx {
                b[*k] -= gamma[j];
            }
        }
        b
    }

    /// Circulation shed by each strip (strength of its trailing legs).
    pub fn strip_circulation(&self, gamma: &[f64]) -> Vec<f64> {
        self.mesh
            .strips
            .iter()
            .map(|s| if self.rings { s.panels.last().map(|&j| gamma[j]).unwrap_or(0.0) } else { s.panels.iter().map(|&j| gamma[j]).sum() })
            .collect()
    }

    /// Total induced velocity at a point (all elements).
    pub fn induced(&self, gamma: &[f64], p: Vec3) -> Vec3 {
        let mut v = Vec3::ZERO;
        for (j, g) in gamma.iter().enumerate() {
            if *g != 0.0 {
                v += self.element_velocity(j, p, false) * *g;
            }
        }
        v
    }
}

/// Factored VLM system with unit solutions.
pub struct VlmSystem<'a> {
    pub lat: Lattice<'a>,
    pub cog: Vec3,
    /// unit solutions: u, v, w, p, q, r
    pub unit: [Vec<f64>; 6],
    /// linear-solve diagnostics (empty when off)
    pub report: SolveReport,
}

impl<'a> VlmSystem<'a> {
    pub fn build(mesh: &'a Mesh, rings: bool, cog: Vec3, ground: Option<f64>, prog: &Progress) -> Result<Self, PlaneError> {
        Self::build_with(mesh, rings, cog, ground, prog, Diagnostics::Off)
    }

    pub fn build_with(mesh: &'a Mesh, rings: bool, cog: Vec3, ground: Option<f64>, prog: &Progress, diag: Diagnostics) -> Result<Self, PlaneError> {
        if mesh.thick {
            return Err(PlaneError::Input("VLM needs a thin mesh".into()));
        }
        let lat = Lattice::new(mesh, rings, ground);
        let n = mesh.panels.len();
        if n == 0 {
            return Err(PlaneError::Geometry("empty mesh".into()));
        }
        debug_assert!(mesh.panels.iter().all(|p| p.kind == PanelKind::Mid));
        prog.stage(format!("Building {n} x {n} vortex influence matrix"));
        let mut a = DenseMatrix::zeros(n, n);
        let done = std::sync::atomic::AtomicUsize::new(0);
        let cancelled = std::sync::atomic::AtomicBool::new(false);
        let sub = prog.sub(0.0, 0.7);
        a.fill_columns_par(|j, col| {
            if cancelled.load(std::sync::atomic::Ordering::Relaxed) {
                return;
            }
            for (i, v) in col.iter_mut().enumerate() {
                let pi = &mesh.panels[i];
                *v = lat.element_velocity(j, pi.ctrl, false).dot(pi.normal);
            }
            let d = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
            if d % 64 == 0 {
                sub.set(d as f64 / n as f64);
                if sub.is_cancelled() {
                    cancelled.store(true, std::sync::atomic::Ordering::Relaxed);
                }
            }
        });
        prog.check()?;
        prog.stage(format!("Factorising ({n} unknowns)"));
        let lu = Lu::factor_with(a, prog)?;
        prog.set(0.9);
        prog.check()?;
        let mut rhs = DenseMatrix::zeros(n, 6);
        for i in 0..n {
            let p = &mesh.panels[i];
            let r = p.ctrl - cog;
            let e = [Vec3::X, Vec3::Y, Vec3::Z];
            for k in 0..3 {
                rhs.set(i, k, -e[k].dot(p.normal));
                rhs.set(i, 3 + k, -r.cross(e[k]).dot(p.normal));
            }
        }
        let b0 = if diag != Diagnostics::Off { rhs.clone() } else { DenseMatrix::default() };
        lu.solve_many(&mut rhs);
        let report = if diag != Diagnostics::Off {
            prog.stage("Checking the solve");
            // rows re-assembled in f64 exactly as in the fill above
            let op = FnRows {
                n,
                block: 1,
                fill: |i: usize, out: &mut [f64]| {
                    let pi = &mesh.panels[i];
                    for (j, o) in out.iter_mut().enumerate() {
                        *o = lat.element_velocity(j, pi.ctrl, false).dot(pi.normal);
                    }
                },
            };
            let te: Vec<usize> = mesh.strips.iter().filter_map(|s| s.panels.last().copied()).collect();
            diagnose(diag, &lu, &op, &b0, &mut rhs, &te, 64)
        } else {
            SolveReport::default()
        };
        let unit = std::array::from_fn(|k| rhs.col(k).to_vec());
        prog.set(1.0);
        Ok(VlmSystem { lat, cog, unit, report })
    }

    /// Circulations for freestream `vinf` (relative air velocity) and body rates `omega`.
    pub fn solve(&self, vinf: Vec3, omega: Vec3) -> Vec<f64> {
        let c = [vinf.x, vinf.y, vinf.z, omega.x, omega.y, omega.z];
        let n = self.unit[0].len();
        let mut g = vec![0.0; n];
        for k in 0..6 {
            if c[k] != 0.0 {
                for (gi, u) in g.iter_mut().zip(&self.unit[k]) {
                    *gi += c[k] * u;
                }
            }
        }
        g
    }

    /// Aerodynamic loads (in N / (rho) i.e. per unit density) for given circulations.
    pub fn loads(&self, gamma: &[f64], vinf: Vec3, omega: Vec3) -> crate::analysis::RawLoads {
        let mesh = self.lat.mesh;
        let bound = self.lat.bound_circulation(gamma);
        let sg = self.lat.strip_circulation(gamma);
        let wdir = vinf.normalized();
        let mut panel_force = vec![Vec3::ZERO; mesh.panels.len()];
        let mut force = Vec3::ZERO;
        let mut moment = Vec3::ZERO;
        for (j, p) in mesh.panels.iter().enumerate() {
            let mid = (p.va + p.vb) * 0.5;
            let vloc = vinf - omega.cross(mid - self.cog);
            let f = vloc.cross(p.vb - p.va) * bound[j];
            panel_force[j] = f;
            force += f;
            moment += (mid - self.cog).cross(f);
        }
        // Trefftz-plane induced velocities at the strip trailing edges
        let tre: Vec<Vec3> = mesh
            .strips
            .par_iter()
            .map(|s| {
                let m = (s.te_a + s.te_b) * 0.5;
                trefftz_velocity(mesh, &sg, m, wdir, self.lat.ground)
            })
            .collect();
        crate::analysis::RawLoads::from_parts(mesh, panel_force, force, moment, sg, tre, vinf)
    }
}

/// Relaxed (free) wake: polylines shed from the trailing-element corners.
#[derive(Clone, Debug, Default)]
pub struct FreeWake {
    pub lines: Vec<Vec<Vec3>>,
    /// per element: (line of corner A, line of corner B) for trailing elements
    pub legs: Vec<Option<(usize, usize)>>,
}

impl<'a> VlmSystem<'a> {
    /// Velocity of element j without its trailing legs.
    fn element_core(&self, j: usize, p: Vec3) -> Vec3 {
        let lat = &self.lat;
        let pj = &lat.mesh.panels[j];
        let (va, vb) = (pj.va, pj.vb);
        let rc = lat.rc;
        if !lat.rings {
            return segment(p, va, vb, rc);
        }
        match lat.next[j] {
            Some(k) => {
                let (na, nb) = (lat.mesh.panels[k].va, lat.mesh.panels[k].vb);
                segment(p, va, vb, rc) + segment(p, vb, nb, rc) + segment(p, nb, na, rc) + segment(p, na, va, rc)
            }
            None => {
                let [la, ta, tb, lb] = pj.nodes;
                let ca = ta + (ta - la) * 0.25;
                let cb = tb + (tb - lb) * 0.25;
                segment(p, va, vb, rc) + segment(p, vb, cb, rc) + segment(p, ca, va, rc)
            }
        }
    }

    fn corners(&self, j: usize) -> (Vec3, Vec3) {
        let pj = &self.lat.mesh.panels[j];
        if !self.lat.rings {
            return (pj.va, pj.vb);
        }
        let [la, ta, tb, lb] = pj.nodes;
        (ta + (ta - la) * 0.25, tb + (tb - lb) * 0.25)
    }

    fn polyline(&self, line: &[Vec3], dir: Vec3, p: Vec3) -> Vec3 {
        let rc = self.lat.rc;
        let mut v = Vec3::ZERO;
        for w in line.windows(2) {
            v += segment(p, w[0], w[1], rc);
        }
        if let Some(last) = line.last() {
            v += semi_infinite(p, *last, dir, rc);
        }
        v
    }

    fn with_ground(&self, p: Vec3, f: impl Fn(Vec3) -> Vec3) -> Vec3 {
        match self.lat.ground {
            None => f(p),
            Some(h) => {
                let m = f(Vec3::new(p.x, p.y, -p.z - 2.0 * h));
                f(p) + Vec3::new(m.x, m.y, -m.z)
            }
        }
    }

    /// Solve with a relaxed wake for one freestream (no superposition).
    /// `length` = wake length (m), `nseg` = segments per wake line.
    pub fn solve_free_wake(&self, vinf: Vec3, iterations: usize, length: f64, nseg: usize, prog: &Progress) -> Result<(Vec<f64>, FreeWake), PlaneError> {
        let mesh = self.lat.mesh;
        let n = mesh.panels.len();
        let dir = vinf.normalized();
        let nseg = nseg.max(2);
        let tol = 1e-9 * length.max(1e-9);
        let mut starts: Vec<Vec3> = Vec::new();
        let mut legs = vec![None; n];
        for j in 0..n {
            if self.lat.next[j].is_none() {
                let (a, b) = self.corners(j);
                let mut idx = |p: Vec3| -> usize {
                    if let Some(i) = starts.iter().position(|q| (*q - p).norm() < tol) {
                        i
                    } else {
                        starts.push(p);
                        starts.len() - 1
                    }
                };
                let ia = idx(a);
                let ib = idx(b);
                legs[j] = Some((ia, ib));
            }
        }
        // geometric spacing along the wake (fine near the trailing edge)
        let r: f64 = 1.15;
        let tot: f64 = (0..nseg).map(|k| r.powi(k as i32)).sum();
        let seg_len: Vec<f64> = (0..nseg).map(|k| length * r.powi(k as i32) / tot).collect();
        let mut lines: Vec<Vec<Vec3>> = starts
            .iter()
            .map(|s| {
                let mut l = vec![*s];
                for d in &seg_len {
                    let last = *l.last().unwrap();
                    l.push(last + dir * *d);
                }
                l
            })
            .collect();
        let mut gamma = vec![0.0; n];
        for it in 0..=iterations {
            prog.check()?;
            prog.stage(format!("Free wake iteration {}/{}", it, iterations));
            let mut a = DenseMatrix::zeros(n, n);
            a.fill_columns_par(|j, col| {
                for (i, v) in col.iter_mut().enumerate() {
                    let p = &mesh.panels[i];
                    let vel = self.with_ground(p.ctrl, |x| {
                        let mut w = self.element_core(j, x);
                        if let Some((la, lb)) = legs[j] {
                            w += self.polyline(&lines[lb], dir, x) - self.polyline(&lines[la], dir, x);
                        }
                        w
                    });
                    *v = vel.dot(p.normal);
                }
            });
            let lu = Lu::factor_with(a, prog)?;
            let mut rhs: Vec<f64> = mesh.panels.iter().map(|p| -vinf.dot(p.normal)).collect();
            lu.solve(&mut rhs);
            gamma = rhs;
            if it == iterations {
                break;
            }
            // net circulation of each wake line
            let mut lg = vec![0.0; lines.len()];
            for j in 0..n {
                if let Some((la, lb)) = legs[j] {
                    lg[lb] += gamma[j];
                    lg[la] -= gamma[j];
                }
            }
            let induced = |x: Vec3| -> Vec3 {
                self.with_ground(x, |y| {
                    let mut v = Vec3::ZERO;
                    for j in 0..n {
                        if gamma[j] != 0.0 {
                            v += self.element_core(j, y) * gamma[j];
                        }
                    }
                    for (l, g) in lines.iter().zip(&lg) {
                        if *g != 0.0 {
                            v += self.polyline(l, dir, y) * *g;
                        }
                    }
                    v
                })
            };
            let new_lines: Vec<Vec<Vec3>> = lines
                .par_iter()
                .map(|l| {
                    let mut out = vec![l[0]];
                    for (k, d) in seg_len.iter().enumerate() {
                        let v = vinf + induced(l[k]);
                        let step = *out.last().unwrap() + v.normalized() * *d;
                        out.push(l[k + 1] * 0.5 + step * 0.5);
                    }
                    out
                })
                .collect();
            lines = new_lines;
            prog.set((it + 1) as f64 / (iterations + 1) as f64);
        }
        Ok((gamma, FreeWake { lines, legs }))
    }
}

/// Half the Trefftz-plane induced velocity (= the velocity the trailing
/// vortex sheet induces at the lifting line) at point `m`, for strip
/// circulations `sg`, with the wake aligned with the wind `w`.
pub fn trefftz_velocity(mesh: &Mesh, sg: &[f64], m: Vec3, w: Vec3, ground: Option<f64>) -> Vec3 {
    let proj = |p: Vec3| p - w * p.dot(w);
    let pm = proj(m);
    let mut v = Vec3::ZERO;
    let add = |v: &mut Vec3, at: Vec3, g: f64, sign: f64| {
        let r = pm - proj(at);
        let r2 = r.norm2();
        if r2 > 1e-16 {
            // infinite line along w with circulation g; half of it
            *v += w.cross(r) * (sign * g / (2.0 * PI * r2) * 0.5);
        }
    };
    let mut segs: Vec<(Vec3, Vec3, f64)> = mesh.strips.iter().zip(sg).map(|(s, g)| (s.te_a, s.te_b, *g)).collect();
    segs.extend(carry_segments(mesh, sg));
    for (a, b, g) in segs {
        if g == 0.0 {
            continue;
        }
        add(&mut v, b, g, 1.0);
        add(&mut v, a, g, -1.0);
        if let Some(h) = ground {
            let mir = |p: Vec3| Vec3::new(p.x, p.y, -p.z - 2.0 * h);
            add(&mut v, mir(b), g, -1.0);
            add(&mut v, mir(a), g, 1.0);
        }
    }
    v
}

/// Carry-through segments across fuselages: (start, end, circulation).
pub fn carry_segments(mesh: &Mesh, sg: &[f64]) -> Vec<(Vec3, Vec3, f64)> {
    mesh.carry.iter().map(|&(l, r)| (mesh.strips[l].te_b, mesh.strips[r].te_a, 0.5 * (sg[l] + sg[r]))).collect()
}
