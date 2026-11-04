//! Discontinuous P1 Galerkin source/doublet triangles (internal Dirichlet).
//!
//! Each triangle owns three doublet values; coincident vertices are deliberately
//! not welded, especially across the trailing edge. With barycentric test and
//! trial functions the equation is (K - M/2 + W T) mu = -S sigma. T takes the
//! upper-minus-lower trace at the two TE endpoints, and W convects that linear
//! trace downstream. Rows are divided by the target triangle's area.
//!
//! Source densities are constant, sigma = -n.Vrel. The inner source and linear
//! doublet integrals are analytic (solid angle and edge logarithms); outer
//! integrals use symmetric degree-five quadrature, subdivided for near pairs.
//! The self doublet term is the exact P1 mass matrix, not a collocation jump.

use crate::analysis::{trefftz_force, RawLoads};
use crate::mesh::{Mesh, Panel, PanelKind};
use crate::panel::{build_wake, combine_unit};
use crate::PlaneError;
use foiley_core::linalg::{diagnose, DenseMatrix, Diagnostics, FnRows, Lu, SolveReport};
use foiley_core::{Progress, Vec3};
use rayon::prelude::*;
use std::f64::consts::PI;
use std::sync::atomic::{AtomicUsize, Ordering};

const FOUR_PI: f64 = 4.0 * PI;
const RULE: [([f64; 3], f64); 7] = [
    ([1.0 / 3.0; 3], 0.225),
    ([0.059715871789770, 0.470142064105115, 0.470142064105115], 0.132394152788506),
    ([0.470142064105115, 0.059715871789770, 0.470142064105115], 0.132394152788506),
    ([0.470142064105115, 0.470142064105115, 0.059715871789770], 0.132394152788506),
    ([0.797426985353087, 0.101286507323456, 0.101286507323456], 0.125939180544827),
    ([0.101286507323456, 0.797426985353087, 0.101286507323456], 0.125939180544827),
    ([0.101286507323456, 0.101286507323456, 0.797426985353087], 0.125939180544827),
];

/// Image of `p` in the ground plane z = -h.
fn mirror(p: Vec3, h: f64) -> Vec3 {
    Vec3::new(p.x, p.y, -p.z - 2.0 * h)
}

#[derive(Clone)]
struct Triangle {
    nodes: [Vec3; 3],
    normal: Vec3,
    grad: [Vec3; 3],
    centre: Vec3,
    area: f64,
    size: f64,
}

impl Triangle {
    /// Same acceptance as the mesher (area > 1e-14): slivers only get large,
    /// finite basis gradients.
    fn new(p: &Panel) -> Result<Self, PlaneError> {
        let [a, b, c, _] = p.nodes;
        if !p.is_tri() || !p.nodes.iter().all(|v| v.is_finite()) || !p.normal.is_finite() || !p.area.is_finite() || p.area <= 1e-14 {
            return Err(PlaneError::Geometry("linear doublets require triangles of positive area".into()));
        }
        let (normal, area2) = (p.normal, 2.0 * p.area);
        Ok(Self {
            nodes: [a, b, c],
            normal,
            centre: p.centre,
            area: p.area,
            size: p.size,
            grad: [normal.cross(c - b) / area2, normal.cross(a - c) / area2, normal.cross(b - a) / area2],
        })
    }

    fn basis(&self, p: Vec3) -> [f64; 3] {
        let r = p - self.nodes[0];
        let b = self.grad[1].dot(r);
        let c = self.grad[2].dot(r);
        [1.0 - b - c, b, c]
    }

    /// Principal-value potentials. Boundary jumps are inserted in assembly.
    fn potentials(&self, p: Vec3) -> (f64, [f64; 3]) {
        let r = self.nodes.map(|v| p - v);
        if (p - self.centre).norm2() > 64.0 * self.size * self.size {
            // Smooth far field: quadrature avoids cancellation between large
            // barycentric coordinates in the edge formula (no monopole shortcut).
            let mut s = 0.0;
            let mut d = [0.0; 3];
            for (b, w) in RULE {
                let q = self.nodes[0] * b[0] + self.nodes[1] * b[1] + self.nodes[2] * b[2];
                let v = p - q;
                let len = v.norm();
                let weight = w * self.area / FOUR_PI;
                s -= weight / len;
                let dd = weight * v.dot(self.normal) / len.powi(3);
                for k in 0..3 {
                    d[k] += b[k] * dd;
                }
            }
            return (s, d);
        }
        let h = r[0].dot(self.normal);
        let lens = r.map(Vec3::norm);
        let omega = if h.abs() <= 1e-14 * self.size {
            0.0
        } else {
            2.0 * r[0].dot(r[1].cross(r[2])).atan2(lens[0] * lens[1] * lens[2] + r[0].dot(r[1]) * lens[2] + r[1].dot(r[2]) * lens[0] + r[2].dot(r[0]) * lens[1])
        };
        let mut edge = Vec3::ZERO;
        let mut boundary = 0.0;
        for k in 0..3 {
            let j = (k + 1) % 3;
            let e = self.nodes[j] - self.nodes[k];
            let len = e.norm();
            let outward = (e / len).cross(self.normal);
            let sum = lens[k] + lens[j];
            let log = (2.0 * len / (sum - len).max(1e-30 * self.size)).ln_1p();
            edge += outward * log;
            boundary -= r[k].dot(outward) * log;
        }
        let b = self.basis(p);
        let d = std::array::from_fn(|k| (b[k] * omega - h * self.grad[k].dot(edge)) / FOUR_PI);
        (-(boundary - h * omega) / FOUR_PI, d)
    }

    /// Direct plus ground image. On its own plane (`own`) the direct doublet
    /// principal value is zero, so only the source and the image remain.
    fn image_potentials(&self, p: Vec3, ground: Option<f64>, own: bool) -> (f64, [f64; 3]) {
        let (mut s, mut d) = self.potentials(p);
        if own {
            d = [0.0; 3];
        }
        if let Some(h) = ground {
            let (si, di) = self.potentials(mirror(p, h));
            s += si;
            for k in 0..3 {
                d[k] += di[k];
            }
        }
        (s, d)
    }
}

#[derive(Clone, Copy)]
struct Sample {
    p: Vec3,
    basis: [f64; 3],
    weight: f64,
}

fn quadrature(t: &Triangle, depth: usize) -> Vec<Sample> {
    fn visit(t: &Triangle, v: [Vec3; 3], depth: usize, fraction: f64, out: &mut Vec<Sample>) {
        if depth == 0 {
            for (b, w) in RULE {
                let p = v[0] * b[0] + v[1] * b[1] + v[2] * b[2];
                out.push(Sample { p, basis: t.basis(p), weight: fraction * w });
            }
        } else {
            let [a, b, c] = v;
            let (ab, bc, ca) = ((a + b) * 0.5, (b + c) * 0.5, (c + a) * 0.5);
            for sub in [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]] {
                visit(t, sub, depth - 1, fraction * 0.25, out);
            }
        }
    }
    let mut out = Vec::new();
    visit(t, t.nodes, depth, 1.0, &mut out);
    out
}

struct LinearWake {
    strip: usize,
    /// DOFs at endpoints A and B on each side of the TE.
    upper: [usize; 2],
    lower: [usize; 2],
    /// Triangle, with the endpoint (A=0 or B=1) associated with each vertex.
    panels: Vec<(Triangle, [usize; 3])>,
}

fn trace_dofs(mesh: &Mesh, ids: &[usize], a: Vec3, b: Vec3) -> Result<[usize; 2], PlaneError> {
    // Only the triangle containing the *whole edge* owns the trace. The other
    // half of the quad may share a vertex, but has an independent P1 value.
    for &id in ids {
        let p = &mesh.panels[id];
        let find = |x: Vec3| (0..3).find(|&k| (p.nodes[k] - x).norm() < 1e-8 * p.size);
        if let (Some(i), Some(j)) = (find(a), find(b)) {
            if i != j {
                return Ok([3 * id + i, 3 * id + j]);
            }
        }
    }
    // No owner, e.g. a collapsed TE triangle dropped by the mesher: take the
    // nearest vertex of the remaining TE triangles for each endpoint.
    let nearest = |x: Vec3| {
        ids.iter().flat_map(|&id| (0..3).map(move |k| (3 * id + k, (mesh.panels[id].nodes[k] - x).norm2()))).min_by(|p, q| p.1.total_cmp(&q.1)).map(|p| p.0)
    };
    match (nearest(a), nearest(b)) {
        (Some(i), Some(j)) => Ok([i, j]),
        _ => Err(PlaneError::Geometry("no triangle on the trailing-edge segment".into())),
    }
}

/// Factored P1 Galerkin system. Unit solutions have three entries per triangle,
/// in its node order, for translations u,v,w and rotations p,q,r.
pub struct LinearPanelSystem<'a> {
    pub mesh: &'a Mesh,
    pub cog: Vec3,
    pub unit_mu: [Vec<f64>; 6],
    pub ground: Option<f64>,
    triangles: Vec<Triangle>,
    /// Seven-point samples of each triangle, for loads.
    coarse: Vec<Vec<Sample>>,
    wake: Vec<LinearWake>,
    /// linear-solve diagnostics (empty when off)
    pub report: SolveReport,
}

impl<'a> LinearPanelSystem<'a> {
    pub fn build(
        mesh: &'a Mesh,
        cog: Vec3,
        wake_length: f64,
        wake_panels: usize,
        wake_ratio: f64,
        ground: Option<f64>,
        prog: &Progress,
    ) -> Result<Self, PlaneError> {
        Self::build_order(mesh, cog, wake_length, wake_panels, wake_ratio, ground, prog, 2, Diagnostics::Off)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn build_with(mesh: &'a Mesh, cog: Vec3, wake_length: f64, wake_panels: usize, wake_ratio: f64, ground: Option<f64>, prog: &Progress, diag: Diagnostics) -> Result<Self, PlaneError> {
        Self::build_order(mesh, cog, wake_length, wake_panels, wake_ratio, ground, prog, 2, diag)
    }

    fn build_order(
        mesh: &'a Mesh,
        cog: Vec3,
        wake_length: f64,
        wake_panels: usize,
        wake_ratio: f64,
        ground: Option<f64>,
        prog: &Progress,
        depth: usize,
        diag: Diagnostics,
    ) -> Result<Self, PlaneError> {
        prog.check()?;
        if !mesh.thick || mesh.panels.is_empty() {
            return Err(PlaneError::Input("linear doublets need a non-empty thick triangle mesh".into()));
        }
        if !cog.is_finite() || !wake_length.is_finite() || wake_length <= 0.0 || wake_panels == 0 || !wake_ratio.is_finite() || wake_ratio <= 0.0 {
            return Err(PlaneError::Input("invalid linear-panel wake settings or centre of gravity".into()));
        }
        let triangles = mesh.panels.iter().map(Triangle::new).collect::<Result<Vec<_>, _>>()?;
        if let Some(h) = ground {
            if !h.is_finite() || triangles.iter().flat_map(|t| t.nodes).any(|p| p.z <= -h) {
                return Err(PlaneError::Input("the mesh must lie above the ground plane".into()));
            }
        }
        let mut wake = Vec::new();
        for w in build_wake(mesh, wake_length, wake_panels, wake_ratio) {
            let s = &mesh.strips[w.strip];
            let upper = trace_dofs(mesh, &w.upper_te, s.te_a, s.te_b)?;
            let lower = trace_dofs(mesh, &w.lower_te, s.te_a, s.te_b)?;
            let mut panels = Vec::new();
            for p in w.panels {
                let [a, a2, b2, b] = p.nodes;
                for ([p0, p1, p2], side) in [([a, a2, b2], [0, 0, 1]), ([a, b2, b], [0, 1, 1])] {
                    // A collapsed wake triangle carries no doublet; skip it.
                    if let Ok(t) = Triangle::new(&Panel::tri(p0, p1, p2, PanelKind::Mid, s.wing, Some(w.strip))) {
                        panels.push((t, side));
                    }
                }
            }
            wake.push(LinearWake { strip: w.strip, upper, lower, panels });
        }
        let nt = triangles.len();
        let n = 3 * nt;
        prog.stage(format!("Building {n} x {n} Galerkin linear-doublet matrix"));
        let coarse: Vec<_> = triangles.iter().map(|t| quadrature(t, 0)).collect();
        let fine: Vec<_> = triangles.iter().map(|t| quadrature(t, depth)).collect();
        // Refined samples on target i only when the source (or its image) is near.
        let tested = |i: usize, source: &Triangle| {
            let (c, r) = (triangles[i].centre, 1.5 * (triangles[i].size + source.size));
            if (c - source.centre).norm() < r || ground.is_some_and(|h| (mirror(c, h) - source.centre).norm() < r) {
                &fine[i]
            } else {
                &coarse[i]
            }
        };
        let axes = [Vec3::X, Vec3::Y, Vec3::Z];
        let sigma: Vec<[f64; 6]> = triangles
            .iter()
            .map(|t| {
                std::array::from_fn(|k| {
                    let v = if k < 3 { axes[k] } else { (t.centre - cog).cross(axes[k - 3]) };
                    -t.normal.dot(v)
                })
            })
            .collect();
        // Tested wake influence on each target: [wake][row][endpoint field].
        let done = AtomicUsize::new(0);
        let wake_rows: Vec<Vec<[[f64; 2]; 3]>> = (0..nt)
            .into_par_iter()
            .map(|i| -> Result<_, PlaneError> {
                prog.check()?;
                let rows = wake
                    .iter()
                    .map(|w| {
                        let mut rows = [[0.0; 2]; 3];
                        for (t, side) in &w.panels {
                            for sample in tested(i, t) {
                                let d = t.image_potentials(sample.p, ground, false).1;
                                for (r, row) in rows.iter_mut().enumerate() {
                                    let wt = sample.weight * sample.basis[r];
                                    for k in 0..3 {
                                        row[side[k]] += wt * d[k];
                                    }
                                }
                            }
                        }
                        rows
                    })
                    .collect();
                let count = done.fetch_add(1, Ordering::Relaxed) + 1;
                prog.set(0.1 * count as f64 / nt as f64);
                Ok(rows)
            })
            .collect::<Result<_, _>>()?;
        // TE-trace columns of each triangle: (local vertex, wake, endpoint, sign).
        let mut trace = vec![Vec::new(); nt];
        for (iw, w) in wake.iter().enumerate() {
            for e in 0..2 {
                trace[w.upper[e] / 3].push((w.upper[e] % 3, iw, e, 1.0));
                trace[w.lower[e] / 3].push((w.lower[e] % 3, iw, e, -1.0));
            }
        }
        // Each worker owns the three columns of source j and loops over all
        // targets; the six RHS columns are reduced from per-worker buffers.
        let mut a = DenseMatrix::zeros(n, n);
        let done = AtomicUsize::new(0);
        let data = a
            .data
            .par_chunks_mut(3 * n)
            .enumerate()
            .try_fold(
                || vec![0.0; 6 * n],
                |mut rhs, (j, cols)| -> Result<Vec<f64>, PlaneError> {
                    prog.check()?;
                    let source = &triangles[j];
                    for i in 0..nt {
                        if i % 32 == 0 {
                            prog.check()?;
                        }
                        for sample in tested(i, source) {
                            let (s, d) = source.image_potentials(sample.p, ground, i == j);
                            for r in 0..3 {
                                let w = sample.weight * sample.basis[r];
                                for k in 0..3 {
                                    cols[k * n + 3 * i + r] += w * d[k];
                                }
                                for k in 0..6 {
                                    rhs[k * n + 3 * i + r] -= w * s * sigma[j][k];
                                }
                            }
                        }
                    }
                    for r in 0..3 {
                        for k in 0..3 {
                            cols[k * n + 3 * j + r] -= if r == k { 1.0 / 12.0 } else { 1.0 / 24.0 };
                        }
                    }
                    for &(k, iw, e, sign) in &trace[j] {
                        for (col, rows) in cols[k * n..(k + 1) * n].chunks_mut(3).zip(&wake_rows) {
                            for r in 0..3 {
                                col[r] += sign * rows[iw][r][e];
                            }
                        }
                    }
                    let count = done.fetch_add(1, Ordering::Relaxed) + 1;
                    prog.set(0.1 + 0.7 * count as f64 / nt as f64);
                    Ok(rhs)
                },
            )
            .try_reduce(
                || vec![0.0; 6 * n],
                |mut x, y| {
                    x.iter_mut().zip(y).for_each(|(x, y)| *x += y);
                    Ok(x)
                },
            )?;
        prog.check()?;
        let mut rhs = DenseMatrix { rows: n, cols: 6, data };
        prog.stage(format!("Factorising Galerkin system ({n} unknowns)"));
        let lu = Lu::factor(a)?;
        prog.check()?;
        let b0 = if diag != Diagnostics::Off { rhs.clone() } else { DenseMatrix::default() };
        lu.solve_many(&mut rhs);
        let report = if diag != Diagnostics::Off {
            prog.stage("Checking the solve");
            // the three Galerkin rows of target triangle i, re-assembled with
            // the same operations, per entry in the same order, as above
            let op = FnRows {
                n,
                block: 3,
                fill: |i: usize, out: &mut [f64]| {
                    out.iter_mut().for_each(|v| *v = 0.0);
                    for (j, source) in triangles.iter().enumerate() {
                        for sample in tested(i, source) {
                            let (_, d) = source.image_potentials(sample.p, ground, i == j);
                            for r in 0..3 {
                                let w = sample.weight * sample.basis[r];
                                for k in 0..3 {
                                    out[r * n + 3 * j + k] += w * d[k];
                                }
                            }
                        }
                        if i == j {
                            for r in 0..3 {
                                for k in 0..3 {
                                    out[r * n + 3 * j + k] -= if r == k { 1.0 / 12.0 } else { 1.0 / 24.0 };
                                }
                            }
                        }
                        for &(k, iw, e, sign) in &trace[j] {
                            for r in 0..3 {
                                out[r * n + 3 * j + k] += sign * wake_rows[i][iw][r][e];
                            }
                        }
                    }
                },
            };
            let te: Vec<usize> = wake.iter().flat_map(|w| w.upper.iter().chain(&w.lower).copied()).collect();
            diagnose(diag, &lu, &op, &b0, &mut rhs, &te, 64)
        } else {
            SolveReport::default()
        };
        drop(wake_rows);
        if rhs.data.iter().any(|v| !v.is_finite()) {
            return Err(PlaneError::Solver("non-finite linear-doublet solution".into()));
        }
        let unit_mu = std::array::from_fn(|k| rhs.col(k).to_vec());
        prog.set(1.0);
        Ok(Self { mesh, cog, unit_mu, ground, triangles, coarse, wake, report })
    }

    pub fn solve(&self, vinf: Vec3, omega: Vec3) -> Vec<f64> {
        combine_unit(&self.unit_mu, vinf, omega)
    }

    /// Exact in-element gradient, with no averaging across the TE.
    pub fn mu_gradient(&self, mu: &[f64]) -> Vec<Vec3> {
        assert_eq!(mu.len(), 3 * self.triangles.len());
        self.triangles.iter().enumerate().map(|(i, t)| t.grad[1] * (mu[3 * i + 1] - mu[3 * i]) + t.grad[2] * (mu[3 * i + 2] - mu[3 * i])).collect()
    }

    pub fn loads(&self, vinf: Vec3, omega: Vec3) -> RawLoads {
        let mu = self.solve(vinf, omega);
        let grad = self.mu_gradient(&mu);
        let v2 = vinf.norm2().max(1e-30);
        let mut panel_force = Vec::with_capacity(self.triangles.len());
        let mut cp = Vec::with_capacity(self.triangles.len());
        let mut moment = Vec3::ZERO;
        for (i, (t, q)) in self.triangles.iter().zip(&self.coarse).enumerate() {
            let mut force = Vec3::ZERO;
            let mut c = 0.0;
            for sample in q {
                let v = vinf - omega.cross(sample.p - self.cog);
                let vt = v - t.normal * v.dot(t.normal) + grad[i];
                let pressure = (1.0 - vt.norm2() / v2) * sample.weight;
                let df = t.normal * (-0.5 * v2 * t.area * pressure);
                c += pressure;
                force += df;
                moment += (sample.p - self.cog).cross(df);
            }
            cp.push(c);
            panel_force.push(force);
        }
        let mut ends = vec![[0.0; 2]; self.mesh.strips.len()];
        for w in &self.wake {
            ends[w.strip] = std::array::from_fn(|k| mu[w.upper[k]] - mu[w.lower[k]]);
        }
        let sg: Vec<_> = ends.iter().map(|g| (g[0] + g[1]) * 0.5).collect();
        let wdir = vinf.normalized();
        let trefftz = |m: Vec3| linear_trefftz(self.mesh, &ends, m, wdir, self.ground);
        let tre = self.mesh.strips.iter().map(|s| trefftz((s.te_a + s.te_b) * 0.5)).collect();
        let mut loads = RawLoads::from_parts(self.mesh, panel_force, Vec3::ZERO, moment, sg.clone(), tre, vinf);
        // Carry-through drag from the same linear model (with ground) as the strips.
        let carry: Vec<_> = carry_segments(self.mesh, &ends).into_iter().map(|(a, b, g)| (a, b, (g[0] + g[1]) * 0.5)).collect();
        loads.induced_force = trefftz_force(self.mesh, &sg, &loads.strip_wi, wdir, &carry, trefftz);
        for (i, s) in self.mesh.strips.iter().enumerate() {
            let f = (vinf - omega.cross(s.c4 - self.cog)).cross(s.te_b - s.te_a) * sg[i];
            loads.strip_force[i] = f;
            loads.force += f;
        }
        loads.wake_force = loads.force;
        loads.cp = cp;
        loads.mu = mu;
        loads
    }
}

/// Fuselage carry-through, linear between the inner TE traces of the adjacent
/// strips so that no vortex is shed at the junctions.
fn carry_segments(mesh: &Mesh, ends: &[[f64; 2]]) -> Vec<(Vec3, Vec3, [f64; 2])> {
    mesh.carry.iter().map(|&(l, r)| (mesh.strips[l].te_b, mesh.strips[r].te_a, [ends[l][1], ends[r][0]])).collect()
}

/// Half-Trefftz velocity of a piecewise linear circulation: endpoint vortices
/// plus the distributed streamwise vorticity -dGamma/ds. The in-sheet integral
/// is a Cauchy principal value. Constant circulation reduces to the VLM kernel.
fn linear_trefftz(mesh: &Mesh, ends: &[[f64; 2]], p: Vec3, w: Vec3, ground: Option<f64>) -> Vec3 {
    fn segment(a: Vec3, b: Vec3, g: [f64; 2], p: Vec3, w: Vec3) -> Vec3 {
        let proj = |v: Vec3| v - w * w.dot(v);
        let ra = proj(p - a);
        let rb = proj(p - b);
        let d = proj(b - a);
        let len = d.norm();
        if len < 1e-14 {
            return Vec3::ZERO;
        }
        let kernel = |r: Vec3| if r.norm2() > 1e-28 { w.cross(r) / (FOUR_PI * r.norm2()) } else { Vec3::ZERO };
        let e = d / len;
        let u = ra.dot(e);
        let hvec = ra - e * u;
        let h = hvec.norm();
        let mut integral = e * (0.5 / len * (ra.norm2().max(1e-28) / rb.norm2().max(1e-28)).ln());
        if h > 1e-12 * len {
            integral += hvec * ((u / h).atan() - ((u - len) / h).atan()) / (len * h);
        }
        kernel(rb) * g[1] - kernel(ra) * g[0] - w.cross(integral) * ((g[1] - g[0]) / FOUR_PI)
    }
    let mut v = Vec3::ZERO;
    let mut segments: Vec<_> = mesh.strips.iter().zip(ends).map(|(s, g)| (s.te_a, s.te_b, *g)).collect();
    segments.extend(carry_segments(mesh, ends));
    for (a, b, g) in segments {
        v += segment(a, b, g, p, w);
        if let Some(h) = ground {
            v -= segment(mirror(a, h), mirror(b, h), g, p, w);
        }
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Plane, Wing, WingKind, WingSection};
    use crate::mesh::{build, MeshOptions};
    use crate::FoilDb;
    use foiley_foil::Foil;

    fn triangle() -> Triangle {
        Triangle::new(&Panel::tri(Vec3::new(-0.2, 0.1, 0.3), Vec3::new(1.3, 0.2, 0.4), Vec3::new(0.1, 0.9, -0.2), PanelKind::Fuse, 0, None)).unwrap()
    }

    #[test]
    fn analytic_kernels_match_independent_surface_integration() {
        let t = triangle();
        for p in [t.centre + t.normal * 0.2, t.centre - t.normal * 0.4, Vec3::new(0.9, -0.2, 1.4), Vec3::new(20.0, 5.0, -12.0)] {
            let (s, d) = t.potentials(p);
            let mut ss = 0.0;
            let mut dd = [0.0; 3];
            for q in quadrature(&t, 5) {
                let r = p - q.p;
                ss -= q.weight * t.area / (FOUR_PI * r.norm());
                for k in 0..3 {
                    dd[k] += q.weight * t.area * q.basis[k] * r.dot(t.normal) / (FOUR_PI * r.norm().powi(3));
                }
            }
            assert!((s - ss).abs() < 2e-9, "source {s} vs {ss}");
            for k in 0..3 {
                assert!((d[k] - dd[k]).abs() < 2e-9, "basis {k}: {} vs {}", d[k], dd[k]);
            }
        }
        let b = [0.2, 0.3, 0.5];
        let p = t.nodes[0] * b[0] + t.nodes[1] * b[1] + t.nodes[2] * b[2];
        for sign in [-1.0, 1.0] {
            let (_, d) = t.potentials(p + t.normal * (sign * 1e-9));
            for k in 0..3 {
                assert!((d[k] - sign * 0.5 * b[k]).abs() < 1e-7);
            }
        }
        assert!(t.potentials(p).1.iter().all(|d| d.abs() < 1e-12));
    }

    #[test]
    fn quadrature_mass_matrix_and_linear_gradient() {
        let t = triangle();
        for depth in [0, 1, 2] {
            let q = quadrature(&t, depth);
            for i in 0..3 {
                for j in 0..3 {
                    let mass: f64 = q.iter().map(|s| s.weight * s.basis[i] * s.basis[j]).sum();
                    assert!((mass - if i == j { 1.0 / 6.0 } else { 1.0 / 12.0 }).abs() < 1e-14);
                }
            }
        }
        let v = Vec3::new(0.6, -0.3, 0.1);
        let mu = t.nodes.map(|p| 2.0 + v.dot(p));
        let g = t.grad[1] * (mu[1] - mu[0]) + t.grad[2] * (mu[2] - mu[0]);
        assert!((g - (v - t.normal * v.dot(t.normal))).norm() < 1e-13);
    }

    fn wing_mesh() -> Mesh {
        let mut db = FoilDb::new();
        let foil = Foil::naca("0012").unwrap();
        db.add_foil(&foil);
        let mut a = WingSection::new(0.0, 1.0, 0.0, &foil.name);
        a.nx = 8;
        a.ny = 8;
        let b = WingSection::new(4.0, 1.0, 0.0, &foil.name);
        let plane = Plane {
            name: String::new(),
            description: String::new(),
            wings: vec![Wing::new("wing", WingKind::Main, vec![a, b])],
            fuses: vec![],
            point_masses: vec![],
        };
        build(&plane, &db, &MeshOptions { thick: true, triangulate: true, ..Default::default() }).unwrap()
    }

    #[test]
    fn quadrature_and_wake_refinement_preserve_te_lift() {
        let m = wing_mesh();
        let v = crate::analysis::wind_dir(5.0, 0.0);
        let mut results = Vec::new();
        for (depth, nw) in [(2, 5), (3, 5), (3, 10)] {
            let sys = LinearPanelSystem::build_order(&m, Vec3::ZERO, 100.0, nw, 1.1, None, &Progress::none(), depth, Diagnostics::Off).unwrap();
            let r = sys.loads(v, Vec3::ZERO);
            let mut te_jump: f64 = 0.0;
            for s in &m.strips {
                te_jump = te_jump.max((r.cp[s.te_up] - r.cp[s.te_low]).abs());
            }
            println!("depth={depth} nw={nw} Fz={} TE dCp={te_jump}", r.force.z);
            results.push(r);
        }
        assert!((results[0].force.z / results[1].force.z - 1.0).abs() < 0.001);
        assert!((results[1].force.z / results[2].force.z - 1.0).abs() < 1e-7);
        // Only actual TE trace values enter the wake, never the opposite vertex.
        let sys = LinearPanelSystem::build(&m, Vec3::ZERO, 100.0, 5, 1.1, None, &Progress::none()).unwrap();
        for w in &sys.wake {
            let s = &m.strips[w.strip];
            for (k, endpoint) in [s.te_a, s.te_b].into_iter().enumerate() {
                for id in [w.upper[k], w.lower[k]] {
                    assert!((m.panels[id / 3].nodes[id % 3] - endpoint).norm() < 1e-10);
                }
            }
        }
        let gamma = vec![0.7; m.strips.len()];
        let ends = vec![[0.7; 2]; m.strips.len()];
        for h in [None, Some(0.3)] {
            let p = Vec3::new(0.8, 0.2, 0.1);
            assert!((linear_trefftz(&m, &ends, p, v, h) - crate::vlm::trefftz_velocity(&m, &gamma, p, v, h)).norm() < 1e-12);
        }
        // Reordering vertices changes local DOF numbers, never the TE trace.
        let mut permuted = m.clone();
        for p in &mut permuted.panels {
            p.nodes = [p.nodes[1], p.nodes[2], p.nodes[0], p.nodes[0]];
        }
        let reordered = LinearPanelSystem::build(&permuted, Vec3::ZERO, 100.0, 5, 1.1, None, &Progress::none()).unwrap();
        let (a, b) = (sys.loads(v, Vec3::ZERO), reordered.loads(v, Vec3::ZERO));
        assert!((a.force - b.force).norm() < 1e-9);
        assert!((a.moment - b.moment).norm() < 1e-9);
        assert!(a.cp.iter().zip(b.cp).all(|(a, b)| (a - b).abs() < 1e-8));
    }

    #[test]
    fn linear_trefftz_matches_refined_constant_sheets() {
        let mut m = wing_mesh();
        // An artificial carry-through gap: its circulation must continue
        // linearly between the adjacent traces, with no junction vortices.
        m.carry = vec![(1, m.strips.len() - 2)];
        let p = Vec3::new(0.7, 0.2, 0.3);
        let w = crate::analysis::wind_dir(5.0, 3.0);
        let ends: Vec<_> = m.strips.iter().map(|s| [0.6 + 0.2 * s.te_a.y, 0.6 + 0.2 * s.te_b.y]).collect();
        let (l, r) = m.carry[0];
        let mut sheets: Vec<_> = m.strips.iter().zip(&ends).map(|(s, g)| (s.te_a, s.te_b, *g)).collect();
        sheets.push((m.strips[l].te_b, m.strips[r].te_a, [ends[l][1], ends[r][0]]));
        assert!((sheets.last().unwrap().1 - sheets.last().unwrap().0).norm() > 1.0);
        let kernel = |at: Vec3| {
            let r = p - at;
            let r = r - w * w.dot(r);
            w.cross(r) / (FOUR_PI * r.norm2())
        };
        for ground in [None, Some(0.6)] {
            let mut refined = Vec3::ZERO;
            for &(ta, tb, g) in &sheets {
                for i in 0..512 {
                    let a = ta.lerp(tb, i as f64 / 512.0);
                    let b = ta.lerp(tb, (i + 1) as f64 / 512.0);
                    let gamma = g[0] + (g[1] - g[0]) * (i as f64 + 0.5) / 512.0;
                    refined += (kernel(b) - kernel(a)) * gamma;
                    if let Some(h) = ground {
                        refined -= (kernel(mirror(b, h)) - kernel(mirror(a, h))) * gamma;
                    }
                }
            }
            assert!((linear_trefftz(&m, &ends, p, w, ground) - refined).norm() < 1e-7);
        }
    }

    #[test]
    fn cancellation_and_degenerate_geometry_are_errors() {
        let mut m = wing_mesh();
        let prog = Progress::new();
        prog.cancel();
        assert!(matches!(LinearPanelSystem::build(&m, Vec3::ZERO, 100.0, 5, 1.1, None, &prog), Err(PlaneError::Cancelled)));
        let p = &m.panels[0];
        m.panels[0] = Panel::tri(p.nodes[0], p.nodes[0], p.nodes[2], p.kind, p.wing, p.strip);
        assert!(matches!(LinearPanelSystem::build(&m, Vec3::ZERO, 100.0, 5, 1.1, None, &Progress::none()), Err(PlaneError::Geometry(_))));
    }

    #[test]
    fn mesher_slivers_and_missing_te_owner_are_accepted() {
        let m = wing_mesh();
        let v = crate::analysis::wind_dir(5.0, 0.0);
        let solve = |m: &Mesh| LinearPanelSystem::build(m, Vec3::ZERO, 100.0, 5, 1.1, Some(2.0), &Progress::none()).unwrap().loads(v, Vec3::ZERO);
        let base = solve(&m);
        // Split the longest non-TE edge at a point 1e-13 inside the triangle:
        // the sliver has area > 1e-14, so the mesher would keep it.
        let te: Vec<_> = m.strips.iter().flat_map(|s| [Some(s.te_up), Some(s.te_low), s.te_up2, s.te_low2]).flatten().collect();
        let len = |&(i, k): &(usize, usize)| (m.panels[i].nodes[(k + 1) % 3] - m.panels[i].nodes[k]).norm();
        let (id, k) = (0..m.panels.len()).filter(|i| !te.contains(i)).flat_map(|i| (0..3).map(move |k| (i, k))).max_by(|x, y| len(x).total_cmp(&len(y))).unwrap();
        let p = m.panels[id].clone();
        let (a, b, c) = (p.nodes[k], p.nodes[(k + 1) % 3], p.nodes[(k + 2) % 3]);
        let (mid, e) = ((a + b) * 0.5, (b - a).normalized());
        let s = mid + (c - mid - e * e.dot(c - mid)).normalized() * (8e-14 / (b - a).norm());
        let mut sliver = m.clone();
        sliver.panels[id] = Panel::tri(a, s, c, p.kind, p.wing, p.strip);
        sliver.panels.push(Panel::tri(s, b, c, p.kind, p.wing, p.strip));
        sliver.panels.push(Panel::tri(a, b, s, p.kind, p.wing, p.strip));
        sliver.neighbours.extend([vec![], vec![]]);
        let t = sliver.panels.last().unwrap();
        assert!(t.area > 1e-14 && 2.0 * t.area < 1e-12 * t.size * t.size, "area {} size {}", t.area, t.size);
        let r = solve(&sliver);
        println!("sliver Fz {} vs {}", r.force.z, base.force.z);
        assert!(r.force.is_finite() && r.moment.is_finite() && r.cp.iter().all(|c| c.is_finite()));
        assert!((r.force.z / base.force.z - 1.0).abs() < 1e-5);
        // A collapsed TE triangle dropped by the mesher leaves no owner of the
        // whole TE segment: the nearest vertices of its neighbour carry the trace.
        let mut dropped = m.clone();
        let st = &mut dropped.strips[0];
        st.te_up = st.te_up2.take().unwrap();
        let r = solve(&dropped);
        println!("missing TE owner Fz {} vs {}", r.force.z, base.force.z);
        assert!(r.force.is_finite() && (r.force.z / base.force.z - 1.0).abs() < 0.015);
    }

    #[test]
    fn sphere_solution_converges_to_analytic_potential_and_pressure() {
        let mut faces = Vec::new();
        for x in [-1.0, 1.0] {
            for y in [-1.0, 1.0] {
                for z in [-1.0, 1.0] {
                    let (a, b, c) = (Vec3::X * x, Vec3::Y * y, Vec3::Z * z);
                    faces.push(if x * y * z > 0.0 { [a, b, c] } else { [a, c, b] });
                }
            }
        }
        let mut errors = Vec::new();
        for level in 1..=4 {
            faces = faces
                .into_iter()
                .flat_map(|[a, b, c]| {
                    let (ab, bc, ca) = ((a + b).normalized(), (b + c).normalized(), (c + a).normalized());
                    [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]
                })
                .collect();
            if level < 3 {
                continue;
            }
            let panels: Vec<_> = faces.iter().map(|&[a, b, c]| Panel::tri(a, b, c, PanelKind::Fuse, usize::MAX, None)).collect();
            let n = panels.len();
            let m = Mesh { panels, strips: vec![], wings: vec![], flaps: vec![], thick: true, neighbours: vec![vec![]; n], n_wing_panels: 0, carry: vec![] };
            let sys = LinearPanelSystem::build(&m, Vec3::ZERO, 100.0, 5, 1.1, None, &Progress::none()).unwrap();
            let loads = sys.loads(Vec3::X, Vec3::ZERO);
            let potential =
                (m.panels.iter().enumerate().map(|(i, p)| (0..3).map(|k| (loads.mu[3 * i + k] - 0.5 * p.nodes[k].x).powi(2)).sum::<f64>()).sum::<f64>()
                    / (3 * n) as f64)
                    .sqrt();
            let pressure = (m.panels.iter().zip(&loads.cp).map(|(p, c)| (c - (1.0 - 2.25 * (1.0 - p.centre.normalized().x.powi(2)))).powi(2)).sum::<f64>()
                / n as f64)
                .sqrt();
            println!("sphere {n}: potential rms {potential} pressure rms {pressure}");
            assert!(loads.panel_force.iter().copied().fold(Vec3::ZERO, |a, b| a + b).norm() < 1e-8);
            errors.push((potential, pressure));
        }
        assert!(errors[1].0 < 0.5 * errors[0].0 && errors[1].0 < 0.015);
        assert!(errors[1].1 < 0.65 * errors[0].1 && errors[1].1 < 0.05);
    }
}
