//! Constant-strength source/doublet panel method on closed (thick) meshes.
//!
//! Formulation (Morino, internal Dirichlet): the perturbation potential
//! inside every body is zero, so on each panel the doublet strength equals
//! the outer perturbation potential, mu = phi, and the source strength is
//! sigma = -n . V_rel.  The trailing-edge condition is enforced by wake
//! doublet panels of strength mu_upper_TE - mu_lower_TE, shed along +x.
//! Surface velocities come from the tangential gradient of mu; pressure
//! forces are integrated over the panels and the induced drag is taken in
//! the Trefftz plane.
//!
//! Kernels: unit source  phi = -1/(4 pi) int dS / r,
//!          unit doublet phi = 1/(4 pi) int n.(P-Q) / |P-Q|^3 dS
//! (+1/2 just above a panel, -1/2 just below), evaluated with the exact
//! flat-quadrilateral formulas of NASA CR-4023 (VSAERO) and point
//! approximations in the far field.

use crate::mesh::{Mesh, Panel, PanelKind};
use crate::PlaneError;
use foiley_core::linalg::{diagnose, DenseMatrix, Diagnostics, FnRows, Lu, SolveReport};
use foiley_core::{Progress, Vec3};
use rayon::prelude::*;
use std::f64::consts::PI;

const FOUR_PI: f64 = 4.0 * PI;
/// far-field threshold in panel sizes
const RFF: f64 = 6.0;

/// (source potential, doublet potential) of a unit-strength flat quad at `p`.
pub fn quad_potentials(pan: &Panel, p: Vec3, want_src: bool) -> (f64, f64) {
    let r = p - pan.centre;
    let d = r.norm();
    if d > RFF * pan.size {
        let pn = r.dot(pan.normal);
        return (-pan.area / (FOUR_PI * d), pan.area * pn / (FOUR_PI * d * d * d));
    }
    let n = pan.normal;
    let (l, m) = (pan.l, pan.m);
    let mut pn = r.dot(n);
    let eps = 1e-10 * pan.size.max(1e-12);
    if pn.abs() < eps {
        pn = if pn < 0.0 { -eps } else { eps };
    }
    let mut cjk_sum = 0.0;
    let mut src_sum = 0.0;
    for k in 0..4 {
        let n1 = pan.nodes[k];
        let n2 = pan.nodes[(k + 1) % 4];
        let s = n2 - n1;
        let sl = s.norm();
        if sl < 1e-14 {
            continue;
        }
        let a = p - n1;
        let b = p - n2;
        let (an, bn) = (a.norm(), b.norm());
        let sm = s.dot(m);
        let sll = s.dot(l);
        let am = a.dot(m);
        let al = a.dot(l);
        let all = am * sll - al * sm;
        let pa = pn * pn * sll + all * am;
        let pb = pa - all * sm;
        let rnum = sm * pn * (bn * pa - an * pb);
        let dnom = pa * pb + pn * pn * an * bn * sm * sm;
        let cjk = rnum.atan2(dnom);
        cjk_sum += cjk;
        if want_src {
            let den = an + bn - sl;
            if den > 1e-14 * sl {
                let gl = ((an + bn + sl) / den).ln() / sl;
                src_sum += all * gl - pn * cjk;
            }
        }
    }
    (-src_sum / FOUR_PI, cjk_sum / FOUR_PI)
}

/// Velocity induced at `p` by a unit-strength flat source panel
/// (gradient of the source potential), NASA CR-4023 eq. 40.
pub fn quad_source_velocity(pan: &Panel, p: Vec3) -> Vec3 {
    let r = p - pan.centre;
    let d = r.norm();
    if d > RFF * pan.size {
        return r * (pan.area / (FOUR_PI * d * d * d));
    }
    let n = pan.normal;
    let (l, m) = (pan.l, pan.m);
    let mut pn = r.dot(n);
    let eps = 1e-10 * pan.size.max(1e-12);
    if pn.abs() < eps {
        pn = if pn < 0.0 { -eps } else { eps };
    }
    let mut v = Vec3::ZERO;
    for k in 0..4 {
        let n1 = pan.nodes[k];
        let n2 = pan.nodes[(k + 1) % 4];
        let s = n2 - n1;
        let sl = s.norm();
        if sl < 1e-14 {
            continue;
        }
        let a = p - n1;
        let b = p - n2;
        let (an, bn) = (a.norm(), b.norm());
        let sm = s.dot(m);
        let sll = s.dot(l);
        let am = a.dot(m);
        let al = a.dot(l);
        let all = am * sll - al * sm;
        let pa = pn * pn * sll + all * am;
        let pb = pa - all * sm;
        let cjk = (sm * pn * (bn * pa - an * pb)).atan2(pa * pb + pn * pn * an * bn * sm * sm);
        let den = an + bn - sl;
        let gl = if den > 1e-14 * sl { ((an + bn + sl) / den).ln() / sl } else { 0.0 };
        v += n * cjk + l * (sm * gl) - m * (sll * gl);
    }
    v * (1.0 / FOUR_PI)
}

/// Wake panels shed from one strip's trailing edge.
#[derive(Clone, Debug)]
pub struct WakeColumn {
    pub strip: usize,
    /// panels whose mean doublet strength defines the upper / lower TE value
    /// (one quad, or the two triangles of a split quad)
    pub upper_te: Vec<usize>,
    pub lower_te: Vec<usize>,
    pub panels: Vec<Panel>,
}

/// Doublets for one operating point from the unit solutions u, v, w, p, q, r.
pub(crate) fn combine_unit(unit: &[Vec<f64>; 6], vinf: Vec3, omega: Vec3) -> Vec<f64> {
    let c = [vinf.x, vinf.y, vinf.z, omega.x, omega.y, omega.z];
    let mut mu = vec![0.0; unit[0].len()];
    for k in 0..6 {
        if c[k] != 0.0 {
            for (m, u) in mu.iter_mut().zip(&unit[k]) {
                *m += c[k] * u;
            }
        }
    }
    mu
}

fn mean_mu(mu: &[f64], ids: &[usize]) -> f64 {
    ids.iter().map(|&i| mu[i]).sum::<f64>() / ids.len().max(1) as f64
}

pub(crate) fn build_wake(mesh: &Mesh, length: f64, n: usize, ratio: f64) -> Vec<WakeColumn> {
    let n = n.max(1);
    let series: f64 = (0..n).map(|k| ratio.powi(k as i32)).sum();
    let xend = mesh.strips.iter().map(|s| s.te_a.x.max(s.te_b.x)).fold(f64::MIN, f64::max) + length;
    mesh.strips
        .iter()
        .enumerate()
        .filter(|(_, s)| s.panels.len() >= 2)
        .map(|(is, s)| {
            let mut lower_te = vec![s.te_low];
            lower_te.extend(s.te_low2);
            let mut upper_te = vec![s.te_up];
            upper_te.extend(s.te_up2);
            let mut panels = Vec::with_capacity(n);
            let (mut a, mut b) = (s.te_a, s.te_b);
            let (mut la, mut lb) = ((xend - a.x) / series, (xend - b.x) / series);
            for _ in 0..n {
                let (a2, b2) = (a + Vec3::X * la, b + Vec3::X * lb);
                let mut p = mesh_panel([a, a2, b2, b], is);
                p.kind = PanelKind::Mid;
                panels.push(p);
                a = a2;
                b = b2;
                la *= ratio;
                lb *= ratio;
            }
            WakeColumn { strip: is, upper_te, lower_te, panels }
        })
        .collect()
}

fn mesh_panel(nodes: [Vec3; 4], strip: usize) -> Panel {
    // reuse the mesh panel constructor through a tiny thin mesh panel
    let [la, ta, tb, lb] = nodes;
    let nn = (tb - la).cross(lb - ta);
    let normal = nn.normalized();
    let centre = (la + ta + tb + lb) * 0.25;
    let mut m = ((lb + tb) * 0.5 - centre).normalized();
    m = (m - normal * normal.dot(m)).normalized();
    let l = m.cross(normal);
    Panel {
        nodes,
        normal,
        area: 0.5 * nn.norm(),
        centre,
        l,
        m,
        size: (tb - la).norm().max((lb - ta).norm()),
        ctrl: centre,
        va: la,
        vb: lb,
        kind: PanelKind::Mid,
        wing: usize::MAX,
        strip: Some(strip),
        trailing: false,
        leading: false,
        flap: None,
    }
}

/// Factored panel system with unit solutions.
pub struct PanelSystem<'a> {
    pub mesh: &'a Mesh,
    pub cog: Vec3,
    pub wake: Vec<WakeColumn>,
    /// unit doublet solutions: u, v, w, p, q, r (one value per panel)
    pub unit_mu: [Vec<f64>; 6],
    pub ground: Option<f64>,
    /// linear-solve diagnostics (empty when off)
    pub report: SolveReport,
}


impl<'a> PanelSystem<'a> {
    pub fn build(mesh: &'a Mesh, cog: Vec3, wake_length: f64, wake_panels: usize, wake_ratio: f64, ground: Option<f64>, prog: &Progress) -> Result<Self, PlaneError> {
        Self::build_with(mesh, cog, wake_length, wake_panels, wake_ratio, ground, prog, Diagnostics::Off)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn build_with(mesh: &'a Mesh, cog: Vec3, wake_length: f64, wake_panels: usize, wake_ratio: f64, ground: Option<f64>, prog: &Progress, diag: Diagnostics) -> Result<Self, PlaneError> {
        if !mesh.thick {
            return Err(PlaneError::Input("the panel method needs a thick mesh".into()));
        }
        let n = mesh.panels.len();
        if n == 0 {
            return Err(PlaneError::Geometry("empty mesh".into()));
        }
        let wake = build_wake(mesh, wake_length, wake_panels, wake_ratio);
        let pot = |pan: &Panel, p: Vec3, src: bool| -> (f64, f64) {
            let (s, d) = quad_potentials(pan, p, src);
            match ground {
                None => (s, d),
                Some(h) => {
                    let pm = Vec3::new(p.x, p.y, -p.z - 2.0 * h);
                    // image body: mirrored geometry; for a doublet the
                    // mirrored normal flips the z part -> evaluate the real
                    // panel at the mirrored point (same value by symmetry)
                    let (s2, d2) = quad_potentials(pan, pm, src);
                    (s + s2, d + d2)
                }
            }
        };
        prog.stage(format!("Building {n} x {n} doublet influence matrix"));
        let mut a = DenseMatrix::zeros(n, n);
        let done = std::sync::atomic::AtomicUsize::new(0);
        let stop = std::sync::atomic::AtomicBool::new(false);
        let sub = prog.sub(0.0, 0.45);
        a.fill_columns_par(|j, col| {
            if stop.load(std::sync::atomic::Ordering::Relaxed) {
                return;
            }
            let pj = &mesh.panels[j];
            for (i, v) in col.iter_mut().enumerate() {
                *v = if i == j {
                    let (_, d2) = match ground {
                        Some(h) => {
                            let p = mesh.panels[i].centre;
                            quad_potentials(pj, Vec3::new(p.x, p.y, -p.z - 2.0 * h), false)
                        }
                        None => (0.0, 0.0),
                    };
                    -0.5 + d2
                } else {
                    pot(pj, mesh.panels[i].centre, false).1
                };
            }
            let d = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
            if d % 128 == 0 {
                sub.set(d as f64 / n as f64);
                if sub.is_cancelled() {
                    stop.store(true, std::sync::atomic::Ordering::Relaxed);
                }
            }
        });
        prog.check()?;
        // wake contributions
        prog.stage("Adding wake influence");
        let wake_cols: Vec<(Vec<usize>, Vec<usize>, Vec<f64>)> = wake
            .par_iter()
            .map(|w| {
                let col: Vec<f64> = (0..n).map(|i| w.panels.iter().map(|wp| pot(wp, mesh.panels[i].centre, false).1).sum()).collect();
                (w.upper_te.clone(), w.lower_te.clone(), col)
            })
            .collect();
        for (up, lo, col) in wake_cols {
            let (ku, kl) = (1.0 / up.len() as f64, 1.0 / lo.len() as f64);
            for (i, v) in col.iter().enumerate() {
                for &u in &up {
                    a.add(i, u, *v * ku);
                }
                for &l in &lo {
                    a.add(i, l, -*v * kl);
                }
            }
        }
        prog.set(0.55);
        prog.check()?;
        // source unit strengths and RHS = - S sigma
        prog.stage("Building source terms");
        let e = [Vec3::X, Vec3::Y, Vec3::Z];
        let sig: Vec<[f64; 6]> = mesh
            .panels
            .iter()
            .map(|p| {
                let r = p.centre - cog;
                std::array::from_fn(|k| if k < 3 { -p.normal.dot(e[k]) } else { -p.normal.dot(r.cross(e[k - 3])) })
            })
            .collect();
        let rhs_rows: Vec<[f64; 6]> = (0..n)
            .into_par_iter()
            .map(|i| {
                let pi = mesh.panels[i].centre;
                let mut acc = [0.0; 6];
                for (j, pj) in mesh.panels.iter().enumerate() {
                    let s = pot(pj, pi, true).0;
                    for k in 0..6 {
                        acc[k] -= s * sig[j][k];
                    }
                }
                acc
            })
            .collect();
        prog.set(0.7);
        prog.check()?;
        prog.stage(format!("Factorising ({n} unknowns)"));
        let lu = Lu::factor_with(a, prog)?;
        prog.set(0.95);
        let mut rhs = DenseMatrix::zeros(n, 6);
        for (i, r) in rhs_rows.iter().enumerate() {
            for k in 0..6 {
                rhs.set(i, k, r[k]);
            }
        }
        let b0 = if diag != Diagnostics::Off { rhs.clone() } else { DenseMatrix::default() };
        lu.solve_many(&mut rhs);
        let report = if diag != Diagnostics::Off {
            prog.stage("Checking the solve");
            // rows re-assembled in f64 with the same operations, in the same
            // order, as the column fill and wake pass above
            let op = FnRows {
                n,
                block: 1,
                fill: |i: usize, out: &mut [f64]| {
                    let ci = mesh.panels[i].centre;
                    for (j, o) in out.iter_mut().enumerate() {
                        let pj = &mesh.panels[j];
                        *o = if i == j {
                            let (_, d2) = match ground {
                                Some(h) => quad_potentials(pj, Vec3::new(ci.x, ci.y, -ci.z - 2.0 * h), false),
                                None => (0.0, 0.0),
                            };
                            -0.5 + d2
                        } else {
                            pot(pj, ci, false).1
                        };
                    }
                    for w in &wake {
                        let v: f64 = w.panels.iter().map(|wp| pot(wp, ci, false).1).sum();
                        let (ku, kl) = (1.0 / w.upper_te.len() as f64, 1.0 / w.lower_te.len() as f64);
                        for &u in &w.upper_te {
                            out[u] += v * ku;
                        }
                        for &l in &w.lower_te {
                            out[l] += -v * kl;
                        }
                    }
                },
            };
            let te: Vec<usize> = wake.iter().flat_map(|w| w.upper_te.iter().chain(&w.lower_te).copied()).collect();
            diagnose(diag, &lu, &op, &b0, &mut rhs, &te, 64)
        } else {
            SolveReport::default()
        };
        let unit_mu = std::array::from_fn(|k| rhs.col(k).to_vec());
        prog.set(1.0);
        Ok(PanelSystem { mesh, cog, wake, unit_mu, ground, report })
    }

    pub fn solve(&self, vinf: Vec3, omega: Vec3) -> Vec<f64> {
        combine_unit(&self.unit_mu, vinf, omega)
    }

    /// Tangential gradient of mu on each panel (least squares on neighbours).
    pub fn mu_gradient(&self, mu: &[f64]) -> Vec<Vec3> {
        let mesh = self.mesh;
        (0..mesh.panels.len())
            .into_par_iter()
            .map(|i| {
                let p = &mesh.panels[i];
                let nb = &mesh.neighbours[i];
                let (mut sll, mut slm, mut smm, mut sbl, mut sbm) = (0.0, 0.0, 0.0, 0.0, 0.0);
                let mut cnt = 0;
                for &j in nb {
                    let q = &mesh.panels[j];
                    let d = q.centre - p.centre;
                    let (dl, dm) = (d.dot(p.l), d.dot(p.m));
                    let w = 1.0 / (dl * dl + dm * dm).max(1e-30);
                    let dmu = mu[j] - mu[i];
                    sll += w * dl * dl;
                    slm += w * dl * dm;
                    smm += w * dm * dm;
                    sbl += w * dl * dmu;
                    sbm += w * dm * dmu;
                    cnt += 1;
                }
                if cnt == 0 {
                    return Vec3::ZERO;
                }
                let det = sll * smm - slm * slm;
                if det.abs() > 1e-12 * (sll * smm).max(1e-300) {
                    let gl = (sbl * smm - sbm * slm) / det;
                    let gm = (sll * sbm - slm * sbl) / det;
                    p.l * gl + p.m * gm
                } else if sll >= smm && sll > 0.0 {
                    p.l * (sbl / sll)
                } else if smm > 0.0 {
                    p.m * (sbm / smm)
                } else {
                    Vec3::ZERO
                }
            })
            .collect()
    }

    /// Loads per unit density for one operating point.
    pub fn loads(&self, vinf: Vec3, omega: Vec3) -> crate::analysis::RawLoads {
        let mesh = self.mesh;
        let mu = self.solve(vinf, omega);
        let grad = self.mu_gradient(&mu);
        let v2 = vinf.norm2().max(1e-30);
        let mut panel_force = vec![Vec3::ZERO; mesh.panels.len()];
        let mut cp = vec![0.0; mesh.panels.len()];
        let mut force = Vec3::ZERO;
        let mut moment = Vec3::ZERO;
        for (i, p) in mesh.panels.iter().enumerate() {
            let vrel = vinf - omega.cross(p.centre - self.cog);
            let vt = vrel - p.normal * p.normal.dot(vrel) + grad[i];
            let c = 1.0 - vt.norm2() / v2;
            cp[i] = c;
            let f = p.normal * (-c * 0.5 * v2 * p.area);
            panel_force[i] = f;
            force += f;
            moment += (p.centre - self.cog).cross(f);
        }
        let mut sg = vec![0.0; mesh.strips.len()];
        for w in &self.wake {
            sg[w.strip] = mean_mu(&mu, &w.upper_te) - mean_mu(&mu, &w.lower_te);
        }
        let wdir = vinf.normalized();
        let _ = force;
        let tre: Vec<Vec3> = mesh.strips.par_iter().map(|s| crate::vlm::trefftz_velocity(mesh, &sg, (s.te_a + s.te_b) * 0.5, wdir, self.ground)).collect();
        let mut r = crate::analysis::RawLoads::from_parts(mesh, panel_force, Vec3::ZERO, moment, sg.clone(), tre, vinf);
        r.pressure_force = force;
        // Forces from the far field (Kutta-Joukowski on the shed circulation,
        // as flow5 does): they converge much faster with chordwise panel
        // count than integrated pressures.  Moments keep the pressures.
        let mut total = Vec3::ZERO;
        for (i, s) in mesh.strips.iter().enumerate() {
            let vloc = vinf - omega.cross(s.c4 - self.cog);
            let f = vloc.cross(s.te_b - s.te_a) * sg[i];
            r.strip_force[i] = f;
            total += f;
        }
        r.force = total;
        r.wake_force = total;
        r.cp = cp;
        r.mu = mu;
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cube(n: usize) -> Vec<Panel> {
        // unit cube centred at origin, outward normals
        let mut v = Vec::new();
        let h = 0.5;
        let faces: [(Vec3, Vec3, Vec3); 6] = [
            (Vec3::X, Vec3::Y, Vec3::Z),
            (-Vec3::X, Vec3::Z, Vec3::Y),
            (Vec3::Y, Vec3::Z, Vec3::X),
            (-Vec3::Y, Vec3::X, Vec3::Z),
            (Vec3::Z, Vec3::X, Vec3::Y),
            (-Vec3::Z, Vec3::Y, Vec3::X),
        ];
        for (nrm, u, w) in faces {
            for i in 0..n {
                for j in 0..n {
                    let f = |a: usize, b: usize| nrm * h + u * (-h + a as f64 / n as f64) + w * (-h + b as f64 / n as f64);
                    let p = mesh_panel([f(i, j), f(i + 1, j), f(i + 1, j + 1), f(i, j + 1)], 0);
                    let p = if p.normal.dot(nrm) < 0.0 { mesh_panel([f(i, j), f(i, j + 1), f(i + 1, j + 1), f(i + 1, j)], 0) } else { p };
                    v.push(p);
                }
            }
        }
        v
    }

    #[test]
    fn source_velocity_is_gradient() {
        let c = cube(1);
        let p = &c[0];
        for q in [p.centre + Vec3::new(0.3, 0.2, 0.1), p.centre + Vec3::new(1.2, -0.4, 0.5), p.centre + Vec3::new(-0.2, 0.1, -0.3)] {
            let v = quad_source_velocity(p, q);
            let h = 1e-6;
            let f = |d: Vec3| quad_potentials(p, q + d, true).0;
            let g = Vec3::new((f(Vec3::X * h) - f(-Vec3::X * h)) / (2.0 * h), (f(Vec3::Y * h) - f(-Vec3::Y * h)) / (2.0 * h), (f(Vec3::Z * h) - f(-Vec3::Z * h)) / (2.0 * h));
            assert!((v - g).norm() < 1e-6 * (1.0 + g.norm()), "{v:?} vs {g:?}");
        }
    }

    #[test]
    fn doublet_solid_angle() {
        let c = cube(3);
        let inside: f64 = c.iter().map(|p| quad_potentials(p, Vec3::new(0.1, -0.2, 0.15), false).1).sum();
        assert!((inside + 1.0).abs() < 1e-9, "inside {inside}");
        let outside: f64 = c.iter().map(|p| quad_potentials(p, Vec3::new(0.9, 0.3, -0.7), false).1).sum();
        assert!(outside.abs() < 1e-9, "outside {outside}");
        let p = &c[4];
        let above = quad_potentials(p, p.centre + p.normal * 1e-7, false).1;
        assert!((above - 0.5).abs() < 1e-4, "above {above}");
    }

    #[test]
    fn source_far_field_matches() {
        let c = cube(1);
        let p = &c[0];
        // far field: the exact kernel tends to the point-source value
        let qf = p.centre + Vec3::new(6.0, 1.6, 1.2);
        let exact_far = quad_potentials(&Panel { size: 1e9, ..p.clone() }, qf, true).0;
        let far = -p.area / (FOUR_PI * (qf - p.centre).norm());
        assert!((exact_far - far).abs() < 0.002 * far.abs(), "{exact_far} vs {far}");
        let q = p.centre + Vec3::new(1.5, 0.4, 0.3);
        let exact = quad_potentials(&Panel { size: 1e9, ..p.clone() }, q, true).0;
        // numeric integration check of the exact kernel
        let mut num = 0.0;
        let m = 200;
        for i in 0..m {
            for j in 0..m {
                let s = (i as f64 + 0.5) / m as f64;
                let t = (j as f64 + 0.5) / m as f64;
                let [a, b, _c, d] = p.nodes;
                let pt = a + (b - a) * s + (d - a) * t;
                num += -1.0 / (FOUR_PI * (q - pt).norm()) * p.area / (m * m) as f64;
            }
        }
        assert!((exact - num).abs() < 1e-5, "{exact} vs {num}");
    }
}

