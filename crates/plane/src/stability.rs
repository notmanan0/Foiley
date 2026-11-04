//! Flight dynamics: trim, stability derivatives, state matrices, modes and
//! time responses (Etkin & Reid small-perturbation equations in stability
//! axes; level flight, theta0 = phi0 = 0).  The state matrix is always the
//! coupled 8x8 [u w q theta v p r phi]; the 4x4 longitudinal and lateral
//! matrices are its diagonal blocks, and the modes are sorted into the two
//! sets by eigenvector energy.  The reference state may have sideslip: the
//! stability axes follow `wind_dir(alpha, beta)`.

use crate::analysis::{wind_dir, wind_normal, wind_side, Ctx, G};
use crate::evaluator::Evaluator;
use crate::PlaneError;
use foiley_core::linalg::{eigen, Complex};
use foiley_core::{Progress, Vec3, DEG};
use serde::{Deserialize, Serialize};

/// Which loads the stability derivatives and the T7 trim use.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum StabilityModel {
    /// potential-flow forces and moments only (pressure + induced)
    #[default]
    Inviscid,
    /// the reduced loads of the operating point: viscous strip drag (with
    /// the strips' local velocity under rotation) and section Cm, extra
    /// drag, fuselage friction
    FullLoad,
}

/// Speed at which the potential-flow lift at unit speed `raw1` equals the
/// weight (|lift| so the Cm curve of the trim scan stays continuous).
fn weight_speed(ctx: &Ctx, raw1: &crate::analysis::RawLoads, alpha: f64) -> Option<f64> {
    let l1 = (raw1.force + raw1.induced_force).dot(wind_normal(alpha)).abs() * ctx.def.fluid.rho;
    (l1 > 1e-12 && ctx.inertia.mass > 0.0).then(|| (ctx.inertia.mass * G / l1).sqrt())
}

/// Pitching moment coefficient about the CoG: inviscid at unit speed, or
/// (`FullLoad`) the reduced Cm at the speed where lift equals weight, which
/// sets the Reynolds numbers of the viscous strips.
fn cm_at(ctx: &Ctx, ev: &Evaluator, alpha: f64, beta: f64, prog: &Progress) -> Result<f64, PlaneError> {
    let dir = wind_dir(alpha, beta);
    let (raw, _) = ev.loads(dir, 1.0, Vec3::ZERO, prog)?;
    match ctx.def.stability_model {
        StabilityModel::Inviscid => Ok(raw.moment.dot(wind_side(alpha, beta)) / (0.5 * ctx.s_ref * ctx.c_ref)),
        StabilityModel::FullLoad => {
            let u = weight_speed(ctx, &raw, alpha).unwrap_or(ctx.def.velocity);
            let (raw, _) = ev.loads(dir, u, Vec3::ZERO, prog)?;
            Ok(ctx.reduce(ev.mesh(), &raw, ev.sys.cog(), alpha, beta, u, ctx.inertia.mass, None).cm)
        }
    }
}

/// Find the trimmed angle of attack (Cm = 0, with the definition's
/// `stability_model`) and the speed at which lift equals weight.  Returns
/// (alpha deg, speed m/s), or None when there is no trim; `Err` only for a solver failure or cancellation.
pub fn trim(ctx: &Ctx, ev: &Evaluator, beta: f64, prog: &Progress) -> Result<Option<(f64, f64)>, PlaneError> {
    let cm_at = |a: f64| cm_at(ctx, ev, a, beta, prog);
    // coarse scan for a stable crossing (Cm decreasing through zero)
    let mut prev: Option<(f64, f64)> = None;
    let mut bracket = None;
    let mut a = -10.0;
    while a <= 30.0 + 1e-9 {
        let c = cm_at(a)?;
        if let Some((pa, pc)) = prev {
            if pc * c <= 0.0 {
                bracket = Some((pa, pc, a, c));
                if pc > c {
                    break; // stable trim found
                }
            }
        }
        prev = Some((a, c));
        a += 2.0;
    }
    if bracket.is_none() {
        let c0 = cm_at(-30.0)?;
        let c1 = cm_at(-10.0)?;
        if c0 * c1 <= 0.0 {
            bracket = Some((-30.0, c0, -10.0, c1));
        }
    }
    let Some((mut a0, mut c0, mut a1, mut c1)) = bracket else {
        return Ok(None);
    };
    // Illinois regula falsi
    let mut side = 0;
    let mut am = a0;
    for _ in 0..60 {
        am = (a0 * c1 - a1 * c0) / (c1 - c0);
        let cm = cm_at(am)?;
        if cm.abs() < 1e-8 || (a1 - a0).abs() < 1e-7 {
            break;
        }
        if cm * c1 > 0.0 {
            a1 = am;
            c1 = cm;
            if side == -1 {
                c0 *= 0.5;
            }
            side = -1;
        } else {
            a0 = am;
            c0 = cm;
            if side == 1 {
                c1 *= 0.5;
            }
            side = 1;
        }
    }
    let (raw, _) = ev.loads(wind_dir(am, beta), 1.0, Vec3::ZERO, prog)?;
    let l1 = (raw.force + raw.induced_force).dot(wind_normal(am)) * ctx.def.fluid.rho;
    if l1 <= 1e-12 || ctx.inertia.mass <= 0.0 {
        return Ok(None);
    }
    Ok(Some((am, (ctx.inertia.mass * G / l1).sqrt())))
}

/// Dimensional derivatives in stability axes (forces N per m/s or per rad/s).
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Derivatives {
    pub xu: f64,
    pub xw: f64,
    pub xq: f64,
    pub zu: f64,
    pub zw: f64,
    pub zq: f64,
    pub mu: f64,
    pub mw: f64,
    pub mq: f64,
    pub yv: f64,
    pub yp: f64,
    pub yr: f64,
    pub lv: f64,
    pub lp: f64,
    pub lr: f64,
    pub nv: f64,
    pub np: f64,
    pub nr: f64,
    // coupling between the longitudinal and lateral sets: zero (to
    // rounding) only for a mirror-symmetric aircraft in symmetric flight
    pub xv: f64,
    pub xp: f64,
    pub xr: f64,
    pub zv: f64,
    pub zp: f64,
    pub zr: f64,
    pub mv: f64,
    pub mp: f64,
    pub mr: f64,
    pub yu: f64,
    pub yw: f64,
    pub yq: f64,
    pub lu: f64,
    pub lw: f64,
    pub lq: f64,
    pub nu: f64,
    pub nw: f64,
    pub nq: f64,
}

/// Non-dimensional stability derivatives (per radian).
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct NonDim {
    pub cxu: f64,
    pub czu: f64,
    pub cmu: f64,
    pub cxa: f64,
    pub cza: f64,
    pub cma: f64,
    pub cxq: f64,
    pub czq: f64,
    pub cmq: f64,
    pub cyb: f64,
    pub cyp: f64,
    pub cyr: f64,
    pub clb: f64,
    pub clp: f64,
    pub clr: f64,
    pub cnb: f64,
    pub cnp: f64,
    pub cnr: f64,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Mode {
    pub name: String,
    pub eigenvalue: Complex,
    /// undamped natural frequency, Hz
    pub freq_hz: f64,
    /// damping ratio
    pub damping: f64,
    /// time to half (stable) or double (unstable) amplitude, s
    pub t_half_or_double: f64,
    pub stable: bool,
    /// eigenvector of the coupled matrix, states [u w q theta v p r phi]
    pub vector: Vec<Complex>,
    /// classified as longitudinal (`longitudinal_fraction` > 0.5)
    pub longitudinal: bool,
    /// share of the eigenvector's energy in [u w q theta] (velocities over
    /// u0, rates over |lambda|, angles as they are)
    #[serde(default)]
    pub longitudinal_fraction: f64,
}

/// States of the coupled matrix `StabilityResult::a`.
pub const STATES: [&str; 8] = ["u", "w", "q", "theta", "v", "p", "r", "phi"];

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct StabilityResult {
    pub alpha: f64,
    /// sideslip of the reference state, degrees
    #[serde(default)]
    pub beta: f64,
    pub u0: f64,
    pub dim: Derivatives,
    pub nondim: NonDim,
    /// neutral point x, m
    pub x_np: f64,
    /// static margin, fraction of the reference chord (positive = stable)
    pub static_margin: f64,
    /// coupled 8x8 state matrix, states [u w q theta v p r phi] (`STATES`)
    #[serde(default)]
    pub a: [[f64; 8]; 8],
    /// longitudinal block of `a` ([u w q theta])
    pub a_long: [[f64; 4]; 4],
    /// lateral block of `a` ([v p r phi])
    pub a_lat: [[f64; 4]; 4],
    /// modes of the coupled matrix
    pub modes: Vec<Mode>,
    /// moments of inertia in stability axes (Ixx, Iyy, Izz, Izx)
    pub inertia_s: [f64; 4],
    /// control derivatives per unit of the control parameter:
    /// (X, Y, Z, L, M, N) in stability axes
    #[serde(default)]
    pub control: Option<[f64; 6]>,
    /// control input column of the coupled state equations
    #[serde(default)]
    pub b: [f64; 8],
    /// control input columns of the 4x4 blocks
    #[serde(default)]
    pub b_long: [f64; 4],
    #[serde(default)]
    pub b_lat: [f64; 4],
}

/// Force and moment (N, N.m about the CoG, body axes) for air velocity
/// `vinf` and body rates `omega`, under the definition's stability model.
/// The raw loads always come from the evaluator (unit solutions); `FullLoad`
/// adds one `reduce` per call.
fn forces(ctx: &Ctx, ev: &Evaluator, vinf: Vec3, omega: Vec3, prog: &Progress) -> Result<(Vec3, Vec3), PlaneError> {
    let (raw, _) = ev.loads_at(vinf, omega, prog)?;
    let rho = ctx.def.fluid.rho;
    match ctx.def.stability_model {
        StabilityModel::Inviscid => Ok(((raw.force + raw.induced_force) * rho, raw.moment * rho)),
        StabilityModel::FullLoad => {
            let v = vinf.norm();
            let alpha = vinf.z.atan2(vinf.x) / DEG;
            let beta = -(vinf.y / v).clamp(-1.0, 1.0).asin() / DEG;
            let (_, f, m) = ctx.reduce_full(ev.mesh(), &raw, ev.sys.cog(), alpha, beta, v, ctx.inertia.mass, Some(omega));
            Ok((f, m))
        }
    }
}

/// Stability axes of the reference state (alpha, beta): x_s against the
/// relative wind `wind_dir(alpha, beta)`, z_s in the plane of symmetry,
/// y_s = z_s x x_s (`wind_side`).
pub fn axes(alpha: f64, beta: f64) -> (Vec3, Vec3, Vec3) {
    (-wind_dir(alpha, beta), wind_side(alpha, beta), -wind_normal(alpha))
}

/// The inertia tensor in the axes (is, js, ks) and its inverse.
fn inertia_axes(ctx: &Ctx, is: Vec3, js: Vec3, ks: Vec3) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let i = &ctx.inertia;
    // the inertia tensor (products with minus signs)
    let jm = [[i.ixx, -i.ixy, -i.ixz], [-i.ixy, i.iyy, -i.iyz], [-i.ixz, -i.iyz, i.izz]];
    let e = [is, js, ks];
    let j: [[f64; 3]; 3] = std::array::from_fn(|r| {
        std::array::from_fn(|c| {
            let (a, b) = ([e[r].x, e[r].y, e[r].z], [e[c].x, e[c].y, e[c].z]);
            let mut s = 0.0;
            for p in 0..3 {
                for q in 0..3 {
                    s += a[p] * jm[p][q] * b[q];
                }
            }
            s
        })
    });
    let det = j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1]) - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0]) + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
    let cof = |r: usize, c: usize| {
        let (r1, r2) = ((r + 1) % 3, (r + 2) % 3);
        let (c1, c2) = ((c + 1) % 3, (c + 2) % 3);
        j[r1][c1] * j[r2][c2] - j[r1][c2] * j[r2][c1]
    };
    let inv = std::array::from_fn(|r| std::array::from_fn(|c| cof(c, r) / det));
    (j, inv)
}

/// Stability derivatives, the coupled state matrix and modes at a trimmed
/// point (alpha, beta, u0).  A relaxed wake is frozen at the reference
/// state (see `Evaluator::derivative_evaluator`) unless the definition asks
/// for `relaxed_derivatives`.
#[allow(clippy::too_many_arguments)]
pub fn analyse(ctx: &Ctx, ev: &Evaluator, _cog: Vec3, alpha: f64, beta: f64, u0: f64, prog: &Progress) -> Result<StabilityResult, PlaneError> {
    let (is, js, ks) = axes(alpha, beta);
    let v0 = -is * u0;
    let du = 1e-3 * u0.max(1.0);
    let dw = 1e-3 * u0.max(1.0) / ctx.c_ref.max(1e-6);
    let m = ctx.inertia.mass;
    if m <= 0.0 {
        return Err(PlaneError::Input("mass is zero".into()));
    }
    let frozen = ev.derivative_evaluator(-is, u0, prog)?;
    let ev = frozen.as_ref().unwrap_or(ev);
    // translations: the air velocity is V0 - e*delta when the aircraft speeds up along e
    let tr = |e: Vec3| -> Result<(Vec3, Vec3), PlaneError> {
        let (fm, mm) = forces(ctx, ev, v0 - e * du, Vec3::ZERO, prog)?;
        let (fp, mp) = forces(ctx, ev, v0 + e * du, Vec3::ZERO, prog)?;
        Ok(((fm - fp) / (2.0 * du), (mm - mp) / (2.0 * du)))
    };
    let (f_u, m_u) = tr(is)?;
    let (f_v, m_v) = tr(js)?;
    let (f_w, m_w) = tr(ks)?;
    let rot = |e: Vec3| -> Result<(Vec3, Vec3), PlaneError> {
        let (fp, mp) = forces(ctx, ev, v0, e * dw, prog)?;
        let (fm, mm) = forces(ctx, ev, v0, -e * dw, prog)?;
        Ok(((fp - fm) / (2.0 * dw), (mp - mm) / (2.0 * dw)))
    };
    let (f_p, m_p) = rot(is)?;
    let (f_q, m_q) = rot(js)?;
    let (f_r, m_r) = rot(ks)?;
    // X, Y, Z, L, M, N of each perturbation
    let comp = |f: Vec3, mo: Vec3| [f.dot(is), f.dot(js), f.dot(ks), mo.dot(is), mo.dot(js), mo.dot(ks)];
    let [xu, yu, zu, lu, mu, nu] = comp(f_u, m_u);
    let [xv, yv, zv, lv, mv, nv] = comp(f_v, m_v);
    let [xw, yw, zw, lw, mw, nw] = comp(f_w, m_w);
    let [xp, yp, zp, lp, mp, np] = comp(f_p, m_p);
    let [xq, yq, zq, lq, mq, nq] = comp(f_q, m_q);
    let [xr, yr, zr, lr, mr, nr] = comp(f_r, m_r);
    let d = Derivatives { xu, xw, xq, zu, zw, zq, mu, mw, mq, yv, yp, yr, lv, lp, lr, nv, np, nr, xv, xp, xr, zv, zp, zr, mv, mp, mr, yu, yw, yq, lu, lw, lq, nu, nw, nq };

    // inertia in stability axes
    let (j, ji) = inertia_axes(ctx, is, js, ks);
    let (ixx, iyy, izz, izx) = (j[0][0], j[1][1], j[2][2], -j[0][2]);
    if !(ixx > 0.0 && iyy > 0.0 && izz > 0.0) {
        return Err(PlaneError::Input("moments of inertia must be positive (add structural or point masses)".into()));
    }
    let g = G;
    // angular accelerations: J^-1 (L, M, N) for each state [u w q v p r]
    let mom = [[lu, lw, lq, lv, lp, lr], [mu, mw, mq, mv, mp, mr], [nu, nw, nq, nv, np, nr]];
    let ang = |axis: usize, s: usize| (0..3).map(|k| ji[axis][k] * mom[k][s]).sum::<f64>();
    // [u w q v p r] -> column of the 8-state vector [u w q theta v p r phi]
    let col = [0usize, 1, 2, 4, 5, 6];
    let mut a = [[0.0; 8]; 8];
    // longitudinal rows (Etkin: Xq neglected)
    let xs = [d.xu, d.xw, 0.0, d.xv, d.xp, d.xr];
    let zs = [d.zu, d.zw, d.zq + m * u0, d.zv, d.zp, d.zr];
    let ys = [d.yu, d.yw, d.yq, d.yv, d.yp, d.yr];
    for (k, &c) in col.iter().enumerate() {
        a[0][c] = xs[k] / m;
        a[1][c] = zs[k] / m;
        a[2][c] = ang(1, k);
        a[4][c] = ys[k] / m;
        a[5][c] = ang(0, k);
        a[6][c] = ang(2, k);
    }
    a[0][3] = -g;
    a[3][2] = 1.0;
    a[4][6] = d.yr / m - u0;
    a[4][7] = g;
    a[7][5] = 1.0;
    let a_long = std::array::from_fn(|r| std::array::from_fn(|c| a[r][c]));
    let a_lat = std::array::from_fn(|r| std::array::from_fn(|c| a[4 + r][4 + c]));
    // non-dimensional
    let rho = ctx.def.fluid.rho;
    let (s, b, c) = (ctx.s_ref, ctx.b_ref, ctx.c_ref);
    let q = 0.5 * rho * u0 * u0;
    let cw0 = m * g / (q * s);
    let h = 0.5 * rho * u0 * s;
    let nd = NonDim {
        cxu: d.xu / h,
        czu: (d.zu + rho * u0 * s * cw0) / h,
        cmu: d.mu / (h * c),
        cxa: d.xw / h,
        cza: d.zw / h,
        cma: d.mw / (h * c),
        cxq: d.xq / (0.25 * rho * u0 * c * s),
        czq: d.zq / (0.25 * rho * u0 * c * s),
        cmq: d.mq / (0.25 * rho * u0 * c * c * s),
        cyb: d.yv * u0 / (q * s),
        cyp: 2.0 * u0 * d.yp / (q * s * b),
        cyr: 2.0 * u0 * d.yr / (q * s * b),
        clb: d.lv * u0 / (q * s * b),
        clp: d.lp * 2.0 * u0 / b / (q * s * b),
        clr: d.lr * 2.0 * u0 / b / (q * s * b),
        cnb: d.nv * u0 / (q * s * b),
        cnp: d.np * 2.0 * u0 / b / (q * s * b),
        cnr: d.nr * 2.0 * u0 / b / (q * s * b),
    };
    let x_np = if nd.cza.abs() > 1e-12 { ctx.inertia.cog.x + nd.cma / nd.cza * c } else { f64::NAN };
    let static_margin = (x_np - ctx.inertia.cog.x) / c;
    let modes = classify(&a, u0)?;
    Ok(StabilityResult { alpha, beta, u0, dim: d, nondim: nd, x_np, static_margin, a, a_long, a_lat, modes, inertia_s: [ixx, iyy, izz, izx], control: None, b: [0.0; 8], b_long: [0.0; 4], b_lat: [0.0; 4] })
}

/// Control derivatives from a second system built with the control
/// parameter increased by `dt`: forces at the same trimmed state, both
/// through their evaluators (so with the same wake model).
pub fn add_control_derivatives(ctx: &Ctx, base: &Evaluator, plus: &Evaluator, res: &mut StabilityResult, dt: f64, prog: &Progress) -> Result<(), PlaneError> {
    let (is, js, ks) = axes(res.alpha, res.beta);
    let v0 = -is * res.u0;
    let (f0, m0) = forces(ctx, base, v0, Vec3::ZERO, prog)?;
    let (f1, m1) = forces(ctx, plus, v0, Vec3::ZERO, prog)?;
    let (df, dm) = ((f1 - f0) / dt, (m1 - m0) / dt);
    let c = [df.dot(is), df.dot(js), df.dot(ks), dm.dot(is), dm.dot(js), dm.dot(ks)];
    let m = ctx.inertia.mass;
    let (_, ji) = inertia_axes(ctx, is, js, ks);
    let ang = |axis: usize| (0..3).map(|k| ji[axis][k] * c[3 + k]).sum::<f64>();
    res.b = [c[0] / m, c[2] / m, ang(1), 0.0, c[1] / m, ang(0), ang(2), 0.0];
    res.b_long = std::array::from_fn(|k| res.b[k]);
    res.b_lat = std::array::from_fn(|k| res.b[4 + k]);
    res.control = Some(c);
    Ok(())
}

/// Shape of a control input for forced responses.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum ControlInput {
    /// constant from t = 0
    Step,
    /// held for `duration` seconds
    Pulse { duration: f64 },
    /// +1 for `duration`, then -1 for `duration`
    Doublet { duration: f64 },
    /// sinusoid of frequency `hz`
    Sine { hz: f64 },
}

impl ControlInput {
    pub fn value(&self, t: f64) -> f64 {
        match *self {
            ControlInput::Step => 1.0,
            ControlInput::Pulse { duration } => {
                if t < duration {
                    1.0
                } else {
                    0.0
                }
            }
            ControlInput::Doublet { duration } => {
                if t < duration {
                    1.0
                } else if t < 2.0 * duration {
                    -1.0
                } else {
                    0.0
                }
            }
            ControlInput::Sine { hz } => (2.0 * std::f64::consts::PI * hz * t).sin(),
        }
    }
}

/// RK4 response of dx/dt = A x + B u(t), x(0) = 0; `amp` scales the input
/// (units of the control parameter).
pub fn forced_response(a: &[[f64; 4]; 4], b: &[f64; 4], input: ControlInput, amp: f64, t_end: f64, dt: f64) -> Vec<(f64, [f64; 4])> {
    let f = |x: &[f64; 4], t: f64| -> [f64; 4] {
        let u = amp * input.value(t);
        std::array::from_fn(|i| (0..4).map(|j| a[i][j] * x[j]).sum::<f64>() + b[i] * u)
    };
    let mut out = Vec::new();
    let mut x = [0.0; 4];
    let mut t = 0.0;
    out.push((t, x));
    let n = ((t_end / dt).ceil() as usize).min(200_000);
    for _ in 0..n {
        let k1 = f(&x, t);
        let x2: [f64; 4] = std::array::from_fn(|i| x[i] + 0.5 * dt * k1[i]);
        let k2 = f(&x2, t + 0.5 * dt);
        let x3: [f64; 4] = std::array::from_fn(|i| x[i] + 0.5 * dt * k2[i]);
        let k3 = f(&x3, t + 0.5 * dt);
        let x4: [f64; 4] = std::array::from_fn(|i| x[i] + dt * k3[i]);
        let k4 = f(&x4, t + dt);
        x = std::array::from_fn(|i| x[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]));
        t += dt;
        if x.iter().any(|v| !v.is_finite() || v.abs() > 1e10) {
            break;
        }
        out.push((t, x));
    }
    out
}

fn mode_of(name: &str, l: Complex, v: Vec<Complex>, long: bool, fraction: f64) -> Mode {
    let wn = l.abs();
    let damping = if wn > 0.0 { -l.re / wn } else { 0.0 };
    let t = if l.re.abs() > 1e-12 { std::f64::consts::LN_2 / l.re.abs() } else { f64::INFINITY };
    Mode { name: name.into(), eigenvalue: l, freq_hz: wn / (2.0 * std::f64::consts::PI), damping, t_half_or_double: t, stable: l.re < 0.0, vector: v, longitudinal: long, longitudinal_fraction: fraction }
}

/// Share of the energy of eigenvector `x` (eigenvalue `l`) in the
/// longitudinal states [u w q theta]: velocities over `u0`, rates over
/// |lambda| (so that they compare with the angles), angles as they are.
pub fn longitudinal_fraction(x: &[Complex], l: Complex, u0: f64) -> f64 {
    let (u0, lm) = (u0.abs().max(1e-12), l.abs().max(1e-9));
    let scale = [1.0 / u0, 1.0 / u0, 1.0 / lm, 1.0, 1.0 / u0, 1.0 / lm, 1.0 / lm, 1.0];
    let e: Vec<f64> = x.iter().zip(scale).map(|(c, s)| (c.abs() * s).powi(2)).collect();
    let tot: f64 = e.iter().sum();
    if tot > 0.0 {
        e[..4].iter().sum::<f64>() / tot
    } else {
        0.0
    }
}

/// Eigen-decompose the coupled 8x8 matrix, sort the modes into
/// longitudinal and lateral by their eigenvector energy and name them.
pub fn classify(a: &[[f64; 8]; 8], u0: f64) -> Result<Vec<Mode>, PlaneError> {
    let rows: Vec<Vec<f64>> = a.iter().map(|r| r.to_vec()).collect();
    let (vals, vecs) = eigen(&rows)?;
    let frac: Vec<f64> = (0..vals.len()).map(|k| longitudinal_fraction(&vecs[k], vals[k], u0)).collect();
    let (long, lat): (Vec<usize>, Vec<usize>) = (0..vals.len()).partition(|&k| frac[k] > 0.5);
    let mode = |name: &str, k: usize, lg: bool| mode_of(name, vals[k], vecs[k].clone(), lg, frac[k]);
    let mut modes = Vec::new();
    // longitudinal: complex pairs (Im > 0 representatives) sorted by |lambda|
    let mut lp: Vec<usize> = long.iter().copied().filter(|&k| vals[k].im > 1e-9).collect();
    lp.sort_by(|&a, &b| vals[b].abs().total_cmp(&vals[a].abs()));
    let mut lr: Vec<usize> = long.iter().copied().filter(|&k| vals[k].im.abs() <= 1e-9).collect();
    lr.sort_by(|&a, &b| vals[b].abs().total_cmp(&vals[a].abs()));
    match lp.len() {
        2 => {
            modes.push(mode("Short period", lp[0], true));
            modes.push(mode("Phugoid", lp[1], true));
        }
        1 => {
            // the larger-magnitude pair is the short period unless both real roots are larger
            let big_real = lr.first().map(|&k| vals[k].abs() > vals[lp[0]].abs()).unwrap_or(false);
            if big_real {
                for (n, &k) in lr.iter().enumerate() {
                    modes.push(mode(&format!("Short period (real {})", n + 1), k, true));
                }
                modes.push(mode("Phugoid", lp[0], true));
            } else {
                modes.push(mode("Short period", lp[0], true));
                for (n, &k) in lr.iter().enumerate() {
                    modes.push(mode(&format!("Phugoid (real {})", n + 1), k, true));
                }
            }
        }
        _ => {
            for (n, &k) in lp.iter().enumerate() {
                modes.push(mode(&format!("Longitudinal oscillation {}", n + 1), k, true));
            }
            for (n, &k) in lr.iter().enumerate() {
                modes.push(mode(&format!("Longitudinal real {}", n + 1), k, true));
            }
        }
    }
    // lateral
    let tp: Vec<usize> = lat.iter().copied().filter(|&k| vals[k].im > 1e-9).collect();
    let mut tr: Vec<usize> = lat.iter().copied().filter(|&k| vals[k].im.abs() <= 1e-9).collect();
    for &k in &tp {
        modes.push(mode("Dutch roll", k, false));
    }
    tr.sort_by(|&a, &b| vals[a].re.total_cmp(&vals[b].re));
    if !tr.is_empty() {
        let roll = tr[0];
        modes.push(mode("Roll damping", roll, false));
        if let Some(&sp) = tr[1..].iter().min_by(|&&a, &&b| vals[a].re.abs().total_cmp(&vals[b].re.abs())) {
            modes.push(mode("Spiral", sp, false));
        }
    }
    Ok(modes)
}

/// RK4 time response of dx/dt = A x from x0 (free response).
pub fn time_response(a: &[[f64; 4]; 4], x0: [f64; 4], t_end: f64, dt: f64) -> Vec<(f64, [f64; 4])> {
    let f = |x: &[f64; 4]| -> [f64; 4] { std::array::from_fn(|i| (0..4).map(|j| a[i][j] * x[j]).sum()) };
    let mut out = Vec::new();
    let mut x = x0;
    let mut t = 0.0;
    let n = ((t_end / dt).ceil() as usize).min(200_000);
    out.push((t, x));
    for _ in 0..n {
        let k1 = f(&x);
        let x2: [f64; 4] = std::array::from_fn(|i| x[i] + 0.5 * dt * k1[i]);
        let k2 = f(&x2);
        let x3: [f64; 4] = std::array::from_fn(|i| x[i] + 0.5 * dt * k2[i]);
        let k3 = f(&x3);
        let x4: [f64; 4] = std::array::from_fn(|i| x[i] + dt * k3[i]);
        let k4 = f(&x4);
        x = std::array::from_fn(|i| x[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]));
        t += dt;
        if x.iter().any(|v| !v.is_finite() || v.abs() > 1e10) {
            break;
        }
        out.push((t, x));
    }
    out
}

/// Initial state that excites a mode (real part of its eigenvector,
/// normalised so the largest component is `amp`).
/// The 4 components are those of the mode's block (`a_long` for a
/// longitudinal mode, `a_lat` otherwise) of the coupled eigenvector.
pub fn modal_initial_state(m: &Mode, amp: f64) -> [f64; 4] {
    let off = if m.vector.len() == 8 && !m.longitudinal { 4 } else { 0 };
    let mut x: [f64; 4] = std::array::from_fn(|i| m.vector.get(off + i).map(|c| c.re).unwrap_or(0.0));
    let mx = x.iter().fold(0.0f64, |a, v| a.max(v.abs()));
    if mx > 0.0 {
        for v in &mut x {
            *v *= amp / mx;
        }
    }
    x
}
