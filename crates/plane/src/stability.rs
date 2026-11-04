//! Flight dynamics: trim, stability derivatives, state matrices, modes and
//! time responses (Etkin & Reid small-perturbation equations in stability
//! axes; level flight, theta0 = 0).

use crate::analysis::{wind_dir, wind_normal, wind_side, Ctx, G};
use crate::evaluator::Evaluator;
use crate::PlaneError;
use foiley_core::linalg::{eigen, Complex};
use foiley_core::{Progress, Vec3, DEG};
use serde::{Deserialize, Serialize};

/// Pitching moment coefficient (inviscid, about the CoG) at unit speed.
fn cm_at(ctx: &Ctx, ev: &Evaluator, alpha: f64, beta: f64, prog: &Progress) -> Result<f64, PlaneError> {
    let (raw, _) = ev.loads(wind_dir(alpha, beta), 1.0, Vec3::ZERO, prog)?;
    Ok(raw.moment.dot(wind_side(alpha, beta)) / (0.5 * ctx.s_ref * ctx.c_ref))
}

/// Find the trimmed angle of attack (Cm = 0) and the speed at which lift
/// equals weight.  Returns (alpha deg, speed m/s), or None when there is no
/// trim; `Err` only for a solver failure or cancellation.
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
    /// state-vector components of the eigenvector
    pub vector: Vec<Complex>,
    pub longitudinal: bool,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct StabilityResult {
    pub alpha: f64,
    pub u0: f64,
    pub dim: Derivatives,
    pub nondim: NonDim,
    /// neutral point x, m
    pub x_np: f64,
    /// static margin, fraction of the reference chord (positive = stable)
    pub static_margin: f64,
    pub a_long: [[f64; 4]; 4],
    pub a_lat: [[f64; 4]; 4],
    pub modes: Vec<Mode>,
    /// moments of inertia in stability axes (Ixx, Iyy, Izz, Izx)
    pub inertia_s: [f64; 4],
    /// control derivatives per unit of the control parameter:
    /// (X, Y, Z, L, M, N) in stability axes
    #[serde(default)]
    pub control: Option<[f64; 6]>,
    /// control input columns of the state equations
    #[serde(default)]
    pub b_long: [f64; 4],
    #[serde(default)]
    pub b_lat: [f64; 4],
}

fn forces(ctx: &Ctx, ev: &Evaluator, vinf: Vec3, omega: Vec3, prog: &Progress) -> Result<(Vec3, Vec3), PlaneError> {
    let (raw, _) = ev.loads_at(vinf, omega, prog)?;
    let rho = ctx.def.fluid.rho;
    Ok(((raw.force + raw.induced_force) * rho, raw.moment * rho))
}

/// Stability derivatives, state matrices and modes at a trimmed point.
/// A relaxed wake is frozen at the reference state (see
/// `Evaluator::derivative_evaluator`) unless the definition asks for
/// `relaxed_derivatives`.
pub fn analyse(ctx: &Ctx, ev: &Evaluator, _cog: Vec3, alpha: f64, u0: f64, prog: &Progress) -> Result<StabilityResult, PlaneError> {
    let (sa, ca) = (alpha * DEG).sin_cos();
    let is = Vec3::new(-ca, 0.0, -sa);
    let js = Vec3::Y;
    let ks = Vec3::new(sa, 0.0, -ca);
    let v0 = -is * u0;
    let du = 1e-3 * u0.max(1.0);
    let dw = 1e-3 * u0.max(1.0) / ctx.c_ref.max(1e-6);
    let m = ctx.inertia.mass;
    if m <= 0.0 {
        return Err(PlaneError::Input("mass is zero".into()));
    }
    let mut d = Derivatives::default();
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
    d.xu = f_u.dot(is);
    d.zu = f_u.dot(ks);
    d.mu = m_u.dot(js);
    d.xw = f_w.dot(is);
    d.zw = f_w.dot(ks);
    d.mw = m_w.dot(js);
    d.xq = f_q.dot(is);
    d.zq = f_q.dot(ks);
    d.mq = m_q.dot(js);
    d.yv = f_v.dot(js);
    d.lv = m_v.dot(is);
    d.nv = m_v.dot(ks);
    d.yp = f_p.dot(js);
    d.lp = m_p.dot(is);
    d.np = m_p.dot(ks);
    d.yr = f_r.dot(js);
    d.lr = m_r.dot(is);
    d.nr = m_r.dot(ks);

    // inertia in stability axes
    let i = &ctx.inertia;
    let j = |a: Vec3, b: Vec3| -> f64 {
        // a^T J b with J the inertia tensor (products with minus signs)
        let jm = [[i.ixx, -i.ixy, -i.ixz], [-i.ixy, i.iyy, -i.iyz], [-i.ixz, -i.iyz, i.izz]];
        let av = [a.x, a.y, a.z];
        let bv = [b.x, b.y, b.z];
        let mut s = 0.0;
        for r in 0..3 {
            for c in 0..3 {
                s += av[r] * jm[r][c] * bv[c];
            }
        }
        s
    };
    let ixx = j(is, is);
    let iyy = j(js, js);
    let izz = j(ks, ks);
    let izx = -j(is, ks);
    if !(ixx > 0.0 && iyy > 0.0 && izz > 0.0) {
        return Err(PlaneError::Input("moments of inertia must be positive (add structural or point masses)".into()));
    }
    let g = G;
    let a_long = [
        [d.xu / m, d.xw / m, 0.0, -g],
        [d.zu / m, d.zw / m, (d.zq + m * u0) / m, 0.0],
        [d.mu / iyy, d.mw / iyy, d.mq / iyy, 0.0],
        [0.0, 0.0, 1.0, 0.0],
    ];
    let den = ixx * izz - izx * izx;
    let ixp = den / izz;
    let izp = den / ixx;
    let izxp = izx / den;
    let a_lat = [
        [d.yv / m, d.yp / m, d.yr / m - u0, g],
        [d.lv / ixp + izxp * d.nv, d.lp / ixp + izxp * d.np, d.lr / ixp + izxp * d.nr, 0.0],
        [izxp * d.lv + d.nv / izp, izxp * d.lp + d.np / izp, izxp * d.lr + d.nr / izp, 0.0],
        [0.0, 1.0, 0.0, 0.0],
    ];
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
    let modes = classify(&a_long, &a_lat)?;
    Ok(StabilityResult { alpha, u0, dim: d, nondim: nd, x_np, static_margin, a_long, a_lat, modes, inertia_s: [ixx, iyy, izz, izx], control: None, b_long: [0.0; 4], b_lat: [0.0; 4] })
}

/// Control derivatives from a second system built with the control
/// parameter increased by `dt`: forces at the same trimmed state, both
/// through their evaluators (so with the same wake model).
pub fn add_control_derivatives(ctx: &Ctx, base: &Evaluator, plus: &Evaluator, res: &mut StabilityResult, dt: f64, prog: &Progress) -> Result<(), PlaneError> {
    let (sa, ca) = (res.alpha * DEG).sin_cos();
    let is = Vec3::new(-ca, 0.0, -sa);
    let js = Vec3::Y;
    let ks = Vec3::new(sa, 0.0, -ca);
    let v0 = -is * res.u0;
    let (f0, m0) = forces(ctx, base, v0, Vec3::ZERO, prog)?;
    let (f1, m1) = forces(ctx, plus, v0, Vec3::ZERO, prog)?;
    let (df, dm) = ((f1 - f0) / dt, (m1 - m0) / dt);
    let c = [df.dot(is), df.dot(js), df.dot(ks), dm.dot(is), dm.dot(js), dm.dot(ks)];
    let m = ctx.inertia.mass;
    let [ixx, iyy, izz, izx] = res.inertia_s;
    let den = ixx * izz - izx * izx;
    let (ixp, izp, izxp) = (den / izz, den / ixx, izx / den);
    res.b_long = [c[0] / m, c[2] / m, c[4] / iyy, 0.0];
    res.b_lat = [c[1] / m, c[3] / ixp + c[5] * izxp, c[3] * izxp + c[5] / izp, 0.0];
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

fn mode_of(name: &str, l: Complex, v: Vec<Complex>, long: bool) -> Mode {
    let wn = l.abs();
    let damping = if wn > 0.0 { -l.re / wn } else { 0.0 };
    let t = if l.re.abs() > 1e-12 { std::f64::consts::LN_2 / l.re.abs() } else { f64::INFINITY };
    Mode { name: name.into(), eigenvalue: l, freq_hz: wn / (2.0 * std::f64::consts::PI), damping, t_half_or_double: t, stable: l.re < 0.0, vector: v, longitudinal: long }
}

/// Eigen-decompose both 4x4 systems and name the modes.
pub fn classify(a_long: &[[f64; 4]; 4], a_lat: &[[f64; 4]; 4]) -> Result<Vec<Mode>, PlaneError> {
    let rows = |a: &[[f64; 4]; 4]| a.iter().map(|r| r.to_vec()).collect::<Vec<_>>();
    let (vl, xl) = eigen(&rows(a_long))?;
    let (vt, xt) = eigen(&rows(a_lat))?;
    let mut modes = Vec::new();
    // longitudinal: complex pairs (Im > 0 representatives) sorted by |lambda|
    let mut lp: Vec<usize> = (0..4).filter(|&k| vl[k].im > 1e-9).collect();
    lp.sort_by(|&a, &b| vl[b].abs().total_cmp(&vl[a].abs()));
    let mut lr: Vec<usize> = (0..4).filter(|&k| vl[k].im.abs() <= 1e-9).collect();
    lr.sort_by(|&a, &b| vl[b].abs().total_cmp(&vl[a].abs()));
    match lp.len() {
        2 => {
            modes.push(mode_of("Short period", vl[lp[0]], xl[lp[0]].clone(), true));
            modes.push(mode_of("Phugoid", vl[lp[1]], xl[lp[1]].clone(), true));
        }
        1 => {
            // the larger-magnitude pair is the short period unless both real roots are larger
            let big_real = lr.first().map(|&k| vl[k].abs() > vl[lp[0]].abs()).unwrap_or(false);
            if big_real {
                for (n, &k) in lr.iter().enumerate() {
                    modes.push(mode_of(&format!("Short period (real {})", n + 1), vl[k], xl[k].clone(), true));
                }
                modes.push(mode_of("Phugoid", vl[lp[0]], xl[lp[0]].clone(), true));
            } else {
                modes.push(mode_of("Short period", vl[lp[0]], xl[lp[0]].clone(), true));
                for (n, &k) in lr.iter().enumerate() {
                    modes.push(mode_of(&format!("Phugoid (real {})", n + 1), vl[k], xl[k].clone(), true));
                }
            }
        }
        _ => {
            for (n, &k) in lr.iter().enumerate() {
                modes.push(mode_of(&format!("Longitudinal real {}", n + 1), vl[k], xl[k].clone(), true));
            }
        }
    }
    // lateral
    let tp: Vec<usize> = (0..4).filter(|&k| vt[k].im > 1e-9).collect();
    let mut tr: Vec<usize> = (0..4).filter(|&k| vt[k].im.abs() <= 1e-9).collect();
    for &k in &tp {
        modes.push(mode_of("Dutch roll", vt[k], xt[k].clone(), false));
    }
    tr.sort_by(|&a, &b| vt[a].re.total_cmp(&vt[b].re));
    // drop the zero root (heading-like) if present, keep the meaningful ones
    if !tr.is_empty() {
        let roll = tr[0];
        modes.push(mode_of("Roll damping", vt[roll], xt[roll].clone(), false));
        if let Some(&sp) = tr[1..].iter().min_by(|&&a, &&b| vt[a].re.abs().total_cmp(&vt[b].re.abs())) {
            modes.push(mode_of("Spiral", vt[sp], xt[sp].clone(), false));
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
pub fn modal_initial_state(m: &Mode, amp: f64) -> [f64; 4] {
    let mut x: [f64; 4] = std::array::from_fn(|i| m.vector.get(i).map(|c| c.re).unwrap_or(0.0));
    let mx = x.iter().fold(0.0f64, |a, v| a.max(v.abs()));
    if mx > 0.0 {
        for v in &mut x {
            *v *= amp / mx;
        }
    }
    x
}
