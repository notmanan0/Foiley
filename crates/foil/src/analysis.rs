//! 2D analysis: XFoil operating points and polars with progress reporting,
//! cancellation, panic isolation and automatic convergence recovery.
//!
//! Recovery strategy for a point that fails to converge (or produces NaN,
//! or panics inside the engine):
//!
//! 1. re-initialise the boundary layer from the inviscid solution and
//!    retry with twice the iteration limit;
//! 2. march to the target from the last converged point (or from the
//!    nearest easy point) in 2, then 4 sub-steps;
//! 3. give up on this point, record it in `FoilPolar::failed`, rebuild a
//!    clean engine state and carry on with the next point.
//!
//! Sweeps start at the point nearest zero lift and march outwards in both
//! directions, so every point warm-starts from a converged neighbour on the
//! attached-flow branch.  Long sweeps are split into independent chunks that
//! run on separate threads (each chunk again marching away from zero).

use crate::foil::{Foil, Paneling};
use crate::polar::{FoilPoint, FoilPolar};
use crate::FoilError;
use foiley::sweep::{polar_point, solve_point};
use foiley::xfoil::{XFoil, DTOR};
use foiley_core::Progress;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::panic::{catch_unwind, AssertUnwindSafe};

/// How Re and Mach vary along a polar (XFoil polar types plus xflr5's
/// type 4 fixed-alpha and type 6 flap-sweep polars).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum PolarType {
    /// Type 1: fixed Re and Mach (fixed chord and speed).
    #[default]
    FixedSpeed,
    /// Type 2: Re*sqrt(CL) and M*sqrt(CL) fixed (fixed wing loading).
    FixedLift,
    /// Type 3: Re*CL fixed (rubber chord).
    RubberChord,
    /// Type 4: fixed alpha, sweep Re.
    FixedAlpha,
    /// Type 6: fixed alpha and Re, sweep trailing-edge flap deflection.
    FlapSweep,
}

impl PolarType {
    pub fn number(self) -> u8 {
        match self {
            PolarType::FixedSpeed => 1,
            PolarType::FixedLift => 2,
            PolarType::RubberChord => 3,
            PolarType::FixedAlpha => 4,
            PolarType::FlapSweep => 6,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            PolarType::FixedSpeed => "Type 1 - fixed speed (Re, M constant)",
            PolarType::FixedLift => "Type 2 - fixed lift (Re*sqrt(CL) constant)",
            PolarType::RubberChord => "Type 3 - rubber chord (Re*CL constant)",
            PolarType::FixedAlpha => "Type 4 - fixed alpha, sweep Re",
            PolarType::FlapSweep => "Type 6 - fixed alpha, sweep flap angle",
        }
    }
    pub const ALL: [PolarType; 5] = [PolarType::FixedSpeed, PolarType::FixedLift, PolarType::RubberChord, PolarType::FixedAlpha, PolarType::FlapSweep];
}

/// What to sweep.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Sweep {
    /// angle of attack, degrees (first, last, step)
    Alpha { from: f64, to: f64, step: f64 },
    /// lift coefficient
    Cl { from: f64, to: f64, step: f64 },
    /// Reynolds number at fixed alpha (type 4)
    Re { from: f64, to: f64, step: f64 },
    /// flap deflection in degrees at fixed alpha (type 6)
    Flap { from: f64, to: f64, step: f64 },
    /// explicit list of angles of attack
    AlphaList(Vec<f64>),
}

impl Sweep {
    /// The spec values, in sweep order.
    pub fn values(&self) -> Vec<f64> {
        let range = |a: f64, b: f64, s: f64| -> Vec<f64> {
            let s = if s == 0.0 { 1.0 } else { s.abs() };
            let n = ((b - a).abs() / s + 0.5).floor() as usize;
            let dir = if b >= a { 1.0 } else { -1.0 };
            (0..=n.min(100_000)).map(|i| a + dir * s * i as f64).collect()
        };
        match self {
            Sweep::Alpha { from, to, step } | Sweep::Cl { from, to, step } | Sweep::Re { from, to, step } | Sweep::Flap { from, to, step } => range(*from, *to, *step),
            Sweep::AlphaList(v) => v.clone(),
        }
    }
    pub fn is_alpha(&self) -> bool {
        matches!(self, Sweep::Alpha { .. } | Sweep::AlphaList(_))
    }
}

/// All 2D analysis parameters.  Stored with every polar so results are
/// reproducible from the polar alone.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct AnalysisSettings {
    pub polar_type: PolarType,
    /// Reynolds number (type 1); Re*sqrt(CL) (type 2); Re*CL (type 3)
    pub re: f64,
    /// Mach number (or M*sqrt(CL) for type 2)
    pub mach: f64,
    /// e^N transition criterion (9 = average wind tunnel, 11-12 = quiet free air, 4-6 = turbulent)
    pub ncrit: f64,
    /// forced transition (trip) location x/c, upper surface (1 = free)
    pub xtr_top: f64,
    /// forced transition location x/c, lower surface (1 = free)
    pub xtr_bot: f64,
    /// viscous (boundary layer) or inviscid analysis
    pub viscous: bool,
    /// fixed alpha for type 4 / type 6 polars, degrees
    pub alpha_fixed: f64,
    /// flap hinge for type 6 polars (x/c, y/t)
    pub flap_hinge: (f64, f64),
    /// Newton iteration limit per point (XFoil ITER, default 100)
    pub max_iter: usize,
    /// viscous acceleration parameter (XFoil VACC, default 0.01)
    pub vaccel: f64,
    /// paneling
    pub paneling: Paneling,
    /// recover unconverged points (re-init BL, sub-stepping)
    pub recover: bool,
    /// start sweeps from the point nearest zero lift and march outwards
    pub from_zero: bool,
    /// run long sweeps in parallel chunks (results then agree with the
    /// sequential sweep to the Newton tolerance, not bit-for-bit)
    pub parallel: bool,
    /// discard converged points with CD below this value (0 = keep all)
    pub min_cd: f64,
    /// moment reference point (x/c, y/c)
    pub cm_ref: (f64, f64),
    /// store Cp and boundary-layer distributions for every point
    pub keep_details: bool,
}

impl Default for AnalysisSettings {
    fn default() -> Self {
        AnalysisSettings {
            polar_type: PolarType::FixedSpeed,
            re: 1.0e6,
            mach: 0.0,
            ncrit: 9.0,
            xtr_top: 1.0,
            xtr_bot: 1.0,
            viscous: true,
            alpha_fixed: 0.0,
            flap_hinge: (0.75, 0.5),
            max_iter: 100,
            vaccel: 0.01,
            paneling: Paneling::default(),
            recover: true,
            from_zero: true,
            parallel: true,
            min_cd: 0.0,
            cm_ref: (0.25, 0.0),
            keep_details: false,
        }
    }
}

impl AnalysisSettings {
    /// Check inputs and explain problems in plain words.
    pub fn validate(&self) -> Result<(), String> {
        if self.viscous && !(self.re > 1.0e3 && self.re < 1.0e10) {
            return Err(format!("Reynolds number {} is out of the supported range 1e3..1e10", self.re));
        }
        if !(0.0..0.95).contains(&self.mach) {
            return Err(format!("Mach {} must be in 0..0.95 (XFoil is a subsonic code)", self.mach));
        }
        if !(0.1..=20.0).contains(&self.ncrit) {
            return Err(format!("Ncrit {} should be between 0.1 and 20 (9 is the usual value)", self.ncrit));
        }
        for (n, v) in [("top", self.xtr_top), ("bottom", self.xtr_bot)] {
            if !(0.0..=1.0).contains(&v) {
                return Err(format!("forced transition {n} x/c = {v} must be between 0 and 1"));
            }
        }
        if self.max_iter == 0 || self.max_iter > 10_000 {
            return Err("iteration limit must be 1..10000".into());
        }
        Ok(())
    }

    /// A descriptive name like "T1_Re1.000_M0.00_N9.0".
    pub fn polar_name(&self) -> String {
        let t = self.polar_type.number();
        let mut s = match self.polar_type {
            PolarType::FixedAlpha => format!("T4_a{:.2}_M{:.2}", self.alpha_fixed, self.mach),
            PolarType::FlapSweep => format!("T6_a{:.2}_Re{:.3}", self.alpha_fixed, self.re / 1e6),
            _ => format!("T{t}_Re{:.3}_M{:.2}", self.re / 1e6, self.mach),
        };
        if self.viscous {
            s.push_str(&format!("_N{:.1}", self.ncrit));
            if self.xtr_top < 1.0 || self.xtr_bot < 1.0 {
                s.push_str(&format!("_trip{:.2}-{:.2}", self.xtr_top, self.xtr_bot));
            }
        } else {
            s.push_str("_inviscid");
        }
        s
    }
}

/// Surface distributions of one solved operating point.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct PointDetail {
    pub point: FoilPoint,
    /// panel node x, y
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// inviscid and viscous Cp at the nodes
    pub cp_inv: Vec<f64>,
    pub cp_visc: Vec<f64>,
    /// boundary layers: upper, lower, wake
    pub top: BlSide,
    pub bot: BlSide,
    pub wake: BlSide,
}

/// Boundary-layer quantities along one side, from the stagnation point.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct BlSide {
    pub x: Vec<f64>,
    /// arc length from stagnation point
    pub s: Vec<f64>,
    /// edge velocity Ue/Vinf
    pub ue: Vec<f64>,
    /// displacement thickness / c
    pub dstar: Vec<f64>,
    /// momentum thickness / c
    pub theta: Vec<f64>,
    /// shape factor H = d*/theta
    pub h: Vec<f64>,
    /// skin friction coefficient
    pub cf: Vec<f64>,
    /// amplification ratio N (laminar) / sqrt(max shear coeff) (turbulent)
    pub n_or_ctau: Vec<f64>,
}

fn extract_detail(xf: &XFoil, p: FoilPoint) -> PointDetail {
    let n = xf.n;
    let mut d = PointDetail { point: p, x: xf.x[1..=n].to_vec(), y: xf.y[1..=n].to_vec(), cp_inv: xf.cpi[1..=n].to_vec(), ..Default::default() };
    d.cp_visc = if xf.lvisc { xf.cpv[1..=n].to_vec() } else { d.cp_inv.clone() };
    if xf.lvisc && xf.lvconv {
        let qinf2 = 0.5 * xf.qinf * xf.qinf;
        let side = |is: usize, from: usize, to: usize| -> BlSide {
            let mut b = BlSide::default();
            for ibl in from..=to {
                let ip = xf.ipan[is][ibl];
                if ip == 0 || ip > n + xf.nw {
                    continue;
                }
                b.x.push(xf.x[ip]);
                b.s.push(xf.xssi[is][ibl]);
                b.ue.push(xf.uedg[is][ibl] / xf.qinf);
                b.dstar.push(xf.dstr[is][ibl]);
                b.theta.push(xf.thet[is][ibl]);
                let th = xf.thet[is][ibl];
                b.h.push(if th > 0.0 { xf.dstr[is][ibl] / th } else { 0.0 });
                b.cf.push(xf.tau[is][ibl] / qinf2);
                b.n_or_ctau.push(xf.ctau[is][ibl]);
            }
            b
        };
        d.top = side(1, 2, xf.iblte[1]);
        d.bot = side(2, 2, xf.iblte[2]);
        if xf.nbl[2] > xf.iblte[2] {
            d.wake = side(2, xf.iblte[2] + 1, xf.nbl[2]);
        }
    }
    d
}

/// Configure a fresh engine for these settings.
fn setup_engine(foil: &Foil, st: &AnalysisSettings, re: f64) -> Result<XFoil, FoilError> {
    let mut xf = foil.engine(&st.paneling)?;
    xf.lvisc = st.viscous;
    xf.reinf1 = re;
    xf.minf1 = st.mach;
    xf.acrit = st.ncrit;
    xf.xstrip[1] = st.xtr_top;
    xf.xstrip[2] = st.xtr_bot;
    xf.itmax = st.max_iter;
    xf.vaccel = st.vaccel;
    xf.xcmref = st.cm_ref.0;
    xf.ycmref = st.cm_ref.1;
    let (ma, re_t) = match st.polar_type {
        PolarType::FixedLift => (2, 2),
        PolarType::RubberChord => (1, 3),
        _ => (1, 1),
    };
    xf.matyp = ma;
    xf.retyp = re_t;
    let mut r = 0.0;
    xf.mrcl(1.0, &mut r);
    xf.comset();
    xf.lvconv = false;
    xf.lblini = false;
    xf.lipan = false;
    Ok(xf)
}

/// Solve one point with panic isolation; returns None on failure/NaN.
fn try_solve(xf: &mut XFoil, lalfa: bool, spec: f64, iters: usize) -> Option<FoilPoint> {
    let r = catch_unwind(AssertUnwindSafe(|| solve_point(xf, lalfa, spec, iters)));
    match r {
        Ok(sp) if sp.converged => {
            let q = sp.point;
            if q.cl.is_finite() && q.cd.is_finite() && q.cm.is_finite() && q.alpha.is_finite() {
                Some(to_foil_point(xf, &sp.point))
            } else {
                None
            }
        }
        _ => None,
    }
}

fn to_foil_point(xf: &XFoil, q: &foiley::polar::PolarPoint) -> FoilPoint {
    let cdf = if xf.lvisc { xf.cdf } else { 0.0 };
    FoilPoint {
        alpha: q.alpha,
        cl: q.cl,
        cd: q.cd,
        cdp: q.cdp,
        cdf,
        cm: q.cm,
        xtr_top: if xf.lvisc { q.xtr_top } else { 1.0 },
        xtr_bot: if xf.lvisc { q.xtr_bot } else { 1.0 },
        cpmin: q.cpmin,
        re: if xf.lvisc { xf.reinf } else { 0.0 },
        mach: xf.minf,
        iterations: if xf.lvisc { xf.viscal_iters as u32 } else { 0 },
        recovered: false,
        recovery: crate::polar::Recovery::None,
        flap: 0.0,
    }
}

/// Spec value in engine units (alpha in radians).
fn eng(lalfa: bool, v: f64) -> f64 {
    if lalfa {
        v * DTOR
    } else {
        v
    }
}

/// Context shared by the points of one marching chunk.
struct Marcher<'a> {
    foil: &'a Foil,
    st: &'a AnalysisSettings,
    re: f64,
    lalfa: bool,
    xf: XFoil,
    last_ok: Option<f64>,
}

impl<'a> Marcher<'a> {
    fn new(foil: &'a Foil, st: &'a AnalysisSettings, re: f64, lalfa: bool) -> Result<Self, FoilError> {
        Ok(Marcher { foil, st, re, lalfa, xf: setup_engine(foil, st, re)?, last_ok: None })
    }

    fn reset_bl(&mut self) {
        self.xf.lblini = false;
        self.xf.lipan = false;
        self.xf.lvconv = false;
    }

    fn rebuild(&mut self) {
        if let Ok(x) = setup_engine(self.foil, self.st, self.re) {
            self.xf = x;
        }
        self.last_ok = None;
    }

    /// Solve `v` (alpha deg or CL), with recovery.  Returns the point and log notes.
    fn solve(&mut self, v: f64, prog: &Progress) -> Option<FoilPoint> {
        let it = self.st.max_iter;
        if let Some(p) = try_solve(&mut self.xf, self.lalfa, eng(self.lalfa, v), it) {
            self.last_ok = Some(v);
            return Some(p);
        }
        if !self.st.viscous || !self.st.recover {
            self.rebuild();
            return None;
        }
        // 1. fresh BL from inviscid solution, more iterations
        self.reset_bl();
        if let Some(mut p) = try_solve(&mut self.xf, self.lalfa, eng(self.lalfa, v), it * 2) {
            p.recovered = true;
            p.recovery = crate::polar::Recovery::ReinitBl;
            prog.debug(format!("recovered {} by re-initialising the boundary layer", self.fmt(v)));
            self.last_ok = Some(v);
            return Some(p);
        }
        // 2. march in sub-steps from the last converged point, or from
        //    neighbouring points that converge from a cold start
        let unit = if self.lalfa { 1.0 } else { 0.1 };
        let mut starts: Vec<f64> = Vec::new();
        if let Some(s) = self.last_ok {
            if (s - v).abs() > 1e-9 {
                starts.push(s);
            }
        }
        // positive k = towards zero lift (usually easier), negative = away
        let dir = if v < 0.0 { -1.0 } else { 1.0 };
        for k in [1.0, -1.0, 2.0, -2.0, 4.0] {
            starts.push(v - k * unit * dir);
        }
        let mut tries = 0;
        for (start, nsub) in starts.iter().flat_map(|s| [(*s, 2usize), (*s, 4)]) {
            tries += 1;
            if tries > 6 {
                break;
            }
            self.rebuild();
            let mut ok = try_solve(&mut self.xf, self.lalfa, eng(self.lalfa, start), it).is_some();
            if !ok {
                continue;
            }
            for k in 1..nsub {
                let vk = start + (v - start) * k as f64 / nsub as f64;
                if try_solve(&mut self.xf, self.lalfa, eng(self.lalfa, vk), it).is_none() {
                    ok = false;
                    break;
                }
            }
            if !ok {
                continue;
            }
            if let Some(mut p) = try_solve(&mut self.xf, self.lalfa, eng(self.lalfa, v), it) {
                p.recovered = true;
                p.recovery = crate::polar::Recovery::March { substeps: nsub as u32, attempts: tries as u32 };
                prog.debug(format!("recovered {} by marching in {nsub} sub-steps from {}", self.fmt(v), self.fmt(start)));
                self.last_ok = Some(v);
                return Some(p);
            }
        }
        // 3. give up; leave a clean state for the next point
        self.rebuild();
        None
    }

    fn fmt(&self, v: f64) -> String {
        if self.lalfa {
            format!("alpha={v:.2}")
        } else {
            format!("CL={v:.3}")
        }
    }
}

/// Order sweep values so marching starts at the value nearest zero and
/// moves outwards: returns chunks (each to be marched in order).
fn march_order(values: &[f64], from_zero: bool) -> Vec<Vec<f64>> {
    if values.is_empty() {
        return vec![];
    }
    if !from_zero {
        return vec![values.to_vec()];
    }
    let mut v = values.to_vec();
    v.sort_by(|a, b| a.total_cmp(b));
    v.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    let i0 = (0..v.len()).min_by(|&a, &b| v[a].abs().total_cmp(&v[b].abs())).unwrap();
    let up: Vec<f64> = v[i0..].to_vec();
    let down: Vec<f64> = v[..=i0].iter().rev().cloned().collect();
    let mut out = vec![up];
    if down.len() > 1 {
        out.push(down);
    }
    out
}

/// Split marching chunks further for parallel execution: a chunk
/// [a0, a1, ..., an] becomes [a0..ak], [ak..]: each later piece repeats its
/// first point, which is solved cold, recovering by marching if needed.
fn split_chunks(chunks: Vec<Vec<f64>>, pieces: usize) -> Vec<(Vec<f64>, bool)> {
    let mut out = Vec::new();
    let total: usize = chunks.iter().map(|c| c.len()).sum();
    for c in chunks {
        let k = ((pieces * c.len()) as f64 / total.max(1) as f64).round().max(1.0) as usize;
        let len = c.len().div_ceil(k).max(6);
        let mut first = true;
        let mut i = 0;
        while i < c.len() {
            let j = (i + len).min(c.len());
            out.push((c[i..j].to_vec(), first));
            first = false;
            i = j;
        }
    }
    out
}

/// Run a polar.  Progress counts one unit per operating point.
pub fn analyze(foil: &Foil, st: &AnalysisSettings, sweep: &Sweep, prog: &Progress) -> Result<(FoilPolar, Vec<PointDetail>), FoilError> {
    foil.validate()?;
    st.validate().map_err(FoilError::Geometry)?;
    foiley::set_quiet(true);
    let mut polar = FoilPolar {
        name: st.polar_name(),
        foil: foil.name.clone(),
        polar_type: st.polar_type.number(),
        re: st.re,
        mach: st.mach,
        ncrit: st.ncrit,
        xtr_top: st.xtr_top,
        xtr_bot: st.xtr_bot,
        viscous: st.viscous,
        points: vec![],
        failed: vec![],
        flap_deg: 0.0,
        settings: Some(st.clone()),
        foil_hash: Some(foil.geom_hash()),
        manifest: Some(foiley_core::provenance::RunManifest::new(foiley_core::provenance::hash_json("foiley.settings2d.v1", st), foil.geom_hash())),
        legacy: false,
    };
    let values = sweep.values();
    if values.is_empty() {
        return Ok((polar, vec![]));
    }
    let total = values.len();
    prog.stage(format!("{}: {} ({} points)", foil.name, polar.name, total));
    let done = std::sync::atomic::AtomicUsize::new(0);
    let tick = |prog: &Progress| {
        let d = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
        prog.set(d as f64 / total as f64);
    };

    let mut details = Vec::new();
    match (st.polar_type, sweep) {
        (PolarType::FixedAlpha, Sweep::Re { .. }) => {
            // march Re from high (easy) to low
            let mut res: Vec<f64> = values.clone();
            res.sort_by(|a, b| b.total_cmp(a));
            let mut m = Marcher::new(foil, st, res[0], true)?;
            for re in res {
                prog.check()?;
                m.re = re;
                m.xf.reinf1 = re;
                let mut r = 0.0;
                m.xf.mrcl(1.0, &mut r);
                m.xf.lvconv = false;
                match m.solve(st.alpha_fixed, prog) {
                    Some(p) => {
                        prog.count_ok();
                        if st.keep_details {
                            details.push(extract_detail(&m.xf, p));
                        }
                        polar.insert(p)
                    }
                    None => {
                        prog.count_failed();
                        prog.warn(format!("{}: Re={re:.0} did not converge", foil.name));
                        polar.failed.push(re);
                    }
                }
                tick(prog);
            }
        }
        (PolarType::FlapSweep, Sweep::Flap { .. }) => {
            // each flap angle is a different geometry: independent, parallel
            let results: Vec<Result<(f64, Option<(FoilPoint, Option<PointDetail>)>), FoilError>> = values
                .par_iter()
                .map(|&d| {
                    foiley::set_quiet(true);
                    prog.check()?;
                    let f = foil.with_flap(st.flap_hinge.0, st.flap_hinge.1, d)?;
                    let mut m = Marcher::new(&f, st, st.re, true)?;
                    let mut r = None;
                    // march alpha from 0 to the fixed alpha for robustness
                    let n = (st.alpha_fixed.abs() / 2.0).ceil() as usize;
                    for k in 1..=n {
                        let _ = m.solve(st.alpha_fixed * k as f64 / n as f64 - st.alpha_fixed / n as f64, prog);
                    }
                    if let Some(mut p) = m.solve(st.alpha_fixed, prog) {
                        p.alpha = st.alpha_fixed;
                        let det = if st.keep_details { Some(extract_detail(&m.xf, p)) } else { None };
                        r = Some((p, det));
                    }
                    tick(prog);
                    Ok((d, r))
                })
                .collect();
            for r in results {
                let (d, res) = r?;
                match res {
                    Some((mut p, det)) => {
                        prog.count_ok();
                        p.flap = d;
                        if let Some(mut dd) = det {
                            dd.point = p;
                            details.push(dd);
                        }
                        polar.insert(p);
                    }
                    None => {
                        prog.count_failed();
                        polar.failed.push(d);
                    }
                }
            }
        }
        _ => {
            let lalfa = sweep.is_alpha();
            if matches!(sweep, Sweep::Re { .. } | Sweep::Flap { .. }) {
                return Err(FoilError::Geometry("this sweep kind needs a type 4 (Re) or type 6 (flap) polar".into()));
            }
            let chunks = march_order(&values, st.from_zero);
            let threads = rayon::current_num_threads();
            let pieces = if st.parallel && total >= 16 { (threads / 2).clamp(1, 8) } else { 1 };
            let jobs = if pieces > 1 { split_chunks(chunks, pieces) } else { chunks.into_iter().map(|c| (c, true)).collect() };
            let results: Vec<Result<Vec<(f64, Option<(FoilPoint, Option<PointDetail>)>)>, FoilError>> = jobs
                .par_iter()
                .map(|(chunk, _anchored)| {
                    foiley::set_quiet(true);
                    let mut m = Marcher::new(foil, st, st.re, lalfa)?;
                    let mut out = Vec::with_capacity(chunk.len());
                    for &v in chunk {
                        prog.check()?;
                        let r = m.solve(v, prog).map(|p| {
                            let d = if st.keep_details { Some(extract_detail(&m.xf, p)) } else { None };
                            (p, d)
                        });
                        match &r {
                            Some((p, _)) => {
                                prog.count_ok();
                                prog.debug(format!("{} {}: CL={:.4} CD={:.5} Cm={:.4}", foil.name, m.fmt(v), p.cl, p.cd, p.cm));
                            }
                            None => {
                                prog.count_failed();
                                prog.warn(format!("{} {}: not converged after recovery attempts (skipped)", foil.name, m.fmt(v)));
                            }
                        }
                        out.push((v, r));
                        tick(prog);
                    }
                    Ok(out)
                })
                .collect();
            for r in results {
                for (v, res) in r? {
                    match res {
                        Some((p, d)) => {
                            if p.cd >= st.min_cd || !st.viscous {
                                if let Some(d) = d {
                                    details.push(d);
                                }
                                polar.insert(p);
                            }
                        }
                        None => {
                            if !polar.failed.iter().any(|f| (f - v).abs() < 1e-9) {
                                polar.failed.push(v);
                            }
                        }
                    }
                }
            }
            // a point solved by two chunks may have failed in one: drop it from `failed`
            polar.failed.retain(|f| !polar.points.iter().any(|p| if lalfa { (p.alpha - f).abs() < 1e-6 } else { (p.cl - f).abs() < 1e-6 }));
            polar.failed.sort_by(|a, b| a.total_cmp(b));
        }
    }
    details.sort_by(|a, b| (a.point.alpha + a.point.flap).total_cmp(&(b.point.alpha + b.point.flap)));
    prog.set(1.0);
    let nf = polar.failed.len();
    if nf > 0 {
        prog.warn(format!("{} / {}: {} of {} points did not converge", foil.name, polar.name, nf, total));
    }
    Ok((polar, details))
}

/// Solve a single operating point and return full distributions.
pub fn analyze_point(foil: &Foil, st: &AnalysisSettings, alpha_or_cl: f64, is_alpha: bool, prog: &Progress) -> Result<Option<PointDetail>, FoilError> {
    foil.validate()?;
    st.validate().map_err(FoilError::Geometry)?;
    foiley::set_quiet(true);
    let mut m = Marcher::new(foil, st, st.re, is_alpha)?;
    prog.stage(format!("{}: solving {}", foil.name, m.fmt(alpha_or_cl)));
    // warm up from zero for high angles
    if is_alpha && alpha_or_cl.abs() > 6.0 && st.viscous {
        let n = (alpha_or_cl.abs() / 3.0).floor() as usize;
        for k in 1..=n {
            prog.check()?;
            let _ = m.solve(alpha_or_cl * k as f64 / (n + 1) as f64, prog);
            prog.set(k as f64 / (n + 1) as f64);
        }
    }
    let r = m.solve(alpha_or_cl, prog).map(|p| extract_detail(&m.xf, p));
    prog.set(1.0);
    Ok(r)
}

/// Inviscid helper used in tests and quick previews.
pub fn inviscid_cl(foil: &Foil, alpha_deg: f64) -> Option<f64> {
    let st = AnalysisSettings { viscous: false, ..Default::default() };
    let mut xf = setup_engine(foil, &st, 0.0).ok()?;
    let r = catch_unwind(AssertUnwindSafe(|| solve_point(&mut xf, true, alpha_deg * DTOR, 1))).ok()?;
    Some(polar_point(&xf).cl).filter(|_| r.converged)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn naca0012_polar() {
        let f = Foil::naca("0012").unwrap();
        let st = AnalysisSettings { re: 1e6, ..Default::default() };
        let p = Progress::new();
        let (pol, _) = analyze(&f, &st, &Sweep::Alpha { from: -4.0, to: 10.0, step: 1.0 }, &p).unwrap();
        assert!(pol.points.len() >= 14, "{} points, failed {:?}", pol.points.len(), pol.failed);
        let p0 = pol.at_alpha(0.0).unwrap().0;
        assert!(p0.cl.abs() < 1e-3);
        assert!((p0.cd - 0.0054).abs() < 0.0008, "cd0 = {}", p0.cd);
        let p5 = pol.at_alpha(5.0).unwrap().0;
        assert!((p5.cl - 0.55).abs() < 0.05, "cl5 = {}", p5.cl);
        assert!((p.fraction() - 1.0).abs() < 1e-9);
        assert!((pol.cl_alpha().unwrap() - 0.11).abs() < 0.01);
    }
    #[test]
    fn hard_sweep_recovers_and_cancels() {
        let f = Foil::naca("4412").unwrap();
        let st = AnalysisSettings { re: 2e5, ..Default::default() };
        let p = Progress::new();
        let (pol, _) = analyze(&f, &st, &Sweep::Alpha { from: -6.0, to: 18.0, step: 0.5 }, &p).unwrap();
        assert!(pol.points.len() > 40, "{} ok, failed {:?}", pol.points.len(), pol.failed);
        let c = Progress::new();
        c.cancel();
        assert!(matches!(analyze(&f, &st, &Sweep::Alpha { from: 0.0, to: 5.0, step: 1.0 }, &c), Err(FoilError::Cancelled)));
    }
    #[test]
    fn cl_sweep_and_type2_and_type4() {
        let f = Foil::naca("2412").unwrap();
        let st = AnalysisSettings { re: 5e5, ..Default::default() };
        let (pol, _) = analyze(&f, &st, &Sweep::Cl { from: 0.2, to: 1.0, step: 0.2 }, &Progress::none()).unwrap();
        assert_eq!(pol.points.len(), 5, "{:?}", pol.failed);
        assert!((pol.points[2].cl - 0.6).abs() < 1e-3);
        let st2 = AnalysisSettings { polar_type: PolarType::FixedLift, re: 3e5, ..Default::default() };
        let (pol2, _) = analyze(&f, &st2, &Sweep::Alpha { from: 1.0, to: 6.0, step: 1.0 }, &Progress::none()).unwrap();
        // Re = Re*sqrt(CL) / sqrt(CL)
        let q = pol2.points.last().unwrap();
        assert!((q.re - 3e5 / q.cl.sqrt()).abs() / q.re < 0.01, "{q:?}");
        let st4 = AnalysisSettings { polar_type: PolarType::FixedAlpha, alpha_fixed: 3.0, ..Default::default() };
        let (pol4, _) = analyze(&f, &st4, &Sweep::Re { from: 2e5, to: 1e6, step: 2e5 }, &Progress::none()).unwrap();
        assert_eq!(pol4.points.len(), 5);
        assert!(pol4.points[0].cd > pol4.points[4].cd);
    }
    #[test]
    fn detail_and_inviscid() {
        let f = Foil::naca("0012").unwrap();
        let d = analyze_point(&f, &AnalysisSettings::default(), 4.0, true, &Progress::none()).unwrap().unwrap();
        assert!(d.top.x.len() > 20 && d.bot.x.len() > 20 && !d.wake.x.is_empty());
        assert!(d.point.xtr_top < d.point.xtr_bot);
        let cl = inviscid_cl(&f, 5.0).unwrap();
        assert!((cl - 0.6).abs() < 0.05, "{cl}");
    }
}
