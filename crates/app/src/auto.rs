//! Automatic 2D polar meshes for 3D viscous analyses, sized from what the
//! analysis will actually need (F021).
//!
//! * The required envelope comes from the operating states the definition
//!   expands to: the alpha / speed / sideslip sweeps, the T8 list, the T6
//!   ranges (speed, alpha, beta, mass, flaps, tilts), the T7 gains
//!   (trimmed), and the mass for the speed of T2/T3/T7.
//! * Re per foil is V c / nu over the strips that use it (min/max chord).
//! * The CL / alpha span comes from a cheap inviscid VLM2 pre-pass that
//!   reuses each configuration's unit solutions (RHS combinations only).
//! * Mach, Ncrit and trips come from the selected family
//!   (`PlanePolarDef::family`).
//! * Only the Re stations and flap angles that no existing family covers are
//!   generated.  Deflected geometries are generated as copies named with
//!   `deflected_name`; `Project::add_foil_polar` files their polars under
//!   the base foil with `flap_deg` set (`adopt`).

use foiley_core::hash::Hash32;
use foiley_core::Vec3;
use foiley_foil::batch::{re_list_geometric, BatchCase};
use foiley_foil::{AnalysisSettings, Foil, FoilPolar, Sweep};
use foiley_plane::analysis::{make_system, mesh_options, wind_dir, wind_normal, ControlTarget, Ctx, PolarKind, G};
use foiley_plane::mesh::build;
use foiley_plane::{FamilyId, FoilDb, Method, Plane, PlanePolarDef};
use std::collections::HashMap;

/// Suffix marking a deflected copy of a foil in an automatic batch.
pub const FLAP_TAG: &str = " @flap";

/// Name of the deflected copy of `base` (flap `deg` degrees, + down).
pub fn deflected_name(base: &str, deg: f64) -> String {
    format!("{base}{FLAP_TAG}{deg:+.2}")
}

/// (base foil, flap angle) of a deflected copy's name.
pub fn parse_deflected(name: &str) -> Option<(&str, f64)> {
    let i = name.rfind(FLAP_TAG)?;
    let d: f64 = name[i + FLAP_TAG.len()..].parse().ok()?;
    Some((&name[..i], d))
}

/// File a polar of a deflected copy under its base foil: `foil` becomes
/// the base, `flap_deg` the angle, `foil_hash` the base geometry's hash,
/// and the name gets the angle (so it never merges with an undeflected
/// polar of the same settings).  Other polars are returned unchanged.
pub fn adopt(mut p: FoilPolar, base_hash: impl Fn(&str) -> Option<Hash32>) -> FoilPolar {
    if let Some((base, d)) = parse_deflected(&p.foil).map(|(b, d)| (b.to_string(), d)) {
        if let Some(h) = base_hash(&base) {
            p.foil = base;
            p.flap_deg = d;
            p.foil_hash = Some(h);
            p.name = format!("{} flap{d:+}", p.name);
        }
    }
    p
}

/// What one foil's polars must cover.
#[derive(Clone, Debug, PartialEq)]
pub struct FoilCoverage {
    pub foil: String,
    /// Reynolds range of the strips using the foil (with a 10 % margin)
    pub re: (f64, f64),
    /// nonzero flap angles of strips using the foil, sorted
    pub flaps: Vec<f64>,
}

/// The envelope a 3D analysis needs from the 2D polars.
#[derive(Clone, Debug, PartialEq)]
pub struct Coverage {
    pub foils: Vec<FoilCoverage>,
    /// effective strip angle of attack range, degrees
    pub alpha: (f64, f64),
    /// strip CL range
    pub cl: (f64, f64),
    /// operating states evaluated
    pub states: usize,
}

/// One operating state of the analysis.
#[derive(Clone, Debug)]
struct State {
    alpha: f64,
    beta: f64,
    /// None: speed from the lift = weight balance
    v: Option<f64>,
    mass: f64,
    /// T7: find the trimmed alpha (Cm = 0) first
    trim: bool,
    flaps: Vec<((usize, usize), f64)>,
    tilts: Vec<(usize, f64)>,
}

fn states(def: &PlanePolarDef, mass: f64) -> Vec<State> {
    let st = |alpha: f64, beta: f64, v: Option<f64>| State { alpha, beta, v, mass, trim: false, flaps: vec![], tilts: vec![] };
    let range = def.range.values();
    match def.kind {
        PolarKind::FixedSpeed => range.iter().map(|&a| st(a, def.beta, Some(def.velocity))).collect(),
        PolarKind::FixedLift | PolarKind::Glide => range.iter().map(|&a| st(a, def.beta, None)).collect(),
        PolarKind::FixedAlpha => range.iter().map(|&v| st(def.alpha, def.beta, Some(v))).collect(),
        PolarKind::Sideslip => range.iter().map(|&b| st(def.alpha, b, Some(def.velocity))).collect(),
        PolarKind::List => def.list.iter().map(|&(a, b, v)| st(a, b, Some(v))).collect(),
        PolarKind::Control => range
            .iter()
            .map(|&t| {
                let val = |target: ControlTarget| def.ranges.iter().find(|r| r.target == target).map(|r| r.min + t * (r.max - r.min));
                let mut s = st(val(ControlTarget::Alpha).unwrap_or(def.alpha), val(ControlTarget::Beta).unwrap_or(def.beta), Some(val(ControlTarget::Velocity).unwrap_or(def.velocity)));
                s.mass = val(ControlTarget::Mass).unwrap_or(mass);
                for r in &def.ranges {
                    let x = r.min + t * (r.max - r.min);
                    match r.target {
                        ControlTarget::Flap { wing, index } => s.flaps.push(((wing, index), x)),
                        ControlTarget::Tilt { wing } => s.tilts.push((wing, x)),
                        _ => {}
                    }
                }
                s
            })
            .collect(),
        PolarKind::Stability => range
            .iter()
            .map(|&t| {
                let mut s = st(0.0, def.beta, None);
                s.trim = true;
                for (target, gain) in &def.gains {
                    match *target {
                        ControlTarget::Flap { wing, index } => s.flaps.push(((wing, index), gain * t)),
                        ControlTarget::Tilt { wing } => s.tilts.push((wing, gain * t)),
                        _ => {}
                    }
                }
                s
            })
            .collect(),
    }
}

type Config = (Vec<((usize, usize), f64)>, Vec<(usize, f64)>);

/// The 2D envelope `def` needs on `plane`, from an inviscid VLM2 pre-pass
/// over its operating states (one factorization per flap / tilt setting;
/// every state of a setting is a combination of its unit solutions).
pub fn required(plane: &Plane, db: &FoilDb, def: &PlanePolarDef) -> Result<Coverage, String> {
    let pre = PlanePolarDef { method: Method::Vlm2, kind: PolarKind::FixedSpeed, viscous: false, free_wake: 0, derivatives: false, keep_cp: false, ..def.clone() };
    let ctx = Ctx::new(plane, db, &pre).map_err(|e| e.to_string())?;
    let mass = ctx.inertia.mass;
    let cog = if mass > 0.0 { ctx.inertia.cog } else { Vec3::ZERO };
    let nu = def.fluid.nu();
    let rho = def.fluid.rho;
    // group the states by configuration
    let mut groups: Vec<(Config, Vec<State>)> = Vec::new();
    for s in states(def, mass) {
        let cfg = (s.flaps.clone(), s.tilts.clone());
        match groups.iter_mut().find(|g| g.0 == cfg) {
            Some(g) => g.1.push(s),
            None => groups.push((cfg, vec![s])),
        }
    }
    let mut per_foil: Vec<FoilCoverage> = Vec::new();
    let (mut amin, mut amax, mut cmin, mut cmax) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    let mut n = 0;
    for ((flaps, tilts), sts) in &groups {
        let flaps: HashMap<(usize, usize), f64> = flaps.iter().cloned().collect();
        let tilts: HashMap<usize, f64> = tilts.iter().cloned().collect();
        let mesh = build(plane, db, &mesh_options(&pre, &flaps, &tilts)).map_err(|e| e.to_string())?;
        let sys = make_system(&pre, &mesh, cog, ctx.c_ref, &foiley_core::Progress::none()).map_err(|e| e.to_string())?;
        for s in sts {
            let mut alpha = s.alpha;
            if s.trim {
                // first stable Cm = 0 crossing on a 1 deg scan
                let cm = |a: f64| ctx.reduce(&mesh, &sys.loads(wind_dir(a, s.beta), Vec3::ZERO), cog, a, s.beta, 1.0, s.mass).cm;
                let mut prev = (-10.0, cm(-10.0));
                let mut found = None;
                let mut a = -9.0;
                while a <= 20.0 {
                    let c = cm(a);
                    if prev.1 > 0.0 && c <= 0.0 {
                        found = Some(prev.0 + (a - prev.0) * prev.1 / (prev.1 - c));
                        break;
                    }
                    prev = (a, c);
                    a += 1.0;
                }
                match found {
                    Some(a) => alpha = a,
                    None => continue,
                }
            }
            let v = match s.v {
                Some(v) => v,
                None => {
                    let r = sys.loads(wind_dir(alpha, s.beta), Vec3::ZERO);
                    let l1 = (r.force + r.induced_force).dot(wind_normal(alpha)) * rho;
                    if l1 <= 1e-12 || s.mass <= 0.0 {
                        continue;
                    }
                    (s.mass * G / l1).sqrt()
                }
            };
            if !(v > 0.0 && v.is_finite()) {
                continue;
            }
            n += 1;
            let raw = sys.loads(wind_dir(alpha, s.beta) * v, Vec3::ZERO);
            let opp = ctx.reduce(&mesh, &raw, cog, alpha, s.beta, v, s.mass);
            for (st, r) in mesh.strips.iter().zip(&opp.strips) {
                amin = amin.min(r.alpha_eff);
                amax = amax.max(r.alpha_eff);
                cmin = cmin.min(r.cl);
                cmax = cmax.max(r.cl);
                let re = v * st.chord / nu;
                for f in [&st.foil_a, &st.foil_b] {
                    let i = match per_foil.iter().position(|c| &c.foil == f) {
                        Some(i) => i,
                        None => {
                            per_foil.push(FoilCoverage { foil: f.clone(), re: (re, re), flaps: vec![] });
                            per_foil.len() - 1
                        }
                    };
                    let c = &mut per_foil[i];
                    c.re = (c.re.0.min(re), c.re.1.max(re));
                    if st.flap_deg.abs() > 1e-9 && !c.flaps.iter().any(|d| (d - st.flap_deg).abs() < 1e-9) {
                        c.flaps.push(st.flap_deg);
                    }
                }
            }
        }
    }
    if n == 0 {
        return Err("automatic 2D polars: none of the analysis' operating states could be evaluated (no lift for the weight balance, or no trim)".into());
    }
    for c in &mut per_foil {
        c.re = (0.9 * c.re.0, 1.1 * c.re.1);
        c.flaps.sort_by(|a, b| a.total_cmp(b));
    }
    Ok(Coverage { foils: per_foil, alpha: (amin, amax), cl: (cmin, cmax), states: n })
}

/// One batch of the plan: foils (deflected copies included) x cases.
#[derive(Clone, Debug)]
pub struct AutoBatch {
    pub foils: Vec<Foil>,
    pub cases: Vec<BatchCase>,
    pub settings: AnalysisSettings,
    pub sweep: Sweep,
}

/// The 2D polars to generate before a 3D analysis.
#[derive(Clone, Debug, Default)]
pub struct AutoPlan {
    pub coverage: Option<Coverage>,
    pub batches: Vec<AutoBatch>,
}

impl AutoPlan {
    pub fn is_empty(&self) -> bool {
        self.polar_count() == 0
    }
    pub fn polar_count(&self) -> usize {
        self.batches.iter().map(|b| b.foils.len() * b.cases.len()).sum()
    }
    pub fn describe(&self) -> String {
        let mut names: Vec<&str> = self.batches.iter().flat_map(|b| b.foils.iter().map(|f| f.name.as_str())).collect();
        names.sort();
        names.dedup();
        let (lo, hi) = self.batches.iter().flat_map(|b| b.cases.iter().map(|c| c.re)).fold((f64::MAX, 0.0f64), |(a, b), r| (a.min(r), b.max(r)));
        format!("{} 2D polars for {} (Re {lo:.0}..{hi:.0})", self.polar_count(), names.join(", "))
    }
}

/// (mach, ncrit, xtr_top, xtr_bot) of the family the plan fills for a
/// foil, and that family if it exists already.  The selection's values
/// (defaults where it leaves them open) when it matches nothing yet.
fn family_params(db: &FoilDb, foil: &str, def: &PlanePolarDef) -> Result<((f64, f64, f64, f64), Option<FamilyId>), String> {
    let sel = &def.family;
    let d = AnalysisSettings::default();
    let fresh = (sel.mach.unwrap_or(d.mach), sel.ncrit.unwrap_or(d.ncrit), sel.xtr_top.unwrap_or(d.xtr_top), sel.xtr_bot.unwrap_or(d.xtr_bot));
    let hits = db.families_of(foil).into_iter().filter(|&i| db.family(i).key.flap_deg == 0.0 && sel.matches(&db.family(i).key)).count();
    if hits == 0 {
        return Ok((fresh, None));
    }
    let id = db.resolve(foil, sel)?.expect("a matching family exists");
    let k = db.family(id).key;
    Ok(((k.mach, k.ncrit, k.xtr_top, k.xtr_bot), Some(id)))
}

/// Stations of `want` that no polar of `have` (Re values) covers: a
/// station is covered by a polar within half a station spacing of it.
fn missing(want: &[f64], have: &[f64]) -> Vec<f64> {
    let r = if want.len() >= 2 { (want[want.len() - 1] / want[0]).powf(1.0 / (want.len() - 1) as f64) } else { 1.2 };
    let tol = r.sqrt().max(1.05);
    want.iter().copied().filter(|&s| !have.iter().any(|&h| h >= s / tol && h <= s * tol)).collect()
}

/// The 2D polars `def` on `plane` still needs, as batches: `count` Re
/// stations per foil over its required range, at the selected family's
/// Mach / Ncrit / trips, and the same for every flap angle its deflected
/// strips use; stations and angles already covered are left out.
pub fn plan(plane: &Plane, foils: &[Foil], db: &FoilDb, def: &PlanePolarDef, count: usize) -> Result<AutoPlan, String> {
    if !(def.viscous || def.method == Method::Llt) || count == 0 {
        return Ok(AutoPlan::default());
    }
    let cov = required(plane, db, def)?;
    let (a0, a1) = cov.alpha;
    let sweep = Sweep::Alpha { from: (a0 - 3.0).floor().clamp(-15.0, -2.0), to: (a1 + 3.0).ceil().clamp(6.0, 22.0), step: 0.5 };
    // (family params, stations) -> foils
    let mut groups: Vec<((f64, f64, f64, f64), Vec<f64>, Vec<Foil>)> = Vec::new();
    let mut add = |params: (f64, f64, f64, f64), re: Vec<f64>, f: Foil| {
        if re.is_empty() {
            return;
        }
        match groups.iter_mut().find(|g| g.0 == params && g.1 == re) {
            Some(g) => g.2.push(f),
            None => groups.push((params, re, vec![f])),
        }
    };
    for c in &cov.foils {
        let Some(foil) = foils.iter().find(|f| f.name == c.foil) else { continue };
        let (params, base) = family_params(db, &c.foil, def)?;
        let lo = c.re.0.max(1e4);
        let want = re_list_geometric(lo, c.re.1.max(lo * 1.2), count.max(2));
        let have = |id: Option<FamilyId>| id.map(|i| db.family(i).polars.iter().map(|p| p.re).collect::<Vec<f64>>()).unwrap_or_default();
        add(params, missing(&want, &have(base)), foil.clone());
        let set = base.map(|b| db.flap_families(b)).unwrap_or_default();
        for &d in &c.flaps {
            let existing = set.iter().find(|(f, _)| (f - d).abs() < 0.25).map(|x| x.1);
            let Some((xh, yh)) = foil.hinge else { continue };
            let Ok(mut g) = foil.with_flap(xh, yh, d) else { continue };
            g.name = deflected_name(&foil.name, d);
            add(params, missing(&want, &have(existing)), g);
        }
    }
    let batches = groups
        .into_iter()
        .map(|((mach, ncrit, xt, xb), re, foils)| AutoBatch {
            foils,
            cases: re.into_iter().map(|re| BatchCase { re, mach, ncrit }).collect(),
            settings: AnalysisSettings { mach, ncrit, xtr_top: xt, xtr_bot: xb, ..Default::default() },
            sweep: sweep.clone(),
        })
        .collect();
    Ok(AutoPlan { coverage: Some(cov), batches })
}
