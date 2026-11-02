//! 3D polars: definitions (types 1-8), the drivers, and the reduction of
//! raw loads to coefficients and span distributions.

use crate::foildb::FoilDb;
use crate::geom::{Inertia, Plane};
use crate::mesh::{build, Mesh, MeshOptions};
use crate::panel::PanelSystem;
use crate::stability::{self, StabilityResult};
use crate::vlm::VlmSystem;
use crate::PlaneError;
use foiley_core::atmosphere::Fluid;
use foiley_core::{Progress, Vec3, DEG};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

pub const G: f64 = 9.80665;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum Method {
    /// Nonlinear lifting line (main wing only, uses the 2D polars for lift)
    Llt,
    /// Horseshoe vortex lattice
    Vlm1,
    /// Ring vortex lattice (recommended default)
    #[default]
    Vlm2,
    /// 3D source/doublet panels on thick surfaces, fuselage included
    Panel,
    /// 3D source/doublet panels on triangles (flow5 default; also STL/STEP bodies)
    TriPanel,
}

impl Method {
    pub const ALL: [Method; 5] = [Method::Llt, Method::Vlm1, Method::Vlm2, Method::Panel, Method::TriPanel];
    /// Thick-surface methods (panels on both sides of the wings, fuselage included).
    pub fn is_panel(self) -> bool {
        matches!(self, Method::Panel | Method::TriPanel)
    }
    pub fn label(self) -> &'static str {
        match self {
            Method::Llt => "LLT - nonlinear lifting line (main wing only)",
            Method::Vlm1 => "VLM1 - horseshoe vortices",
            Method::Vlm2 => "VLM2 - ring vortices (recommended)",
            Method::Panel => "3D panels - thick surfaces + fuselage",
            Method::TriPanel => "3D panels on triangles (uniform doublets)",
        }
    }
    pub fn short(self) -> &'static str {
        match self {
            Method::Llt => "LLT",
            Method::Vlm1 => "VLM1",
            Method::Vlm2 => "VLM2",
            Method::Panel => "Panel",
            Method::TriPanel => "TriPanel",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum PolarKind {
    /// T1: fixed speed, sweep alpha
    #[default]
    FixedSpeed,
    /// T2: fixed lift (speed from weight), sweep alpha
    FixedLift,
    /// T3: glide polar (speed from weight and total aerodynamic force)
    Glide,
    /// T4: fixed alpha, sweep speed
    FixedAlpha,
    /// T5: sideslip sweep at fixed alpha and speed
    Sideslip,
    /// T6: control polar, ranges on speed, alpha, beta, mass, CoG, flaps, tilts
    Control,
    /// T7: stability polar: trimmed points vs control parameter, with modes
    Stability,
    /// T8: explicit list of (alpha, beta, speed)
    List,
}

impl PolarKind {
    pub const ALL: [PolarKind; 8] = [PolarKind::FixedSpeed, PolarKind::FixedLift, PolarKind::Glide, PolarKind::FixedAlpha, PolarKind::Sideslip, PolarKind::Control, PolarKind::Stability, PolarKind::List];
    pub fn number(self) -> u8 {
        match self {
            PolarKind::FixedSpeed => 1,
            PolarKind::FixedLift => 2,
            PolarKind::Glide => 3,
            PolarKind::FixedAlpha => 4,
            PolarKind::Sideslip => 5,
            PolarKind::Control => 6,
            PolarKind::Stability => 7,
            PolarKind::List => 8,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            PolarKind::FixedSpeed => "T1 - fixed speed, sweep angle of attack",
            PolarKind::FixedLift => "T2 - fixed lift = weight, sweep angle of attack",
            PolarKind::Glide => "T3 - glide (speed polar), sweep angle of attack",
            PolarKind::FixedAlpha => "T4 - fixed angle of attack, sweep speed",
            PolarKind::Sideslip => "T5 - sideslip sweep",
            PolarKind::Control => "T6 - control polar (ranges vs control parameter)",
            PolarKind::Stability => "T7 - stability (trimmed, modes)",
            PolarKind::List => "T8 - list of operating points",
        }
    }
    /// Name of the swept variable.
    pub fn swept(self) -> &'static str {
        match self {
            PolarKind::FixedSpeed | PolarKind::FixedLift | PolarKind::Glide => "alpha (deg)",
            PolarKind::FixedAlpha => "speed (m/s)",
            PolarKind::Sideslip => "beta (deg)",
            PolarKind::Control | PolarKind::Stability => "control parameter",
            PolarKind::List => "list index",
        }
    }
    pub fn needs_mass(self) -> bool {
        matches!(self, PolarKind::FixedLift | PolarKind::Glide | PolarKind::Stability)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Range {
    pub from: f64,
    pub to: f64,
    pub step: f64,
}

impl Range {
    pub fn values(&self) -> Vec<f64> {
        let s = if self.step == 0.0 { 1.0 } else { self.step.abs() };
        let n = ((self.to - self.from).abs() / s + 0.5).floor() as usize;
        let d = if self.to >= self.from { 1.0 } else { -1.0 };
        (0..=n.min(10_000)).map(|i| self.from + d * s * i as f64).collect()
    }
}

/// Reference dimensions for coefficients.
#[derive(Clone, Copy, Debug, PartialEq, Default, Serialize, Deserialize)]
pub enum RefDims {
    /// main wing planform projected on the xy plane (flow5 default)
    #[default]
    Projected,
    /// main wing developed planform
    Planform,
    Custom { area: f64, span: f64, chord: f64 },
}

/// What a control-polar range acts on.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum ControlTarget {
    Velocity,
    Alpha,
    Beta,
    Mass,
    CogX,
    CogZ,
    /// trailing-edge flap `index` of wing `wing` (degrees, + down)
    Flap { wing: usize, index: usize },
    /// extra tilt of wing `wing` (degrees, + nose up)
    Tilt { wing: usize },
}

impl ControlTarget {
    pub fn label(&self, plane: &Plane) -> String {
        let wn = |w: usize| plane.wings.get(w).map(|w| w.name.clone()).unwrap_or_else(|| format!("wing {w}"));
        match self {
            ControlTarget::Velocity => "speed (m/s)".into(),
            ControlTarget::Alpha => "alpha (deg)".into(),
            ControlTarget::Beta => "beta (deg)".into(),
            ControlTarget::Mass => "mass (kg)".into(),
            ControlTarget::CogX => "CoG x (m)".into(),
            ControlTarget::CogZ => "CoG z (m)".into(),
            ControlTarget::Flap { wing, index } => format!("{} flap {} (deg)", wn(*wing), index + 1),
            ControlTarget::Tilt { wing } => format!("{} tilt (deg)", wn(*wing)),
        }
    }
}

/// A range driven by the control parameter t: value = min + t (max - min).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ControlRange {
    pub target: ControlTarget,
    pub min: f64,
    pub max: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ExtraDrag {
    pub name: String,
    pub area: f64,
    pub cd: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct LltSettings {
    pub stations: usize,
    pub max_iter: usize,
    /// extra damping of the induced-angle iteration (1 = automatic per-station
    /// relaxation; increase if a case oscillates)
    pub relax: f64,
    /// convergence on the induced angle, degrees
    pub tolerance: f64,
}

impl Default for LltSettings {
    fn default() -> Self {
        LltSettings { stations: 40, max_iter: 400, relax: 1.0, tolerance: 0.005 }
    }
}

/// Full definition of a 3D analysis.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct PlanePolarDef {
    pub name: String,
    pub kind: PolarKind,
    pub method: Method,
    /// speed for T1/T5/T7 start, m/s
    pub velocity: f64,
    /// alpha for T4/T5, degrees
    pub alpha: f64,
    /// sideslip for T1-T4, degrees
    pub beta: f64,
    /// main sweep range (alpha, speed, beta or control parameter)
    pub range: Range,
    /// T8 list of (alpha, beta, speed)
    pub list: Vec<(f64, f64, f64)>,
    /// T6 ranges
    pub ranges: Vec<ControlRange>,
    /// T7 control gains: (target, degrees per unit control parameter)
    pub gains: Vec<(ControlTarget, f64)>,
    /// preset flap deflections for all polar types, degrees
    pub flaps: Vec<((usize, usize), f64)>,
    /// add viscous drag interpolated from the 2D polars
    pub viscous: bool,
    /// interpolate viscous drag from the strip CL (xflr5 default) rather than its angle
    pub viscous_from_cl: bool,
    pub fluid: Fluid,
    pub reference: RefDims,
    /// override plane mass properties
    pub inertia: Option<Inertia>,
    /// ground effect: height of the ground below the origin, m
    pub ground_height: Option<f64>,
    pub extra_drag: Vec<ExtraDrag>,
    /// fuselage skin-friction drag (turbulent flat plate x form factor)
    pub fuse_drag: bool,
    /// include fuselage panels (panel method)
    pub include_fuse: bool,
    /// compute stability derivatives and modes for every point (needs mass)
    pub derivatives: bool,
    /// wake length for the panel method, in reference chords
    pub wake_length: f64,
    pub wake_panels: usize,
    pub wake_ratio: f64,
    pub llt: LltSettings,
    /// keep panel Cp in the results (memory heavy for big meshes)
    pub keep_cp: bool,
    /// VLM free-wake relaxation iterations (0 = rigid wake along x)
    pub free_wake: usize,
    /// treat foils without loaded geometry as flat plates
    pub allow_missing_foils: bool,
}

impl Default for PlanePolarDef {
    fn default() -> Self {
        PlanePolarDef {
            name: String::new(),
            kind: PolarKind::FixedSpeed,
            method: Method::Vlm2,
            velocity: 10.0,
            alpha: 0.0,
            beta: 0.0,
            range: Range { from: -4.0, to: 10.0, step: 1.0 },
            list: vec![],
            ranges: vec![],
            gains: vec![],
            flaps: vec![],
            viscous: true,
            viscous_from_cl: true,
            fluid: Fluid::default(),
            reference: RefDims::Projected,
            inertia: None,
            ground_height: None,
            extra_drag: vec![],
            fuse_drag: false,
            include_fuse: true,
            derivatives: false,
            wake_length: 30.0,
            wake_panels: 5,
            wake_ratio: 1.1,
            llt: LltSettings::default(),
            keep_cp: true,
            allow_missing_foils: false,
            free_wake: 0,
        }
    }
}

impl PlanePolarDef {
    pub fn auto_name(&self) -> String {
        let mut s = format!("T{}-{}", self.kind.number(), self.method.short());
        match self.kind {
            PolarKind::FixedSpeed | PolarKind::Sideslip => s.push_str(&format!("-{:.1}m/s", self.velocity)),
            PolarKind::FixedAlpha => s.push_str(&format!("-a{:.1}", self.alpha)),
            _ => {}
        }
        if self.kind == PolarKind::Sideslip {
            s.push_str(&format!("-a{:.1}", self.alpha));
        }
        if !self.viscous {
            s.push_str("-inviscid");
        }
        if self.ground_height.is_some() {
            s.push_str("-ground");
        }
        s
    }

    pub fn validate(&self) -> Result<(), String> {
        if !(self.fluid.rho > 0.0 && self.fluid.mu > 0.0) {
            return Err("fluid density and viscosity must be positive".into());
        }
        if matches!(self.kind, PolarKind::FixedSpeed | PolarKind::Sideslip | PolarKind::Control) && !(self.velocity > 0.0) && !self.ranges.iter().any(|r| r.target == ControlTarget::Velocity) {
            return Err("speed must be positive".into());
        }
        if self.kind == PolarKind::FixedAlpha && self.range.values().iter().any(|v| *v <= 0.0) {
            return Err("speed range of a T4 polar must be positive".into());
        }
        if self.kind == PolarKind::List && self.list.is_empty() {
            return Err("T8 polar needs a list of (alpha, beta, speed)".into());
        }
        if self.method == Method::Llt && !matches!(self.kind, PolarKind::FixedSpeed | PolarKind::FixedLift | PolarKind::Glide | PolarKind::FixedAlpha) {
            return Err("LLT supports polar types 1-4 only; use VLM or panels".into());
        }
        if self.wake_length <= 0.0 || self.wake_panels == 0 {
            return Err("wake length and panel count must be positive".into());
        }
        Ok(())
    }
}

/// Raw loads per unit density from a potential solver.
#[derive(Clone, Debug, Default)]
pub struct RawLoads {
    pub panel_force: Vec<Vec3>,
    /// inviscid force without induced drag (near field), per unit density
    pub force: Vec3,
    pub moment: Vec3,
    pub strip_gamma: Vec<f64>,
    /// induced velocity at each strip (half the Trefftz value)
    pub strip_wi: Vec<Vec3>,
    pub strip_force: Vec<Vec3>,
    /// induced-drag force from the Trefftz plane, per unit density
    pub induced_force: Vec3,
    pub cp: Vec<f64>,
    pub mu: Vec<f64>,
}

impl RawLoads {
    pub fn from_parts(mesh: &Mesh, panel_force: Vec<Vec3>, force: Vec3, moment: Vec3, sg: Vec<f64>, tre: Vec<Vec3>, vinf: Vec3) -> RawLoads {
        let mut strip_force = vec![Vec3::ZERO; mesh.strips.len()];
        for (p, f) in mesh.panels.iter().zip(&panel_force) {
            if let Some(s) = p.strip {
                strip_force[s] += *f;
            }
        }
        let wdir = vinf.normalized();
        let mut ind = Vec3::ZERO;
        for (s, (g, w)) in mesh.strips.iter().zip(sg.iter().zip(&tre)) {
            let f = w.cross(s.te_b - s.te_a) * *g;
            ind += wdir * f.dot(wdir);
        }
        // circulation carried across fuselages
        for (a, b, g) in crate::vlm::carry_segments(mesh, &sg) {
            let w = crate::vlm::trefftz_velocity(mesh, &sg, (a + b) * 0.5, wdir, None);
            let f = w.cross(b - a) * g;
            ind += wdir * f.dot(wdir);
        }
        let v2 = vinf.norm2().max(1e-30);
        let cp = mesh.panels.iter().zip(&panel_force).map(|(p, f)| -f.dot(p.normal) / (0.5 * v2 * p.area.max(1e-30))).collect();
        RawLoads { panel_force, force, moment, strip_gamma: sg, strip_wi: tre, strip_force, induced_force: ind, cp, mu: vec![] }
    }
}

/// Results for one spanwise strip.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct StripResult {
    pub wing: usize,
    pub span_pos: f64,
    pub chord: f64,
    pub cl: f64,
    pub cdi: f64,
    pub cdv: f64,
    /// induced angle, degrees (negative = downwash)
    pub ai: f64,
    pub re: f64,
    pub xtr_top: f64,
    pub xtr_bot: f64,
    /// circulation, m^2/s
    pub gamma: f64,
    /// lift per unit span, N/m
    pub lift_per_span: f64,
    pub twist: f64,
    /// effective angle of attack seen by the section, degrees
    pub alpha_eff: f64,
    pub cm_c4: f64,
    /// bending moment at this station from the outboard loads, N.m
    pub bending: f64,
    /// viscous interpolation outside the 2D polar envelope
    pub out_of_polar: bool,
}

/// One 3D operating point.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct PlaneOpp {
    pub alpha: f64,
    pub beta: f64,
    pub velocity: f64,
    pub ctrl: f64,
    pub cl: f64,
    pub cd: f64,
    pub cdi: f64,
    pub cdv: f64,
    pub cy: f64,
    pub cm: f64,
    /// rolling moment coefficient
    pub cr: f64,
    /// yawing moment coefficient
    pub cn: f64,
    pub lift: f64,
    pub drag: f64,
    pub x_cp: f64,
    /// span efficiency (Oswald) CL^2 / (pi AR CDi)
    pub efficiency: f64,
    pub mass: f64,
    pub cog: Vec3,
    pub strips: Vec<StripResult>,
    /// panel pressure coefficients (or loading for VLM), may be empty
    pub cp: Vec<f32>,
    pub warnings: Vec<String>,
    pub stability: Option<StabilityResult>,
    /// root bending moment of the main wing, N.m
    pub root_bending: f64,
    /// relaxed wake lines (free-wake VLM), for display
    #[serde(default)]
    pub wake: Vec<Vec<[f32; 3]>>,
}

impl PlaneOpp {
    pub fn ld(&self) -> f64 {
        if self.cd.abs() > 1e-12 {
            self.cl / self.cd
        } else {
            0.0
        }
    }
    /// sink rate for a gliding point, m/s
    pub fn vz(&self) -> f64 {
        let g = (self.cd / self.cl.max(1e-12)).atan();
        self.velocity * g.sin()
    }
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct PlanePolar {
    pub name: String,
    pub plane: String,
    pub def: Option<PlanePolarDef>,
    pub points: Vec<PlaneOpp>,
    pub failed: Vec<String>,
    /// reference dimensions used (area, span, chord)
    pub ref_dims: (f64, f64, f64),
    /// number of panels / unknowns
    pub unknowns: usize,
}

impl PlanePolar {
    /// Neutral point estimate x_NP = dX_cp*CL / dCL (first to last point).
    pub fn neutral_point(&self) -> Option<f64> {
        let a = self.points.first()?;
        let b = self.points.last()?;
        let d = b.cl - a.cl;
        if d.abs() < 1e-6 {
            return None;
        }
        Some((b.x_cp * b.cl - a.x_cp * a.cl) / d)
    }
    pub fn to_csv(&self) -> String {
        let mut s = String::from("alpha,beta,V,ctrl,CL,CD,CDi,CDv,CY,Cm,Cl_roll,Cn,L/D,Lift_N,Drag_N,Xcp,efficiency,Vz,root_bending\n");
        for p in &self.points {
            s.push_str(&format!(
                "{:.3},{:.3},{:.4},{:.4},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.4},{:.4},{:.5},{:.5},{:.4},{:.5},{:.4}\n",
                p.alpha, p.beta, p.velocity, p.ctrl, p.cl, p.cd, p.cdi, p.cdv, p.cy, p.cm, p.cr, p.cn, p.ld(), p.lift, p.drag, p.x_cp, p.efficiency, p.vz(), p.root_bending
            ));
        }
        s
    }
}

/// Air-relative velocity direction for alpha, beta (degrees);
/// beta > 0 = wind from the right (AVL / flow5 convention).
pub fn wind_dir(alpha: f64, beta: f64) -> Vec3 {
    let (a, b) = (alpha * DEG, -beta * DEG);
    Vec3::new(a.cos() * b.cos(), b.sin(), a.sin() * b.cos())
}
pub fn wind_normal(alpha: f64) -> Vec3 {
    let a = alpha * DEG;
    Vec3::new(-a.sin(), 0.0, a.cos())
}
pub fn wind_side(alpha: f64, beta: f64) -> Vec3 {
    let (a, b) = (alpha * DEG, -beta * DEG);
    Vec3::new(-a.cos() * b.sin(), b.cos(), -a.sin() * b.sin())
}

/// A factored potential-flow system.
pub enum System<'a> {
    Vlm(VlmSystem<'a>),
    Panel(PanelSystem<'a>),
}

impl<'a> System<'a> {
    pub fn loads(&self, vinf: Vec3, omega: Vec3) -> RawLoads {
        match self {
            System::Vlm(s) => {
                let g = s.solve(vinf, omega);
                s.loads(&g, vinf, omega)
            }
            System::Panel(s) => s.loads(vinf, omega),
        }
    }
    pub fn mesh(&self) -> &Mesh {
        match self {
            System::Vlm(s) => s.lat.mesh,
            System::Panel(s) => s.mesh,
        }
    }
    pub fn cog(&self) -> Vec3 {
        match self {
            System::Vlm(s) => s.cog,
            System::Panel(s) => s.cog,
        }
    }
}

/// Everything needed to reduce loads to an operating point.
pub struct Ctx<'a> {
    pub plane: &'a Plane,
    pub db: &'a FoilDb,
    pub def: &'a PlanePolarDef,
    pub s_ref: f64,
    pub b_ref: f64,
    pub c_ref: f64,
    pub inertia: Inertia,
}

impl<'a> Ctx<'a> {
    pub fn new(plane: &'a Plane, db: &'a FoilDb, def: &'a PlanePolarDef) -> Result<Ctx<'a>, PlaneError> {
        let (s, b, c) = match def.reference {
            RefDims::Projected => plane.reference(true),
            RefDims::Planform => plane.reference(false),
            RefDims::Custom { area, span, chord } => (area, span, chord),
        };
        if !(s > 0.0 && b > 0.0 && c > 0.0) {
            return Err(PlaneError::Input(format!("reference dimensions must be positive (area {s}, span {b}, chord {c})")));
        }
        let inertia = def.inertia.unwrap_or_else(|| plane.inertia(db));
        Ok(Ctx { plane, db, def, s_ref: s, b_ref: b, c_ref: c, inertia })
    }

    /// Viscous drag, transition and section Cm for a strip: interpolated
    /// between the strip's two end foils.  Returns None if no polar exists.
    pub fn strip_viscous(&self, foil_a: &str, foil_b: &str, tau: f64, re: f64, cl: f64, alpha_eff: f64) -> Option<(f64, f64, f64, f64, bool)> {
        let get = |f: &str| if self.def.viscous_from_cl { self.db.at_re_cl(f, re, cl) } else { self.db.at_re_alpha(f, re, alpha_eff) };
        let (pa, oa) = get(foil_a)?;
        let (pb, ob) = if foil_b == foil_a { (pa, oa) } else { get(foil_b)? };
        let l = |x: f64, y: f64| x + tau * (y - x);
        Some((l(pa.cd, pb.cd), l(pa.xtr_top, pb.xtr_top), l(pa.xtr_bot, pb.xtr_bot), l(pa.cm, pb.cm), oa || ob))
    }

    /// Reduce raw loads to an operating point.
    pub fn reduce(&self, mesh: &Mesh, raw: &RawLoads, cog: Vec3, alpha: f64, beta: f64, v: f64, mass: f64) -> PlaneOpp {
        let rho = self.def.fluid.rho;
        let nu = self.def.fluid.nu();
        let q = 0.5 * rho * v * v;
        let (wd, wn, ws) = (wind_dir(alpha, beta), wind_normal(alpha), wind_side(alpha, beta));
        let qs = q * self.s_ref;
        let mut warnings = Vec::new();
        let mut strips = Vec::with_capacity(mesh.strips.len());
        let mut f_visc = Vec3::ZERO;
        let mut m_visc = Vec3::ZERO;
        let mut missing_polars: Vec<String> = Vec::new();
        let mut n_out = 0;
        for (i, s) in mesh.strips.iter().enumerate() {
            let f = raw.strip_force[i] * rho;
            let area = s.area.max(1e-30);
            let cl = f.dot(s.normal) / (q * area);
            let wi = raw.strip_wi.get(i).cloned().unwrap_or_default();
            let ai = (wi.dot(s.normal)).atan2(v) / DEG;
            let gamma = raw.strip_gamma.get(i).cloned().unwrap_or(0.0);
            let fi = wi.cross(s.te_b - s.te_a) * (gamma * rho);
            let cdi = fi.dot(wd) / (q * area);
            let re = v * s.chord / nu;
            let alpha_eff = alpha + s.twist + ai;
            let mut sr = StripResult {
                wing: s.wing,
                span_pos: s.span_pos,
                chord: s.chord,
                cl,
                cdi,
                ai,
                re,
                xtr_top: 1.0,
                xtr_bot: 1.0,
                gamma: gamma * v.signum(),
                lift_per_span: f.dot(wn) / s.width.max(1e-30),
                twist: s.twist,
                alpha_eff,
                ..Default::default()
            };
            if self.def.viscous {
                match self.strip_viscous(&s.foil_a, &s.foil_b, s.tau, re, cl, alpha_eff) {
                    Some((cd, xt, xb, cm, out)) => {
                        sr.cdv = cd;
                        sr.xtr_top = xt;
                        sr.xtr_bot = xb;
                        sr.cm_c4 = cm;
                        sr.out_of_polar = out;
                        if out {
                            n_out += 1;
                        }
                        let dfv = wd * (q * area * cd);
                        f_visc += dfv;
                        m_visc += (s.c4 - cog).cross(dfv);
                    }
                    None => {
                        for f in [&s.foil_a, &s.foil_b] {
                            if !missing_polars.contains(f) {
                                missing_polars.push(f.clone());
                            }
                        }
                    }
                }
            }
            strips.push(sr);
        }
        if !missing_polars.is_empty() {
            warnings.push(format!("no 2D polars for foil(s) {}: viscous drag not counted for those strips", missing_polars.join(", ")));
        }
        if n_out > 0 {
            warnings.push(format!("{n_out} strip(s) outside the 2D polar envelope (Re or CL out of range): viscous drag extrapolated"));
        }
        // bending moments (per wing, from the tip inwards, about the strip centre)
        for w in &mesh.wings {
            let r = w.strips.clone();
            let idx: Vec<usize> = r.collect();
            for (k, &i) in idx.iter().enumerate() {
                let si = &mesh.strips[i];
                let side_left = si.left;
                let mut bm = 0.0;
                for &j in &idx[..] {
                    let sj = &mesh.strips[j];
                    let outboard = if side_left { j < idx[k] || (j == i) } else { j > i || j == i };
                    if sj.left == side_left && outboard && j != i {
                        let arm = (sj.c4 - si.c4).dot(Vec3::Y).abs().max((sj.c4 - si.c4).dot(Vec3::Z).abs());
                        bm += strips[j].lift_per_span * sj.width * arm;
                    }
                }
                strips[i].bending = bm;
            }
        }
        // extra drags
        let mut d_extra = 0.0;
        for e in &self.def.extra_drag {
            d_extra += q * e.area * e.cd;
        }
        if self.def.fuse_drag {
            for f in &self.plane.fuses {
                let len = f.length();
                let re = v * len / nu;
                if re > 1e3 {
                    let cf = 0.455 / re.log10().powf(2.58);
                    let d = f.fineness_inverse();
                    let ff = 1.0 + 1.5 * d.powf(1.5) + 50.0 * d.powi(3);
                    d_extra += q * cf * ff * f.wetted_area();
                }
            }
        }
        let f_inv = raw.force * rho;
        let f_ind = raw.induced_force * rho;
        let f_tot = f_inv + f_ind + f_visc + wd * d_extra;
        let m_tot = raw.moment * rho + m_visc;
        let lift = f_tot.dot(wn);
        let cl = lift / qs;
        let cdi = f_ind.dot(wd) / qs;
        let cdv = (f_visc.dot(wd) + d_extra) / qs;
        let cd = cdi + cdv;
        let cm = m_tot.dot(ws) / (qs * self.c_ref);
        let cr = (raw.moment * rho).dot(wd) / (qs * self.b_ref);
        let cn = m_tot.dot(wn) / (qs * self.b_ref);
        let cy = f_tot.dot(ws) / qs;
        let x_cp = if lift.abs() > 1e-9 { cog.x - (raw.moment * rho).dot(Vec3::Y) / f_inv.dot(wn) } else { cog.x };
        let ar = self.b_ref * self.b_ref / self.s_ref;
        let efficiency = if cdi.abs() > 1e-12 { cl * cl / (std::f64::consts::PI * ar * cdi) } else { 0.0 };
        // root bending of the main wing = max over its strips
        let root_bending = mesh.wings.iter().enumerate().find(|(_, w)| w.kind == crate::geom::WingKind::Main).map(|(_, w)| w.strips.clone().map(|i| strips[i].bending.abs()).fold(0.0, f64::max)).unwrap_or(0.0);
        PlaneOpp {
            alpha,
            beta,
            velocity: v,
            ctrl: 0.0,
            cl,
            cd,
            cdi,
            cdv,
            cy,
            cm,
            cr,
            cn,
            lift,
            drag: cd * qs,
            x_cp,
            efficiency,
            mass,
            cog,
            strips,
            cp: if self.def.keep_cp { raw.cp.iter().map(|c| *c as f32).collect() } else { vec![] },
            warnings,
            stability: None,
            root_bending,
            wake: vec![],
        }
    }
}

/// Build the mesh options for a definition and control value.
pub fn mesh_options(def: &PlanePolarDef, extra_flaps: &HashMap<(usize, usize), f64>, tilts: &HashMap<usize, f64>) -> MeshOptions {
    let mut flaps: HashMap<(usize, usize), f64> = def.flaps.iter().cloned().collect();
    for (k, v) in extra_flaps {
        *flaps.entry(*k).or_insert(0.0) += v;
    }
    MeshOptions { thick: def.method.is_panel(), include_fuse: def.include_fuse, flaps, tilts: tilts.clone(), allow_missing_foils: def.allow_missing_foils, triangulate: def.method == Method::TriPanel, fuse_in_thin: false }
}

/// Build and factor the solver for a mesh.
pub fn make_system<'a>(def: &PlanePolarDef, mesh: &'a Mesh, cog: Vec3, c_ref: f64, prog: &Progress) -> Result<System<'a>, PlaneError> {
    Ok(match def.method {
        Method::Vlm1 | Method::Vlm2 => System::Vlm(VlmSystem::build(mesh, def.method == Method::Vlm2, cog, def.ground_height, prog)?),
        Method::Panel | Method::TriPanel => System::Panel(PanelSystem::build(mesh, cog, def.wake_length * c_ref, def.wake_panels, def.wake_ratio, def.ground_height, prog)?),
        Method::Llt => return Err(PlaneError::Input("LLT has no panel system".into())),
    })
}

/// Speed at which lift equals weight for the loads at unit speed.
fn balance_speed(ctx: &Ctx, raw1: &RawLoads, alpha: f64, mass: f64) -> Option<f64> {
    let l1 = (raw1.force + raw1.induced_force).dot(wind_normal(alpha)) * ctx.def.fluid.rho;
    if l1 <= 1e-12 || mass <= 0.0 {
        return None;
    }
    Some((mass * G / l1).sqrt())
}

/// Run a 3D polar.  Points are streamed to `on_point` as they are computed.
pub fn run_polar(plane: &Plane, db: &FoilDb, def: &PlanePolarDef, prog: &Progress, on_point: &(dyn Fn(&PlaneOpp) + Sync)) -> Result<PlanePolar, PlaneError> {
    def.validate().map_err(PlaneError::Input)?;
    plane.validate().map_err(PlaneError::Geometry)?;
    let ctx = Ctx::new(plane, db, def)?;
    let mass = ctx.inertia.mass;
    if def.kind.needs_mass() && mass <= 0.0 {
        return Err(PlaneError::Input(format!("{} needs the plane mass: set structural/point masses or an inertia override", def.kind.label())));
    }
    let cog = if mass > 0.0 { ctx.inertia.cog } else { Vec3::ZERO };
    let name = if def.name.is_empty() { def.auto_name() } else { def.name.clone() };
    let mut polar = PlanePolar { name: name.clone(), plane: plane.name.clone(), def: Some(def.clone()), ref_dims: (ctx.s_ref, ctx.b_ref, ctx.c_ref), ..Default::default() };
    prog.stage(format!("{}: {}", plane.name, name));
    prog.info(format!("{} - {} - {}", plane.name, def.kind.label(), def.method.label()));
    if mass > 0.0 {
        prog.info(format!("mass {:.3} kg, CoG ({:.4}, {:.4}, {:.4}) m", mass, cog.x, cog.y, cog.z));
    }
    if def.method == Method::Llt {
        let pts = crate::llt::run(&ctx, prog, on_point, &mut polar.failed)?;
        polar.points = pts;
        polar.unknowns = def.llt.stations;
        return Ok(polar);
    }
    match def.kind {
        PolarKind::Control => run_control(&ctx, &mut polar, prog, on_point)?,
        PolarKind::Stability => run_stability(&ctx, &mut polar, prog, on_point)?,
        _ => {
            prog.stage("Meshing");
            let mesh = build(plane, db, &mesh_options(def, &HashMap::new(), &HashMap::new()))?;
            polar.unknowns = mesh.panels.len();
            prog.info(format!("mesh: {} panels, {} strips", mesh.panels.len(), mesh.strips.len()));
            let sys = make_system(def, &mesh, cog, ctx.c_ref, &prog.sub(0.0, 0.5))?;
            let ops: Vec<(f64, f64, f64)> = match def.kind {
                PolarKind::FixedSpeed | PolarKind::FixedLift | PolarKind::Glide => def.range.values().into_iter().map(|a| (a, def.beta, def.velocity)).collect(),
                PolarKind::FixedAlpha => def.range.values().into_iter().map(|v| (def.alpha, def.beta, v)).collect(),
                PolarKind::Sideslip => def.range.values().into_iter().map(|b| (def.alpha, b, def.velocity)).collect(),
                PolarKind::List => def.list.clone(),
                _ => unreachable!(),
            };
            let rest = prog.sub(0.5, 1.0);
            let n = ops.len();
            let results: Vec<Result<Option<PlaneOpp>, PlaneError>> = {
                use rayon::prelude::*;
                ops.par_iter()
                    .enumerate()
                    .map(|(k, &(a, b, v))| {
                        rest.check()?;
                        let r = solve_point(&ctx, &sys, &mesh, cog, mass, a, b, v);
                        rest.part(k, n).set(1.0);
                        match &r {
                            Ok(Some(p)) => {
                                rest.count_ok();
                                rest.debug(format!("alpha={:.2} beta={:.2} V={:.2}: CL={:.4} CD={:.5} Cm={:.4}", p.alpha, p.beta, p.velocity, p.cl, p.cd, p.cm));
                                on_point(p);
                            }
                            Ok(None) => rest.count_failed(),
                            Err(_) => {}
                        }
                        r
                    })
                    .collect()
            };
            for (k, r) in results.into_iter().enumerate() {
                match r? {
                    Some(p) => polar.points.push(p),
                    None => polar.failed.push(format!("alpha={:.2} beta={:.2} V={:.2}: no positive lift for the speed balance", ops[k].0, ops[k].1, ops[k].2)),
                }
            }
        }
    }
    let key = |p: &PlaneOpp| match def.kind {
        PolarKind::FixedAlpha => p.velocity,
        PolarKind::Sideslip => p.beta,
        PolarKind::Control | PolarKind::Stability => p.ctrl,
        PolarKind::List => 0.0,
        _ => p.alpha,
    };
    if def.kind != PolarKind::List {
        polar.points.sort_by(|a, b| key(a).total_cmp(&key(b)));
    }
    for f in &polar.failed {
        prog.warn(f.clone());
    }
    prog.set(1.0);
    Ok(polar)
}

#[allow(clippy::too_many_arguments)]
fn solve_point(ctx: &Ctx, sys: &System, mesh: &Mesh, cog: Vec3, mass: f64, alpha: f64, beta: f64, v_spec: f64) -> Result<Option<PlaneOpp>, PlaneError> {
    let def = ctx.def;
    let wd = wind_dir(alpha, beta);
    let v = match def.kind {
        PolarKind::FixedLift | PolarKind::Glide => {
            let raw1 = sys.loads(wd, Vec3::ZERO);
            match balance_speed(ctx, &raw1, alpha, mass) {
                Some(v) => v,
                None => return Ok(None),
            }
        }
        _ => v_spec,
    };
    let (raw, wake) = match (sys, def.free_wake) {
        (System::Vlm(s), it) if it > 0 => {
            let (g, fw) = s.solve_free_wake(wd * v, it, def.wake_length * ctx.c_ref, 24, &Progress::none())?;
            (s.loads(&g, wd * v, Vec3::ZERO), fw.lines.iter().map(|l| l.iter().map(|p| [p.x as f32, p.y as f32, p.z as f32]).collect()).collect())
        }
        _ => (sys.loads(wd * v, Vec3::ZERO), vec![]),
    };
    let mut opp = ctx.reduce(mesh, &raw, cog, alpha, beta, v, mass);
    opp.wake = wake;
    if def.kind == PolarKind::Glide {
        // speed from the total aerodynamic force (viscous drag depends on Re)
        let mut vv = v;
        for _ in 0..30 {
            let ca = (opp.cl * opp.cl + opp.cd * opp.cd).sqrt();
            if ca <= 1e-9 {
                return Ok(None);
            }
            let vn = (2.0 * mass * G / (def.fluid.rho * ctx.s_ref * ca)).sqrt();
            let done = (vn - vv).abs() < 0.005;
            vv = vn;
            let raw = sys.loads(wd * vv, Vec3::ZERO);
            opp = ctx.reduce(mesh, &raw, cog, alpha, beta, vv, mass);
            if done {
                break;
            }
        }
    }
    if def.derivatives && mass > 0.0 {
        opp.stability = stability::analyse(ctx, sys, cog, alpha, opp.velocity).ok();
    }
    Ok(Some(opp))
}

fn control_value(def: &PlanePolarDef, t: f64, target: ControlTarget) -> Option<f64> {
    def.ranges.iter().find(|r| r.target == target).map(|r| r.min + t * (r.max - r.min))
}

fn run_control(ctx: &Ctx, polar: &mut PlanePolar, prog: &Progress, on_point: &(dyn Fn(&PlaneOpp) + Sync)) -> Result<(), PlaneError> {
    let def = ctx.def;
    let ts = def.range.values();
    let n = ts.len();
    for (k, &t) in ts.iter().enumerate() {
        prog.check()?;
        let sub = prog.part(k, n);
        sub.stage(format!("Control parameter {t:.3} ({}/{n})", k + 1));
        let v = control_value(def, t, ControlTarget::Velocity).unwrap_or(def.velocity);
        let a = control_value(def, t, ControlTarget::Alpha).unwrap_or(def.alpha);
        let b = control_value(def, t, ControlTarget::Beta).unwrap_or(def.beta);
        let mut inertia = ctx.inertia;
        if let Some(m) = control_value(def, t, ControlTarget::Mass) {
            inertia.mass = m;
        }
        if let Some(x) = control_value(def, t, ControlTarget::CogX) {
            inertia.cog.x = x;
        }
        if let Some(z) = control_value(def, t, ControlTarget::CogZ) {
            inertia.cog.z = z;
        }
        let mut flaps = HashMap::new();
        let mut tilts = HashMap::new();
        for r in &def.ranges {
            let val = r.min + t * (r.max - r.min);
            match r.target {
                ControlTarget::Flap { wing, index } => {
                    flaps.insert((wing, index), val);
                }
                ControlTarget::Tilt { wing } => {
                    tilts.insert(wing, val);
                }
                _ => {}
            }
        }
        let mesh = build(ctx.plane, ctx.db, &mesh_options(def, &flaps, &tilts))?;
        polar.unknowns = mesh.panels.len();
        let sys = make_system(def, &mesh, inertia.cog, ctx.c_ref, &sub.sub(0.0, 0.9))?;
        let raw = sys.loads(wind_dir(a, b) * v, Vec3::ZERO);
        let mut opp = ctx.reduce(&mesh, &raw, inertia.cog, a, b, v, inertia.mass);
        opp.ctrl = t;
        if def.derivatives && inertia.mass > 0.0 {
            let c2 = Ctx { inertia, ..Ctx::new(ctx.plane, ctx.db, def)? };
            opp.stability = stability::analyse(&c2, &sys, inertia.cog, a, v).ok();
        }
        sub.set(1.0);
        prog.count_ok();
        on_point(&opp);
        polar.points.push(opp);
    }
    Ok(())
}

fn run_stability(ctx: &Ctx, polar: &mut PlanePolar, prog: &Progress, on_point: &(dyn Fn(&PlaneOpp) + Sync)) -> Result<(), PlaneError> {
    let def = ctx.def;
    let ts = def.range.values();
    let n = ts.len();
    let cog = ctx.inertia.cog;
    let mass = ctx.inertia.mass;
    for (k, &t) in ts.iter().enumerate() {
        prog.check()?;
        let sub = prog.part(k, n);
        sub.stage(format!("Stability: control {t:.3} ({}/{n}) - trimming", k + 1));
        let mut flaps = HashMap::new();
        let mut tilts = HashMap::new();
        for (target, gain) in &def.gains {
            match target {
                ControlTarget::Flap { wing, index } => {
                    flaps.insert((*wing, *index), gain * t);
                }
                ControlTarget::Tilt { wing } => {
                    tilts.insert(*wing, gain * t);
                }
                _ => {}
            }
        }
        let mesh = build(ctx.plane, ctx.db, &mesh_options(def, &flaps, &tilts))?;
        polar.unknowns = mesh.panels.len();
        let sys = make_system(def, &mesh, cog, ctx.c_ref, &sub.sub(0.0, 0.7))?;
        match stability::trim(ctx, &sys, cog, def.beta) {
            Some((alpha, u0)) => {
                let raw = sys.loads(wind_dir(alpha, def.beta) * u0, Vec3::ZERO);
                let mut opp = ctx.reduce(&mesh, &raw, cog, alpha, def.beta, u0, mass);
                opp.ctrl = t;
                sub.stage(format!("Stability: control {t:.3} - derivatives and modes"));
                match stability::analyse(ctx, &sys, cog, alpha, u0) {
                    Ok(mut s) => {
                        // control derivatives: rebuild with the control parameter nudged
                        let gmax = def.gains.iter().map(|g| g.1.abs()).fold(0.0, f64::max);
                        if gmax > 0.0 {
                            let dt = 0.5 / gmax;
                            let mut f2 = HashMap::new();
                            let mut t2 = HashMap::new();
                            for (target, gain) in &def.gains {
                                match target {
                                    ControlTarget::Flap { wing, index } => {
                                        f2.insert((*wing, *index), gain * (t + dt));
                                    }
                                    ControlTarget::Tilt { wing } => {
                                        t2.insert(*wing, gain * (t + dt));
                                    }
                                    _ => {}
                                }
                            }
                            if let Ok(m2) = build(ctx.plane, ctx.db, &mesh_options(def, &f2, &t2)) {
                                if let Ok(s2) = make_system(def, &m2, cog, ctx.c_ref, &Progress::none()) {
                                    stability::add_control_derivatives(ctx, &sys, &s2, &mut s, dt);
                                }
                            }
                        }
                        opp.stability = Some(s)
                    }
                    Err(e) => opp.warnings.push(format!("stability analysis failed: {e}")),
                }
                prog.count_ok();
                prog.info(format!("ctrl {t:.3}: trimmed at alpha={alpha:.3} deg, V={u0:.3} m/s, CL={:.4}", opp.cl));
                on_point(&opp);
                polar.points.push(opp);
            }
            None => {
                prog.count_failed();
                polar.failed.push(format!("ctrl {t:.3}: no trim found (Cm does not cross zero between -30 and 30 deg, or no positive lift at trim)"));
            }
        }
        sub.set(1.0);
    }
    let _ = mass;
    Ok(())
}
