//! 3D polars: definitions (types 1-8), the drivers, and the reduction of
//! raw loads to coefficients and span distributions.

use crate::evaluator::{Evaluator, WakeReport};
use crate::foildb::{Envelope, FamilyId, FamilySelector, FoilDb};
use crate::geom::{Inertia, Plane};
use crate::linear_panel::LinearPanelSystem;
use crate::mesh::{build, Mesh, MeshOptions, Strip};
use crate::panel::PanelSystem;
use crate::stability::{self, StabilityModel, StabilityResult};
use crate::vlm::VlmSystem;
use crate::PlaneError;
use foiley_core::atmosphere::Fluid;
use foiley_core::linalg::{Diagnostics, SolveReport};
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
    /// Discontinuous linear-doublet triangles with Galerkin testing.
    TriLinear,
}

impl Method {
    pub const ALL: [Method; 6] = [Method::Llt, Method::Vlm1, Method::Vlm2, Method::Panel, Method::TriPanel, Method::TriLinear];
    /// Thick-surface methods (panels on both sides of the wings, fuselage included).
    pub fn is_panel(self) -> bool {
        matches!(self, Method::Panel | Method::TriPanel | Method::TriLinear)
    }
    /// Number of algebraic unknowns for a panel mesh (LLT uses its station count).
    pub fn unknowns(self, panels: usize) -> usize {
        panels * if self == Method::TriLinear { 3 } else { 1 }
    }
    pub fn label(self) -> &'static str {
        match self {
            Method::Llt => "LLT - nonlinear lifting line (main wing only)",
            Method::Vlm1 => "VLM1 - horseshoe vortices",
            Method::Vlm2 => "VLM2 - ring vortices (recommended)",
            Method::Panel => "3D panels - thick surfaces + fuselage",
            Method::TriPanel => "3D panels on triangles (uniform doublets)",
            Method::TriLinear => "3D panels on triangles (linear doublets, Galerkin)",
        }
    }
    /// The one compatibility table between methods, polar types and optional
    /// settings.  `PlanePolarDef::validate`, the GUI (which greys out what is
    /// unsupported) and the CLI all go through it.
    pub fn supports(self, kind: PolarKind, setting: Setting) -> bool {
        use PolarKind::*;
        let vlm = matches!(self, Method::Vlm1 | Method::Vlm2);
        let llt = self == Method::Llt;
        match setting {
            Setting::PolarType => !llt || matches!(kind, FixedSpeed | FixedLift | Glide | FixedAlpha),
            // every VLM driver goes through the wake-aware evaluator
            Setting::FreeWake => vlm,
            // the lifting line is a symmetric main-wing model in free air:
            // no sideslip, ground, control surfaces, fuselage or derivatives
            Setting::Sideslip | Setting::GroundEffect | Setting::Flaps | Setting::Tilts | Setting::FuseDrag | Setting::Derivatives => !llt,
        }
    }
    pub fn short(self) -> &'static str {
        match self {
            Method::Llt => "LLT",
            Method::Vlm1 => "VLM1",
            Method::Vlm2 => "VLM2",
            Method::Panel => "Panel",
            Method::TriPanel => "TriPanel",
            Method::TriLinear => "TriLinear",
        }
    }
}

/// Optional parts of a polar definition that not every method honours.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Setting {
    /// the polar type itself
    PolarType,
    /// nonzero sideslip (beta, T5 sweep, T6 beta range, T8 betas)
    Sideslip,
    GroundEffect,
    /// flap deflections (presets, T6 ranges, T7 gains)
    Flaps,
    /// wing tilts (T6 ranges, T7 gains)
    Tilts,
    FuseDrag,
    Derivatives,
    /// free-wake relaxation
    FreeWake,
}

impl Setting {
    pub const ALL: [Setting; 8] = [Setting::PolarType, Setting::Sideslip, Setting::GroundEffect, Setting::Flaps, Setting::Tilts, Setting::FuseDrag, Setting::Derivatives, Setting::FreeWake];
    pub fn label(self) -> &'static str {
        match self {
            Setting::PolarType => "this polar type",
            Setting::Sideslip => "sideslip",
            Setting::GroundEffect => "ground effect",
            Setting::Flaps => "flap deflections",
            Setting::Tilts => "wing tilts",
            Setting::FuseDrag => "fuselage friction drag",
            Setting::Derivatives => "stability derivatives",
            Setting::FreeWake => "free-wake relaxation",
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
    /// loads behind the stability derivatives and the T7 trim: inviscid
    /// (default) or the full reduced loads (viscous strip drag and section
    /// Cm, extra drag, fuselage friction)
    pub stability_model: StabilityModel,
    /// wake length for the panel method, in reference chords
    pub wake_length: f64,
    pub wake_panels: usize,
    pub wake_ratio: f64,
    pub llt: LltSettings,
    /// keep panel Cp in the results (memory heavy for big meshes)
    pub keep_cp: bool,
    /// keep the span distributions (strips) in the results; turn off for
    /// large sweeps, where they dominate memory and project size
    pub keep_strips: bool,
    /// VLM free-wake relaxation: the iteration cap (0 = rigid wake along x)
    pub free_wake: usize,
    /// free wake: converged when the max node displacement / c_ref and
    /// ||dGamma||inf / ||Gamma||inf of a step are both below this
    pub wake_tol: f64,
    /// free wake: segments per trailing line
    pub wake_nseg: usize,
    /// free wake: under-relaxation of the node update, in (0, 1]
    pub wake_relax: f64,
    /// free wake: fail a point whose wake did not converge (default: warn)
    pub strict_wake: bool,
    /// free wake: re-relax the wake for each of the 12 derivative
    /// perturbations (12 relaxations per point) instead of freezing the
    /// converged wake of the reference state and superposing unit solutions
    pub relaxed_derivatives: bool,
    /// linear-solve and force-balance checks (off by default)
    pub diagnostics: Diagnostics,
    /// treat foils without loaded geometry as flat plates
    pub allow_missing_foils: bool,
    /// fail a point if any strip's viscous interpolation is clamped to the
    /// 2D polar envelope (or has no polar), instead of warning
    pub strict_envelope: bool,
    /// which 2D polar family of each foil the viscous interpolation (and
    /// the lifting line) uses; must select exactly one per foil
    pub family: FamilySelector,
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
            stability_model: StabilityModel::Inviscid,
            wake_length: 30.0,
            wake_panels: 5,
            wake_ratio: 1.1,
            llt: LltSettings::default(),
            keep_cp: true,
            keep_strips: true,
            allow_missing_foils: false,
            strict_envelope: false,
            family: FamilySelector::default(),
            free_wake: 0,
            wake_tol: 1e-3,
            wake_nseg: 24,
            wake_relax: 0.5,
            strict_wake: false,
            relaxed_derivatives: false,
            diagnostics: Diagnostics::Off,
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

    /// Whether this definition asks for `setting` (see `Method::supports`).
    pub fn uses(&self, setting: Setting) -> bool {
        let range_on = |f: &dyn Fn(&ControlTarget) -> bool| self.ranges.iter().any(|r| f(&r.target) && (r.min != 0.0 || r.max != 0.0));
        let gain_on = |f: &dyn Fn(&ControlTarget) -> bool| self.gains.iter().any(|(t, g)| f(t) && *g != 0.0);
        match setting {
            Setting::PolarType => true,
            Setting::Sideslip => match self.kind {
                PolarKind::Sideslip => true,
                PolarKind::List => self.list.iter().any(|p| p.1 != 0.0),
                PolarKind::Control => self.beta != 0.0 || range_on(&|t| *t == ControlTarget::Beta),
                _ => self.beta != 0.0,
            },
            Setting::GroundEffect => self.ground_height.is_some(),
            Setting::Flaps => {
                let fl = |t: &ControlTarget| matches!(t, ControlTarget::Flap { .. });
                self.flaps.iter().any(|f| f.1 != 0.0) || (self.kind == PolarKind::Control && range_on(&fl)) || (self.kind == PolarKind::Stability && gain_on(&fl))
            }
            Setting::Tilts => {
                let tl = |t: &ControlTarget| matches!(t, ControlTarget::Tilt { .. });
                (self.kind == PolarKind::Control && range_on(&tl)) || (self.kind == PolarKind::Stability && gain_on(&tl))
            }
            Setting::FuseDrag => self.fuse_drag,
            Setting::Derivatives => self.derivatives,
            Setting::FreeWake => self.free_wake > 0,
        }
    }

    /// Settings this definition uses that its method does not support.
    pub fn unsupported(&self) -> Vec<Setting> {
        Setting::ALL.into_iter().filter(|s| self.uses(*s) && !self.method.supports(self.kind, *s)).collect()
    }

    /// Reset every setting the method does not support to its default (the
    /// GUI calls this when the method or polar type changes, so a greyed-out
    /// option cannot keep a value that would fail validation).
    pub fn clear_unsupported(&mut self) {
        let (m, k) = (self.method, self.kind);
        let ok = |s| m.supports(k, s);
        if !ok(Setting::Sideslip) {
            self.beta = 0.0;
            self.list.iter_mut().for_each(|p| p.1 = 0.0);
            self.ranges.retain(|r| r.target != ControlTarget::Beta);
        }
        if !ok(Setting::GroundEffect) {
            self.ground_height = None;
        }
        if !ok(Setting::Flaps) {
            self.flaps.clear();
            self.ranges.retain(|r| !matches!(r.target, ControlTarget::Flap { .. }));
            self.gains.retain(|g| !matches!(g.0, ControlTarget::Flap { .. }));
        }
        if !ok(Setting::Tilts) {
            self.ranges.retain(|r| !matches!(r.target, ControlTarget::Tilt { .. }));
            self.gains.retain(|g| !matches!(g.0, ControlTarget::Tilt { .. }));
        }
        if !ok(Setting::FuseDrag) {
            self.fuse_drag = false;
        }
        if !ok(Setting::Derivatives) {
            self.derivatives = false;
        }
        if !ok(Setting::FreeWake) {
            self.free_wake = 0;
        }
    }

    /// Check the definition, including every operating point it expands to.
    pub fn validate(&self) -> Result<(), String> {
        let fin = |name: &str, v: f64| if v.is_finite() { Ok(()) } else { Err(format!("{name} must be a finite number (got {v})")) };
        fin("speed", self.velocity)?;
        fin("alpha", self.alpha)?;
        fin("sideslip", self.beta)?;
        fin("range start", self.range.from)?;
        fin("range end", self.range.to)?;
        fin("range step", self.range.step)?;
        fin("fluid density", self.fluid.rho)?;
        fin("fluid viscosity", self.fluid.mu)?;
        fin("speed of sound", self.fluid.a)?;
        fin("wake length", self.wake_length)?;
        fin("wake ratio", self.wake_ratio)?;
        fin("LLT relaxation", self.llt.relax)?;
        fin("LLT tolerance", self.llt.tolerance)?;
        if let RefDims::Custom { area, span, chord } = self.reference {
            for (n, v) in [("reference area", area), ("reference span", span), ("reference chord", chord)] {
                fin(n, v)?;
                if v <= 0.0 {
                    return Err(format!("{n} must be positive"));
                }
            }
        }
        if let Some(i) = &self.inertia {
            for (n, v) in [("mass", i.mass), ("CoG x", i.cog.x), ("CoG y", i.cog.y), ("CoG z", i.cog.z), ("Ixx", i.ixx), ("Iyy", i.iyy), ("Izz", i.izz), ("Ixz", i.ixz)] {
                fin(n, v)?;
            }
            if i.mass < 0.0 {
                return Err("mass must not be negative".into());
            }
        }
        if let Some(h) = self.ground_height {
            fin("ground height", h)?;
        }
        for e in &self.extra_drag {
            fin(&format!("extra drag '{}' area", e.name), e.area)?;
            fin(&format!("extra drag '{}' CD", e.name), e.cd)?;
        }
        for (_, v) in &self.flaps {
            fin("flap deflection", *v)?;
        }
        for (_, g) in &self.gains {
            fin("control gain", *g)?;
        }
        for r in &self.ranges {
            fin("control range minimum", r.min)?;
            fin("control range maximum", r.max)?;
        }
        if !(self.fluid.rho > 0.0 && self.fluid.mu > 0.0) {
            return Err("fluid density and viscosity must be positive".into());
        }
        if !(self.wake_length > 0.0) || self.wake_panels == 0 || !(self.wake_ratio > 0.0) {
            return Err("wake length, ratio and panel count must be positive".into());
        }
        if self.free_wake > 0 {
            fin("free-wake tolerance", self.wake_tol)?;
            fin("free-wake relaxation", self.wake_relax)?;
            if self.wake_tol < 0.0 {
                return Err("free-wake tolerance must not be negative".into());
            }
            if !(self.wake_relax > 0.0 && self.wake_relax <= 1.0) {
                return Err(format!("free-wake relaxation must be in (0, 1] (got {})", self.wake_relax));
            }
            if self.wake_nseg < 2 {
                return Err(format!("free wake needs at least 2 segments per line (got {})", self.wake_nseg));
            }
        }
        if let Some(s) = self.unsupported().first() {
            return Err(if *s == Setting::PolarType {
                format!("{} does not support {}", self.method.short(), self.kind.label())
            } else {
                format!("{} does not support {} in {}", self.method.short(), s.label(), self.kind.label())
            });
        }
        let v_range = self.ranges.iter().find(|r| r.target == ControlTarget::Velocity);
        match self.kind {
            PolarKind::FixedSpeed | PolarKind::Sideslip if !(self.velocity > 0.0) => return Err("speed must be positive".into()),
            PolarKind::Control => match v_range {
                Some(r) if !(r.min > 0.0 && r.max > 0.0) => return Err(format!("T6 speed range must stay positive (from {} to {} m/s)", r.min, r.max)),
                None if !(self.velocity > 0.0) => return Err("speed must be positive".into()),
                _ => {}
            },
            PolarKind::FixedAlpha if self.range.values().iter().any(|v| *v <= 0.0) => return Err("speed range of a T4 polar must be positive".into()),
            PolarKind::List => {
                if self.list.is_empty() {
                    return Err("T8 polar needs a list of (alpha, beta, speed)".into());
                }
                for (k, &(a, b, v)) in self.list.iter().enumerate() {
                    if !(a.is_finite() && b.is_finite() && v.is_finite()) {
                        return Err(format!("T8 point {} has a non-finite value ({a}, {b}, {v})", k + 1));
                    }
                    if v <= 0.0 {
                        return Err(format!("T8 point {} needs a positive speed (got {v} m/s)", k + 1));
                    }
                }
            }
            _ => {}
        }
        if let Some(r) = self.ranges.iter().find(|r| r.target == ControlTarget::Mass) {
            if r.min < 0.0 || r.max < 0.0 {
                return Err("T6 mass range must not be negative".into());
            }
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
    /// total force from the near field (integrated pressures, or
    /// Kutta-Joukowski on the bound vortices for VLM), per unit density
    pub pressure_force: Vec3,
    /// total force from Kutta-Joukowski on the shed (wake) circulation, per unit density
    pub wake_force: Vec3,
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
        let ind = trefftz_force(mesh, &sg, &tre, wdir, &crate::vlm::carry_segments(mesh, &sg), |m| crate::vlm::trefftz_velocity(mesh, &sg, m, wdir, None));
        let v2 = vinf.norm2().max(1e-30);
        let cp = mesh.panels.iter().zip(&panel_force).map(|(p, f)| -f.dot(p.normal) / (0.5 * v2 * p.area.max(1e-30))).collect();
        let pressure_force = panel_force.iter().fold(Vec3::ZERO, |a, f| a + *f);
        let wake_force = mesh.strips.iter().zip(&sg).fold(Vec3::ZERO, |a, (s, g)| a + vinf.cross(s.te_b - s.te_a) * *g);
        RawLoads { panel_force, force, moment, strip_gamma: sg, strip_wi: tre, strip_force, induced_force: ind, cp, mu: vec![], pressure_force, wake_force }
    }
}

/// Induced-drag force from the half-Trefftz velocities `tre` at the strips and
/// `carry_w` at the centre of each fuselage carry-through segment (a, b, gamma).
pub(crate) fn trefftz_force(mesh: &Mesh, sg: &[f64], tre: &[Vec3], wdir: Vec3, carry: &[(Vec3, Vec3, f64)], carry_w: impl Fn(Vec3) -> Vec3) -> Vec3 {
    let mut ind = Vec3::ZERO;
    for (s, (g, w)) in mesh.strips.iter().zip(sg.iter().zip(tre)) {
        let f = w.cross(s.te_b - s.te_a) * *g;
        ind += wdir * f.dot(wdir);
    }
    // circulation carried across fuselages
    for &(a, b, g) in carry {
        let f = carry_w((a + b) * 0.5).cross(b - a) * g;
        ind += wdir * f.dot(wdir);
    }
    ind
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
    /// viscous interpolation outside the 2D polar envelope (clamped)
    pub out_of_polar: bool,
    /// `foildb::Envelope` of the viscous interpolation, as a u8
    #[serde(default)]
    pub envelope: u8,
    /// deflected strip with no 2D polar family at its flap angle: its
    /// section data came from the nearest family's angle (F003)
    #[serde(default)]
    pub flap_clamped: bool,
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
    /// side-force coefficient; None when the method has no lateral model (LLT)
    pub cy: Option<f64>,
    pub cm: f64,
    /// rolling moment coefficient; None when the method has no lateral model
    pub cr: Option<f64>,
    /// yawing moment coefficient; None when the method has no lateral model
    pub cn: Option<f64>,
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
    /// why `stability` is missing when derivatives were asked for
    #[serde(default)]
    pub stability_error: Option<String>,
    /// root bending moment of the main wing, N.m
    pub root_bending: f64,
    /// relaxed wake lines (free-wake VLM), for display
    #[serde(default)]
    pub wake: Vec<Vec<[f32; 3]>>,
    /// T3 speed balance: |resultant - weight| / weight at the reported speed
    #[serde(default)]
    pub balance_residual: Option<f64>,
    /// numerical diagnostics (when the definition asks for them)
    #[serde(default)]
    pub diagnostics: Option<PointDiagnostics>,
}

/// Numerical checks of one operating point (`PlanePolarDef.diagnostics`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct PointDiagnostics {
    /// min|u_ii| / max|u_ii| of the LU factors (a crude warning sign)
    #[serde(alias = "pivot_ratio")]
    pub diag_ratio: f64,
    /// Hager-Higham estimate of the 1-norm condition number
    pub cond1_estimate: Option<f64>,
    /// ||b - Ax||_inf / (||A||_inf ||x||_inf + ||b||_inf) on the checked rows
    /// (operator re-assembled in f64); None for a relaxed-wake solve, which
    /// re-factors per step and is not checked
    pub backward_error: Option<f64>,
    pub residual_inf: Option<f64>,
    /// rows checked (all of them for `Full`)
    pub rows_checked: usize,
    pub unknowns: usize,
    /// `Full`: iterative refinement kept for at least one unit solution
    pub refined: bool,
    /// near-field minus Trefftz-plane drag coefficient
    pub cd_near_minus_trefftz: f64,
    /// pressure (near-field) minus wake (Kutta-Joukowski on the shed
    /// circulation) lift coefficient
    pub cl_pressure_minus_wake: f64,
    /// relaxed-wake convergence; set on every relaxed-wake point, also with
    /// diagnostics off (the other fields are then zero)
    #[serde(default)]
    pub wake: Option<WakeReport>,
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
    /// provenance of the run (F016); None for polars from older versions
    #[serde(default)]
    pub manifest: Option<foiley_core::provenance::RunManifest>,
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
                "{:.3},{:.3},{:.4},{:.4},{:.6},{:.6},{:.6},{:.6},{},{:.6},{},{},{:.4},{:.4},{:.5},{:.5},{:.4},{:.5},{:.4}\n",
                p.alpha, p.beta, p.velocity, p.ctrl, p.cl, p.cd, p.cdi, p.cdv, opt6(p.cy), p.cm, opt6(p.cr), opt6(p.cn), p.ld(), p.lift, p.drag, p.x_cp, p.efficiency, p.vz(), p.root_bending
            ));
        }
        s
    }
}

/// A CSV field for an optional coefficient (empty when absent).
fn opt6(v: Option<f64>) -> String {
    v.map(|v| format!("{v:.6}")).unwrap_or_default()
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
    TriLinear(LinearPanelSystem<'a>),
}

impl<'a> System<'a> {
    pub fn loads(&self, vinf: Vec3, omega: Vec3) -> RawLoads {
        match self {
            System::Vlm(s) => {
                let g = s.solve(vinf, omega);
                s.loads(&g, vinf, omega)
            }
            System::Panel(s) => s.loads(vinf, omega),
            System::TriLinear(s) => s.loads(vinf, omega),
        }
    }
    pub fn mesh(&self) -> &Mesh {
        match self {
            System::Vlm(s) => s.lat.mesh,
            System::Panel(s) => s.mesh,
            System::TriLinear(s) => s.mesh,
        }
    }
    /// Linear-solve diagnostics of the factored system.
    pub fn report(&self) -> &SolveReport {
        match self {
            System::Vlm(s) => &s.report,
            System::Panel(s) => &s.report,
            System::TriLinear(s) => &s.report,
        }
    }
    /// The six unit solutions (u, v, w, p, q, r).
    pub fn unit(&self) -> &[Vec<f64>; 6] {
        match self {
            System::Vlm(s) => &s.unit,
            System::Panel(s) => &s.unit_mu,
            System::TriLinear(s) => &s.unit_mu,
        }
    }
    /// Diagnostics of a point solved with `vinf`, `omega` (`relaxed`: the
    /// circulation came from a relaxed-wake re-solve, not the unit solutions).
    pub fn point_diagnostics(&self, raw: &RawLoads, vinf: Vec3, omega: Vec3, s_ref: f64, relaxed: bool) -> PointDiagnostics {
        let rep = self.report();
        let c = [vinf.x, vinf.y, vinf.z, omega.x, omega.y, omega.z];
        let x_inf = (0..self.unit()[0].len()).map(|i| c.iter().zip(self.unit()).map(|(c, u)| c * u[i]).sum::<f64>().abs()).fold(0.0, f64::max);
        let qs = 0.5 * vinf.norm2() * s_ref;
        let a = vinf.normalized();
        let alpha = a.z.atan2(a.x) / DEG;
        let (wd, wn) = (a, wind_normal(alpha));
        let dp = raw.pressure_force - raw.wake_force;
        PointDiagnostics {
            diag_ratio: rep.diag_ratio,
            cond1_estimate: rep.cond1_estimate,
            backward_error: if relaxed { None } else { rep.backward_error(&c, x_inf) },
            residual_inf: if relaxed { None } else { rep.residual_inf(&c) },
            rows_checked: rep.residuals.rows.len(),
            unknowns: rep.n,
            refined: rep.refined.iter().any(|r| *r),
            cd_near_minus_trefftz: (raw.pressure_force - raw.induced_force).dot(wd) / qs,
            cl_pressure_minus_wake: dp.dot(wn) / qs,
            wake: None,
        }
    }
    pub fn cog(&self) -> Vec3 {
        match self {
            System::Vlm(s) => s.cog,
            System::Panel(s) => s.cog,
            System::TriLinear(s) => s.cog,
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
    /// polar family of every wing section's foil: `[wing][section][left]`
    /// (resolved once, so the per-strip lookup is two array indexes)
    pub section_family: Vec<Vec<[Option<FamilyId>; 2]>>,
    /// the same per foil name (setup code only)
    pub foil_family: Vec<(String, Option<FamilyId>)>,
    /// per resolved family (indexed by `FamilyId`): it and the families that
    /// differ from it only in flap angle, sorted by angle
    pub flap_sets: Vec<Vec<(f64, FamilyId)>>,
}

/// Section data a strip's viscous interpolation produced.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StripViscous {
    pub cd: f64,
    pub xtr_top: f64,
    pub xtr_bot: f64,
    pub cm: f64,
    pub envelope: Envelope,
    /// deflected strip without a family at its flap angle (see `StripResult`)
    pub flap_clamped: bool,
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
        // the 2D polar family of every foil, once; ambiguity is an error
        let uses_polars = def.viscous || def.method == Method::Llt;
        let mut foil_family: Vec<(String, Option<FamilyId>)> = Vec::new();
        let mut errors: Vec<String> = Vec::new();
        for f in plane.foil_names() {
            let id = if uses_polars {
                db.resolve(&f, &def.family).unwrap_or_else(|e| {
                    errors.push(e);
                    None
                })
            } else {
                None
            };
            foil_family.push((f, id));
        }
        if !errors.is_empty() {
            return Err(PlaneError::Input(errors.join("
")));
        }
        let fam = |name: &str| foil_family.iter().find(|(n, _)| n == name).and_then(|(_, id)| *id);
        let section_family = plane.wings.iter().map(|w| (0..w.sections.len()).map(|i| [fam(w.foil(i, false)), fam(w.foil(i, true))]).collect()).collect();
        let mut flap_sets = vec![Vec::new(); db.families().len()];
        for id in foil_family.iter().filter_map(|f| f.1) {
            flap_sets[id.0 as usize] = db.flap_families(id);
        }
        Ok(Ctx { plane, db, def, s_ref: s, b_ref: b, c_ref: c, inertia, section_family, foil_family, flap_sets })
    }

    /// The resolved polar family of a foil (setup code; not for hot loops).
    pub fn family_of(&self, foil: &str) -> Option<FamilyId> {
        self.foil_family.iter().find(|(n, _)| n == foil).and_then(|(_, id)| *id)
    }

    /// Polar families of a strip's two end foils (no string lookups).
    pub fn strip_families(&self, s: &Strip) -> (Option<FamilyId>, Option<FamilyId>) {
        let at = |i: usize| self.section_family.get(s.wing).and_then(|w| w.get(i)).and_then(|f| f[s.left as usize]);
        (at(s.surf), at(s.surf + 1))
    }

    /// Section data of one foil at a flap angle: linear between the two
    /// nearest families of its flap set, then in Re.  The flag is set when
    /// the angle lies outside the deflected families (none, or beyond them).
    fn section_at(&self, base: FamilyId, flap: f64, re: f64, cl: f64, alpha_eff: f64) -> Option<(foiley_foil::FoilPoint, Envelope, bool)> {
        let get = |f: FamilyId| if self.def.viscous_from_cl { self.db.at_re_cl(f, re, cl) } else { self.db.at_re_alpha(f, re, alpha_eff) };
        if flap.abs() <= 1e-9 {
            let (p, e) = get(base)?;
            return Some((p, e, false));
        }
        let set = self.flap_sets.get(base.0 as usize).map(|v| v.as_slice()).unwrap_or(&[]);
        if set.len() < 2 {
            let (p, e) = get(base)?;
            return Some((p, e, true));
        }
        let (first, last) = (set[0], set[set.len() - 1]);
        if flap <= first.0 {
            let (p, e) = get(first.1)?;
            return Some((p, e, flap < first.0 - 1e-6));
        }
        if flap >= last.0 {
            let (p, e) = get(last.1)?;
            return Some((p, e, flap > last.0 + 1e-6));
        }
        let j = set.partition_point(|x| x.0 < flap);
        let (lo, hi) = (set[j - 1], set[j]);
        let (pa, ea) = get(lo.1)?;
        let (pb, eb) = get(hi.1)?;
        let t = (flap - lo.0) / (hi.0 - lo.0);
        Some((foiley_foil::polar::lerp_point(&pa, &pb, t), ea.combine(eb), false))
    }

    /// Viscous drag, transition and section Cm for a strip: interpolated
    /// between the strip's two end foils (each at the strip's flap angle).
    /// Returns None if no polar exists.
    pub fn strip_viscous(&self, s: &Strip, re: f64, cl: f64, alpha_eff: f64) -> Option<StripViscous> {
        let (fa, fb) = self.strip_families(s);
        let (fa, fb) = (fa?, fb?);
        let (pa, oa, ca) = self.section_at(fa, s.flap_deg, re, cl, alpha_eff)?;
        let (pb, ob, cb) = if fb == fa { (pa, oa, ca) } else { self.section_at(fb, s.flap_deg, re, cl, alpha_eff)? };
        let tau = s.tau;
        let l = |x: f64, y: f64| x + tau * (y - x);
        Some(StripViscous { cd: l(pa.cd, pb.cd), xtr_top: l(pa.xtr_top, pb.xtr_top), xtr_bot: l(pa.xtr_bot, pb.xtr_bot), cm: l(pa.cm, pb.cm), envelope: oa.combine(ob), flap_clamped: ca || cb })
    }

    /// Reduce raw loads to an operating point.  `omega` (body rates, rad/s)
    /// makes every strip use its local air velocity V - omega x (c4 - cog)
    /// for Re, drag direction and dynamic pressure.
    #[allow(clippy::too_many_arguments)]
    pub fn reduce(&self, mesh: &Mesh, raw: &RawLoads, cog: Vec3, alpha: f64, beta: f64, v: f64, mass: f64, omega: Option<Vec3>) -> PlaneOpp {
        self.reduce_full(mesh, raw, cog, alpha, beta, v, mass, omega).0
    }

    /// [`Ctx::reduce`], also returning the total force and moment (about
    /// `cog`, N and N.m, body axes).  With `StabilityModel::FullLoad` the
    /// moment (and Cm) include the sections' viscous Cm.
    #[allow(clippy::too_many_arguments)]
    pub fn reduce_full(&self, mesh: &Mesh, raw: &RawLoads, cog: Vec3, alpha: f64, beta: f64, v: f64, mass: f64, omega: Option<Vec3>) -> (PlaneOpp, Vec3, Vec3) {
        let rho = self.def.fluid.rho;
        let section_cm = self.def.stability_model == StabilityModel::FullLoad;
        let nu = self.def.fluid.nu();
        let q = 0.5 * rho * v * v;
        let (wd, wn, ws) = (wind_dir(alpha, beta), wind_normal(alpha), wind_side(alpha, beta));
        let qs = q * self.s_ref;
        let mut warnings = Vec::new();
        let mut strips = Vec::with_capacity(mesh.strips.len());
        let mut f_visc = Vec3::ZERO;
        let mut m_visc = Vec3::ZERO;
        let mut missing_polars: Vec<String> = Vec::new();
        let mut n_env = [0usize; 5];
        let mut m_sec = Vec3::ZERO;
        let mut n_flap = 0usize;
        for (i, s) in mesh.strips.iter().enumerate() {
            // local air velocity of the strip (the freestream without rotation)
            let (v, wd, q, da) = match omega {
                Some(w) => {
                    let vloc = wd * v - w.cross(s.c4 - cog);
                    let vl = vloc.norm().max(1e-30);
                    let da = ((vloc.dot(s.normal) / vl).clamp(-1.0, 1.0).asin() - (wd.dot(s.normal)).clamp(-1.0, 1.0).asin()) / DEG;
                    (vl, vloc / vl, 0.5 * rho * vl * vl, da)
                }
                None => (v, wd, q, 0.0),
            };
            let f = raw.strip_force[i] * rho;
            let area = s.area.max(1e-30);
            let cl = f.dot(s.normal) / (q * area);
            let wi = raw.strip_wi.get(i).cloned().unwrap_or_default();
            let ai = (wi.dot(s.normal)).atan2(v) / DEG;
            let gamma = raw.strip_gamma.get(i).cloned().unwrap_or(0.0);
            let fi = wi.cross(s.te_b - s.te_a) * (gamma * rho);
            let cdi = fi.dot(wd) / (q * area);
            let re = v * s.chord / nu;
            let alpha_eff = alpha + s.twist + ai + da;
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
                match self.strip_viscous(s, re, cl, alpha_eff) {
                    Some(StripViscous { cd, xtr_top: xt, xtr_bot: xb, cm, envelope: env, flap_clamped }) => {
                        sr.cdv = cd;
                        sr.xtr_top = xt;
                        sr.xtr_bot = xb;
                        sr.cm_c4 = cm;
                        sr.flap_clamped = flap_clamped;
                        n_flap += flap_clamped as usize;
                        sr.out_of_polar = env.is_clamped();
                        sr.envelope = env as u8;
                        n_env[env as usize] += 1;
                        let dfv = wd * (q * area * cd);
                        f_visc += dfv;
                        m_visc += (s.c4 - cog).cross(dfv);
                        if section_cm {
                            // nose-up couple about n x chord
                            let chord_dir = ((s.te_a + s.te_b) * 0.5 - s.le).normalized();
                            m_sec += s.normal.cross(chord_dir) * (q * area * s.chord * cm);
                        }
                    }
                    None => {
                        sr.envelope = Envelope::Missing as u8;
                        n_env[Envelope::Missing as usize] += 1;
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
        if let Some(w) = envelope_warning(&n_env, "strip") {
            warnings.push(w);
        }
        if n_flap > 0 {
            warnings.push(format!("{n_flap} deflected strip(s) have no 2D polar family at their flap angle: section data from the nearest angle (add deflected polars)"));
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
        let m_tot = raw.moment * rho + m_visc + m_sec;
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
        let opp = PlaneOpp {
            alpha,
            beta,
            velocity: v,
            ctrl: 0.0,
            cl,
            cd,
            cdi,
            cdv,
            cy: Some(cy),
            cm,
            cr: Some(cr),
            cn: Some(cn),
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
            stability_error: None,
            root_bending,
            wake: vec![],
            balance_residual: None,
            diagnostics: None,
        };
        (opp, f_tot, m_tot)
    }
}

/// Build the mesh options for a definition and control value.
pub fn mesh_options(def: &PlanePolarDef, extra_flaps: &HashMap<(usize, usize), f64>, tilts: &HashMap<usize, f64>) -> MeshOptions {
    let mut flaps: HashMap<(usize, usize), f64> = def.flaps.iter().cloned().collect();
    for (k, v) in extra_flaps {
        *flaps.entry(*k).or_insert(0.0) += v;
    }
    MeshOptions { thick: def.method.is_panel(), include_fuse: def.include_fuse, flaps, tilts: tilts.clone(), allow_missing_foils: def.allow_missing_foils, triangulate: matches!(def.method, Method::TriPanel | Method::TriLinear), fuse_in_thin: false, flap_normals: matches!(def.method, Method::Vlm1 | Method::Vlm2) }
}

/// Build and factor the solver for a mesh.
pub fn make_system<'a>(def: &PlanePolarDef, mesh: &'a Mesh, cog: Vec3, c_ref: f64, prog: &Progress) -> Result<System<'a>, PlaneError> {
    Ok(match def.method {
        Method::Vlm1 | Method::Vlm2 => System::Vlm(VlmSystem::build_with(mesh, def.method == Method::Vlm2, cog, def.ground_height, prog, def.diagnostics)?),
        Method::Panel | Method::TriPanel => System::Panel(PanelSystem::build_with(mesh, cog, def.wake_length * c_ref, def.wake_panels, def.wake_ratio, def.ground_height, prog, def.diagnostics)?),
        Method::TriLinear => System::TriLinear(LinearPanelSystem::build_with(mesh, cog, def.wake_length * c_ref, def.wake_panels, def.wake_ratio, def.ground_height, prog, def.diagnostics)?),
        Method::Llt => return Err(PlaneError::Input("LLT has no panel system".into())),
    })
}

/// "viscous drag clamped to the 2D polar envelope on N strip(s): a in Re,
/// b in CL, c in both", or None when nothing was clamped.
pub fn envelope_warning(counts: &[usize; 5], what: &str) -> Option<String> {
    let n: usize = counts[1..4].iter().sum();
    if n == 0 {
        return None;
    }
    let parts: Vec<String> = [(1, "in Re"), (2, "in CL"), (3, "in both")].iter().filter(|(k, _)| counts[*k] > 0).map(|(k, l)| format!("{} {l}", counts[*k])).collect();
    Some(format!("viscous drag clamped to the 2D polar envelope on {n} {what}(s): {}", parts.join(", ")))
}

/// With `strict_envelope`, why a point fails (any clamped or polar-less strip).
pub fn envelope_failure(def: &PlanePolarDef, o: &PlaneOpp) -> Option<String> {
    if !def.strict_envelope || !def.viscous {
        return None;
    }
    let mut counts = [0usize; 5];
    for s in &o.strips {
        counts[(s.envelope as usize).min(4)] += 1;
    }
    let bad: usize = counts[1..].iter().sum();
    let flap = o.strips.iter().filter(|s| s.flap_clamped).count();
    if bad == 0 && flap == 0 {
        return None;
    }
    let mut parts: Vec<String> = Envelope::KINDS[1..].iter().filter(|e| counts[**e as usize] > 0).map(|e| format!("{} {}", counts[*e as usize], e.label())).collect();
    if flap > 0 {
        parts.push(format!("{flap} without a polar family at their flap angle"));
    }
    Some(format!("strict envelope: {} strip(s) outside the 2D polars ({})", o.strips.iter().filter(|s| s.envelope != 0 || s.flap_clamped).count(), parts.join(", ")))
}

/// Run the stability analysis for a point, recording a failure on the point
/// (`stability_error` plus a warning) instead of dropping it.  Only a
/// cancellation propagates.
fn attach_stability(opp: &mut PlaneOpp, res: Result<StabilityResult, PlaneError>) -> Result<(), PlaneError> {
    match res {
        Ok(s) => opp.stability = Some(s),
        Err(PlaneError::Cancelled) => return Err(PlaneError::Cancelled),
        Err(e) => {
            let msg = format!("stability analysis failed: {e}");
            opp.warnings.push(msg.clone());
            opp.stability_error = Some(msg);
        }
    }
    Ok(())
}

/// Stability needs a mass; say so on the point when it is missing.
fn stability_needs_mass(opp: &mut PlaneOpp) {
    let msg = "stability derivatives not computed: the plane has no mass (set masses or an inertia override)".to_string();
    opp.warnings.push(msg.clone());
    opp.stability_error = Some(msg);
}

/// Why an operating point is not reportable, if any coefficient is not finite.
pub fn nonfinite_reason(o: &PlaneOpp) -> Option<String> {
    let vals = [("alpha", o.alpha), ("V", o.velocity), ("CL", o.cl), ("CD", o.cd), ("CY", o.cy.unwrap_or(0.0)), ("Cm", o.cm), ("Cl", o.cr.unwrap_or(0.0)), ("Cn", o.cn.unwrap_or(0.0)), ("lift", o.lift), ("drag", o.drag)];
    let bad: Vec<&str> = vals.iter().filter(|(_, v)| !v.is_finite()).map(|(n, _)| *n).collect();
    if bad.is_empty() {
        None
    } else {
        Some(format!("non-finite {}", bad.join(", ")))
    }
}

/// Tolerance on |resultant - weight| / weight for the T3 speed balance.
pub const GLIDE_TOL: f64 = 1e-4;
/// Passes of the T3 speed balance before the point is declared unconverged.
pub const GLIDE_MAX_ITER: usize = 50;

/// Speed at which lift equals weight for the loads at unit speed.
fn balance_speed(ctx: &Ctx, raw1: &RawLoads, alpha: f64, mass: f64) -> Option<f64> {
    let l1 = (raw1.force + raw1.induced_force).dot(wind_normal(alpha)) * ctx.def.fluid.rho;
    if l1 <= 1e-12 || mass <= 0.0 {
        return None;
    }
    Some((mass * G / l1).sqrt())
}

/// With `strict_wake`, why a point fails (its relaxed wake did not converge).
pub fn wake_failure(def: &PlanePolarDef, o: &PlaneOpp) -> Option<String> {
    if !def.strict_wake {
        return None;
    }
    o.diagnostics.as_ref().and_then(|d| d.wake.as_ref()).and_then(|w| w.warning()).map(|w| format!("strict wake: {w}"))
}

/// Why a computed point is not reportable (non-finite, strict envelope, strict wake).
fn point_failure(def: &PlanePolarDef, o: &PlaneOpp) -> Option<String> {
    nonfinite_reason(o).or_else(|| envelope_failure(def, o)).or_else(|| wake_failure(def, o))
}

/// Loads of one state through the evaluator, reduced to an operating point,
/// with its diagnostics, relaxed wake lines and wake report.
#[allow(clippy::too_many_arguments)]
fn evaluate(ctx: &Ctx, ev: &Evaluator, cog: Vec3, alpha: f64, beta: f64, v: f64, mass: f64, prog: &Progress) -> Result<PlaneOpp, PlaneError> {
    let wd = wind_dir(alpha, beta);
    let (raw, wrep) = ev.loads(wd, v, Vec3::ZERO, prog)?;
    let mut opp = ctx.reduce(ev.mesh(), &raw, cog, alpha, beta, v, mass, None);
    if ctx.def.diagnostics != Diagnostics::Off {
        opp.diagnostics = Some(ev.sys.point_diagnostics(&raw, wd * v, Vec3::ZERO, ctx.s_ref, wrep.relaxed));
    }
    if wrep.relaxed {
        opp.wake = ev.wake_lines(wd, v, Vec3::ZERO);
        if let Some(w) = wrep.warning() {
            opp.warnings.push(w);
        }
        opp.diagnostics.get_or_insert_with(Default::default).wake = Some(wrep);
    }
    Ok(opp)
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
    let mut polar = PlanePolar { name: name.clone(), plane: plane.name.clone(), def: Some(def.clone()), ref_dims: (ctx.s_ref, ctx.b_ref, ctx.c_ref), manifest: Some(crate::provenance::run_manifest(&ctx)), ..Default::default() };
    prog.stage(format!("{}: {}", plane.name, name));
    prog.info(format!("{} - {} - {}", plane.name, def.kind.label(), def.method.label()));
    if mass > 0.0 {
        prog.info(format!("mass {:.3} kg, CoG ({:.4}, {:.4}, {:.4}) m", mass, cog.x, cog.y, cog.z));
    }
    if def.method == Method::Llt {
        let pts = crate::llt::run(&ctx, prog, on_point, &mut polar.failed)?;
        polar.points = pts;
        polar.unknowns = def.llt.stations;
        drop_strips(&mut polar);
        return Ok(polar);
    }
    match def.kind {
        PolarKind::Control => run_control(&ctx, &mut polar, prog, on_point)?,
        PolarKind::Stability => run_stability(&ctx, &mut polar, prog, on_point)?,
        _ => {
            prog.stage("Meshing");
            let mesh = build(plane, db, &mesh_options(def, &HashMap::new(), &HashMap::new()))?;
            polar.unknowns = def.method.unknowns(mesh.panels.len());
            prog.info(format!("mesh: {} panels, {} strips", mesh.panels.len(), mesh.strips.len()));
            let ev = Evaluator::new(make_system(def, &mesh, cog, ctx.c_ref, &prog.sub(0.0, 0.5))?, def, ctx.c_ref);
            let ops: Vec<(f64, f64, f64)> = match def.kind {
                PolarKind::FixedSpeed | PolarKind::FixedLift | PolarKind::Glide => def.range.values().into_iter().map(|a| (a, def.beta, def.velocity)).collect(),
                PolarKind::FixedAlpha => def.range.values().into_iter().map(|v| (def.alpha, def.beta, v)).collect(),
                PolarKind::Sideslip => def.range.values().into_iter().map(|b| (def.alpha, b, def.velocity)).collect(),
                PolarKind::List => def.list.clone(),
                _ => unreachable!(),
            };
            let rest = prog.sub(0.5, 1.0);
            let n = ops.len();
            let results: Vec<Result<Result<PlaneOpp, String>, PlaneError>> = {
                use rayon::prelude::*;
                ops.par_iter()
                    .enumerate()
                    .map(|(k, &(a, b, v))| {
                        rest.check()?;
                        let part = rest.part(k, n);
                        let r = solve_point(&ctx, &ev, cog, mass, a, b, v, &part);
                        if !matches!(r, Err(PlaneError::Cancelled)) {
                            part.set(1.0);
                        }
                        match &r {
                            Ok(Ok(p)) => {
                                rest.count_ok();
                                rest.debug(format!("alpha={:.2} beta={:.2} V={:.2}: CL={:.4} CD={:.5} Cm={:.4}", p.alpha, p.beta, p.velocity, p.cl, p.cd, p.cm));
                                on_point(p);
                            }
                            Ok(Err(_)) => rest.count_failed(),
                            Err(_) => {}
                        }
                        r
                    })
                    .collect()
            };
            for (k, r) in results.into_iter().enumerate() {
                match r? {
                    Ok(p) => polar.points.push(p),
                    Err(why) => polar.failed.push(format!("alpha={:.2} beta={:.2} V={:.2}: {why}", ops[k].0, ops[k].1, ops[k].2)),
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
    drop_strips(&mut polar);
    prog.set(1.0);
    Ok(polar)
}

/// Remove the span distributions when the definition does not keep them.
fn drop_strips(polar: &mut PlanePolar) {
    if polar.def.as_ref().is_some_and(|d| !d.keep_strips) {
        for p in &mut polar.points {
            p.strips = Vec::new();
        }
    }
}

#[allow(clippy::too_many_arguments)]
/// Solve one operating point.  `Ok(Err(reason))` is a point that failed
/// (it goes to the polar's `failed` list); `Err` aborts the whole polar.
fn solve_point(ctx: &Ctx, ev: &Evaluator, cog: Vec3, mass: f64, alpha: f64, beta: f64, v_spec: f64, prog: &Progress) -> Result<Result<PlaneOpp, String>, PlaneError> {
    let def = ctx.def;
    let v = match def.kind {
        PolarKind::FixedLift | PolarKind::Glide => {
            // lift = weight from the loads at unit speed (a relaxed wake does
            // not depend on the speed, so the balance holds at the result)
            let (raw1, _) = ev.loads(wind_dir(alpha, beta), 1.0, Vec3::ZERO, prog)?;
            match balance_speed(ctx, &raw1, alpha, mass) {
                Some(v) => v,
                None => return Ok(Err("no positive lift for the speed balance".into())),
            }
        }
        _ => v_spec,
    };
    let mut opp = evaluate(ctx, ev, cog, alpha, beta, v, mass, prog)?;
    if def.kind == PolarKind::Glide {
        // speed from the total aerodynamic force (viscous drag depends on Re)
        let weight = mass * G;
        let residual = |o: &PlaneOpp| ((o.cl * o.cl + o.cd * o.cd).sqrt() * 0.5 * def.fluid.rho * o.velocity * o.velocity * ctx.s_ref - weight).abs() / weight;
        let mut res = residual(&opp);
        let mut pass = 0;
        while res >= GLIDE_TOL {
            if pass == GLIDE_MAX_ITER {
                return Ok(Err(format!("glide speed balance did not converge in {GLIDE_MAX_ITER} passes (|resultant - weight|/weight = {res:.2e})")));
            }
            pass += 1;
            let ca = (opp.cl * opp.cl + opp.cd * opp.cd).sqrt();
            if ca <= 1e-9 {
                return Ok(Err("no aerodynamic force for the glide speed balance".into()));
            }
            let vn = (2.0 * weight / (def.fluid.rho * ctx.s_ref * ca)).sqrt();
            opp = evaluate(ctx, ev, cog, alpha, beta, vn, mass, prog)?;
            res = residual(&opp);
        }
        opp.balance_residual = Some(res);
    }
    if let Some(why) = point_failure(def, &opp) {
        return Ok(Err(why));
    }
    if def.derivatives {
        if mass > 0.0 {
            let res = stability::analyse(ctx, ev, cog, alpha, beta, opp.velocity, prog);
            attach_stability(&mut opp, res)?;
        } else {
            stability_needs_mass(&mut opp);
        }
    }
    Ok(Ok(opp))
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
        let (cx, cz) = (control_value(def, t, ControlTarget::CogX), control_value(def, t, ControlTarget::CogZ));
        if let Some(x) = cx {
            inertia.cog.x = x;
        }
        if let Some(z) = cz {
            inertia.cog.z = z;
        }
        // moments about the CoG as in T1 (the origin for a massless plane)
        let cog = if inertia.mass > 0.0 || cx.is_some() || cz.is_some() { inertia.cog } else { Vec3::ZERO };
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
        polar.unknowns = def.method.unknowns(mesh.panels.len());
        let ev = Evaluator::new(make_system(def, &mesh, cog, ctx.c_ref, &sub.sub(0.0, 0.9))?, def, ctx.c_ref);
        let mut opp = evaluate(ctx, &ev, cog, a, b, v, inertia.mass, &sub)?;
        opp.ctrl = t;
        if let Some(why) = point_failure(def, &opp) {
            sub.set(1.0);
            prog.count_failed();
            polar.failed.push(format!("ctrl {t:.3}: {why}"));
            continue;
        }
        if def.derivatives {
            if inertia.mass > 0.0 {
                let c2 = Ctx { inertia, ..Ctx::new(ctx.plane, ctx.db, def)? };
                attach_stability(&mut opp, stability::analyse(&c2, &ev, cog, a, b, v, &sub))?;
            } else {
                stability_needs_mass(&mut opp);
            }
        }
        sub.set(1.0);
        prog.count_ok();
        on_point(&opp);
        polar.points.push(opp);
    }
    Ok(())
}

/// Flap and tilt deflections of the T7 gains at control parameter `t`.
fn gain_deflections(def: &PlanePolarDef, t: f64) -> (HashMap<(usize, usize), f64>, HashMap<usize, f64>) {
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
    (flaps, tilts)
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
        let (flaps, tilts) = gain_deflections(def, t);
        let mesh = build(ctx.plane, ctx.db, &mesh_options(def, &flaps, &tilts))?;
        polar.unknowns = def.method.unknowns(mesh.panels.len());
        let ev = Evaluator::new(make_system(def, &mesh, cog, ctx.c_ref, &sub.sub(0.0, 0.7))?, def, ctx.c_ref);
        match stability::trim(ctx, &ev, def.beta, &sub)? {
            Some((alpha, u0)) => {
                let mut opp = evaluate(ctx, &ev, cog, alpha, def.beta, u0, mass, &sub)?;
                opp.ctrl = t;
                if let Some(why) = point_failure(def, &opp) {
                    prog.count_failed();
                    polar.failed.push(format!("ctrl {t:.3}: {why}"));
                    sub.set(1.0);
                    continue;
                }
                sub.stage(format!("Stability: control {t:.3} - derivatives and modes"));
                match stability::analyse(ctx, &ev, cog, alpha, def.beta, u0, &sub) {
                    Ok(mut s) => {
                        // control derivatives: rebuild with the control parameter nudged
                        let gmax = def.gains.iter().map(|g| g.1.abs()).fold(0.0, f64::max);
                        if gmax > 0.0 {
                            let dt = 0.5 / gmax;
                            let (f2, t2) = gain_deflections(def, t + dt);
                            let r = (|| -> Result<(), PlaneError> {
                                let m2 = build(ctx.plane, ctx.db, &mesh_options(def, &f2, &t2))?;
                                let ev2 = Evaluator::new(make_system(def, &m2, cog, ctx.c_ref, &Progress::none())?, def, ctx.c_ref);
                                stability::add_control_derivatives(ctx, &ev, &ev2, &mut s, dt, &sub)
                            })();
                            match r {
                                Ok(()) => {}
                                Err(PlaneError::Cancelled) => return Err(PlaneError::Cancelled),
                                Err(e) => opp.warnings.push(format!("control derivatives not computed: {e}")),
                            }
                        }
                        opp.stability = Some(s)
                    }
                    Err(e) => attach_stability(&mut opp, Err(e))?,
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
    Ok(())
}
