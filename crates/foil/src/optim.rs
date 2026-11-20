//! Airfoil shape optimisation (xflr5 Optim2d single-point and Xoptfoil2-style
//! multi-point operating-point table with structural spar thickness constraint).
//!
//! The design variables change the shape of the base foil.  As in 1.1 the
//! foil is first moved and scaled to unit chord with its leading edge at the
//! origin (the flap hinge and the analysis work in x/c), but it is not
//! rotated, so the design angle of attack is the one the base foil is
//! analysed at:
//! either the amplitudes of Hicks-Henne bump functions added along the local
//! surface normal, or changes of the weights of a CST (class-shape
//! transformation) representation added to each surface, which move the
//! thickness, camber and leading-edge radius globally without waviness.
//! A trailing-edge flap angle can be added as a further variable.
//!
//! Two evaluation modes are supported:
//! * **Single-point mode** (`operating_points` empty): each candidate is
//!   analysed at one angle of attack and Reynolds number, and one or more
//!   objectives (`Cl`, `Cd`, `Cl/Cd`, `Cpmin`, `Cm`) are minimised, maximised
//!   or driven to a target value.
//! * **Multi-point operating-point table** (`operating_points` non-empty):
//!   each candidate is analysed across every enabled `(Re, target CL)` row in
//!   XFoil CL-specified mode and scored against each row's [`DragGoal`]
//!   (`Minimise`, `AtMostOriginal`, `AtMost`, `Target`).  Optional
//!   [`SparConstraint`] enforces a minimum `t/c` at a chordwise spar station
//!   before calling XFoil.
//!
//! Two search algorithms are available, following xflr5:
//! * [`Algorithm::Mopso`]: multi-objective particle swarm with a Pareto
//!   archive of non-dominated designs (plain PSO with a single objective);
//! * [`Algorithm::Genetic`]: roulette selection, blend crossover and
//!   Gaussian mutation, with elitism, on a scalarised fitness.
//!
//! Candidates of a generation are analysed in parallel (rayon); the random
//! numbers are drawn only in the serial parts, so a run is reproducible for
//! a given `seed`.

use crate::analysis::{analyze_point, analyze_points, AnalysisSettings};
use crate::foil::Foil;
use crate::FoilError;
use foiley_core::Progress;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Score assigned to unconverged / invalid candidates.
const LARGE: f64 = 1.0e10;

/// Aerodynamic quantity used as an objective.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Objective {
    Cl,
    Cd,
    ClCd,
    Cpmin,
    Cm,
}

impl Objective {
    pub fn label(&self) -> &'static str {
        match self {
            Objective::Cl => "Cl",
            Objective::Cd => "Cd",
            Objective::ClCd => "Cl/Cd",
            Objective::Cpmin => "Cpmin",
            Objective::Cm => "Cm",
        }
    }
    fn value(&self, cl: f64, cd: f64, cm: f64, cpmin: f64) -> f64 {
        match self {
            Objective::Cl => cl,
            Objective::Cd => cd,
            Objective::ClCd => cl / cd,
            Objective::Cpmin => cpmin,
            Objective::Cm => cm,
        }
    }
}

/// What to do with an objective (xflr5 type -1 / 0 / +1).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum Goal {
    Minimise,
    /// Reach `value`; anything within `tolerance` counts as satisfied.
    Target { value: f64, tolerance: f64 },
    Maximise,
}

impl Goal {
    /// Score to minimise for an objective value.
    fn score(&self, v: f64) -> f64 {
        match *self {
            Goal::Minimise => v,
            Goal::Maximise => -v,
            Goal::Target { value, tolerance } => {
                let d = (v - value).abs();
                if d <= tolerance.abs() {
                    0.0
                } else {
                    d
                }
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct ObjectiveSpec {
    pub objective: Objective,
    pub enabled: bool,
    pub goal: Goal,
}

impl Default for ObjectiveSpec {
    fn default() -> Self {
        ObjectiveSpec { objective: Objective::ClCd, enabled: true, goal: Goal::Maximise }
    }
}

/// Drag goal for one row of the multi-point operating-point table
/// (Xoptfoil2 style: `min-drag`, keep at or below original, cap, or target).
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DragGoal {
    /// Minimise profile drag coefficient CD at this operating point.
    #[default]
    Minimise,
    /// Keep CD at or below the original (base) airfoil's CD at this point.
    AtMostOriginal,
    /// Keep CD at or below a fixed upper bound.
    AtMost { cd: f64 },
    /// Drive CD towards a target value.
    Target { cd: f64 },
}

impl DragGoal {
    pub fn label(&self) -> &'static str {
        match self {
            DragGoal::Minimise => "minimise",
            DragGoal::AtMostOriginal => "<= original",
            DragGoal::AtMost { .. } => "<= max CD",
            DragGoal::Target { .. } => "target CD",
        }
    }

    pub fn is_cap(&self) -> bool {
        matches!(self, DragGoal::AtMostOriginal | DragGoal::AtMost { .. })
    }
}

/// One row of the multi-point operating-point table.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct OperatingPoint {
    /// Optional human-readable label (e.g. "Slow flight", "Cruise", "Fast glide").
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub label: String,
    /// Include this operating point in the optimisation.
    pub enabled: bool,
    /// Chord Reynolds number for this operating point.
    pub re: f64,
    /// Target lift coefficient CL (solved in XFoil CL-specified mode).
    pub cl: f64,
    /// Drag objective or cap at this operating point.
    pub goal: DragGoal,
    /// Relative weighting (must be positive, default 1.0).
    pub weight: f64,
}

impl Default for OperatingPoint {
    fn default() -> Self {
        OperatingPoint {
            label: String::new(),
            enabled: true,
            re: 2.0e5,
            cl: 0.5,
            goal: DragGoal::Minimise,
            weight: 1.0,
        }
    }
}

impl OperatingPoint {
    pub fn display_label(&self, index: usize) -> String {
        if self.label.trim().is_empty() {
            format!("OP {} (CL {:.2})", index + 1, self.cl)
        } else {
            self.label.clone()
        }
    }
}

/// Minimum airfoil thickness required at a chordwise spar station.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct SparConstraint {
    /// Chordwise spar location x/c in (0, 1), typically 0.25..0.35.
    pub x: f64,
    /// Minimum thickness-to-chord ratio t/c at `x` (e.g. 0.115 for 11.5%).
    pub min_thickness: f64,
}

impl Default for SparConstraint {
    fn default() -> Self {
        SparConstraint { x: 0.30, min_thickness: 0.11 }
    }
}

/// Evaluation of one airfoil at one operating-point row.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OpPointEval {
    pub label: String,
    pub re: f64,
    pub target_cl: f64,
    pub goal: DragGoal,
    pub weight: f64,
    pub converged: bool,
    pub alpha: f64,
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    pub cpmin: f64,
    pub xtr_top: f64,
    pub xtr_bot: f64,
}

impl OpPointEval {
    /// Whether this evaluated row satisfies its drag goal compared to `orig`.
    /// `Minimise` has no hard cap, so it always satisfies feasibility once converged
    /// (use `cd < orig.cd` separately to colour improvements vs trade-offs).
    pub fn meets_goal(&self, orig: Option<&OpPointEval>) -> Option<bool> {
        if !self.converged {
            return Some(false);
        }
        match self.goal {
            DragGoal::AtMostOriginal => orig.filter(|o| o.converged).map(|o| self.cd <= o.cd + 1e-7),
            DragGoal::AtMost { cd } => Some(self.cd <= cd + 1e-7),
            DragGoal::Target { cd } => Some((self.cd - cd).abs() <= 1e-4),
            DragGoal::Minimise => Some(true),
        }
    }
}

/// Evaluation of the original (base) airfoil at the same operating points and
/// spar station, for original-versus-candidate comparison.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct BaselineEval {
    pub foil_name: String,
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    pub cpmin: f64,
    pub converged: bool,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub op_points: Vec<OpPointEval>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub spar_thickness: Option<f64>,
    pub max_thickness: f64,
}

/// How the genes change the shape.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Parameterization {
    /// `bumps_per_side` Hicks-Henne bumps per surface.
    #[default]
    HicksHenne,
    /// `cst_order + 1` CST weight changes per surface.
    Cst,
}

impl Parameterization {
    pub fn label(&self) -> &'static str {
        match self {
            Parameterization::HicksHenne => "Hicks-Henne",
            Parameterization::Cst => "CST",
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum Algorithm {
    #[default]
    Mopso,
    Genetic,
}

impl Algorithm {
    pub fn label(&self) -> &'static str {
        match self {
            Algorithm::Mopso => "MOPSO",
            Algorithm::Genetic => "GA",
        }
    }
}

/// All optimiser parameters.  Defaults follow xflr5's Optim2d.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct OptimSettings {
    pub algorithm: Algorithm,
    /// Hicks-Henne bumps or CST weights.
    pub parameterization: Parameterization,
    /// CST order: `cst_order + 1` weights per side.
    pub cst_order: usize,
    /// Largest change of a CST weight (genes range +-this).
    pub cst_max_delta: f64,
    /// Hicks-Henne functions per side (design variables = 2 x this).
    pub bumps_per_side: usize,
    /// HH t1: bump j (0-based) peaks at x = t1 (j+1) / (2n+1).
    pub t1: f64,
    /// HH t2: bump width exponent (larger = narrower).
    pub t2: f64,
    /// Maximum bump amplitude, chord fraction (genes range +-this).
    pub max_amplitude: f64,
    /// TE flap angle range (min, max) in degrees, positive down; `None` = no flap variable.
    pub flap: Option<(f64, f64)>,
    /// Flap hinge (x/c, y relative 0 = lower .. 1 = upper).
    pub hinge: (f64, f64),
    /// Angle of attack of the design point, degrees (single-point mode).
    pub alpha: f64,
    pub analysis: AnalysisSettings,
    pub objectives: Vec<ObjectiveSpec>,
    /// Multi-point operating-point table (Xoptfoil2 style).  When non-empty,
    /// each candidate is evaluated at every enabled row's `(Re, target CL)`
    /// and scored by its [`DragGoal`] instead of the single `alpha` / `objectives`.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub operating_points: Vec<OperatingPoint>,
    /// Minimum thickness-to-chord ratio at a chordwise spar station.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub spar: Option<SparConstraint>,
    /// Swarm / population size.
    pub population: usize,
    /// Number of iterations / generations after the initial one.
    pub iterations: usize,
    // PSO
    pub inertia: f64,
    pub cognitive: f64,
    pub social: f64,
    pub regen_probability: f64,
    pub archive_size: usize,
    // GA
    pub crossover: f64,
    pub mutation_probability: f64,
    /// Mutation standard deviation as a fraction of the gene range half-width.
    pub mutation_sigma: f64,
    pub seed: u64,
}

impl OptimSettings {
    /// xflr5's objective table (Cl max, Cd min enabled; Cl/Cd, Cpmin, Cm off).
    pub fn default_objectives() -> Vec<ObjectiveSpec> {
        vec![
            ObjectiveSpec { objective: Objective::Cl, enabled: true, goal: Goal::Maximise },
            ObjectiveSpec { objective: Objective::Cd, enabled: true, goal: Goal::Minimise },
            ObjectiveSpec { objective: Objective::ClCd, enabled: false, goal: Goal::Maximise },
            ObjectiveSpec { objective: Objective::Cpmin, enabled: false, goal: Goal::Maximise },
            ObjectiveSpec { objective: Objective::Cm, enabled: false, goal: Goal::Target { value: 0.0, tolerance: 1.0e-3 } },
        ]
    }

    /// Default 3-point glider table (slow flight, cruise, fast glide).
    pub fn default_operating_points() -> Vec<OperatingPoint> {
        vec![
            OperatingPoint { label: "Slow flight".into(), enabled: true, re: 1.2e5, cl: 0.90, goal: DragGoal::AtMostOriginal, weight: 1.0 },
            OperatingPoint { label: "Cruise".into(), enabled: true, re: 1.8e5, cl: 0.50, goal: DragGoal::Minimise, weight: 1.0 },
            OperatingPoint { label: "Fast glide".into(), enabled: true, re: 2.8e5, cl: 0.20, goal: DragGoal::Minimise, weight: 1.0 },
        ]
    }

    /// True when the multi-point operating-point table is active.
    pub fn is_multipoint(&self) -> bool {
        !self.operating_points.is_empty()
    }

    /// Enabled rows of the operating-point table.
    pub fn enabled_operating_points(&self) -> Vec<OperatingPoint> {
        self.operating_points.iter().filter(|o| o.enabled).cloned().collect()
    }

    /// Shape variables per surface.
    pub fn shape_variables_per_side(&self) -> usize {
        match self.parameterization {
            Parameterization::HicksHenne => self.bumps_per_side,
            Parameterization::Cst => self.cst_order.saturating_add(1),
        }
    }

    /// Check the shape-variable, operating-point and spar settings (the base
    /// analysis settings are checked separately): finite ranges, and counts
    /// small enough to allocate.
    pub fn validate(&self) -> Result<(), String> {
        const MAX_VARIABLES: usize = 100_000;
        if self.shape_variables_per_side() > MAX_VARIABLES {
            return Err(format!("{} shape variables per side; at most {MAX_VARIABLES}", self.shape_variables_per_side()));
        }
        let evaluations = self.population.max(2).checked_mul(self.iterations.saturating_add(1));
        if evaluations.is_none_or(|e| e > 1_000_000_000) {
            return Err(format!("population {} x {} iterations is too many evaluations", self.population, self.iterations));
        }
        match self.parameterization {
            Parameterization::HicksHenne => {
                if !self.max_amplitude.is_finite() || !self.t1.is_finite() || !self.t2.is_finite() {
                    return Err("the Hicks-Henne amplitude and shape parameters must be finite".into());
                }
            }
            Parameterization::Cst => {
                if self.cst_order > 30 {
                    return Err(format!("CST order {}; at most 30", self.cst_order));
                }
                if !self.cst_max_delta.is_finite() {
                    return Err("the largest CST weight change must be finite".into());
                }
            }
        }
        if let Some((a, b)) = self.flap {
            if !a.is_finite() || !b.is_finite() {
                return Err("the flap angle range must be finite".into());
            }
        }
        if let Some(sp) = self.spar {
            if !sp.x.is_finite() || !(0.01..=0.99).contains(&sp.x) {
                return Err(format!("spar station x/c = {} must be between 0.01 and 0.99", sp.x));
            }
            if !sp.min_thickness.is_finite() || !(0.001..=0.5).contains(&sp.min_thickness) {
                return Err(format!("spar minimum thickness t/c = {} must be between 0.1% and 50%", sp.min_thickness));
            }
        }
        for (i, op) in self.operating_points.iter().enumerate() {
            if !op.enabled {
                continue;
            }
            let name = op.display_label(i);
            if !op.re.is_finite() || op.re <= 1.0e3 || op.re >= 1.0e10 {
                return Err(format!("{name}: Reynolds number {} is out of the supported range 1e3..1e10", op.re));
            }
            if !op.cl.is_finite() || !(-3.0..=3.5).contains(&op.cl) {
                return Err(format!("{name}: target CL {} must be between -3 and 3.5", op.cl));
            }
            if !op.weight.is_finite() || op.weight <= 0.0 {
                return Err(format!("{name}: weight {} must be positive", op.weight));
            }
            match op.goal {
                DragGoal::AtMost { cd } | DragGoal::Target { cd } if !cd.is_finite() || cd <= 0.0 => {
                    return Err(format!("{name}: drag target {cd} must be positive"));
                }
                _ => {}
            }
        }
        Ok(())
    }

    /// Number of design variables.
    pub fn dimension(&self) -> usize {
        self.shape_variables_per_side().saturating_mul(2).saturating_add(usize::from(self.flap.is_some()))
    }

    /// (min, max) of each design variable; empty for settings that fail
    /// [`OptimSettings::validate`].
    pub fn bounds(&self) -> Vec<(f64, f64)> {
        if self.validate().is_err() {
            return Vec::new();
        }
        let a = match self.parameterization {
            Parameterization::HicksHenne => self.max_amplitude.abs(),
            Parameterization::Cst => self.cst_max_delta.abs(),
        };
        let mut b = vec![(-a, a); 2 * self.shape_variables_per_side()];
        if let Some((lo, hi)) = self.flap {
            b.push((lo.min(hi), lo.max(hi)));
        }
        b
    }

    /// Genes that reproduce the base foil (zero bumps, flap nearest 0).
    pub fn neutral_genes(&self) -> Vec<f64> {
        self.bounds().iter().map(|&(lo, hi)| 0.0f64.clamp(lo, hi)).collect()
    }

    fn enabled(&self) -> Vec<ObjectiveSpec> {
        self.objectives.iter().filter(|o| o.enabled).cloned().collect()
    }
}

impl Default for OptimSettings {
    fn default() -> Self {
        OptimSettings {
            algorithm: Algorithm::Mopso,
            parameterization: Parameterization::HicksHenne,
            cst_order: 5,
            cst_max_delta: 0.05,
            bumps_per_side: 3,
            t1: 1.7,
            t2: 1.0,
            max_amplitude: 0.01,
            flap: None,
            hinge: (0.7, 0.5),
            alpha: 0.0,
            analysis: AnalysisSettings { re: 1.0e6, mach: 0.0, ncrit: 9.0, xtr_top: 1.0, xtr_bot: 1.0, ..Default::default() },
            objectives: Self::default_objectives(),
            operating_points: Vec::new(),
            spar: None,
            population: 11,
            iterations: 30,
            inertia: 0.3,
            cognitive: 0.7,
            social: 0.7,
            regen_probability: 0.07,
            archive_size: 11,
            crossover: 0.5,
            mutation_probability: 0.15,
            mutation_sigma: 0.2,
            seed: 1,
        }
    }
}

/// One evaluated design.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Candidate {
    pub genes: Vec<f64>,
    pub foil: Foil,
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    pub cpmin: f64,
    /// Value of each enabled objective (same order as the enabled specs, or
    /// the objective rows in multi-point mode).
    pub objectives: Vec<f64>,
    /// Per enabled objective, the score being minimised (-value for
    /// maximise, distance outside the tolerance band for targets).
    pub scores: Vec<f64>,
    /// Converged analysis on a valid geometry satisfying geometric constraints.
    pub feasible: bool,
    /// Per-row evaluation in multi-point mode (empty in single-point mode).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub op_points: Vec<OpPointEval>,
    /// Thickness-to-chord ratio at `spar.x` when a spar constraint is active.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub spar_thickness: Option<f64>,
    /// Total normalised drag-cap constraint violation (0.0 when all caps met).
    #[serde(default)]
    pub violation: f64,
}

impl Candidate {
    /// True when the candidate converged, met geometric constraints and
    /// satisfied every drag cap (`<= original`, `<= max CD`).
    pub fn constraints_met(&self) -> bool {
        self.feasible && self.violation == 0.0 && self.violation.is_finite()
    }

    /// Constrained dominance (Deb): feasible with zero violation beats any
    /// violation; smaller violation beats larger; equal zero violation uses
    /// standard Pareto dominance on `scores`.
    fn dominates(&self, b: &Candidate) -> bool {
        if !self.feasible {
            return false;
        }
        if !b.feasible {
            return true;
        }
        if self.violation > 0.0 || b.violation > 0.0 {
            return self.violation < b.violation;
        }
        let mut better = false;
        for (x, y) in self.scores.iter().zip(&b.scores) {
            if x > y {
                return false;
            }
            if x < y {
                better = true;
            }
        }
        better
    }
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct OptimResult {
    /// Best compromise: lowest scalarised score (the only objective when
    /// a single one is enabled).
    pub best: Option<Candidate>,
    /// Non-dominated archive.
    pub pareto: Vec<Candidate>,
    /// (iteration, metric of `best`), iteration 0 = initial population.
    /// Single-point: first enabled objective value; multi-point: weighted CD
    /// ratio relative to original (1.0 = original, < 1.0 = lower drag).
    pub history: Vec<(usize, f64)>,
    pub evaluations: usize,
    /// Evaluation of the original (base) foil for original-vs-candidate comparison.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub baseline: Option<BaselineEval>,
}

// ---------------------------------------------------------------- RNG

/// SplitMix64: small, fast, deterministic.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        let mut r = Rng(seed ^ 0x9E37_79B9_7F4A_7C15);
        r.next_u64();
        r
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in [0, 1).
    fn f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.f64()
    }
    fn below(&mut self, n: usize) -> usize {
        if n == 0 {
            0
        } else {
            (self.next_u64() % n as u64) as usize
        }
    }
    /// Standard normal (Box-Muller).
    fn normal(&mut self) -> f64 {
        let u1 = self.f64().max(1e-300);
        let u2 = self.f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

// ---------------------------------------------------------------- geometry

/// Hicks-Henne bump (xflr5 `HH`): sin(pi x^(ln 0.5 / ln t1))^t2, peak at x = t1.
fn hicks_henne(x: f64, t1: f64, t2: f64) -> f64 {
    if x <= 0.0 || x >= 1.0 {
        return 0.0;
    }
    let t1 = t1.clamp(1e-3, 1.0 - 1e-3);
    (std::f64::consts::PI * x.powf(0.5f64.ln() / t1.ln())).sin().abs().powf(t2)
}

/// Linear interpolation of y(x) on a polyline whose x increases.
fn interp(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let n = xs.len();
    for i in 1..n {
        let (x0, x1) = (xs[i - 1], xs[i]);
        if (x0 <= x && x <= x1) || (x1 <= x && x <= x0) {
            let d = x1 - x0;
            return if d.abs() < 1e-15 { ys[i] } else { ys[i - 1] + (ys[i] - ys[i - 1]) * (x - x0) / d };
        }
    }
    if (x - xs[0]).abs() < (x - xs[n - 1]).abs() {
        ys[0]
    } else {
        ys[n - 1]
    }
}

/// True if the lower surface crosses above the upper one anywhere.
fn self_intersects(f: &Foil) -> bool {
    let ile = f.le_index();
    let n = f.n();
    if ile == 0 || ile + 1 >= n {
        return true;
    }
    let ux: Vec<f64> = f.x[..=ile].iter().rev().cloned().collect();
    let uy: Vec<f64> = f.y[..=ile].iter().rev().cloned().collect();
    let lx = &f.x[ile..];
    let ly = &f.y[ile..];
    let (x0, x1) = (f.x[ile], f.x[0].max(f.x[n - 1]));
    (1..200).any(|k| {
        let x = x0 + (x1 - x0) * k as f64 / 200.0;
        interp(&ux, &uy, x) - interp(lx, ly, x) < -1e-9
    })
}

/// Move and scale `base` to unit chord with its leading edge at the origin,
/// keeping its chord angle (so zero genes keep the base foil's angle of
/// attack).  Both candidate generation and baseline evaluation share this
/// frame so dimensional / shifted foils compare identically.
fn unit_chord_base(base: &Foil) -> Result<(Foil, Foil, usize, f64, f64), FoilError> {
    base.validate()?;
    let b = base.normalized();
    let n = b.n();
    let ile = b.le_index();
    let (xl, yl) = (base.x[ile], base.y[ile]);
    let (ang, ch) = {
        let (xt, yt) = (0.5 * (base.x[0] + base.x[n - 1]), 0.5 * (base.y[0] + base.y[n - 1]));
        ((yt - yl).atan2(xt - xl), ((xt - xl).powi(2) + (yt - yl).powi(2)).sqrt().max(1e-12))
    };
    let (sa, ca) = ang.sin_cos();
    let mut f = base.clone();
    f.recipe = None;
    for i in 0..n {
        f.x[i] = (base.x[i] - xl) / ch;
        f.y[i] = (base.y[i] - yl) / ch;
    }
    Ok((f, b, ile, sa, ca))
}

/// Build the foil described by `genes`: [upper shape variables..., lower
/// shape variables..., (flap angle)].  The changes are applied to the
/// normalised base foil and vanish at the leading and trailing edges:
/// Hicks-Henne bumps along the outward surface normal, CST weight changes
/// as `dy = x^0.5 (1-x) sum dA_i B_i(x)` on each surface.
pub fn apply_genes(base: &Foil, s: &OptimSettings, genes: &[f64]) -> Result<Foil, FoilError> {
    s.validate().map_err(FoilError::Geometry)?;
    let nb = s.shape_variables_per_side();
    if genes.len() != s.dimension() {
        return Err(FoilError::Geometry(format!("{} genes given, {} expected", genes.len(), s.dimension())));
    }
    if genes.iter().any(|g| !g.is_finite()) {
        return Err(FoilError::Geometry("non-finite design variable".into()));
    }
    let (mut f, b, ile, sa, ca) = unit_chord_base(base)?;
    let n = b.n();
    let mut displace = |i: usize, dxn: f64, dyn_: f64| {
        f.x[i] += ca * dxn - sa * dyn_;
        f.y[i] += sa * dxn + ca * dyn_;
    };
    let any_shape = nb > 0 && genes[..2 * nb].iter().any(|&g| g != 0.0);
    if any_shape {
        for i in 1..n - 1 {
            if i == ile {
                continue;
            }
            let off = if i < ile { 0 } else { nb };
            match s.parameterization {
                Parameterization::HicksHenne => {
                    // outward normal of a counter-clockwise contour
                    let tx = b.x[i + 1] - b.x[i - 1];
                    let ty = b.y[i + 1] - b.y[i - 1];
                    let tl = (tx * tx + ty * ty).sqrt();
                    if tl < 1e-15 {
                        continue;
                    }
                    let (nx, ny) = (ty / tl, -tx / tl);
                    let mut h = 0.0;
                    for j in 0..nb {
                        let t1 = s.t1 * (j + 1) as f64 / (2 * nb + 1) as f64;
                        h += hicks_henne(b.x[i], t1, s.t2) * genes[off + j];
                    }
                    displace(i, nx * h, ny * h);
                }
                Parameterization::Cst => {
                    let basis = crate::gen::cst::basis(nb - 1, 0.5, 1.0, b.x[i]);
                    displace(i, 0.0, basis.iter().zip(&genes[off..off + nb]).map(|(p, g)| p * g).sum::<f64>());
                }
            }
        }
    }
    f.validate()?;
    if self_intersects(&f) {
        return Err(FoilError::Geometry("upper and lower surfaces cross".into()));
    }
    if s.flap.is_some() {
        let deg = genes[2 * nb];
        if deg.abs() > 1e-9 {
            f = f.with_flap(s.hinge.0, s.hinge.1, deg)?;
        }
    }
    f.note = match s.parameterization {
        Parameterization::HicksHenne => "optimised (Hicks-Henne)".into(),
        Parameterization::Cst => format!("optimised (CST order {})", s.cst_order),
    };
    Ok(f)
}

// ---------------------------------------------------------------- evaluation

enum EvalMode {
    Single {
        specs: Vec<ObjectiveSpec>,
        scales: Vec<f64>,
    },
    Multi {
        ops: Vec<OperatingPoint>,
        base_cd: Vec<f64>,
        /// Indices into `ops` that form Pareto score dimensions (`Minimise`/`Target`).
        /// When empty (all rows are drag caps), all rows collapse into one weighted CD ratio.
        obj_rows: Vec<usize>,
        /// Normalised scale per score dimension so `sum(score / scale)` is O(1).
        scales: Vec<f64>,
        /// Sum of cap weights (for normalising drag-cap violations).
        cap_weight: f64,
    },
}

struct Ctx<'a> {
    base: &'a Foil,
    s: &'a OptimSettings,
    mode: EvalMode,
    prog: &'a Progress,
    done: AtomicUsize,
    total: usize,
}

impl Ctx<'_> {
    fn score_dim(&self) -> usize {
        match &self.mode {
            EvalMode::Single { specs, .. } => specs.len(),
            EvalMode::Multi { scales, .. } => scales.len(),
        }
    }

    fn infeasible(&self, genes: Vec<f64>, foil: Foil, spar_thickness: Option<f64>, op_points: Vec<OpPointEval>) -> Candidate {
        let k = self.score_dim();
        Candidate {
            genes,
            foil,
            cl: 0.0,
            cd: 0.0,
            cm: 0.0,
            cpmin: 0.0,
            objectives: vec![f64::NAN; k],
            scores: vec![LARGE; k],
            feasible: false,
            op_points,
            spar_thickness,
            violation: LARGE,
        }
    }

    fn evaluate(&self, genes: Vec<f64>) -> Result<Candidate, FoilError> {
        foiley::set_quiet(true);
        self.prog.check()?;
        let c = self.evaluate_inner(genes)?;
        self.record_eval(&c)?;
        Ok(c)
    }

    fn record_eval(&self, c: &Candidate) -> Result<(), FoilError> {
        if c.feasible {
            self.prog.count_ok();
        } else {
            self.prog.count_failed();
        }
        let d = self.done.fetch_add(1, Ordering::Relaxed) + 1;
        self.prog.set(d as f64 / self.total.max(1) as f64);
        self.prog.check()?;
        Ok(())
    }

    fn score_multi_evals(&self, genes: Vec<f64>, foil: Foil, spar_thickness: Option<f64>, op_points: Vec<OpPointEval>) -> Candidate {
        let EvalMode::Multi { ops, base_cd, obj_rows, .. } = &self.mode else {
            unreachable!();
        };
        let all_ok = op_points.len() == ops.len()
            && op_points.iter().zip(ops).all(|(ev, op)| {
                ev.converged
                    && ev.cd > 0.0
                    && [ev.alpha, ev.cl, ev.cd, ev.cm, ev.cpmin].iter().all(|v| v.is_finite())
                    && (ev.cl - op.cl).abs() < 0.015
                    && ev.alpha.abs() < 30.0
            });
        if !all_ok {
            return self.infeasible(genes, foil, spar_thickness, op_points);
        }
        let mut violation = 0.0;
        let mut any_cap_violated = false;
        for (i, (ev, op)) in op_points.iter().zip(ops).enumerate() {
            let limit = match op.goal {
                DragGoal::AtMostOriginal => Some(base_cd[i]),
                DragGoal::AtMost { cd } => Some(cd),
                _ => None,
            };
            if let Some(lim) = limit {
                if ev.cd > lim + 1e-9 {
                    any_cap_violated = true;
                    let rel = (ev.cd - lim) / lim.max(1e-6);
                    // Keep violation strictly positive even for tiny positive weights.
                    let weighted = (op.weight * rel).max(rel * 1e-12);
                    violation += if weighted > 0.0 && weighted.is_finite() { weighted } else { rel.max(f64::MIN_POSITIVE) };
                }
            }
        }
        if any_cap_violated && (violation <= 0.0 || !violation.is_finite()) {
            violation = f64::MIN_POSITIVE;
        }
        let (objectives, scores) = if !obj_rows.is_empty() {
            let mut objs = Vec::with_capacity(obj_rows.len());
            let mut scs = Vec::with_capacity(obj_rows.len());
            for &i in obj_rows {
                let ev = &op_points[i];
                objs.push(ev.cd);
                scs.push(match ops[i].goal {
                    DragGoal::Minimise => ev.cd,
                    DragGoal::Target { cd } => (ev.cd - cd).abs(),
                    _ => ev.cd,
                });
            }
            (objs, scs)
        } else {
            let wsum: f64 = ops.iter().map(|o| o.weight).sum::<f64>().max(1e-12);
            let ratio: f64 = op_points.iter().enumerate().map(|(i, ev)| ops[i].weight * ev.cd / base_cd[i].max(1e-6)).sum::<f64>() / wsum;
            (vec![ratio], vec![ratio])
        };
        let rep_i = obj_rows.first().copied().unwrap_or(0);
        let rep = &op_points[rep_i];
        Candidate {
            genes,
            foil,
            cl: rep.cl,
            cd: rep.cd,
            cm: rep.cm,
            cpmin: rep.cpmin,
            objectives,
            scores,
            feasible: true,
            op_points,
            spar_thickness,
            violation,
        }
    }

    fn evaluate_inner(&self, genes: Vec<f64>) -> Result<Candidate, FoilError> {
        let foil = match apply_genes(self.base, self.s, &genes) {
            Ok(f) => f,
            Err(FoilError::Cancelled) => return Err(FoilError::Cancelled),
            Err(_) => return Ok(self.infeasible(genes, self.base.clone(), None, Vec::new())),
        };
        let spar_thickness = self.s.spar.map(|sp| foil.thickness_at(sp.x));
        if let (Some(sp), Some(t)) = (self.s.spar, spar_thickness) {
            if t + 1e-6 < sp.min_thickness {
                return Ok(self.infeasible(genes, foil, spar_thickness, Vec::new()));
            }
        }
        match &self.mode {
            EvalMode::Single { specs, .. } => {
                let pt = match analyze_point(&foil, &self.s.analysis, self.s.alpha, true, &Progress::none()) {
                    Ok(Some(d)) => d.point,
                    Ok(None) => return Ok(self.infeasible(genes, foil, spar_thickness, Vec::new())),
                    Err(FoilError::Cancelled) => return Err(FoilError::Cancelled),
                    Err(_) => return Ok(self.infeasible(genes, foil, spar_thickness, Vec::new())),
                };
                let vals: Vec<f64> = specs.iter().map(|o| o.objective.value(pt.cl, pt.cd, pt.cm, pt.cpmin)).collect();
                let ok = pt.cd > 0.0 && [pt.cl, pt.cd, pt.cm, pt.cpmin].iter().all(|v| v.is_finite()) && vals.iter().all(|v| v.is_finite());
                if !ok {
                    return Ok(self.infeasible(genes, foil, spar_thickness, Vec::new()));
                }
                let scores = specs.iter().zip(&vals).map(|(o, &v)| o.goal.score(v)).collect();
                Ok(Candidate {
                    genes,
                    foil,
                    cl: pt.cl,
                    cd: pt.cd,
                    cm: pt.cm,
                    cpmin: pt.cpmin,
                    objectives: vals,
                    scores,
                    feasible: true,
                    op_points: Vec::new(),
                    spar_thickness,
                    violation: 0.0,
                })
            }
            EvalMode::Multi { ops, .. } => {
                let specs: Vec<(f64, f64, bool)> = ops.iter().map(|op| (op.re, op.cl, false)).collect();
                // Zero-span sub-view shares only cancellation with the optimiser progress handle.
                let cancel_view = self.prog.sub(0.0, 0.0);
                let pts = match analyze_points(&foil, &self.s.analysis, &specs, &cancel_view) {
                    Ok(v) => v,
                    Err(FoilError::Cancelled) => return Err(FoilError::Cancelled),
                    Err(_) => return Ok(self.infeasible(genes, foil, spar_thickness, Vec::new())),
                };
                let op_points: Vec<OpPointEval> = ops
                    .iter()
                    .enumerate()
                    .zip(pts)
                    .map(|((i, op), pt)| make_op_eval(i, op, pt))
                    .collect();
                Ok(self.score_multi_evals(genes, foil, spar_thickness, op_points))
            }
        }
    }

    fn evaluate_all(&self, genes: Vec<Vec<f64>>) -> Result<Vec<Candidate>, FoilError> {
        genes.into_par_iter().map(|g| self.evaluate(g)).collect()
    }

    /// Scalarised score (lower is better) among feasible candidates.
    fn scalar(&self, c: &Candidate) -> f64 {
        if !c.feasible {
            return LARGE;
        }
        match &self.mode {
            EvalMode::Single { scales, .. } => c.scores.iter().zip(scales).map(|(s, k)| s / k).sum(),
            EvalMode::Multi { scales, cap_weight, .. } => {
                let obj: f64 = c.scores.iter().zip(scales).map(|(s, k)| s / k).sum();
                if !c.constraints_met() {
                    1.0e4 + 1.0e3 * (c.violation / cap_weight.max(1e-12)) + obj
                } else {
                    obj
                }
            }
        }
    }

    /// Strict feasibility-first ordering (smaller is better):
    /// constraint-satisfying < cap-violating (by smaller violation) < infeasible.
    fn compare_candidate(&self, a: &Candidate, b: &Candidate) -> std::cmp::Ordering {
        match (a.feasible, b.feasible) {
            (true, false) => return std::cmp::Ordering::Less,
            (false, true) => return std::cmp::Ordering::Greater,
            (false, false) => return std::cmp::Ordering::Equal,
            (true, true) => {}
        }
        match (a.constraints_met(), b.constraints_met()) {
            (true, false) => std::cmp::Ordering::Less,
            (false, true) => std::cmp::Ordering::Greater,
            (false, false) => a.violation.total_cmp(&b.violation),
            (true, true) => self.scalar(a).total_cmp(&self.scalar(b)),
        }
    }

    /// Metric recorded in `OptimResult::history` for `best`.
    fn history_metric(&self, best: Option<&Candidate>) -> f64 {
        let Some(b) = best else { return f64::NAN };
        if !b.constraints_met() {
            return f64::NAN;
        }
        match &self.mode {
            EvalMode::Single { .. } => b.objectives.first().copied().unwrap_or(f64::NAN),
            EvalMode::Multi { scales, .. } => b.scores.iter().zip(scales).map(|(s, k)| s / k).sum(),
        }
    }

    /// All enabled objectives satisfied (only possible with target goals).
    fn satisfied(&self, c: &Candidate) -> bool {
        if !c.constraints_met() {
            return false;
        }
        match &self.mode {
            EvalMode::Single { specs, .. } => specs.iter().zip(&c.scores).all(|(o, &sc)| matches!(o.goal, Goal::Target { .. }) && sc == 0.0),
            EvalMode::Multi { ops, .. } => !ops.is_empty() && ops.iter().all(|o| matches!(o.goal, DragGoal::Target { .. })) && c.scores.iter().all(|&sc| sc <= 1e-5),
        }
    }
}

fn make_op_eval(index: usize, op: &OperatingPoint, pt: Option<crate::polar::FoilPoint>) -> OpPointEval {
    let label = op.display_label(index);
    match pt {
        Some(p) => OpPointEval {
            label,
            re: op.re,
            target_cl: op.cl,
            goal: op.goal,
            weight: op.weight,
            converged: true,
            alpha: p.alpha,
            cl: p.cl,
            cd: p.cd,
            cm: p.cm,
            cpmin: p.cpmin,
            xtr_top: p.xtr_top,
            xtr_bot: p.xtr_bot,
        },
        None => OpPointEval {
            label,
            re: op.re,
            target_cl: op.cl,
            goal: op.goal,
            weight: op.weight,
            converged: false,
            alpha: f64::NAN,
            cl: f64::NAN,
            cd: f64::NAN,
            cm: f64::NAN,
            cpmin: f64::NAN,
            xtr_top: f64::NAN,
            xtr_bot: f64::NAN,
        },
    }
}

/// Insert a candidate into a non-dominated archive.
fn archive_insert(archive: &mut Vec<Candidate>, c: &Candidate) -> bool {
    if !c.feasible {
        return false;
    }
    if archive.iter().any(|a| a.dominates(c) || (a.violation == c.violation && a.scores == c.scores)) {
        return false;
    }
    archive.retain(|a| !c.dominates(a));
    archive.push(c.clone());
    true
}

/// Trim the archive to `max` members, removing the most crowded ones
/// (smallest distance to a neighbour in normalised score space).
fn archive_trim(archive: &mut Vec<Candidate>, max: usize) {
    let max = max.max(1);
    while archive.len() > max {
        let k = archive[0].scores.len();
        let mut lo = vec![f64::MAX; k];
        let mut hi = vec![f64::MIN; k];
        for a in archive.iter() {
            for j in 0..k {
                lo[j] = lo[j].min(a.scores[j]);
                hi[j] = hi[j].max(a.scores[j]);
            }
        }
        let norm = |a: &Candidate, b: &Candidate| -> f64 {
            (0..k).map(|j| ((a.scores[j] - b.scores[j]) / (hi[j] - lo[j]).max(1e-12)).powi(2)).sum::<f64>()
        };
        let mut worst = 0;
        let mut dmin = f64::MAX;
        for i in 0..archive.len() {
            // keep the extremes of each objective
            if (0..k).any(|j| archive[i].scores[j] == lo[j]) {
                continue;
            }
            let d = (0..archive.len()).filter(|&o| o != i).map(|o| norm(&archive[i], &archive[o])).fold(f64::MAX, f64::min);
            if d < dmin {
                dmin = d;
                worst = i;
            }
        }
        if dmin == f64::MAX {
            worst = archive.len() - 1;
        }
        archive.remove(worst);
    }
}

fn clamp_genes(g: &mut [f64], bounds: &[(f64, f64)]) {
    for (v, &(lo, hi)) in g.iter_mut().zip(bounds) {
        *v = v.clamp(lo, hi);
    }
}

fn random_genes(rng: &mut Rng, bounds: &[(f64, f64)]) -> Vec<f64> {
    bounds.iter().map(|&(lo, hi)| rng.range(lo, hi)).collect()
}

/// Draw random genes, preferring shapes that satisfy the spar thickness
/// constraint when one is active (keeps the RNG sequence identical when
/// `s.spar` is `None`).
fn random_feasible_genes(rng: &mut Rng, base: &Foil, s: &OptimSettings, bounds: &[(f64, f64)]) -> Vec<f64> {
    let mut g = random_genes(rng, bounds);
    let Some(sp) = s.spar else { return g };
    for _ in 0..16 {
        if let Ok(f) = apply_genes(base, s, &g) {
            if f.thickness_at(sp.x) + 1e-6 >= sp.min_thickness {
                return g;
            }
        }
        g = random_genes(rng, bounds);
    }
    g
}

/// Perturb `center` within `scale * half_width` of each gene bound, screening
/// against an active spar constraint.  Used in multi-point mode so constrained
/// high-dimensional CST / HH searches explore smoothly around feasible shapes.
fn local_feasible_genes(rng: &mut Rng, base: &Foil, s: &OptimSettings, bounds: &[(f64, f64)], center: &[f64], scale: f64) -> Vec<f64> {
    let draw = |rng: &mut Rng| -> Vec<f64> {
        center
            .iter()
            .zip(bounds)
            .map(|(&c, &(lo, hi))| {
                let half = 0.5 * (hi - lo) * scale;
                (c + rng.range(-half, half)).clamp(lo, hi)
            })
            .collect()
    };
    let mut g = draw(rng);
    let Some(sp) = s.spar else { return g };
    for _ in 0..16 {
        if let Ok(f) = apply_genes(base, s, &g) {
            if f.thickness_at(sp.x) + 1e-6 >= sp.min_thickness {
                return g;
            }
        }
        g = draw(rng);
    }
    g
}

struct Particle {
    pos: Vec<f64>,
    vel: Vec<f64>,
    /// personal bests (position, scores, violation)
    best: Vec<(Vec<f64>, Vec<f64>, f64)>,
    cand: Option<Candidate>,
}

impl Particle {
    fn random(rng: &mut Rng, base: &Foil, s: &OptimSettings, bounds: &[(f64, f64)], nbest: usize) -> Particle {
        let pos = random_feasible_genes(rng, base, s, bounds);
        let vel = bounds.iter().map(|&(lo, hi)| {
            let dv = (hi - lo) / 2.0;
            rng.range(-dv / 2.0, dv / 2.0)
        });
        let vel = vel.collect();
        Particle { best: vec![(pos.clone(), Vec::new(), f64::INFINITY); nbest], pos, vel, cand: None }
    }

    fn around(rng: &mut Rng, base: &Foil, s: &OptimSettings, bounds: &[(f64, f64)], center: &[f64], scale: f64, nbest: usize) -> Particle {
        let pos = local_feasible_genes(rng, base, s, bounds, center, scale);
        let vel = bounds
            .iter()
            .map(|&(lo, hi)| {
                let dv = (hi - lo) * 0.25 * scale;
                rng.range(-dv, dv)
            })
            .collect();
        Particle { best: vec![(pos.clone(), Vec::new(), f64::INFINITY); nbest], pos, vel, cand: None }
    }

    /// Constrained xflr5 `updateBest`: replace the first personal best that
    /// the current feasible position improves on (lower violation, or better in
    /// any objective once constraints are met).
    fn update_best(&mut self) {
        let Some(c) = &self.cand else { return };
        if !c.feasible {
            return;
        }
        for b in self.best.iter_mut() {
            let better = if b.1.is_empty() {
                true
            } else if c.violation > 0.0 || b.2 > 0.0 {
                c.violation < b.2
            } else {
                b.1.iter().zip(&c.scores).any(|(o, n)| o > n)
            };
            if better {
                b.0 = self.pos.clone();
                b.1 = c.scores.clone();
                b.2 = c.violation;
                return;
            }
        }
    }
}

/// Run the optimisation.  `on_iteration` is called after the initial
/// population (iteration 0) and after every iteration with a snapshot.
pub fn optimise(base: &Foil, s: &OptimSettings, prog: &Progress, on_iteration: &(dyn Fn(usize, &OptimResult) + Sync)) -> Result<OptimResult, FoilError> {
    base.validate()?;
    s.validate().map_err(FoilError::Geometry)?;
    s.analysis.validate().map_err(FoilError::Geometry)?;
    if s.dimension() == 0 {
        return Err(FoilError::Geometry("no design variables (bumps per side = 0 and no flap)".into()));
    }
    let multipoint = s.is_multipoint();
    let ops = if multipoint {
        if !s.analysis.viscous {
            return Err(FoilError::Geometry("operating-point drag optimisation requires viscous analysis".into()));
        }
        let ops = s.enabled_operating_points();
        if ops.is_empty() {
            return Err(FoilError::Geometry("no operating point is enabled".into()));
        }
        ops
    } else {
        Vec::new()
    };
    let single_specs = if multipoint {
        Vec::new()
    } else {
        let sp = s.enabled();
        if sp.is_empty() {
            return Err(FoilError::Geometry("no optimisation objective is enabled".into()));
        }
        sp
    };

    let pop = s.population.max(2);
    let total = pop * (s.iterations + 1); // bounded by validate()
    let bounds = s.bounds();
    let mut rng = Rng::new(s.seed);
    let algo = s.algorithm.label();
    // Solve baseline and candidates in the exact same unit-chord frame so
    // dimensional or translated base foils produce consistent baseline drag
    // and `AtMostOriginal` caps.
    let (base_frame, _, _, _, _) = unit_chord_base(base)?;
    let base_max_t = base_frame.props().thickness;
    let base_spar_t = s.spar.map(|sp| base_frame.thickness_at(sp.x));
    let neutral_is_base = s.flap.is_none_or(|(lo, hi)| lo.min(hi) <= 0.0 && 0.0 <= lo.max(hi));

    if let (Some(sp), Some(t)) = (s.spar, base_spar_t) {
        if t + 1e-6 < sp.min_thickness {
            prog.warn(format!(
                "original foil '{}' is thinner ({:.2}%) than the spar minimum ({:.2}%) at x/c {:.2}",
                base.name,
                t * 100.0,
                sp.min_thickness * 100.0,
                sp.x
            ));
        }
    }

    // Evaluate baseline (original foil) first for multi-point mode so we know
    // original drag per row (required for `AtMostOriginal` and scaling).
    let (mut ctx, baseline, cached_base_evals) = if multipoint {
        prog.stage(format!("{algo}: evaluating original foil across {} operating points", ops.len()));
        prog.info(format!("{algo} multi-point optimisation of {}: {} variables, {} operating point(s), population {pop}, {} iterations", base.name, s.dimension(), ops.len(), s.iterations));
        if let (Some(sp), Some(t)) = (s.spar, base_spar_t) {
            prog.info(format!("spar constraint: t/c >= {:.2}% at x/c = {:.2} (original {:.2}%)", sp.min_thickness * 100.0, sp.x, t * 100.0));
        }
        let specs: Vec<(f64, f64, bool)> = ops.iter().map(|op| (op.re, op.cl, false)).collect();
        let base_pts = analyze_points(&base_frame, &s.analysis, &specs, prog)?;
        let op_evals: Vec<OpPointEval> = ops
            .iter()
            .enumerate()
            .zip(base_pts)
            .map(|((i, op), pt)| make_op_eval(i, op, pt))
            .collect();
        for (i, (op, ev)) in ops.iter().zip(&op_evals).enumerate() {
            let valid = ev.converged && ev.cd > 0.0 && ev.cd.is_finite() && (ev.cl - op.cl).abs() < 0.015;
            if !valid {
                let lbl = op.display_label(i);
                if matches!(op.goal, DragGoal::AtMostOriginal) {
                    return Err(FoilError::Geometry(format!(
                        "original foil '{}' did not converge at '{}' (Re {:.0}, CL {:.3}), which has a '<= original' drag goal",
                        base.name, lbl, op.re, op.cl
                    )));
                }
                prog.warn(format!("original foil '{}' did not converge at '{}' (Re {:.0}, CL {:.3})", base.name, lbl, op.re, op.cl));
            }
        }
        let base_cd: Vec<f64> = ops
            .iter()
            .zip(&op_evals)
            .map(|(op, ev)| {
                if ev.converged && ev.cd > 1e-6 && ev.cd.is_finite() {
                    ev.cd
                } else if let DragGoal::AtMost { cd } | DragGoal::Target { cd } = op.goal {
                    cd.max(1e-4)
                } else {
                    0.01
                }
            })
            .collect();
        let obj_rows: Vec<usize> = ops.iter().enumerate().filter(|(_, o)| !o.goal.is_cap()).map(|(i, _)| i).collect();
        let scales: Vec<f64> = if !obj_rows.is_empty() {
            let wsum: f64 = obj_rows.iter().map(|&i| ops[i].weight).sum::<f64>().max(1e-12);
            obj_rows.iter().map(|&i| base_cd[i] * wsum / ops[i].weight).collect()
        } else {
            vec![1.0]
        };
        let cap_weight: f64 = ops.iter().filter(|o| o.goal.is_cap()).map(|o| o.weight).sum::<f64>().max(1.0);
        let rep_i = obj_rows.first().copied().unwrap_or(0);
        let rep = &op_evals[rep_i];
        let baseline = BaselineEval {
            foil_name: base.name.clone(),
            cl: rep.cl,
            cd: rep.cd,
            cm: rep.cm,
            cpmin: rep.cpmin,
            converged: op_evals.iter().all(|e| e.converged),
            op_points: op_evals.clone(),
            spar_thickness: base_spar_t,
            max_thickness: base_max_t,
        };
        let ctx = Ctx {
            base,
            s,
            mode: EvalMode::Multi { ops, base_cd, obj_rows, scales, cap_weight },
            prog,
            done: AtomicUsize::new(0),
            total,
        };
        (ctx, Some(baseline), Some(op_evals))
    } else {
        prog.stage(format!("{algo}: evaluating initial population"));
        prog.info(format!("{algo} optimisation of {}: {} variables, {} objective(s), population {pop}, {} iterations", base.name, s.dimension(), single_specs.len(), s.iterations));
        let ctx = Ctx {
            base,
            s,
            mode: EvalMode::Single { scales: vec![1.0; single_specs.len()], specs: single_specs.clone() },
            prog,
            done: AtomicUsize::new(0),
            total,
        };
        (ctx, None, None)
    };

    // initial population: the base foil plus random designs (in multi-point mode,
    // half explore locally around the base foil so constrained runs improve promptly)
    let neutral = s.neutral_genes();
    let mut init = vec![neutral.clone()];
    let local_count = if multipoint { (pop - 1) / 2 } else { 0 };
    for k in 0..(pop - 1) {
        if k < local_count {
            init.push(local_feasible_genes(&mut rng, base, s, &bounds, &neutral, 0.35));
        } else {
            init.push(random_feasible_genes(&mut rng, base, s, &bounds));
        }
    }
    let cands = if let (Some(base_evals), true) = (cached_base_evals, neutral_is_base) {
        // Reuse the already-solved unit-chord baseline operating points for
        // cands[0] when neutral_genes() reproduces base_frame (saves 1 candidate
        // evaluation and guarantees exact agreement with AtMostOriginal).
        let c0 = match apply_genes(base, s, &init[0]) {
            Ok(foil) => {
                let spar_ok = s.spar.is_none_or(|sp| base_spar_t.unwrap_or(0.0) + 1e-6 >= sp.min_thickness);
                if spar_ok {
                    ctx.score_multi_evals(init[0].clone(), foil, base_spar_t, base_evals)
                } else {
                    ctx.infeasible(init[0].clone(), foil, base_spar_t, base_evals)
                }
            }
            Err(_) => ctx.infeasible(init[0].clone(), base_frame.clone(), base_spar_t, base_evals),
        };
        ctx.record_eval(&c0)?;
        let mut rest = ctx.evaluate_all(init[1..].to_vec())?;
        let mut all = Vec::with_capacity(pop);
        all.push(c0);
        all.append(&mut rest);
        all
    } else {
        ctx.evaluate_all(init.clone())?
    };

    // Single-point normalisation scales and baseline from cands[0] (or unit-chord base foil).
    let baseline = if let Some(b) = baseline {
        b
    } else {
        let c0 = &cands[0];
        if let EvalMode::Single { specs, scales } = &mut ctx.mode {
            *scales = specs
                .iter()
                .enumerate()
                .map(|(j, o)| {
                    let v = if c0.feasible { c0.objectives[j].abs() } else { 0.0 };
                    match o.goal {
                        Goal::Target { value, tolerance } => v.max(value.abs()).max(tolerance.abs()).max(1e-6),
                        _ => if v > 1e-6 { v } else { 1.0 },
                    }
                })
                .collect();
        }
        if !c0.feasible {
            prog.warn(format!("the base foil did not converge at alpha = {}", s.alpha));
        }
        if neutral_is_base && c0.feasible {
            BaselineEval {
                foil_name: base.name.clone(),
                cl: c0.cl,
                cd: c0.cd,
                cm: c0.cm,
                cpmin: c0.cpmin,
                converged: true,
                op_points: Vec::new(),
                spar_thickness: base_spar_t,
                max_thickness: base_max_t,
            }
        } else {
            let pt = analyze_point(&base_frame, &s.analysis, s.alpha, true, &Progress::none()).ok().flatten().map(|d| d.point);
            match pt {
                Some(p) => BaselineEval {
                    foil_name: base.name.clone(),
                    cl: p.cl,
                    cd: p.cd,
                    cm: p.cm,
                    cpmin: p.cpmin,
                    converged: true,
                    op_points: Vec::new(),
                    spar_thickness: base_spar_t,
                    max_thickness: base_max_t,
                },
                None => BaselineEval {
                    foil_name: base.name.clone(),
                    spar_thickness: base_spar_t,
                    max_thickness: base_max_t,
                    ..Default::default()
                },
            }
        }
    };

    let mut result = OptimResult { baseline: Some(baseline), ..Default::default() };
    let mut archive: Vec<Candidate> = Vec::new();
    // Track the best constraint-satisfying compromise independently of archive
    // crowding so crowding trim can never discard it and infeasible/cap-violating
    // candidates are never delivered as `result.best` or `result.pareto`.
    let mut best_feasible: Option<Candidate> = None;
    let update_best = |best_feasible: &mut Option<Candidate>, c: &Candidate, ctx: &Ctx| {
        if c.constraints_met() && best_feasible.as_ref().is_none_or(|b| ctx.scalar(c) < ctx.scalar(b)) {
            *best_feasible = Some(c.clone());
        }
    };
    for c in &cands {
        archive_insert(&mut archive, c);
        update_best(&mut best_feasible, c, &ctx);
    }
    archive_trim(&mut archive, s.archive_size);

    let report = |it: usize, archive: &[Candidate], best_feasible: Option<&Candidate>, result: &mut OptimResult, ctx: &Ctx| -> bool {
        let delivered_pareto: Vec<Candidate> = archive.iter().filter(|a| a.constraints_met()).cloned().collect();
        result.best = best_feasible.cloned();
        result.pareto = delivered_pareto;
        result.evaluations = ctx.done.load(Ordering::Relaxed);
        let v = ctx.history_metric(result.best.as_ref());
        result.history.push((it, v));
        match &ctx.mode {
            EvalMode::Single { specs, .. } => {
                ctx.prog.stage(format!("{algo} iteration {it}/{} - best {} {}", s.iterations, specs[0].objective.label(), fmt_val(v)));
            }
            EvalMode::Multi { .. } => {
                let status = match &result.best {
                    Some(_) => format!("weighted CD {:.3}x orig", v),
                    None if archive.iter().any(|a| a.feasible) => "seeking feasible drag caps".into(),
                    None => "no converged design".into(),
                };
                ctx.prog.stage(format!("{algo} iteration {it}/{} - {status}", s.iterations));
            }
        }
        on_iteration(it, result);
        best_feasible.is_some_and(|b| ctx.satisfied(b)) || result.pareto.iter().any(|a| ctx.satisfied(a))
    };

    let mut done = report(0, &archive, best_feasible.as_ref(), &mut result, &ctx);

    match s.algorithm {
        Algorithm::Mopso => {
            let nbest = ctx.score_dim().clamp(1, 3);
            let mut swarm: Vec<Particle> = init
                .into_iter()
                .zip(cands)
                .map(|(g, c)| {
                    let mut p = Particle::random(&mut rng, base, s, &bounds, nbest);
                    p.pos = g;
                    p.best = vec![(p.pos.clone(), Vec::new(), f64::INFINITY); nbest];
                    p.cand = Some(c);
                    p.update_best();
                    p
                })
                .collect();
            let mut it = 1;
            while it <= s.iterations && !done {
                prog.check()?;
                // move (serial: all random draws happen here)
                for p in swarm.iter_mut() {
                    let feasible = p.cand.as_ref().map(|c| c.feasible).unwrap_or(false);
                    let in_front = p.cand.as_ref().map(|c| archive.iter().any(|a| a.genes == c.genes)).unwrap_or(false);
                    let regen = (rng.f64() < s.regen_probability && !in_front) || !feasible || archive.is_empty();
                    if regen {
                        let np = if multipoint && !archive.is_empty() && rng.f64() < 0.6 {
                            let anchor = &archive[rng.below(archive.len())].genes;
                            Particle::around(&mut rng, base, s, &bounds, anchor, 0.30, nbest)
                        } else {
                            Particle::random(&mut rng, base, s, &bounds, nbest)
                        };
                        p.pos = np.pos;
                        p.vel = np.vel;
                        p.best = np.best;
                        continue;
                    }
                    let gbest = &archive[rng.below(archive.len())].genes;
                    let valid: Vec<usize> = (0..p.best.len()).filter(|&i| !p.best[i].1.is_empty()).collect();
                    let pbest = if valid.is_empty() { p.pos.clone() } else { p.best[valid[rng.below(valid.len())]].0.clone() };
                    for j in 0..p.pos.len() {
                        let (r1, r2) = (rng.f64(), rng.f64());
                        let mut v = s.inertia * p.vel[j] + s.cognitive * r1 * (pbest[j] - p.pos[j]) + s.social * r2 * (gbest[j] - p.pos[j]);
                        let np = p.pos[j] + v;
                        if np < bounds[j].0 || np > bounds[j].1 {
                            v = -v;
                        }
                        p.vel[j] = v;
                        p.pos[j] += v;
                    }
                    clamp_genes(&mut p.pos, &bounds);
                }
                // evaluate (parallel)
                let cands = ctx.evaluate_all(swarm.iter().map(|p| p.pos.clone()).collect())?;
                for (p, c) in swarm.iter_mut().zip(cands) {
                    update_best(&mut best_feasible, &c, &ctx);
                    p.cand = Some(c);
                    p.update_best();
                    archive_insert(&mut archive, p.cand.as_ref().unwrap());
                }
                archive_trim(&mut archive, s.archive_size);
                done = report(it, &archive, best_feasible.as_ref(), &mut result, &ctx);
                it += 1;
            }
        }
        Algorithm::Genetic => {
            let mut popn: Vec<Candidate> = cands;
            let mut elite: Candidate = popn.iter().min_by(|a, b| ctx.compare_candidate(a, b)).cloned().unwrap();
            let mut it = 1;
            while it <= s.iterations && !done {
                prog.check()?;
                // Feasibility-first rank-based selection: constraint-satisfying
                // designs always rank ahead of cap violators without a large
                // penalty swamping selection pressure among feasible candidates.
                let mut feas_idx: Vec<usize> = (0..popn.len()).filter(|&i| popn[i].feasible).collect();
                feas_idx.sort_by(|&a, &b| ctx.compare_candidate(&popn[a], &popn[b]));
                let m = feas_idx.len();
                let mut w = vec![0.0; popn.len()];
                for (rank, &idx) in feas_idx.iter().enumerate() {
                    w[idx] = (m - rank) as f64;
                }
                let wsum: f64 = w.iter().sum();
                let mut genes: Vec<Vec<f64>> = (0..pop)
                    .map(|_| {
                        if wsum <= 0.0 {
                            return random_feasible_genes(&mut rng, base, s, &bounds);
                        }
                        let mut r = rng.f64() * wsum;
                        for (i, wi) in w.iter().enumerate() {
                            r -= wi;
                            if r <= 0.0 && *wi > 0.0 {
                                return popn[i].genes.clone();
                            }
                        }
                        popn[w.iter().rposition(|&x| x > 0.0).unwrap_or(0)].genes.clone()
                    })
                    .collect();
                // blend crossover on random pairs
                let alpha = 0.5;
                let mut order: Vec<usize> = (0..pop).collect();
                for i in (1..pop).rev() {
                    order.swap(i, rng.below(i + 1));
                }
                for pair in order.chunks(2) {
                    if pair.len() < 2 || rng.f64() >= s.crossover {
                        continue;
                    }
                    let (a, b) = (genes[pair[0]].clone(), genes[pair[1]].clone());
                    for child in pair {
                        let g = &mut genes[*child];
                        for j in 0..g.len() {
                            let f = rng.range(-alpha, 1.0 + alpha);
                            g[j] = f * a[j] + (1.0 - f) * b[j];
                        }
                    }
                }
                // Gaussian mutation
                for g in genes.iter_mut() {
                    for j in 0..g.len() {
                        if rng.f64() < s.mutation_probability {
                            let half = 0.5 * (bounds[j].1 - bounds[j].0);
                            g[j] += rng.normal() * s.mutation_sigma * half;
                        }
                    }
                    clamp_genes(g, &bounds);
                }
                popn = ctx.evaluate_all(genes)?;
                // elitism: the best design so far (feasibility-first) replaces the worst child
                let (iw, _) = popn.iter().enumerate().max_by(|a, b| ctx.compare_candidate(a.1, b.1)).unwrap();
                let (ib, cb) = popn.iter().enumerate().min_by(|a, b| ctx.compare_candidate(a.1, b.1)).unwrap();
                if ctx.compare_candidate(cb, &elite) == std::cmp::Ordering::Less {
                    elite = popn[ib].clone();
                } else {
                    popn[iw] = elite.clone();
                }
                for c in &popn {
                    archive_insert(&mut archive, c);
                    update_best(&mut best_feasible, c, &ctx);
                }
                archive_trim(&mut archive, s.archive_size);
                done = report(it, &archive, best_feasible.as_ref(), &mut result, &ctx);
                it += 1;
            }
        }
    }
    prog.set(1.0);
    if let Some(b) = &result.best {
        match &ctx.mode {
            EvalMode::Single { specs, .. } => {
                let parts: Vec<String> = specs.iter().zip(&b.objectives).map(|(o, v)| format!("{} {}", o.objective.label(), fmt_val(*v))).collect();
                prog.info(format!("{algo} finished after {} evaluations: {}", result.evaluations, parts.join(", ")));
            }
            EvalMode::Multi { .. } => {
                let parts: Vec<String> = b.op_points.iter().map(|e| format!("{} CD {:.5}", e.label, e.cd)).collect();
                prog.info(format!("{algo} finished after {} evaluations: {}", result.evaluations, parts.join(", ")));
            }
        }
    } else if multipoint {
        prog.warn(format!("{algo} finished without any design satisfying all constraints"));
    } else {
        prog.warn(format!("{algo} finished without any converged design"));
    }
    Ok(result)
}

fn fmt_val(v: f64) -> String {
    if !v.is_finite() {
        "-".into()
    } else if v.abs() >= 10.0 {
        format!("{v:.1}")
    } else {
        format!("{v:.4}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{Duration, Instant};

    fn ld_settings() -> OptimSettings {
        let mut s = OptimSettings { population: 8, iterations: 4, alpha: 4.0, seed: 7, ..Default::default() };
        s.analysis.re = 5.0e5;
        s.objectives = vec![ObjectiveSpec { objective: Objective::ClCd, enabled: true, goal: Goal::Maximise }];
        s
    }

    #[test]
    fn optim_zero_genes_is_base() {
        let base = Foil::naca("0012").unwrap();
        let s = OptimSettings { flap: Some((-5.0, 10.0)), ..Default::default() };
        let f = apply_genes(&base, &s, &s.neutral_genes()).unwrap();
        assert_eq!(f.n(), base.n());
        for i in 0..f.n() {
            assert!((f.x[i] - base.x[i]).abs() < 1e-9 && (f.y[i] - base.y[i]).abs() < 1e-9, "point {i}");
        }
        // a positive upper bump thickens, TE unchanged
        let mut g = s.neutral_genes();
        g[1] = 0.01;
        let f = apply_genes(&base, &s, &g).unwrap();
        assert!(f.props().thickness > base.props().thickness);
        assert!((f.y[0] - base.y[0]).abs() < 1e-12);
        // huge negative bumps on both sides cross over: rejected
        let s2 = OptimSettings { max_amplitude: 0.2, ..Default::default() };
        assert!(apply_genes(&base, &s2, &[-0.2; 6]).is_err());
    }

    #[test]
    fn cst_genes_change_the_shape_smoothly() {
        let base = Foil::naca("2412").unwrap();
        let s = OptimSettings { parameterization: Parameterization::Cst, cst_order: 4, flap: Some((-5.0, 5.0)), ..Default::default() };
        assert_eq!(s.dimension(), 11);
        let f = apply_genes(&base, &s, &s.neutral_genes()).unwrap();
        for i in 0..f.n() {
            assert_eq!((f.x[i], f.y[i]), (base.x[i], base.y[i]), "point {i}");
        }
        let b = base.normalized();
        // raising A0 of the upper surface enlarges the nose radius and thickens
        let mut g = s.neutral_genes();
        g[0] = 0.04;
        let f = apply_genes(&base, &s, &g).unwrap();
        let (p, q) = (b.props(), f.props());
        assert!(q.le_radius > p.le_radius && q.thickness > p.thickness, "{p:?} {q:?}");
        assert_eq!((f.y[0], f.y[f.n() - 1]), (b.y[0], b.y[b.n() - 1]));
        assert!(f.note.contains("CST"));
        // absurd orders and non-finite ranges are input errors, not panics
        let base0 = Foil::naca("0012").unwrap();
        for order in [usize::MAX, usize::MAX / 2, 31] {
            let bad = OptimSettings { parameterization: Parameterization::Cst, cst_order: order, ..Default::default() };
            assert!(matches!(optimise(&base0, &bad, &Progress::none(), &|_, _| {}), Err(FoilError::Geometry(_))));
            assert!(apply_genes(&base0, &bad, &[]).is_err());
            assert!(bad.bounds().is_empty() && bad.neutral_genes().is_empty());
        }
        // a NaN range is rejected even with a correctly sized gene vector
        let nan = OptimSettings { parameterization: Parameterization::Cst, cst_order: 2, cst_max_delta: f64::NAN, ..Default::default() };
        assert!(matches!(optimise(&base0, &nan, &Progress::none(), &|_, _| {}), Err(FoilError::Geometry(_))));
        assert!(apply_genes(&base0, &nan, &[0.0; 6]).is_err());
        // zero genes keep a cambered generated foil's angle (its LE is not
        // its leftmost point), and any base is brought to unit chord with
        // the LE at the origin, as in 1.1: shifted or dimensional
        // coordinates give the same candidate
        let tilted = Foil::naca("23112").unwrap();
        let z = OptimSettings { parameterization: Parameterization::Cst, cst_order: 3, ..Default::default() };
        let g = apply_genes(&tilted, &z, &z.neutral_genes()).unwrap();
        let chord_angle = |f: &Foil| {
            let (i, n) = (f.le_index(), f.n());
            (0.5 * (f.y[0] + f.y[n - 1]) - f.y[i]).atan2(0.5 * (f.x[0] + f.x[n - 1]) - f.x[i])
        };
        assert!((chord_angle(&g) - chord_angle(&tilted)).abs() < 1e-14);
        let (gn, tn) = (g.normalized(), tilted.normalized());
        assert!(gn.x.iter().zip(&tn.x).chain(gn.y.iter().zip(&tn.y)).all(|(a, b)| (a - b).abs() < 1e-12), "same shape");
        let mut moved = tilted.clone();
        for v in moved.x.iter_mut() {
            *v = 100.0 * *v + 300.0;
        }
        for v in moved.y.iter_mut() {
            *v *= 100.0;
        }
        let mut genes = z.neutral_genes();
        genes[1] = 0.02;
        let (a, b) = (apply_genes(&tilted, &z, &genes).unwrap(), apply_genes(&moved, &z, &genes).unwrap());
        assert!(a.x.iter().zip(&b.x).chain(a.y.iter().zip(&b.y)).all(|(p, q)| (p - q).abs() < 1e-12));
        let flapped = OptimSettings { flap: Some((-5.0, 10.0)), ..z.clone() };
        let mut fg = flapped.neutral_genes();
        fg[z.dimension()] = 8.0;
        apply_genes(&moved, &flapped, &fg).unwrap();
        // 1.1 accepted more than 50 bumps; huge workloads are refused
        assert!(OptimSettings { bumps_per_side: 51, ..Default::default() }.validate().is_ok());
        assert!(OptimSettings { iterations: usize::MAX, ..Default::default() }.validate().is_err());
        // settings written before CST existed still load as Hicks-Henne
        let old: OptimSettings = serde_json::from_str(r#"{"bumps_per_side": 4}"#).unwrap();
        assert_eq!((old.parameterization, old.bumps_per_side, old.dimension()), (Parameterization::HicksHenne, 4, 8));
    }

    #[test]
    fn optim_mopso_improves_ld() {
        let base = Foil::naca("0012").unwrap();
        let s = ld_settings();
        let calls = AtomicUsize::new(0);
        let prog = Progress::new();
        let r = optimise(&base, &s, &prog, &|_, _| {
            calls.fetch_add(1, Ordering::Relaxed);
        })
        .unwrap();
        let st = AnalysisSettings { re: 5.0e5, ..Default::default() };
        let p0 = analyze_point(&base, &st, 4.0, true, &Progress::none()).unwrap().unwrap().point;
        let best = r.best.expect("a best design");
        eprintln!("base L/D {:.2}, best {:.2}, evals {}, history {:?}", p0.ld(), best.objectives[0], r.evaluations, r.history);
        assert!(best.objectives[0] >= p0.ld() - 1e-6);
        assert_eq!(r.evaluations, 8 * 5);
        assert_eq!(calls.load(Ordering::Relaxed), 5);
        assert!((prog.fraction() - 1.0).abs() < 1e-9);
        assert!(r.history.windows(2).all(|w| w[1].1 >= w[0].1 - 1e-9));
    }

    #[test]
    fn optim_cancel_is_prompt() {
        let base = Foil::naca("0012").unwrap();
        let s = OptimSettings { iterations: 100, ..ld_settings() };
        let prog = Progress::new();
        prog.cancel();
        let t = Instant::now();
        assert_eq!(optimise(&base, &s, &prog, &|_, _| {}), Err(FoilError::Cancelled));
        assert!(t.elapsed() < Duration::from_secs(2));
        // cancelled from the callback after the first generation
        let prog = Progress::new();
        let p2 = prog.clone();
        let r = optimise(&base, &s, &prog, &move |_, _| p2.cancel());
        assert_eq!(r, Err(FoilError::Cancelled));
    }

    #[test]
    fn optim_genetic_runs() {
        let base = Foil::naca("0012").unwrap();
        let mut s = ld_settings();
        s.algorithm = Algorithm::Genetic;
        s.objectives = vec![
            ObjectiveSpec { objective: Objective::Cl, enabled: true, goal: Goal::Target { value: 0.5, tolerance: 0.005 } },
            ObjectiveSpec { objective: Objective::Cd, enabled: true, goal: Goal::Minimise },
        ];
        let r = optimise(&base, &s, &Progress::none(), &|_, _| {}).unwrap();
        let best = r.best.expect("a best design");
        assert!(best.feasible && !r.pareto.is_empty());
        assert!(r.pareto.iter().all(|c| c.objectives.len() == 2));
        eprintln!("GA best Cl {:.4} Cd {:.5}, pareto {}", best.cl, best.cd, r.pareto.len());
    }

    #[test]
    fn multipoint_glider_section_reduces_drag_keeps_slow_and_spar() {
        let base = Foil::naca("2412").unwrap();
        let s = OptimSettings {
            algorithm: Algorithm::Mopso,
            parameterization: Parameterization::Cst,
            cst_order: 4,
            cst_max_delta: 0.04,
            operating_points: vec![
                OperatingPoint { label: "Slow flight".into(), enabled: true, re: 1.15e5, cl: 0.90, goal: DragGoal::AtMostOriginal, weight: 1.0 },
                OperatingPoint { label: "Cruise".into(), enabled: true, re: 1.75e5, cl: 0.50, goal: DragGoal::Minimise, weight: 1.0 },
                OperatingPoint { label: "Fast glide".into(), enabled: true, re: 2.70e5, cl: 0.20, goal: DragGoal::Minimise, weight: 1.0 },
            ],
            spar: Some(SparConstraint { x: 0.30, min_thickness: 0.115 }),
            population: 10,
            iterations: 8,
            archive_size: 2, // also checks crowding trim never worsens MOPSO history
            seed: 42,
            ..Default::default()
        };
        let r = optimise(&base, &s, &Progress::none(), &|_, _| {}).unwrap();
        let base_eval = r.baseline.as_ref().expect("baseline evaluation");
        assert_eq!(base_eval.op_points.len(), 3);
        assert!(base_eval.op_points.iter().all(|p| p.converged));
        let base_spar = base_eval.spar_thickness.expect("baseline spar thickness");
        assert!((base_spar - 0.12).abs() < 2e-3, "base spar t/c = {base_spar}");

        let best = r.best.as_ref().expect("best candidate");
        assert!(best.constraints_met(), "violation {}", best.violation);
        assert_eq!(best.op_points.len(), 3);
        let best_spar = best.spar_thickness.expect("best spar thickness");
        assert!(best_spar + 1e-6 >= 0.115, "spar violated: {best_spar}");
        // Slow flight drag must stay at or below original
        assert!(
            best.op_points[0].cd <= base_eval.op_points[0].cd + 1e-7,
            "slow cd {} > orig {}",
            best.op_points[0].cd,
            base_eval.op_points[0].cd
        );
        // Weighted cruise + fast drag ratio must strictly improve and never worsen across iterations
        assert_eq!(r.history.len(), 9);
        assert!(r.history.last().unwrap().1 < 0.99, "expected >1% weighted drag reduction, got {:?}", r.history);
        assert!(r.history.windows(2).all(|w| w[1].1 <= w[0].1 + 1e-12), "history regressed: {:?}", r.history);
        assert!(
            best.op_points[1].cd < base_eval.op_points[1].cd || best.op_points[2].cd < base_eval.op_points[2].cd,
            "expected drag reduction in cruise or fast glide"
        );
        eprintln!(
            "multi-point best: slow {:.5} (orig {:.5}), cruise {:.5} (orig {:.5}), fast {:.5} (orig {:.5}), spar {:.2}%, ratio {:.4}",
            best.op_points[0].cd,
            base_eval.op_points[0].cd,
            best.op_points[1].cd,
            base_eval.op_points[1].cd,
            best.op_points[2].cd,
            base_eval.op_points[2].cd,
            best_spar * 100.0,
            r.history.last().unwrap().1
        );
    }

    #[test]
    fn scaled_base_foil_matches_unit_chord_baseline() {
        let unit = Foil::naca("2412").unwrap();
        let mut scaled = unit.clone();
        for x in &mut scaled.x {
            *x = 2.5 * *x + 12.0;
        }
        for y in &mut scaled.y {
            *y = 2.5 * *y - 4.0;
        }
        let s = OptimSettings {
            algorithm: Algorithm::Mopso,
            parameterization: Parameterization::Cst,
            cst_order: 3,
            cst_max_delta: 0.02,
            operating_points: vec![
                OperatingPoint { label: "Slow".into(), enabled: true, re: 1.15e5, cl: 0.90, goal: DragGoal::AtMostOriginal, weight: 1.0 },
                OperatingPoint { label: "Cruise".into(), enabled: true, re: 1.75e5, cl: 0.50, goal: DragGoal::Minimise, weight: 1.0 },
            ],
            population: 3,
            iterations: 0,
            seed: 5,
            ..Default::default()
        };
        let ru = optimise(&unit, &s, &Progress::none(), &|_, _| {}).unwrap();
        let rs = optimise(&scaled, &s, &Progress::none(), &|_, _| {}).unwrap();
        let bu = ru.baseline.as_ref().unwrap();
        let bs = rs.baseline.as_ref().unwrap();
        for (a, b) in bu.op_points.iter().zip(&bs.op_points) {
            assert!((a.cd - b.cd).abs() < 1e-7, "scaled baseline CD {} != unit {}", b.cd, a.cd);
        }
        // Re-evaluating the returned zero-gene foil must match the cached baseline CD.
        let best_s = rs.best.as_ref().unwrap();
        let check = analyze_points(&best_s.foil, &s.analysis, &[(1.15e5, 0.90, false)], &Progress::none()).unwrap();
        assert!((check[0].unwrap().cd - best_s.op_points[0].cd).abs() < 1e-7);
    }

    #[test]
    fn impossible_drag_cap_and_tiny_weight_deliver_no_best() {
        let base = Foil::naca("2412").unwrap();
        for algo in [Algorithm::Mopso, Algorithm::Genetic] {
            let s = OptimSettings {
                algorithm: algo,
                parameterization: Parameterization::Cst,
                cst_order: 2,
                cst_max_delta: 0.02,
                operating_points: vec![
                    OperatingPoint { label: "Impossible".into(), enabled: true, re: 1.5e5, cl: 0.50, goal: DragGoal::AtMost { cd: 1.0e-7 }, weight: 1.0e-300 },
                    OperatingPoint { label: "Cruise".into(), enabled: true, re: 2.0e5, cl: 0.40, goal: DragGoal::Minimise, weight: 1.0 },
                ],
                population: 4,
                iterations: 2,
                seed: 1,
                ..Default::default()
            };
            let r = optimise(&base, &s, &Progress::none(), &|_, _| {}).unwrap();
            assert!(r.best.is_none(), "{algo:?} must not return a cap-violating best design");
            assert!(r.pareto.is_empty(), "{algo:?} must not return cap-violating Pareto designs");
            assert!(r.baseline.is_some());
        }
    }

    #[test]
    fn multipoint_validation_and_ga_respect_spar_and_caps() {
        let base = Foil::naca("2412").unwrap();
        // Validation rejects empty enabled operating points, bad spar, bad weight, inviscid
        let mut bad = OptimSettings {
            operating_points: vec![OperatingPoint { enabled: false, ..Default::default() }],
            ..Default::default()
        };
        assert!(matches!(optimise(&base, &bad, &Progress::none(), &|_, _| {}), Err(FoilError::Geometry(_))));
        bad.operating_points[0].enabled = true;
        bad.analysis.viscous = false;
        assert!(matches!(optimise(&base, &bad, &Progress::none(), &|_, _| {}), Err(FoilError::Geometry(_))));
        bad.analysis.viscous = true;
        bad.spar = Some(SparConstraint { x: 1.5, min_thickness: 0.12 });
        assert!(bad.validate().is_err());
        bad.spar = Some(SparConstraint { x: 0.30, min_thickness: 0.115 });
        bad.operating_points[0].weight = 0.0;
        assert!(bad.validate().is_err());

        // GA multi-point respects spar thickness and drag cap
        let ga = OptimSettings {
            algorithm: Algorithm::Genetic,
            parameterization: Parameterization::Cst,
            cst_order: 3,
            cst_max_delta: 0.03,
            operating_points: vec![
                OperatingPoint { label: "Slow".into(), enabled: true, re: 1.3e5, cl: 0.85, goal: DragGoal::AtMostOriginal, weight: 1.0 },
                OperatingPoint { label: "Cruise".into(), enabled: true, re: 2.0e5, cl: 0.45, goal: DragGoal::Minimise, weight: 1.0 },
            ],
            spar: Some(SparConstraint { x: 0.30, min_thickness: 0.118 }),
            population: 6,
            iterations: 2,
            seed: 3,
            ..Default::default()
        };
        let r = optimise(&base, &ga, &Progress::none(), &|_, _| {}).unwrap();
        let best = r.best.as_ref().expect("GA best");
        let base_eval = r.baseline.as_ref().unwrap();
        assert!(best.constraints_met());
        assert!(best.spar_thickness.unwrap() + 1e-6 >= 0.118);
        assert!(best.op_points[0].cd <= base_eval.op_points[0].cd + 1e-7);
        assert!(best.op_points[1].cd <= base_eval.op_points[1].cd + 1e-7);

        // JSON roundtrip preserves operating_points and spar
        let json = serde_json::to_string(&ga).unwrap();
        let back: OptimSettings = serde_json::from_str(&json).unwrap();
        assert_eq!(back.operating_points, ga.operating_points);
        assert_eq!(back.spar, ga.spar);
    }
}
