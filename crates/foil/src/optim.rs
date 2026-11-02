//! Airfoil shape optimisation (equivalent of xflr5's Optim2d).
//!
//! The design variables are the amplitudes of Hicks-Henne bump functions
//! added to the upper and lower surfaces of the (normalised) base foil along
//! the local surface normal, plus optionally a trailing-edge flap angle.
//! Each candidate is analysed at one operating point with the XFoil engine,
//! and one or more objectives (Cl, Cd, Cl/Cd, Cpmin, Cm) are minimised,
//! maximised or driven to a target value.
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

use crate::analysis::{analyze_point, AnalysisSettings};
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
    /// Angle of attack of the design point, degrees.
    pub alpha: f64,
    pub analysis: AnalysisSettings,
    pub objectives: Vec<ObjectiveSpec>,
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

    /// Number of design variables.
    pub fn dimension(&self) -> usize {
        2 * self.bumps_per_side + usize::from(self.flap.is_some())
    }

    /// (min, max) of each design variable.
    pub fn bounds(&self) -> Vec<(f64, f64)> {
        let a = self.max_amplitude.abs();
        let mut b = vec![(-a, a); 2 * self.bumps_per_side];
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
            bumps_per_side: 3,
            t1: 1.7,
            t2: 1.0,
            max_amplitude: 0.01,
            flap: None,
            hinge: (0.7, 0.5),
            alpha: 0.0,
            analysis: AnalysisSettings { re: 1.0e6, mach: 0.0, ncrit: 9.0, xtr_top: 1.0, xtr_bot: 1.0, ..Default::default() },
            objectives: Self::default_objectives(),
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
    /// Value of each enabled objective (same order as the enabled specs).
    pub objectives: Vec<f64>,
    /// Per enabled objective, the score being minimised (-value for
    /// maximise, distance outside the tolerance band for targets).
    pub scores: Vec<f64>,
    /// Converged analysis on a valid geometry.
    pub feasible: bool,
}

impl Candidate {
    /// `a` dominates `b`: no worse in every score, better in one.
    fn dominates(&self, b: &Candidate) -> bool {
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
    /// (iteration, first enabled objective value of `best`), iteration 0 = initial population.
    pub history: Vec<(usize, f64)>,
    pub evaluations: usize,
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

/// Build the foil described by `genes`: [upper bumps..., lower bumps...,
/// (flap angle)].  Bumps are applied to the normalised base foil along the
/// outward surface normal; the trailing edge and leading edge do not move.
pub fn apply_genes(base: &Foil, s: &OptimSettings, genes: &[f64]) -> Result<Foil, FoilError> {
    let nb = s.bumps_per_side;
    if genes.len() != s.dimension() {
        return Err(FoilError::Geometry(format!("{} genes given, {} expected", genes.len(), s.dimension())));
    }
    if genes.iter().any(|g| !g.is_finite()) {
        return Err(FoilError::Geometry("non-finite design variable".into()));
    }
    base.validate()?;
    let b = base.normalized();
    let n = b.n();
    let ile = b.le_index();
    let mut f = b.clone();
    if nb > 0 {
        for i in 1..n - 1 {
            if i == ile {
                continue;
            }
            // outward normal of a counter-clockwise contour
            let tx = b.x[i + 1] - b.x[i - 1];
            let ty = b.y[i + 1] - b.y[i - 1];
            let tl = (tx * tx + ty * ty).sqrt();
            if tl < 1e-15 {
                continue;
            }
            let (nx, ny) = (ty / tl, -tx / tl);
            let off = if i < ile { 0 } else { nb };
            let mut h = 0.0;
            for j in 0..nb {
                let t1 = s.t1 * (j + 1) as f64 / (2 * nb + 1) as f64;
                h += hicks_henne(b.x[i], t1, s.t2) * genes[off + j];
            }
            f.x[i] += nx * h;
            f.y[i] += ny * h;
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
    f.note = "optimised (Hicks-Henne)".into();
    Ok(f)
}

// ---------------------------------------------------------------- evaluation

struct Ctx<'a> {
    base: &'a Foil,
    s: &'a OptimSettings,
    specs: Vec<ObjectiveSpec>,
    /// per-objective normalisation for the scalarised score
    scales: Vec<f64>,
    prog: &'a Progress,
    done: AtomicUsize,
    total: usize,
}

impl Ctx<'_> {
    fn infeasible(&self, genes: Vec<f64>, foil: Foil) -> Candidate {
        let k = self.specs.len();
        Candidate { genes, foil, cl: 0.0, cd: 0.0, cm: 0.0, cpmin: 0.0, objectives: vec![f64::NAN; k], scores: vec![LARGE; k], feasible: false }
    }

    fn evaluate(&self, genes: Vec<f64>) -> Result<Candidate, FoilError> {
        foiley::set_quiet(true);
        self.prog.check()?;
        let r = self.evaluate_inner(genes);
        let c = r?;
        if c.feasible {
            self.prog.count_ok();
        } else {
            self.prog.count_failed();
        }
        let d = self.done.fetch_add(1, Ordering::Relaxed) + 1;
        self.prog.set(d as f64 / self.total.max(1) as f64);
        self.prog.check()?;
        Ok(c)
    }

    fn evaluate_inner(&self, genes: Vec<f64>) -> Result<Candidate, FoilError> {
        let foil = match apply_genes(self.base, self.s, &genes) {
            Ok(f) => f,
            Err(FoilError::Cancelled) => return Err(FoilError::Cancelled),
            Err(_) => return Ok(self.infeasible(genes, self.base.clone())),
        };
        let pt = match analyze_point(&foil, &self.s.analysis, self.s.alpha, true, &Progress::none()) {
            Ok(Some(d)) => d.point,
            Ok(None) => return Ok(self.infeasible(genes, foil)),
            Err(FoilError::Cancelled) => return Err(FoilError::Cancelled),
            Err(_) => return Ok(self.infeasible(genes, foil)),
        };
        let vals: Vec<f64> = self.specs.iter().map(|o| o.objective.value(pt.cl, pt.cd, pt.cm, pt.cpmin)).collect();
        let ok = pt.cd > 0.0 && [pt.cl, pt.cd, pt.cm, pt.cpmin].iter().all(|v| v.is_finite()) && vals.iter().all(|v| v.is_finite());
        if !ok {
            return Ok(self.infeasible(genes, foil));
        }
        let scores = self.specs.iter().zip(&vals).map(|(o, &v)| o.goal.score(v)).collect();
        Ok(Candidate { genes, foil, cl: pt.cl, cd: pt.cd, cm: pt.cm, cpmin: pt.cpmin, objectives: vals, scores, feasible: true })
    }

    fn evaluate_all(&self, genes: Vec<Vec<f64>>) -> Result<Vec<Candidate>, FoilError> {
        genes.into_par_iter().map(|g| self.evaluate(g)).collect()
    }

    /// Scalarised score (lower is better) used to pick the best compromise
    /// and as the GA fitness.
    fn scalar(&self, c: &Candidate) -> f64 {
        if !c.feasible {
            return LARGE;
        }
        c.scores.iter().zip(&self.scales).map(|(s, k)| s / k).sum()
    }

    /// All enabled objectives satisfied (only possible with target goals).
    fn satisfied(&self, c: &Candidate) -> bool {
        c.feasible && self.specs.iter().zip(&c.scores).all(|(o, &sc)| matches!(o.goal, Goal::Target { .. }) && sc == 0.0)
    }
}

/// Insert a candidate into a non-dominated archive.
fn archive_insert(archive: &mut Vec<Candidate>, c: &Candidate) -> bool {
    if !c.feasible {
        return false;
    }
    if archive.iter().any(|a| a.dominates(c) || a.scores == c.scores) {
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

struct Particle {
    pos: Vec<f64>,
    vel: Vec<f64>,
    /// personal bests (position, scores)
    best: Vec<(Vec<f64>, Vec<f64>)>,
    cand: Option<Candidate>,
}

impl Particle {
    fn random(rng: &mut Rng, bounds: &[(f64, f64)], nbest: usize) -> Particle {
        let pos = random_genes(rng, bounds);
        let vel = bounds.iter().map(|&(lo, hi)| {
            let dv = (hi - lo) / 2.0;
            rng.range(-dv / 2.0, dv / 2.0)
        });
        let vel = vel.collect();
        Particle { best: vec![(pos.clone(), Vec::new()); nbest], pos, vel, cand: None }
    }

    /// xflr5 `updateBest`: replace the first personal best that the current
    /// position improves on in any objective.
    fn update_best(&mut self) {
        let Some(c) = &self.cand else { return };
        for b in self.best.iter_mut() {
            if b.1.is_empty() || b.1.iter().zip(&c.scores).any(|(o, n)| o > n) {
                b.0 = self.pos.clone();
                b.1 = c.scores.clone();
                return;
            }
        }
    }
}

/// Run the optimisation.  `on_iteration` is called after the initial
/// population (iteration 0) and after every iteration with a snapshot.
pub fn optimise(base: &Foil, s: &OptimSettings, prog: &Progress, on_iteration: &(dyn Fn(usize, &OptimResult) + Sync)) -> Result<OptimResult, FoilError> {
    base.validate()?;
    s.analysis.validate().map_err(FoilError::Geometry)?;
    let specs = s.enabled();
    if specs.is_empty() {
        return Err(FoilError::Geometry("no optimisation objective is enabled".into()));
    }
    if s.dimension() == 0 {
        return Err(FoilError::Geometry("no design variables (bumps per side = 0 and no flap)".into()));
    }
    let pop = s.population.max(2);
    let total = pop * (s.iterations + 1);
    let bounds = s.bounds();
    let mut rng = Rng::new(s.seed);
    let mut ctx = Ctx { base, s, specs: specs.clone(), scales: vec![1.0; specs.len()], prog, done: AtomicUsize::new(0), total };
    let algo = s.algorithm.label();
    prog.stage(format!("{algo}: evaluating initial population"));
    prog.info(format!("{algo} optimisation of {}: {} variables, {} objective(s), population {pop}, {} iterations", base.name, s.dimension(), specs.len(), s.iterations));

    // initial population: the base foil plus random designs
    let mut init = vec![s.neutral_genes()];
    while init.len() < pop {
        init.push(random_genes(&mut rng, &bounds));
    }
    let cands = ctx.evaluate_all(init.clone())?;
    // normalisation scales from the base design
    let c0 = &cands[0];
    ctx.scales = specs
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
    if !c0.feasible {
        prog.warn(format!("the base foil did not converge at alpha = {}", s.alpha));
    }

    let mut result = OptimResult::default();
    let mut archive: Vec<Candidate> = Vec::new();
    for c in &cands {
        archive_insert(&mut archive, c);
    }
    archive_trim(&mut archive, s.archive_size);

    let report = |it: usize, archive: &[Candidate], extra: Option<&Candidate>, result: &mut OptimResult, ctx: &Ctx| -> bool {
        let mut best: Option<&Candidate> = archive.iter().min_by(|a, b| ctx.scalar(a).total_cmp(&ctx.scalar(b)));
        if let Some(e) = extra {
            if best.map(|b| ctx.scalar(e) < ctx.scalar(b)).unwrap_or(e.feasible) {
                best = Some(e);
            }
        }
        result.best = best.cloned();
        result.pareto = archive.to_vec();
        result.evaluations = ctx.done.load(Ordering::Relaxed);
        let v = result.best.as_ref().map(|b| b.objectives[0]).unwrap_or(f64::NAN);
        result.history.push((it, v));
        ctx.prog.stage(format!("{algo} iteration {it}/{} - best {} {}", s.iterations, specs[0].objective.label(), fmt_val(v)));
        on_iteration(it, result);
        archive.iter().any(|a| ctx.satisfied(a))
    };

    let mut done = report(0, &archive, None, &mut result, &ctx);

    match s.algorithm {
        Algorithm::Mopso => {
            let nbest = specs.len().min(3);
            let mut swarm: Vec<Particle> = init
                .into_iter()
                .zip(cands)
                .map(|(g, c)| {
                    let mut p = Particle::random(&mut rng, &bounds, nbest);
                    p.pos = g;
                    p.best = vec![(p.pos.clone(), Vec::new()); nbest];
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
                        let np = Particle::random(&mut rng, &bounds, nbest);
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
                    p.cand = Some(c);
                    p.update_best();
                    archive_insert(&mut archive, p.cand.as_ref().unwrap());
                }
                archive_trim(&mut archive, s.archive_size);
                done = report(it, &archive, None, &mut result, &ctx);
                it += 1;
            }
        }
        Algorithm::Genetic => {
            let mut popn: Vec<Candidate> = cands;
            let mut elite: Candidate = popn.iter().min_by(|a, b| ctx.scalar(a).total_cmp(&ctx.scalar(b))).cloned().unwrap();
            let mut it = 1;
            while it <= s.iterations && !done {
                prog.check()?;
                // roulette selection on (worst feasible - score)
                let sc: Vec<f64> = popn.iter().map(|c| ctx.scalar(c)).collect();
                let feas: Vec<f64> = sc.iter().cloned().filter(|v| *v < LARGE).collect();
                let worst = feas.iter().cloned().fold(f64::MIN, f64::max);
                let bestv = feas.iter().cloned().fold(f64::MAX, f64::min);
                let margin = 0.05 * (worst - bestv).max(0.0) + 1e-12;
                let w: Vec<f64> = sc.iter().map(|&v| if v < LARGE { worst - v + margin } else { 0.0 }).collect();
                let wsum: f64 = w.iter().sum();
                let mut genes: Vec<Vec<f64>> = (0..pop)
                    .map(|_| {
                        if wsum <= 0.0 {
                            return random_genes(&mut rng, &bounds);
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
                // elitism: the best design so far replaces the worst child
                let (iw, _) = popn.iter().enumerate().max_by(|a, b| ctx.scalar(a.1).total_cmp(&ctx.scalar(b.1))).unwrap();
                let (ib, cb) = popn.iter().enumerate().min_by(|a, b| ctx.scalar(a.1).total_cmp(&ctx.scalar(b.1))).unwrap();
                if ctx.scalar(cb) < ctx.scalar(&elite) {
                    elite = popn[ib].clone();
                } else {
                    popn[iw] = elite.clone();
                }
                for c in &popn {
                    archive_insert(&mut archive, c);
                }
                archive_trim(&mut archive, s.archive_size);
                done = report(it, &archive, Some(&elite), &mut result, &ctx) || ctx.satisfied(&elite);
                it += 1;
            }
        }
    }
    prog.set(1.0);
    if let Some(b) = &result.best {
        let parts: Vec<String> = specs.iter().zip(&b.objectives).map(|(o, v)| format!("{} {}", o.objective.label(), fmt_val(*v))).collect();
        prog.info(format!("{algo} finished after {} evaluations: {}", result.evaluations, parts.join(", ")));
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
}
