//! One wake-aware evaluator for every 3D driver (F007).
//!
//! Every driver (T1-T8 speed balances, the glide loop, T6 controls, T7 trim
//! and the stability derivatives) gets its raw loads from
//! [`Evaluator::loads`], so the wake model of the definition is honoured
//! everywhere.  Direct `System::loads` calls in the drivers are forbidden
//! (`crates/plane/tests/f007.rs::no_direct_sys_loads`).
//!
//! A relaxed wake depends on the freestream direction and on the body rates
//! divided by the speed, but not on the speed itself (every velocity in the
//! relaxation scales with it), so solutions are cached per direction and
//! scaled: a speed balance re-uses the wake relaxed at unit speed.

use crate::analysis::{Method, PlanePolarDef, RawLoads, System};
use crate::mesh::Mesh;
use crate::vlm::{FreeWake, WakeSettings};
use crate::PlaneError;
use foiley_core::{Progress, Vec3};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

/// How the trailing wake is modelled.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WakeModel {
    /// straight wake along x, superposition of the six unit solutions
    Rigid,
    /// trailing lines relaxed along the local flow (VLM only)
    Relaxed { tol: f64, max_iter: usize, nseg: usize, relax: f64 },
}

impl WakeModel {
    /// The wake model a definition asks for (`free_wake` = iteration cap).
    pub fn from_def(def: &PlanePolarDef) -> WakeModel {
        if def.free_wake > 0 && matches!(def.method, Method::Vlm1 | Method::Vlm2) {
            WakeModel::Relaxed { tol: def.wake_tol, max_iter: def.free_wake, nseg: def.wake_nseg, relax: def.wake_relax }
        } else {
            WakeModel::Rigid
        }
    }
}

/// Convergence of the wake behind one evaluation.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct WakeReport {
    /// false for a rigid wake (the other fields are then zero)
    pub relaxed: bool,
    /// relaxation steps taken
    pub iterations: usize,
    pub converged: bool,
    /// max wake-node displacement of the last step / reference chord
    pub displacement: f64,
    /// ||Gamma_k - Gamma_(k-1)||inf / ||Gamma_k||inf of the last two solves
    pub circulation: f64,
}

impl WakeReport {
    /// The warning for an unconverged relaxed wake.
    pub fn warning(&self) -> Option<String> {
        (self.relaxed && !self.converged).then(|| {
            format!(
                "relaxed wake not converged after {} iteration(s): node displacement {:.2e} c_ref, circulation change {:.2e}",
                self.iterations, self.displacement, self.circulation
            )
        })
    }
}

/// A relaxed solution at `speed`, valid for any speed with the same omega/speed.
struct Relaxed {
    speed: f64,
    gamma: Vec<f64>,
    wake: FreeWake,
    report: WakeReport,
}

/// Raw loads of a factored system under a wake model.
pub struct Evaluator<'a> {
    pub sys: System<'a>,
    pub wake: WakeModel,
    /// re-relax the wake for every derivative perturbation instead of
    /// freezing it at the reference state (12 relaxations per point)
    pub relaxed_derivatives: bool,
    /// wake length, m
    length: f64,
    c_ref: f64,
    cache: Mutex<HashMap<[u64; 6], Arc<Relaxed>>>,
}

impl<'a> Evaluator<'a> {
    pub fn new(sys: System<'a>, def: &PlanePolarDef, c_ref: f64) -> Evaluator<'a> {
        Evaluator { sys, wake: WakeModel::from_def(def), relaxed_derivatives: def.relaxed_derivatives, length: def.wake_length * c_ref, c_ref, cache: Mutex::new(HashMap::new()) }
    }

    /// A rigid-wake evaluator around `sys`.
    pub fn rigid(sys: System<'a>) -> Evaluator<'a> {
        Evaluator { sys, wake: WakeModel::Rigid, relaxed_derivatives: false, length: 1.0, c_ref: 1.0, cache: Mutex::new(HashMap::new()) }
    }

    pub fn mesh(&self) -> &Mesh {
        self.sys.mesh()
    }

    pub fn is_relaxed(&self) -> bool {
        matches!((self.wake, &self.sys), (WakeModel::Relaxed { .. }, System::Vlm(_)))
    }

    /// Raw loads for air arriving along the unit vector `dir` at `speed`,
    /// with body rates `omega` (rad/s).
    pub fn loads(&self, dir: Vec3, speed: f64, omega: Vec3, prog: &Progress) -> Result<(RawLoads, WakeReport), PlaneError> {
        match (&self.sys, self.is_relaxed()) {
            (System::Vlm(s), true) => {
                let e = self.relaxed(dir, speed, omega, prog)?;
                let k = speed / e.speed;
                let g: Vec<f64> = e.gamma.iter().map(|g| g * k).collect();
                Ok((s.loads(&g, dir * speed, omega), e.report.clone()))
            }
            _ => Ok((self.sys.loads(dir * speed, omega), WakeReport::default())),
        }
    }

    /// [`Evaluator::loads`] for a relative air velocity vector.
    pub fn loads_at(&self, vinf: Vec3, omega: Vec3, prog: &Progress) -> Result<(RawLoads, WakeReport), PlaneError> {
        if !self.is_relaxed() {
            return Ok((self.sys.loads(vinf, omega), WakeReport::default()));
        }
        let v = vinf.norm();
        self.loads(vinf / v, v, omega, prog)
    }

    /// Relaxed wake lines behind the state (empty for a rigid wake), for display.
    pub fn wake_lines(&self, dir: Vec3, speed: f64, omega: Vec3) -> Vec<Vec<[f32; 3]>> {
        if !self.is_relaxed() {
            return vec![];
        }
        let cache = self.cache.lock().unwrap_or_else(|e| e.into_inner());
        match cache.get(&key(dir, speed, omega)) {
            Some(e) => e.wake.lines.iter().map(|l| l.iter().map(|p| [p.x as f32, p.y as f32, p.z as f32]).collect()).collect(),
            None => vec![],
        }
    }

    /// The evaluator stability derivatives use about the reference state
    /// (`dir`, `speed`, no rotation): for a relaxed wake, unless
    /// `relaxed_derivatives` is set, a system with the converged wake
    /// frozen, factored once, whose unit solutions superpose.  `None`
    /// means: use `self`.
    pub fn derivative_evaluator(&self, dir: Vec3, speed: f64, prog: &Progress) -> Result<Option<Evaluator<'a>>, PlaneError> {
        match &self.sys {
            System::Vlm(s) if self.is_relaxed() && !self.relaxed_derivatives => {
                let e = self.relaxed(dir, speed, Vec3::ZERO, prog)?;
                Ok(Some(Evaluator::rigid(System::Vlm(s.frozen(&e.wake, prog)?))))
            }
            _ => Ok(None),
        }
    }

    fn relaxed(&self, dir: Vec3, speed: f64, omega: Vec3, prog: &Progress) -> Result<Arc<Relaxed>, PlaneError> {
        if !(speed > 0.0 && speed.is_finite()) {
            return Err(PlaneError::Input(format!("relaxed wake needs a positive speed (got {speed} m/s)")));
        }
        let k = key(dir, speed, omega);
        if let Some(e) = self.cache.lock().unwrap_or_else(|e| e.into_inner()).get(&k) {
            return Ok(e.clone());
        }
        let (System::Vlm(s), WakeModel::Relaxed { tol, max_iter, nseg, relax }) = (&self.sys, self.wake) else {
            return Err(PlaneError::Input("relaxed wake needs a VLM system".into()));
        };
        let ws = WakeSettings { max_iter, tol, nseg, relax, length: self.length, c_ref: self.c_ref };
        let (gamma, wake, report) = s.solve_free_wake(dir * speed, omega, &ws, prog)?;
        let e = Arc::new(Relaxed { speed, gamma, wake, report });
        self.cache.lock().unwrap_or_else(|e| e.into_inner()).insert(k, e.clone());
        Ok(e)
    }
}

/// Cache key: the direction and omega / speed, bit for bit.
fn key(dir: Vec3, speed: f64, omega: Vec3) -> [u64; 6] {
    let w = omega / speed;
    [dir.x.to_bits(), dir.y.to_bits(), dir.z.to_bits(), w.x.to_bits(), w.y.to_bits(), w.z.to_bits()]
}
