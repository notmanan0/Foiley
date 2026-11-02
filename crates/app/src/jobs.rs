//! Background job manager.
//!
//! Design rules (each one addresses a failure seen in xflr5/flow5):
//! * jobs own immutable copies of their inputs, so editing or deleting a
//!   foil/plane while it is being analysed is always safe;
//! * every job ends in exactly one terminal state (Done, Failed, Cancelled),
//!   even if the solver panics, so no dialog or caller can wait forever;
//! * cancellation is per job (no global flags);
//! * compute runs on a dedicated pool with one core left for the UI, and
//!   nothing here ever blocks the caller: the UI polls `jobs()` each frame;
//! * partial results (finished points / polars) are published as they come.

use foiley_core::Progress;
use foiley_foil::analysis::{analyze_point, PointDetail};
use foiley_foil::batch::{run_batch, BatchCase};
use foiley_foil::{analyze, AnalysisSettings, Foil, FoilPolar, Sweep};
use foiley_plane::mesh::{build, Mesh, MeshOptions};
use foiley_plane::{run_polar, FoilDb, Plane, PlaneOpp, PlanePolar, PlanePolarDef};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

pub type JobId = u64;

/// What to compute.  Inputs are owned copies.
pub enum JobKind {
    FoilPolar { foil: Foil, settings: AnalysisSettings, sweep: Sweep },
    FoilPoint { foil: Foil, settings: AnalysisSettings, value: f64, is_alpha: bool },
    Batch { foils: Vec<Foil>, cases: Vec<BatchCase>, settings: AnalysisSettings, sweep: Sweep },
    PlanePolar { plane: Plane, db: FoilDb, def: PlanePolarDef },
    Mesh { plane: Plane, db: FoilDb, options: MeshOptions },
    Optimise { foil: Foil, settings: foiley_foil::optim::OptimSettings },
    Inverse { foil: Foil, alpha: f64, mode: foiley_foil::inverse::InverseMode, paneling: foiley_foil::foil::Paneling, target: Option<(Vec<f64>, Option<(usize, usize)>)>, iterations: usize },
}

impl JobKind {
    pub fn kind_label(&self) -> &'static str {
        match self {
            JobKind::FoilPolar { .. } => "2D polar",
            JobKind::FoilPoint { .. } => "2D point",
            JobKind::Batch { .. } => "2D batch",
            JobKind::PlanePolar { .. } => "3D polar",
            JobKind::Mesh { .. } => "mesh",
            JobKind::Optimise { .. } => "optimiser",
            JobKind::Inverse { .. } => "inverse design",
        }
    }
}

pub enum JobOutput {
    FoilPolar { polar: FoilPolar, details: Vec<PointDetail> },
    FoilPoint(Option<PointDetail>),
    Batch(Vec<FoilPolar>),
    PlanePolar(PlanePolar),
    Mesh(Mesh),
    Optimise(foiley_foil::optim::OptimResult),
    /// (new foil if a target was given, speed distribution)
    Inverse(Option<Foil>, foiley_foil::inverse::SpeedDist),
}

#[derive(Clone, Debug, PartialEq)]
pub enum JobStatus {
    Running,
    Done,
    Failed(String),
    Cancelled,
}

impl JobStatus {
    pub fn is_finished(&self) -> bool {
        !matches!(self, JobStatus::Running)
    }
}

/// Results published while a job runs.
#[derive(Default)]
pub struct Partial {
    pub plane_points: Vec<PlaneOpp>,
    pub foil_polars: Vec<FoilPolar>,
    /// optimiser state after each iteration
    pub optim: Option<foiley_foil::optim::OptimResult>,
}

pub struct Job {
    pub id: JobId,
    pub title: String,
    pub kind: &'static str,
    pub progress: Progress,
    status: Mutex<JobStatus>,
    output: Mutex<Option<JobOutput>>,
    pub partial: Mutex<Partial>,
    started: Instant,
    finished: Mutex<Option<Duration>>,
}

impl Job {
    pub fn status(&self) -> JobStatus {
        self.status.lock().unwrap().clone()
    }
    pub fn elapsed(&self) -> Duration {
        self.finished.lock().unwrap().unwrap_or_else(|| self.started.elapsed())
    }
    pub fn cancel(&self) {
        self.progress.cancel();
    }
    /// Take the final output (once).
    pub fn take_output(&self) -> Option<JobOutput> {
        self.output.lock().unwrap().take()
    }
    pub fn has_output(&self) -> bool {
        self.output.lock().unwrap().is_some()
    }
}

type Notify = Arc<dyn Fn() + Send + Sync>;

pub struct JobManager {
    pool: Arc<rayon::ThreadPool>,
    jobs: Vec<Arc<Job>>,
    next: JobId,
    notify: Notify,
    pub threads: usize,
}

impl JobManager {
    /// `threads` = compute threads (0 = all cores but one).
    pub fn new(threads: usize) -> Self {
        let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
        let n = if threads == 0 { cores.saturating_sub(1).max(1) } else { threads };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .thread_name(|i| format!("foiley-worker-{i}"))
            .build()
            .expect("thread pool");
        JobManager { pool: Arc::new(pool), jobs: Vec::new(), next: 1, notify: Arc::new(|| {}), threads: n }
    }

    /// Called (from worker threads) whenever a job changes state or publishes
    /// partial results; the GUI uses it to request a repaint.
    pub fn set_notify(&mut self, f: impl Fn() + Send + Sync + 'static) {
        self.notify = Arc::new(f);
    }

    pub fn jobs(&self) -> &[Arc<Job>] {
        &self.jobs
    }
    pub fn get(&self, id: JobId) -> Option<Arc<Job>> {
        self.jobs.iter().find(|j| j.id == id).cloned()
    }
    pub fn running(&self) -> usize {
        self.jobs.iter().filter(|j| !j.status().is_finished()).count()
    }
    pub fn cancel_all(&self) {
        for j in &self.jobs {
            j.cancel();
        }
    }
    /// Remove finished jobs whose output was collected.
    pub fn clear_finished(&mut self) {
        self.jobs.retain(|j| !j.status().is_finished() || j.has_output());
    }

    /// Start a job; returns immediately.
    pub fn submit(&mut self, title: impl Into<String>, kind: JobKind) -> Arc<Job> {
        let progress = Progress::new();
        let job = Arc::new(Job {
            id: self.next,
            title: title.into(),
            kind: kind.kind_label(),
            progress,
            status: Mutex::new(JobStatus::Running),
            output: Mutex::new(None),
            partial: Mutex::new(Partial::default()),
            started: Instant::now(),
            finished: Mutex::new(None),
        });
        self.next += 1;
        self.jobs.push(job.clone());
        let pool = self.pool.clone();
        let notify = self.notify.clone();
        let j = job.clone();
        std::thread::Builder::new()
            .name(format!("foiley-job-{}", job.id))
            .spawn(move || {
                let res = catch_unwind(AssertUnwindSafe(|| pool.install(|| execute(&j, kind, &notify))));
                let status = match res {
                    Ok(Ok(out)) => {
                        *j.output.lock().unwrap() = Some(out);
                        if j.progress.is_cancelled() {
                            JobStatus::Cancelled
                        } else {
                            JobStatus::Done
                        }
                    }
                    Ok(Err(e)) if j.progress.is_cancelled() => {
                        let _ = e;
                        JobStatus::Cancelled
                    }
                    Ok(Err(e)) => JobStatus::Failed(e),
                    Err(p) => {
                        let msg = p.downcast_ref::<&str>().map(|s| s.to_string()).or_else(|| p.downcast_ref::<String>().cloned()).unwrap_or_else(|| "unknown".into());
                        JobStatus::Failed(format!("internal solver error ({msg}); the rest of the application is unaffected"))
                    }
                };
                match &status {
                    JobStatus::Failed(e) => j.progress.error(e.clone()),
                    JobStatus::Cancelled => j.progress.warn("cancelled"),
                    _ => {}
                }
                *j.finished.lock().unwrap() = Some(j.started.elapsed());
                *j.status.lock().unwrap() = status;
                notify();
            })
            .expect("spawn job thread");
        job
    }

    /// Block until a job finishes (CLI use only; never call from a UI thread).
    pub fn wait(&self, id: JobId, mut tick: impl FnMut(&Job)) -> JobStatus {
        let Some(j) = self.get(id) else { return JobStatus::Failed("no such job".into()) };
        loop {
            tick(&j);
            let s = j.status();
            if s.is_finished() {
                return s;
            }
            std::thread::sleep(Duration::from_millis(100));
        }
    }
}

fn execute(job: &Job, kind: JobKind, notify: &Notify) -> Result<JobOutput, String> {
    let p = &job.progress;
    match kind {
        JobKind::FoilPolar { foil, settings, sweep } => {
            let (polar, details) = analyze(&foil, &settings, &sweep, p).map_err(|e| e.to_string())?;
            p.info(format!("{}: {} points, {} not converged", polar.name, polar.points.len(), polar.failed.len()));
            Ok(JobOutput::FoilPolar { polar, details })
        }
        JobKind::FoilPoint { foil, settings, value, is_alpha } => {
            let d = analyze_point(&foil, &settings, value, is_alpha, p).map_err(|e| e.to_string())?;
            if d.is_none() {
                p.warn(format!("{}: point did not converge", foil.name));
            }
            Ok(JobOutput::FoilPoint(d))
        }
        JobKind::Batch { foils, cases, settings, sweep } => {
            let res = run_batch(&foils, &cases, &settings, &sweep, p, &|r| {
                if let Ok(pol) = &r.polar {
                    job.partial.lock().unwrap().foil_polars.push(pol.clone());
                    notify();
                }
            });
            p.check().map_err(|e| e.to_string())?;
            let mut out = Vec::new();
            let mut errs = 0;
            for r in res {
                match r.polar {
                    Ok(pol) => out.push(pol),
                    Err(_) => errs += 1,
                }
            }
            if errs > 0 {
                p.warn(format!("{errs} polar(s) failed"));
            }
            Ok(JobOutput::Batch(out))
        }
        JobKind::PlanePolar { plane, db, def } => {
            let pol = run_polar(&plane, &db, &def, p, &|o| {
                job.partial.lock().unwrap().plane_points.push(o.clone());
                notify();
            })
            .map_err(|e| e.to_string())?;
            p.info(format!("{}: {} points ({} unknowns), {} failed", pol.name, pol.points.len(), pol.unknowns, pol.failed.len()));
            Ok(JobOutput::PlanePolar(pol))
        }
        JobKind::Optimise { foil, settings } => {
            let r = foiley_foil::optim::optimise(&foil, &settings, p, &|_, res| {
                job.partial.lock().unwrap().optim = Some(res.clone());
                notify();
            })
            .map_err(|e| e.to_string())?;
            p.info(format!("optimiser: {} evaluations, {} Pareto designs", r.evaluations, r.pareto.len()));
            Ok(JobOutput::Optimise(r))
        }
        JobKind::Inverse { foil, alpha, mode, paneling, target, iterations } => {
            p.stage("Inverse design");
            let out = match target {
                None => (None, foiley_foil::inverse::speed_distribution(&foil, alpha, mode, &paneling).map_err(|e| e.to_string())?),
                Some((t, seg)) => {
                    let (f, d) = foiley_foil::inverse::solve(&foil, alpha, mode, &paneling, &t, seg, iterations).map_err(|e| e.to_string())?;
                    (Some(f), d)
                }
            };
            p.set(1.0);
            Ok(JobOutput::Inverse(out.0, out.1))
        }
        JobKind::Mesh { plane, db, options } => {
            p.stage("Meshing");
            let m = build(&plane, &db, &options).map_err(|e| e.to_string())?;
            p.set(1.0);
            Ok(JobOutput::Mesh(m))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn jobs_finish_cancel_and_survive_errors() {
        let mut m = JobManager::new(2);
        let f = Foil::naca("0012").unwrap();
        let a = m.submit("ok", JobKind::FoilPolar { foil: f.clone(), settings: AnalysisSettings::default(), sweep: Sweep::Alpha { from: 0.0, to: 4.0, step: 1.0 } });
        let b = m.submit("bad", JobKind::FoilPolar { foil: f.clone(), settings: AnalysisSettings { re: -1.0, ..Default::default() }, sweep: Sweep::Alpha { from: 0.0, to: 4.0, step: 1.0 } });
        let c = m.submit("cancel", JobKind::FoilPolar { foil: f, settings: AnalysisSettings::default(), sweep: Sweep::Alpha { from: -10.0, to: 15.0, step: 0.05 } });
        c.cancel();
        assert_eq!(m.wait(a.id, |_| {}), JobStatus::Done);
        assert!(matches!(m.wait(b.id, |_| {}), JobStatus::Failed(_)));
        assert_eq!(m.wait(c.id, |_| {}), JobStatus::Cancelled);
        assert!(matches!(a.take_output(), Some(JobOutput::FoilPolar { .. })));
    }

    #[test]
    fn running_batch_job_cancels_promptly() {
        let mut m = JobManager::new(0);
        let foils = vec![Foil::naca("4412").unwrap()];
        let cases = foiley_foil::batch::re_list_geometric(5e4, 1e6, 40).into_iter().flat_map(|re| [7.0, 8.0, 9.0].map(|ncrit| BatchCase { re, mach: 0.0, ncrit })).collect();
        let j = m.submit("batch", JobKind::Batch { foils, cases, settings: AnalysisSettings::default(), sweep: Sweep::Alpha { from: -10.0, to: 12.0, step: 0.5 } });
        std::thread::sleep(Duration::from_millis(1500));
        assert_eq!(j.status(), JobStatus::Running);
        let t = Instant::now();
        let before = j.progress.fraction();
        j.cancel();
        assert_eq!(m.wait(j.id, |_| {}), JobStatus::Cancelled);
        let after = j.progress.fraction();
        println!("batch job stopped {:.2} s after cancel, progress {:.0}% -> {:.0}%", t.elapsed().as_secs_f64(), before * 100.0, after * 100.0);
        // cancellation lands between points; unoptimised builds solve ~10x slower
        let limit = Duration::from_secs(if cfg!(debug_assertions) { 15 } else { 2 });
        assert!(t.elapsed() < limit, "took {:?}", t.elapsed());
        // polars skipped by the cancel must not count as done
        assert!(after < 0.9, "progress jumped to {after}");
    }
}
