//! Batch 2D analysis: polar meshes over foils x (Re, Mach, Ncrit) cases.
//!
//! Every (foil, case) pair is an independent job run on the rayon pool; each
//! job reports into its own slice of the progress range and is isolated from
//! the others (a failure is reported, the rest carry on).

use crate::analysis::{analyze, AnalysisSettings, Sweep};
use crate::foil::Foil;
use crate::polar::FoilPolar;
use crate::FoilError;
use foiley_core::Progress;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// One row of a batch: the varying parameters.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct BatchCase {
    pub re: f64,
    pub mach: f64,
    pub ncrit: f64,
}

/// Re values spaced geometrically (typical for polar meshes feeding 3D analyses).
/// Values are rounded to 1000, so dense lists drop the duplicates that creates.
pub fn re_list_geometric(from: f64, to: f64, count: usize) -> Vec<f64> {
    if count <= 1 || from <= 0.0 || to <= 0.0 {
        return vec![from];
    }
    let r = (to / from).powf(1.0 / (count - 1) as f64);
    let mut v: Vec<f64> = (0..count).map(|i| (from * r.powi(i as i32) / 1000.0).round() * 1000.0).collect();
    v.dedup();
    v
}

#[derive(Clone, Debug)]
pub struct BatchResult {
    pub foil: String,
    pub case: BatchCase,
    pub polar: Result<FoilPolar, FoilError>,
}

/// Run all foils x cases.  `base` provides everything except Re/Mach/Ncrit.
/// `on_done` is called (from worker threads) as each polar completes, so a
/// caller can stream results into a UI.
pub fn run_batch(
    foils: &[Foil],
    cases: &[BatchCase],
    base: &AnalysisSettings,
    sweep: &Sweep,
    prog: &Progress,
    on_done: &(dyn Fn(&BatchResult) + Sync),
) -> Vec<BatchResult> {
    let jobs: Vec<(usize, usize)> = (0..foils.len()).flat_map(|i| (0..cases.len()).map(move |j| (i, j))).collect();
    let n = jobs.len();
    prog.stage(format!("Batch: {} foils x {} cases = {} polars", foils.len(), cases.len(), n));
    let res: Vec<BatchResult> = jobs
        .par_iter()
        .enumerate()
        .map(|(k, &(i, j))| {
            foiley::set_quiet(true);
            let f = &foils[i];
            let c = cases[j];
            let st = AnalysisSettings { re: c.re, mach: c.mach, ncrit: c.ncrit, parallel: false, ..base.clone() };
            let sub = prog.part(k, n).prefixed(&format!("[{} Re={:.0}] ", f.name, c.re));
            let polar = if prog.is_cancelled() { Err(FoilError::Cancelled) } else { analyze(f, &st, sweep, &sub).map(|r| r.0) };
            match &polar {
                Ok(p) => sub.info(format!("done: {} points, {} not converged", p.points.len(), p.failed.len())),
                Err(FoilError::Cancelled) => {}
                Err(e) => sub.error(format!("failed: {e}")),
            }
            // skipped polars must not advance the bar, or a cancelled batch reads 100%
            if !matches!(polar, Err(FoilError::Cancelled)) {
                sub.set(1.0);
            }
            let r = BatchResult { foil: f.name.clone(), case: c, polar };
            on_done(&r);
            r
        })
        .collect();
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mesh_runs() {
        let foils = vec![Foil::naca("0012").unwrap(), Foil::naca("2412").unwrap()];
        let cases: Vec<BatchCase> = re_list_geometric(2e5, 1e6, 3).into_iter().map(|re| BatchCase { re, mach: 0.0, ncrit: 9.0 }).collect();
        let count = std::sync::atomic::AtomicUsize::new(0);
        let p = Progress::new();
        let r = run_batch(&foils, &cases, &AnalysisSettings::default(), &Sweep::Alpha { from: -2.0, to: 8.0, step: 1.0 }, &p, &|_| {
            count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        });
        assert_eq!(r.len(), 6);
        assert_eq!(count.load(std::sync::atomic::Ordering::Relaxed), 6);
        for b in &r {
            match &b.polar {
                Ok(p) => assert!(p.points.len() >= 10, "{} {:?}: {} pts, failed {:?}", b.foil, b.case, p.points.len(), p.failed),
                Err(e) => panic!("{e}"),
            }
        }
        assert!((p.fraction() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn dense_re_lists_have_no_duplicates() {
        let v = re_list_geometric(5e4, 1e6, 500);
        assert!(v.len() > 40 && v.len() < 500);
        assert!(v.windows(2).all(|w| w[1] > w[0]));
        assert_eq!((v[0], *v.last().unwrap()), (5e4, 1e6));
        assert_eq!(re_list_geometric(5e4, 1e6, 40).len(), 40);
    }

    #[test]
    fn low_re_batch_cancels_promptly() {
        let foils = vec![Foil::naca("4412").unwrap()];
        let cases: Vec<BatchCase> = re_list_geometric(5e4, 1e6, 40).into_iter().flat_map(|re| [7.0, 8.0, 9.0].map(|ncrit| BatchCase { re, mach: 0.0, ncrit })).collect();
        let p = Progress::new();
        let started = std::time::Instant::now();
        std::thread::scope(|s| {
            s.spawn(|| {
                std::thread::sleep(std::time::Duration::from_millis(300));
                p.cancel();
            });
            run_batch(&foils, &cases, &AnalysisSettings::default(), &Sweep::Alpha { from: -10.0, to: 12.0, step: 0.5 }, &p, &|_| {});
        });
        let t = started.elapsed();
        println!("batch returned {:.2} s after start (cancelled at 0.3 s)", t.as_secs_f64());
        // cancellation lands between points; unoptimised builds solve ~10x slower,
        // and shared CI runners take 2-3 s for what takes 0.5 s here.
        let limit = if cfg!(debug_assertions) { 15.0 } else if std::env::var_os("CI").is_some() { 6.0 } else { 1.5 };
        assert!(t.as_secs_f64() < limit, "cancellation took {t:?}");
    }
}
