//! Progress reporting and cooperative cancellation.
//!
//! Every long computation takes a [`Progress`].  It is a cheap clonable
//! handle around shared atomics: the worker writes the completed fraction,
//! a short stage description and log lines; the UI thread reads them at its
//! own pace and may request cancellation at any time.  Workers poll
//! [`Progress::check`] between units of work (operating points, matrix row
//! blocks, iterations), so cancellation is honoured within milliseconds and
//! the UI thread never waits on a solver.
//!
//! Sub-tasks get a [`Progress::sub`] view that maps their own 0..1 range
//! onto a window of the parent's range.

use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum LogLevel {
    Debug,
    Info,
    Warn,
    Error,
}

#[derive(Clone, Debug)]
pub struct LogEntry {
    pub level: LogLevel,
    pub at: Duration,
    pub text: String,
}

/// Returned by [`Progress::check`] when the user asked to stop.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Cancelled;

impl std::fmt::Display for Cancelled {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("cancelled by user")
    }
}
impl std::error::Error for Cancelled {}

struct Inner {
    cancel: AtomicBool,
    frac: AtomicU64,
    stage: Mutex<String>,
    log: Mutex<Vec<LogEntry>>,
    echo: AtomicBool,
    min_echo: Mutex<LogLevel>,
    ok: AtomicUsize,
    failed: AtomicUsize,
    started: Instant,
}

/// Thread-safe progress + cancellation handle.  `Progress::none()` is a no-op
/// sink for library callers that do not care.
#[derive(Clone)]
pub struct Progress {
    inner: Option<Arc<Inner>>,
    /// completed fraction of this view's own window (shared by clones)
    local: Arc<AtomicU64>,
    span: f64,
    prefix: Arc<str>,
}

fn atomic_add_f64(a: &AtomicU64, d: f64) {
    let mut cur = a.load(Ordering::Relaxed);
    loop {
        let new = (f64::from_bits(cur) + d).to_bits();
        match a.compare_exchange_weak(cur, new, Ordering::Relaxed, Ordering::Relaxed) {
            Ok(_) => return,
            Err(c) => cur = c,
        }
    }
}

impl Default for Progress {
    fn default() -> Self {
        Progress::none()
    }
}

impl std::fmt::Debug for Progress {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Progress({:.1}%)", self.fraction() * 100.0)
    }
}

impl Progress {
    pub fn new() -> Self {
        Progress {
            inner: Some(Arc::new(Inner {
                cancel: AtomicBool::new(false),
                frac: AtomicU64::new(0f64.to_bits()),
                stage: Mutex::new(String::new()),
                log: Mutex::new(Vec::new()),
                echo: AtomicBool::new(false),
                min_echo: Mutex::new(LogLevel::Info),
                ok: AtomicUsize::new(0),
                failed: AtomicUsize::new(0),
                started: Instant::now(),
            })),
            local: Arc::new(AtomicU64::new(0f64.to_bits())),
            span: 1.0,
            prefix: Arc::from(""),
        }
    }

    /// A handle that records nothing and can never be cancelled.
    pub fn none() -> Self {
        Progress { inner: None, local: Arc::new(AtomicU64::new(0)), span: 1.0, prefix: Arc::from("") }
    }

    /// Echo log lines at or above `min` to stderr as they arrive (CLI use).
    pub fn with_echo(self, on: bool, min: LogLevel) -> Self {
        if let Some(i) = &self.inner {
            i.echo.store(on, Ordering::Relaxed);
            *i.min_echo.lock().unwrap() = min;
        }
        self
    }

    /// A view for a sub-task worth the fraction `b - a` of this view's
    /// window.  Sub-views may run concurrently: each one adds only its own
    /// increments to the total, so parallel parts aggregate correctly.
    /// (`a` and `b` only define the size of the share.)
    pub fn sub(&self, a: f64, b: f64) -> Progress {
        Progress { inner: self.inner.clone(), local: Arc::new(AtomicU64::new(0f64.to_bits())), span: self.span * (b - a).max(0.0), prefix: self.prefix.clone() }
    }

    /// Same as [`sub`](Self::sub) for the `i`-th of `n` equal parts.
    pub fn part(&self, i: usize, n: usize) -> Progress {
        let n = n.max(1) as f64;
        self.sub(i as f64 / n, (i + 1) as f64 / n)
    }

    /// A view whose log lines and stage text are prefixed (e.g. with the polar name).
    pub fn prefixed(&self, p: &str) -> Progress {
        let mut s = self.clone();
        s.prefix = Arc::from(format!("{}{}", self.prefix, p));
        s
    }

    /// Set completion of this view's window (0..1).  Monotonic: a value
    /// lower than the previous one is ignored.
    pub fn set(&self, f: f64) {
        if let Some(i) = &self.inner {
            let f = f.clamp(0.0, 1.0);
            let mut cur = self.local.load(Ordering::Relaxed);
            loop {
                let old = f64::from_bits(cur);
                if f <= old {
                    return;
                }
                match self.local.compare_exchange_weak(cur, f.to_bits(), Ordering::Relaxed, Ordering::Relaxed) {
                    Ok(_) => {
                        atomic_add_f64(&i.frac, (f - old) * self.span);
                        return;
                    }
                    Err(c) => cur = c,
                }
            }
        }
    }

    /// Overall completed fraction (0..1).
    pub fn fraction(&self) -> f64 {
        self.inner.as_ref().map(|i| f64::from_bits(i.frac.load(Ordering::Relaxed)).clamp(0.0, 1.0)).unwrap_or(0.0)
    }

    /// Describe what is being done right now ("Building influence matrix").
    pub fn stage(&self, s: impl Into<String>) {
        if let Some(i) = &self.inner {
            let s: String = s.into();
            let s = if self.prefix.is_empty() { s } else { format!("{}{}", self.prefix, s) };
            *i.stage.lock().unwrap() = s;
        }
    }

    pub fn stage_text(&self) -> String {
        self.inner.as_ref().map(|i| i.stage.lock().unwrap().clone()).unwrap_or_default()
    }

    pub fn log(&self, level: LogLevel, text: impl Into<String>) {
        if let Some(i) = &self.inner {
            let text: String = text.into();
            let text = if self.prefix.is_empty() { text } else { format!("{}{}", self.prefix, text) };
            let at = i.started.elapsed();
            if i.echo.load(Ordering::Relaxed) && level >= *i.min_echo.lock().unwrap() {
                let tag = match level {
                    LogLevel::Debug => "debug",
                    LogLevel::Info => "info",
                    LogLevel::Warn => "warn",
                    LogLevel::Error => "error",
                };
                eprintln!("[{:7.2}s {tag:5}] {text}", at.as_secs_f64());
            }
            i.log.lock().unwrap().push(LogEntry { level, at, text });
        }
    }
    pub fn info(&self, t: impl Into<String>) {
        self.log(LogLevel::Info, t)
    }
    pub fn warn(&self, t: impl Into<String>) {
        self.log(LogLevel::Warn, t)
    }
    pub fn error(&self, t: impl Into<String>) {
        self.log(LogLevel::Error, t)
    }
    pub fn debug(&self, t: impl Into<String>) {
        self.log(LogLevel::Debug, t)
    }

    /// Take all log lines recorded since the last call.
    pub fn drain_log(&self) -> Vec<LogEntry> {
        self.inner.as_ref().map(|i| std::mem::take(&mut *i.log.lock().unwrap())).unwrap_or_default()
    }

    /// Count a successfully computed unit (operating point, polar...).
    pub fn count_ok(&self) {
        if let Some(i) = &self.inner {
            i.ok.fetch_add(1, Ordering::Relaxed);
        }
    }
    /// Count a unit that failed (did not converge, invalid input...).
    pub fn count_failed(&self) {
        if let Some(i) = &self.inner {
            i.failed.fetch_add(1, Ordering::Relaxed);
        }
    }
    /// (successful, failed) unit counts.
    pub fn counts(&self) -> (usize, usize) {
        self.inner.as_ref().map(|i| (i.ok.load(Ordering::Relaxed), i.failed.load(Ordering::Relaxed))).unwrap_or((0, 0))
    }

    pub fn cancel(&self) {
        if let Some(i) = &self.inner {
            i.cancel.store(true, Ordering::Relaxed);
        }
    }
    pub fn is_cancelled(&self) -> bool {
        self.inner.as_ref().map(|i| i.cancel.load(Ordering::Relaxed)).unwrap_or(false)
    }
    /// `Err(Cancelled)` once cancellation was requested; use with `?`.
    pub fn check(&self) -> Result<(), Cancelled> {
        if self.is_cancelled() {
            Err(Cancelled)
        } else {
            Ok(())
        }
    }
    pub fn elapsed(&self) -> Duration {
        self.inner.as_ref().map(|i| i.started.elapsed()).unwrap_or_default()
    }
    /// Rough remaining-time estimate from the completed fraction.
    pub fn eta(&self) -> Option<Duration> {
        let f = self.fraction();
        if !(0.02..1.0).contains(&f) {
            return None;
        }
        let e = self.elapsed().as_secs_f64();
        Some(Duration::from_secs_f64(e * (1.0 - f) / f))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sub_ranges_map() {
        let p = Progress::new();
        let s = p.sub(0.5, 1.0).part(1, 2);
        s.set(0.5);
        assert!((p.fraction() - 0.125).abs() < 1e-12);
        let a = p.part(0, 2);
        a.set(1.0);
        s.set(1.0);
        a.set(0.3); // ignored: monotonic
        assert!((p.fraction() - 0.75).abs() < 1e-12);
        p.cancel();
        assert!(s.check().is_err());
        Progress::none().set(1.0);
        assert!(Progress::none().check().is_ok());
    }
}
