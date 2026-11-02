//! Terminal progress display for long jobs (stderr), never interleaving
//! with machine-readable output on stdout.

use foiley_app::{Job, JobManager, JobStatus};
use foiley_core::progress::LogLevel;
use std::io::{IsTerminal, Write};
use std::time::{Duration, Instant};

pub struct Ui {
    pub quiet: bool,
    pub verbose: bool,
    pub json: bool,
    tty: bool,
    last_line: Instant,
}

fn fmt_dur(d: Duration) -> String {
    let s = d.as_secs_f64();
    if s < 60.0 {
        format!("{s:.1}s")
    } else {
        format!("{}m{:02}s", (s / 60.0) as u64, (s as u64) % 60)
    }
}

impl Ui {
    pub fn new(quiet: bool, verbose: bool, json: bool) -> Ui {
        Ui { quiet, verbose, json, tty: std::io::stderr().is_terminal(), last_line: Instant::now() - Duration::from_secs(10) }
    }

    fn clear(&self) {
        if self.tty {
            eprint!("\r\x1b[2K");
        }
    }

    /// Draw progress / flush log lines for a running job.
    pub fn tick(&mut self, job: &Job) {
        let p = &job.progress;
        let min = if self.verbose { LogLevel::Debug } else { LogLevel::Info };
        for l in p.drain_log() {
            if self.quiet && l.level < LogLevel::Warn {
                continue;
            }
            if l.level < min {
                continue;
            }
            self.clear();
            let tag = match l.level {
                LogLevel::Debug => "  ",
                LogLevel::Info => "  ",
                LogLevel::Warn => "! ",
                LogLevel::Error => "x ",
            };
            if self.json {
                eprintln!("{}", serde_json::json!({"event": "log", "level": format!("{:?}", l.level), "t": l.at.as_secs_f64(), "text": l.text}));
            } else {
                eprintln!("{tag}[{:>7}] {}", fmt_dur(l.at), l.text);
            }
        }
        if self.quiet {
            return;
        }
        let f = p.fraction();
        let (ok, bad) = p.counts();
        if self.json {
            if self.last_line.elapsed() > Duration::from_millis(500) {
                self.last_line = Instant::now();
                eprintln!("{}", serde_json::json!({"event": "progress", "job": job.id, "fraction": f, "stage": p.stage_text(), "ok": ok, "failed": bad, "elapsed": job.elapsed().as_secs_f64()}));
            }
            return;
        }
        let eta = p.eta().map(|d| format!(" ETA {}", fmt_dur(d))).unwrap_or_default();
        let counts = if ok + bad > 0 { format!(" | {ok} ok") + &if bad > 0 { format!(", {bad} failed") } else { String::new() } } else { String::new() };
        if self.tty {
            let w = 28;
            let n = (f * w as f64).round() as usize;
            let mut stage = p.stage_text();
            if stage.chars().count() > 60 {
                stage = stage.chars().take(57).collect::<String>() + "...";
            }
            eprint!("\r\x1b[2K[{}{}] {:5.1}% {} {}{}{}", "#".repeat(n.min(w)), ".".repeat(w - n.min(w)), f * 100.0, fmt_dur(job.elapsed()), stage, counts, eta);
            let _ = std::io::stderr().flush();
        } else if self.last_line.elapsed() > Duration::from_secs(2) {
            self.last_line = Instant::now();
            eprintln!("  {:5.1}% {} {}{}{}", f * 100.0, fmt_dur(job.elapsed()), p.stage_text(), counts, eta);
        }
    }

    /// Run a job to completion, showing progress.  Ctrl-C cancels gracefully.
    pub fn wait(&mut self, mgr: &JobManager, job: &std::sync::Arc<Job>) -> Result<(), String> {
        install_ctrlc(job.clone());
        let st = mgr.wait(job.id, |j| self.tick(j));
        self.tick(job);
        if !self.quiet && !self.json {
            self.clear();
        }
        match st {
            JobStatus::Done => {
                if !self.quiet && !self.json {
                    let (ok, bad) = job.progress.counts();
                    eprintln!("  done in {} ({ok} ok{})", fmt_dur(job.elapsed()), if bad > 0 { format!(", {bad} failed") } else { String::new() });
                }
                Ok(())
            }
            JobStatus::Cancelled => Err("cancelled".into()),
            JobStatus::Failed(e) => Err(e),
            JobStatus::Running => unreachable!(),
        }
    }
}

/// Ctrl-C: first press cancels the running job (results so far are kept
/// where possible), second press exits.
fn install_ctrlc(job: std::sync::Arc<Job>) {
    use std::sync::Mutex;
    static CURRENT: Mutex<Option<std::sync::Arc<Job>>> = Mutex::new(None);
    *CURRENT.lock().unwrap() = Some(job);
    static ONCE: std::sync::Once = std::sync::Once::new();
    ONCE.call_once(|| {
        let _ = ctrlc_hook(Box::new(|| {
            static PRESSED: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
            if PRESSED.fetch_add(1, std::sync::atomic::Ordering::SeqCst) > 0 {
                std::process::exit(130);
            }
            if let Some(j) = CURRENT.lock().unwrap().as_ref() {
                eprintln!("\n  cancelling... (press Ctrl-C again to quit immediately)");
                j.cancel();
            }
        }));
    });
}

#[cfg(windows)]
fn ctrlc_hook(f: Box<dyn Fn() + Send + Sync>) -> Result<(), ()> {
    use std::sync::OnceLock;
    static HANDLER: OnceLock<Box<dyn Fn() + Send + Sync>> = OnceLock::new();
    let _ = HANDLER.set(f);
    #[link(name = "kernel32")]
    extern "system" {
        fn SetConsoleCtrlHandler(h: Option<unsafe extern "system" fn(u32) -> i32>, add: i32) -> i32;
    }
    unsafe extern "system" fn h(ty: u32) -> i32 {
        if ty == 0 {
            if let Some(f) = HANDLER.get() {
                f();
            }
            1
        } else {
            0
        }
    }
    let ok = unsafe { SetConsoleCtrlHandler(Some(h), 1) };
    if ok != 0 {
        Ok(())
    } else {
        Err(())
    }
}

#[cfg(not(windows))]
fn ctrlc_hook(_f: Box<dyn Fn() + Send + Sync>) -> Result<(), ()> {
    // Default SIGINT behaviour (terminate) on other platforms.
    Err(())
}
