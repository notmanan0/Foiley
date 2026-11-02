//! Foiley: a Rust port of Mark Drela's XFoil 6.97 subsonic airfoil
//! analysis code, with multithreading and SIMD-friendly kernels.

pub mod bl;
pub mod blsolv;
pub mod blsys;
pub mod cli;
pub mod extra;
pub mod fmt;
pub mod gdes;
pub mod gdes_geom;
pub mod geom;
pub mod mdes;
pub mod opextra;
pub mod opextra_blplot;
pub mod panel;
pub mod polar;
pub mod qdes;
pub mod simd;
pub mod spline;
pub mod sweep;
pub mod xfoil;

use std::sync::atomic::{AtomicBool, Ordering};

static VERBOSE: AtomicBool = AtomicBool::new(false);

/// Global verbosity for low-level diagnostics (TRCHEK2 failures etc.).
pub fn verbose() -> bool {
    VERBOSE.load(Ordering::Relaxed)
}
pub fn set_verbose(v: bool) {
    VERBOSE.store(v, Ordering::Relaxed);
}

use std::cell::RefCell;
use std::io::{BufWriter, IsTerminal, Stdout, Write};

thread_local! {
    static OUT: RefCell<BufWriter<Stdout>> = RefCell::new(BufWriter::with_capacity(1 << 16, std::io::stdout()));
    static QUIET: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Suppress console output from this thread (library / worker use).
pub fn set_quiet(q: bool) {
    QUIET.with(|c| c.set(q));
}

/// Whether console output is suppressed on this thread.
pub fn is_quiet() -> bool {
    QUIET.with(|c| c.get())
}

/// Silence engine console output until the guard is dropped.
pub struct QuietGuard(bool);
impl Drop for QuietGuard {
    fn drop(&mut self) {
        set_quiet(self.0);
    }
}
pub fn quiet_scope() -> QuietGuard {
    let g = QuietGuard(is_quiet());
    set_quiet(true);
    g
}

#[doc(hidden)]
pub fn out_write(args: std::fmt::Arguments, nl: bool) {
    if QUIET.with(|c| c.get()) {
        return;
    }
    OUT.with(|o| {
        let mut o = o.borrow_mut();
        let _ = o.write_fmt(args);
        if nl {
            let _ = o.write_all(b"\n");
        }
    });
}

/// Flush buffered console output.
pub fn flush_out() {
    OUT.with(|o| {
        let _ = o.borrow_mut().flush();
    });
}

/// Flush only when talking to a human (keeps batch runs fast).
pub fn flush_if_tty() {
    if std::io::stdin().is_terminal() {
        flush_out();
    }
}

/// Buffered println (stdout is block-buffered like gfortran, not line-buffered).
#[macro_export]
macro_rules! outln {
    () => { $crate::out_write(format_args!(""), true) };
    ($($t:tt)*) => { $crate::out_write(format_args!($($t)*), true) };
}
/// Buffered print.
#[macro_export]
macro_rules! out {
    ($($t:tt)*) => { $crate::out_write(format_args!($($t)*), false) };
}

/// Whether to use fine-grained (intra-Newton-iteration) parallelism.
/// Off by default: on hybrid CPUs the fork/join overhead exceeds the gain
/// for these microsecond-scale kernels. Enable with FOILEY_FINE_PAR=1.
pub fn par_fine() -> bool {
    static F: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *F.get_or_init(|| std::env::var("FOILEY_FINE_PAR").map(|v| v == "1").unwrap_or(false))
}
