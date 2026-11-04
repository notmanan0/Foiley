//! F015: raw-bit comparison of Foiley with the instrumented reference XFoil.
//!
//! The reference is the double-precision build in
//! `$FOILEY_REF_ROOT/build-ref-dp`, built by `build_ref.sh` with
//! `bench/xfoil-rawdump.patch`.  Without it a test prints `SKIPPED(ref)`;
//! with `FOILEY_REQUIRE_REF=1` it fails instead.
//!
//! Each Foiley run is a child process (this test binary re-run on
//! `f015_child_entry`), so that the configuration switches, which the
//! engine reads once per process, take effect.

use foiley::rawdump::{compare, parse, Comparison, Point};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const CHILD: &str = "FOILEY_F015_SCRIPT";

/// Child mode: run the command script named by `FOILEY_F015_SCRIPT` through
/// a Foiley session in this process (a no-op in a normal test run).
#[test]
fn f015_child_entry() {
    let Ok(script) = std::env::var(CHILD) else { return };
    let text = std::fs::read_to_string(&script).unwrap();
    foiley::set_quiet(true);
    let mut s = foiley::cli::Session::new(foiley::cli::Input::from_string(&text));
    s.run(None);
    foiley::flush_out();
}

/// `SKIPPED(ref)`, or a failure under `FOILEY_REQUIRE_REF=1`.
fn skip(why: &str) {
    if std::env::var("FOILEY_REQUIRE_REF").is_ok_and(|v| v == "1") {
        panic!("FOILEY_REQUIRE_REF=1: {why}");
    }
    println!("SKIPPED(ref): {why}");
}

fn ref_exe() -> Result<PathBuf, String> {
    let root = std::env::var("FOILEY_REF_ROOT").map(PathBuf::from).unwrap_or_else(|_| Path::new(env!("CARGO_MANIFEST_DIR")).to_path_buf());
    let exe = root.join("build-ref-dp").join(if cfg!(windows) { "xfoil.exe" } else { "xfoil" });
    if exe.exists() {
        Ok(exe)
    } else {
        Err(format!("no reference build at {} (set FOILEY_REF_ROOT, run build_ref.sh)", exe.display()))
    }
}

/// A fresh scratch directory for one run.
fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("foiley-f015-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// XFoil command script: geometry, then OPER commands with the polar file
/// name substituted for `@POL`.
pub fn script(geom: &str, oper: &str) -> String {
    format!("PLOP\nG F\n\n{geom}\nOPER\n{}\n\nQUIT\n", oper.replace("@POL", "pol.txt"))
}

fn read_dump(f: &Path) -> Result<Vec<Point>, String> {
    let t = std::fs::read_to_string(f).map_err(|e| format!("{}: {e}", f.display()))?;
    parse(&t)
}

/// Run the reference on `script` in `dir`; its dump.
fn run_ref(exe: &Path, dir: &Path, script: &str) -> Result<Vec<Point>, String> {
    use std::io::Write;
    let raw = dir.join("ref.raw");
    let mut c = Command::new(exe).current_dir(dir).env("RAWDUMP", &raw).stdin(Stdio::piped()).stdout(Stdio::null()).stderr(Stdio::null()).spawn().map_err(|e| e.to_string())?;
    c.stdin.take().unwrap().write_all(script.as_bytes()).map_err(|e| e.to_string())?;
    c.wait().map_err(|e| e.to_string())?;
    if !raw.exists() {
        return Err(format!("{} wrote no RAWDUMP file: not the instrumented build (rebuild with build_ref.sh)", exe.display()));
    }
    read_dump(&raw)
}

/// Run Foiley (a child of this test binary, or `exe`) on `script` in `dir`
/// with extra environment `env`; its dump.
pub fn run_foiley(exe: Option<&Path>, dir: &Path, script: &str, env: &[(&str, &str)]) -> Result<Vec<Point>, String> {
    let raw = dir.join("foiley.raw");
    let _ = std::fs::remove_file(&raw);
    let mut c = match exe {
        None => {
            let f = dir.join("script.in");
            std::fs::write(&f, script).unwrap();
            let mut c = Command::new(std::env::current_exe().unwrap());
            c.args(["--exact", "f015_child_entry", "--test-threads=1", "--nocapture"]).env(CHILD, &f).stdin(Stdio::null());
            c
        }
        Some(e) => {
            let mut c = Command::new(e);
            c.stdin(Stdio::piped());
            c
        }
    };
    c.current_dir(dir).env("RAWDUMP", &raw).env("FOILEY_JOBS", "1").stdout(Stdio::null()).stderr(Stdio::null());
    for (k, v) in env {
        c.env(k, v);
    }
    let mut ch = c.spawn().map_err(|e| e.to_string())?;
    if let Some(mut i) = ch.stdin.take() {
        use std::io::Write;
        i.write_all(script.as_bytes()).map_err(|e| e.to_string())?;
    }
    let st = ch.wait().map_err(|e| e.to_string())?;
    if !st.success() {
        return Err(format!("Foiley run failed ({st})"));
    }
    read_dump(&raw)
}

#[test]
fn f015_hooks_agree_on_naca0012() {
    assert!(foiley::rawdump::ENABLED, "the root crate's tests are built with the rawdump feature");
    let exe = match ref_exe() {
        Ok(e) => e,
        Err(e) => return skip(&e),
    };
    let s = script("NACA 0012", "VISC 1e6\nITER 100\nPACC\n@POL\n\nASEQ 0 10 1");
    let dir = scratch("naca0012");
    let r = match run_ref(&exe, &dir, &s) {
        Ok(r) => r,
        Err(e) => return skip(&e),
    };
    let f = run_foiley(None, &dir, &s, &[]).unwrap();
    let c: Comparison = compare(&r, &f);
    println!("NACA 0012 Re 1e6 ASEQ 0..10, sequential, vs the double-precision reference: {}", c.summary());
    assert_eq!(r.len(), 11, "reference stored {} points", r.len());
    assert!(r.iter().all(|p| p.nbl[0] > 10 && !p.bl.is_empty()));
    // The same points, BL stations and transition stations...
    assert!(c.structure.is_empty(), "{}", c.summary());
    // ...and, measured on 4 November 2026, the same values to ~1e-10 of each
    // quantity's scale but not the same bits (the inviscid solution already
    // differs in its last ~12 bits).  FOILEY_REQUIRE_BITS=1 demands identity.
    for (n, d) in &c.max_scaled {
        assert!(*d <= MAX_SCALED, "{n}: scaled difference {d:.2e} > {MAX_SCALED:.0e} ({})", c.summary());
    }
    if std::env::var("FOILEY_REQUIRE_BITS").is_ok_and(|v| v == "1") {
        assert!(c.identical(), "FOILEY_REQUIRE_BITS=1: {}", c.summary());
    }
}

/// Largest accepted |Foiley - reference| / max |reference| per quantity.
pub const MAX_SCALED: f64 = 1e-8;

#[test]
fn f015_compare_reports_ulps() {
    // the harness itself: identical dumps, then one CL off by 3 ulp and a
    // missing BL station
    let p = Point { nbl: [3, 3], scalars: [1.0f64.to_bits(); 7], bl: vec![(1, 2, [0.5f64.to_bits(); 4]), (1, 3, [0.25f64.to_bits(); 4])] };
    let a = vec![p.clone(), p.clone()];
    assert!(compare(&a, &a).identical());
    assert_eq!(compare(&a, &a).summary(), "bit-identical (2 points)");
    let mut b = a.clone();
    b[1].scalars[1] += 3;
    b[0].bl[1].2[2] -= 1;
    let c = compare(&a, &b);
    assert!(!c.identical());
    assert_eq!(c.max_ulp.iter().find(|(n, _)| *n == "CL").unwrap().1, 3);
    assert_eq!(c.max_ulp.iter().find(|(n, _)| *n == "DSTR").unwrap().1, 1);
    assert!(c.summary().starts_with("max ulp: CL 3, DSTR 1 (largest scaled difference"), "{}", c.summary());
    b[1].bl.pop();
    assert!(compare(&a, &b).summary().contains("BL stations"));
    assert_eq!(foiley::rawdump::ulp_distance(0.0f64.to_bits(), (-0.0f64).to_bits()), 1);
    assert_eq!(foiley::rawdump::ulp_distance(1.0f64.to_bits(), f64::NAN.to_bits()), u64::MAX);
    // and the text format round-trips the reference's layout
    let txt = "POINT    81   103\nALFA 0000000000000000\nCL 3E03E56503400000\nCD 3F7622EC12FDA8D7\nCDP 3F5298D6A5AC4E98\nCM BDE05EC3BB184000\nXTR1 3FE5FC16310F0C93\nXTR2 3FE5FC16E8CE74C0\nBL     1     2 3FB318FAF20D1643 3F00D2C5201DEEAB 3F12C0FD14337C91 0000000000000000\n";
    let d = parse(txt).unwrap();
    assert_eq!(d.len(), 1);
    assert_eq!(d[0].nbl, [81, 103]);
    assert_eq!(d[0].scalars[1], 0x3E03E56503400000);
    assert_eq!(d[0].bl[0], (1, 2, [0x3FB318FAF20D1643, 0x3F00D2C5201DEEAB, 0x3F12C0FD14337C91, 0]));
}

