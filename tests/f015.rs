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

use foiley::rawdump::{compare, compare_aligned, parse, Comparison, Point};
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
    // `exe` is another build of this same test binary (the x86-64-v2 copy):
    // it is always run on the child entry only, never on its whole test list,
    // or every child would run `f015_matrix` and spawn children of its own.
    let f = dir.join("script.in");
    std::fs::write(&f, script).unwrap();
    let bin = exe.map(Path::to_path_buf).unwrap_or_else(|| std::env::current_exe().unwrap());
    let mut c = Command::new(bin);
    c.args(["--exact", "f015_child_entry", "--test-threads=1", "--nocapture"]).env(CHILD, &f).stdin(Stdio::null());
    c.current_dir(dir).env("RAWDUMP", &raw).env("FOILEY_JOBS", "1").stdout(Stdio::null()).stderr(Stdio::null());
    for (k, v) in env {
        c.env(k, v);
    }
    let mut ch = c.spawn().map_err(|e| e.to_string())?;
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

// ------------------------------------------------------------------ F015b

/// One corpus case: geometry ("NACA xxxx" or a file of tests/data/foils),
/// commands between the geometry and OPER (inverse design), OPER commands.
#[derive(Clone, Debug)]
pub struct Case {
    pub name: String,
    pub geom: String,
    pub pre: String,
    pub oper: String,
}

impl Case {
    fn new(name: &str, geom: &str, pre: &str, oper: &str) -> Case {
        Case { name: name.into(), geom: geom.into(), pre: pre.into(), oper: oper.into() }
    }
    /// The full command script (a geometry file is LOADed and re-paneled).
    pub fn script(&self) -> String {
        let g = if self.geom.starts_with("NACA") { self.geom.clone() } else { format!("LOAD {}\nPANE", self.geom) };
        format!("PLOP\nG F\n\n{g}\n{}OPER\n{}\n\nQUIT\n", self.pre, self.oper.replace("@POL", "pol.txt"))
    }
    /// Copy the geometry file (if any) into `dir`.
    fn stage(&self, dir: &Path) {
        if !self.geom.starts_with("NACA") {
            let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/foils").join(&self.geom);
            std::fs::copy(&src, dir.join(&self.geom)).unwrap_or_else(|e| panic!("{}: {e}", src.display()));
        }
    }
}

/// The 8 cases of bench/run_bench.py and 40 generated ones: NACA 4- and
/// 5-digit, e387, dae11, trips, Mach 0-0.6, Re 5e4-6e6, ASEQ/CSEQ, and
/// inverse-design (MDES) round trips.
pub fn corpus() -> Vec<Case> {
    let v = |re: &str, seq: &str| format!("VISC {re}\nITER 100\nPACC\n@POL\n\n{seq}");
    let vm = |re: &str, m: &str, seq: &str| format!("VISC {re}\nMACH {m}\nITER 100\nPACC\n@POL\n\n{seq}");
    let trip = |vpar: &str, re: &str, seq: &str| format!("VPAR\n{vpar}\n\nVISC {re}\nITER 100\nPACC\n@POL\n\n{seq}");
    let inv = |seq: &str| format!("PACC\n@POL\n\n{seq}");
    let c = Case::new;
    vec![
        // bench/run_bench.py
        c("n0012_re1e6", "NACA 0012", "", &v("1e6", "ASEQ 0 16 0.5")),
        c("n2412_re3e5", "NACA 2412", "", &v("3e5", "ASEQ -5 14 0.25")),
        c("n4415_re6e6_m3", "NACA 4415", "", &vm("6e6", "0.3", "ASEQ -4 12 0.5")),
        c("n23012_cseq", "NACA 23012", "", &v("2e6", "CSEQ 0 1.2 0.1")),
        c("n0012_trip", "NACA 0012", "", &trip("XTR 0.1 0.2\nN 7", "5e5", "ASEQ 0 10 1")),
        c("e387_re2e5", "e387.dat", "", &v("2e5", "ASEQ -2 12 0.5")),
        c("dae11_re5e5", "dae11.dat", "", &v("5e5", "ASEQ 0 10 0.5")),
        c("n0012_inv", "NACA 0012", "", &inv("ASEQ -5 15 1")),
        // NACA 4-digit
        c("g01_n0006", "NACA 0006", "", &v("2e5", "ASEQ 0 6 1")),
        c("g02_n0015", "NACA 0015", "", &v("1e6", "ASEQ -2 12 1")),
        c("g03_n0018", "NACA 0018", "", &v("3e6", "ASEQ 0 14 2")),
        c("g04_n4412", "NACA 4412", "", &v("5e5", "ASEQ -4 12 1")),
        c("g05_n6409", "NACA 6409", "", &v("2e5", "ASEQ 0 10 1")),
        c("g06_n2415", "NACA 2415", "", &v("6e6", "ASEQ -2 14 2")),
        c("g07_n1408", "NACA 1408", "", &v("1e5", "ASEQ 0 8 1")),
        c("g08_n4421", "NACA 4421", "", &v("2e6", "ASEQ 0 12 2")),
        c("g09_n2408", "NACA 2408", "", &v("5e4", "ASEQ 2 8 1")),
        c("g10_n0009", "NACA 0009", "", &v("3e5", "ASEQ 0 8 0.5")),
        // NACA 5-digit
        c("g11_n23015", "NACA 23015", "", &v("1e6", "ASEQ -2 12 1")),
        c("g12_n22012", "NACA 22012", "", &v("3e6", "ASEQ 0 10 1")),
        c("g13_n21012", "NACA 21012", "", &v("5e5", "ASEQ 0 8 1")),
        c("g14_n24012", "NACA 24012", "", &v("2e6", "ASEQ 0 12 2")),
        c("g15_n25012_cseq", "NACA 25012", "", &v("1e6", "CSEQ 0.2 1.0 0.1")),
        // e387, dae11
        c("g16_e387_re1e5", "e387.dat", "", &v("1e5", "ASEQ 0 8 1")),
        c("g17_e387_cseq", "e387.dat", "", &v("5e5", "CSEQ 0.2 1.0 0.1")),
        c("g18_dae11_re2e5", "dae11.dat", "", &v("2e5", "ASEQ 0 8 1")),
        c("g19_dae11_re1e6", "dae11.dat", "", &v("1e6", "ASEQ 2 10 1")),
        // trips and Ncrit
        c("g20_n2412_trip", "NACA 2412", "", &trip("XTR 0.05 0.05", "1e6", "ASEQ 0 8 1")),
        c("g21_n0012_trip_n5", "NACA 0012", "", &trip("XTR 0.3 1.0\nN 5", "3e5", "ASEQ 0 8 1")),
        c("g22_n4415_n11", "NACA 4415", "", &trip("N 11", "2e6", "ASEQ 0 10 1")),
        c("g23_e387_trip", "e387.dat", "", &trip("XTR 0.2 0.5", "2e5", "ASEQ 0 8 1")),
        // Mach
        c("g24_n0012_m2", "NACA 0012", "", &vm("3e6", "0.2", "ASEQ 0 8 1")),
        c("g25_n2412_m4", "NACA 2412", "", &vm("3e6", "0.4", "ASEQ 0 6 1")),
        c("g26_n0012_m6", "NACA 0012", "", &vm("6e6", "0.6", "ASEQ 0 3 0.5")),
        c("g27_n0012_inv_m5", "NACA 0012", "", "MACH 0.5\nPACC\n@POL\n\nASEQ 0 4 1"),
        c("g28_n4412_m3_cseq", "NACA 4412", "", &vm("1e6", "0.3", "CSEQ 0.4 1.2 0.1")),
        // Reynolds number range
        c("g29_n0012_re5e4", "NACA 0012", "", &v("5e4", "ASEQ 0 6 1")),
        c("g30_n0012_re6e6", "NACA 0012", "", &v("6e6", "ASEQ 0 14 1")),
        c("g31_n2412_re1e5_cseq", "NACA 2412", "", &v("1e5", "CSEQ 0.2 1.0 0.1")),
        // sequences
        c("g32_n0012_down", "NACA 0012", "", &v("1e6", "ASEQ 10 0 -1")),
        c("g33_n2412_cseq_down", "NACA 2412", "", &v("1e6", "CSEQ 1.2 0.0 -0.1")),
        c("g34_n23012_inv_cseq", "NACA 23012", "", &inv("CSEQ 0 1 0.25")),
        c("g35_n6409_inv", "NACA 6409", "", &inv("ASEQ -4 8 2")),
        c("g36_n0012_dense", "NACA 0012", "", "VISC 1e6\nITER 200\nPACC\n@POL\n\nASEQ 0 16 0.25"),
        // inverse-design round trips: MDES of the current Qspec, then PANE
        c("g37_n2412_mdes", "NACA 2412", "MDES\nEXEC\n\nPANE\n", &v("1e6", "ASEQ 0 6 1")),
        c("g38_n0012_mdes_filt", "NACA 0012", "MDES\nFILT\nEXEC\n\nPANE\n", &v("1e6", "ASEQ 0 6 1")),
        c("g39_n4412_mdes", "NACA 4412", "MDES\nEXEC\n\nPANE\n", &inv("ASEQ 0 8 2")),
        c("g40_dae11_mdes", "dae11.dat", "MDES\nFILT\nEXEC\n\nPANE\n", &v("5e5", "ASEQ 2 8 1")),
    ]
}

/// A Foiley configuration of the matrix.
struct Config {
    label: String,
    exe: Option<PathBuf>,
    env: Vec<(&'static str, &'static str)>,
    /// parallel sweeps replay coefficients only: compare the scalars
    coefficients_only: bool,
}

/// Drop the BL arrays (to compare coefficients only).
fn scalars_only(d: &[Point]) -> Vec<Point> {
    d.iter().map(|p| Point { nbl: [0, 0], scalars: p.scalars, bl: vec![] }).collect()
}

/// The x86-64-v2 build of this test binary (`bench/build_v2.sh`), for the
/// SIMD axis: `FOILEY_V2_TEST_BIN`, else target/x86-64-v2/release/deps/f015-*.
fn v2_binary() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("FOILEY_V2_TEST_BIN") {
        return Some(PathBuf::from(p)).filter(|p| p.exists());
    }
    let deps = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/x86-64-v2/release/deps");
    let mut c: Vec<PathBuf> = std::fs::read_dir(deps).ok()?.filter_map(|e| e.ok().map(|e| e.path())).filter(|p| {
        let n = p.file_name().unwrap().to_string_lossy().to_string();
        n.starts_with("f015-") && (n.ends_with(".exe") || !n.contains('.'))
    }).collect();
    c.sort_by_key(|p| std::fs::metadata(p).and_then(|m| m.modified()).ok());
    c.pop()
}

/// Worst comparison of a configuration over the corpus.
#[derive(Default)]
struct Tally {
    cases: usize,
    identical: usize,
    max_ulp: Vec<(&'static str, u64)>,
    max_scaled: Vec<(&'static str, f64)>,
    structure: Vec<String>,
}

impl Tally {
    fn add(&mut self, case: &str, c: &Comparison) {
        self.cases += 1;
        if c.identical() {
            self.identical += 1;
        }
        if self.max_ulp.is_empty() {
            self.max_ulp = c.max_ulp.clone();
            self.max_scaled = c.max_scaled.clone();
        }
        for (m, (_, u)) in self.max_ulp.iter_mut().zip(&c.max_ulp) {
            m.1 = m.1.max(*u);
        }
        for (m, (_, d)) in self.max_scaled.iter_mut().zip(&c.max_scaled) {
            m.1 = m.1.max(*d);
        }
        self.structure.extend(c.structure.iter().map(|s| format!("{case}: {s}")));
    }
    fn label(&self) -> String {
        if self.cases == 0 {
            return "no cases".into();
        }
        if self.identical == self.cases {
            return format!("bit-identical on all {} cases", self.cases);
        }
        let u: Vec<String> = self.max_ulp.iter().filter(|(_, u)| *u > 0).map(|(n, u)| if *u == u64::MAX { format!("{n} NaN") } else { format!("{n} {u}") }).collect();
        let w = self.max_scaled.iter().fold(("", 0.0f64), |a, b| if b.1 > a.1 { *b } else { a });
        let mut s = format!("bit-identical on {}/{} cases; max ulp {}; largest scaled difference {:.1e} ({})", self.identical, self.cases, if u.is_empty() { "0".into() } else { u.join(", ") }, w.1, w.0);
        if !self.structure.is_empty() {
            s.push_str(&format!("; {} structural difference(s), first: {}", self.structure.len(), self.structure[0]));
        }
        s
    }
}

#[test]
fn f015_matrix() {
    assert!(foiley::rawdump::ENABLED);
    let cases = corpus();
    assert_eq!(cases.len(), 48);
    let seq = [("FOILEY_PIPELINE", "0"), ("FOILEY_FINE_PAR", "0")];
    let mut configs = vec![
        Config { label: "sequential (FOILEY_JOBS=1 FOILEY_PIPELINE=0)".into(), exe: None, env: seq.to_vec(), coefficients_only: false },
        Config { label: "FOILEY_PIPELINE=1".into(), exe: None, env: vec![("FOILEY_PIPELINE", "1")], coefficients_only: false },
        Config { label: "FOILEY_FINE_PAR=1".into(), exe: None, env: vec![("FOILEY_PIPELINE", "0"), ("FOILEY_FINE_PAR", "1")], coefficients_only: false },
        Config { label: "FOILEY_JOBS=4 FOILEY_CHUNK_PIPE=1 (parallel sweep, coefficients)".into(), exe: None, env: vec![("FOILEY_JOBS", "4"), ("FOILEY_CHUNK_PIPE", "1")], coefficients_only: true },
    ];
    let v2 = v2_binary();
    if let Some(b) = &v2 {
        configs.push(Config { label: "SIMD: x86-64-v2 build, FOILEY_SIMD=0 (scalar)".into(), exe: Some(b.clone()), env: vec![("FOILEY_PIPELINE", "0"), ("FOILEY_SIMD", "0")], coefficients_only: false });
        configs.push(Config { label: "SIMD: x86-64-v2 build, runtime AVX2".into(), exe: Some(b.clone()), env: vec![("FOILEY_PIPELINE", "0")], coefficients_only: false });
    }
    let reference = ref_exe();
    // run everything, a few cases at a time
    let workers = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(2).clamp(1, 6);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let results: std::sync::Mutex<Vec<(usize, Option<Result<Vec<Point>, String>>, Vec<Result<Vec<Point>, String>>)>> = std::sync::Mutex::new(Vec::new());
    std::thread::scope(|s| {
        for _ in 0..workers {
            s.spawn(|| loop {
                let k = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                let Some(case) = cases.get(k) else { break };
                let dir = scratch(&format!("m{k:02}"));
                case.stage(&dir);
                let sc = case.script();
                let r = reference.as_ref().ok().map(|e| run_ref(e, &dir, &sc));
                let mut fs = Vec::new();
                for (j, cfg) in configs.iter().enumerate() {
                    let d = dir.join(format!("c{j}"));
                    std::fs::create_dir_all(&d).unwrap();
                    case.stage(&d);
                    fs.push(run_foiley(cfg.exe.as_deref(), &d, &sc, &cfg.env));
                }
                results.lock().unwrap().push((k, r, fs));
            });
        }
    });
    let mut results = results.into_inner().unwrap();
    results.sort_by_key(|r| r.0);
    // against Foiley's sequential run, and against the reference
    let mut vs_seq: Vec<Tally> = configs.iter().map(|_| Tally::default()).collect();
    let mut vs_ref: Vec<Tally> = configs.iter().map(|_| Tally::default()).collect();
    let mut ref_failed = Vec::new();
    for (k, r, fs) in &results {
        let name = &cases[*k].name;
        let base = fs[0].as_ref().unwrap_or_else(|e| panic!("{name}: sequential Foiley run failed: {e}"));
        for (j, f) in fs.iter().enumerate() {
            let f = f.as_ref().unwrap_or_else(|e| panic!("{name} [{}]: {e}", configs[j].label));
            let (a, b) = if configs[j].coefficients_only { (scalars_only(base), scalars_only(f)) } else { (base.clone(), f.clone()) };
            if j > 0 {
                vs_seq[j].add(name, &compare_aligned(&a, &b));
            }
            match r {
                Some(Ok(rd)) => {
                    let rd = if configs[j].coefficients_only { scalars_only(rd) } else { rd.clone() };
                    vs_ref[j].add(name, &compare_aligned(&rd, &b));
                }
                Some(Err(e)) if j == 0 => ref_failed.push(format!("{name}: {e}")),
                _ => {}
            }
        }
    }
    println!("F015 raw-bit matrix, {} cases ({}):", cases.len(), foiley::simd::mode());
    for (j, cfg) in configs.iter().enumerate() {
        if j > 0 {
            println!("  {:<62} vs sequential: {}", cfg.label, vs_seq[j].label());
        }
        match &reference {
            Ok(_) if !ref_failed.is_empty() => {}
            Ok(_) => println!("  {:<62} vs reference (dp): {}", cfg.label, vs_ref[j].label()),
            Err(_) => {}
        }
    }
    if v2.is_none() {
        println!("  SIMD: SKIPPED(build): no x86-64-v2 build of this test (run bench/build_v2.sh)");
    }
    for d in &vs_ref[0].structure {
        println!("  sequential vs reference, structure: {d}");
    }
    if let Some(t) = vs_ref.first().filter(|t| t.cases > 0) {
        println!("  sequential vs reference, scaled differences: {}", t.max_scaled.iter().map(|(n, d)| format!("{n} {d:.1e}")).collect::<Vec<_>>().join(", "));
    }
    match &reference {
        Err(e) => skip(e),
        Ok(_) if !ref_failed.is_empty() => skip(&ref_failed.join("; ")),
        Ok(_) => {}
    }
    // every configuration is labelled; the sequential paths of this build
    // keep the structure of the sequential run exactly
    for (j, cfg) in configs.iter().enumerate().skip(1) {
        assert_eq!(vs_seq[j].cases, cases.len(), "{} unlabelled", cfg.label);
        if !cfg.coefficients_only {
            assert!(vs_seq[j].structure.is_empty(), "{}: {}", cfg.label, vs_seq[j].label());
        }
    }
}
