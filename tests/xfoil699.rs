use foiley::cli::{Input, Session};
use foiley::xfoil::{SolverVersion, XFoil};
use std::path::{Path, PathBuf};

fn rows(text: &str) -> Vec<Vec<f64>> {
    text.lines()
        .skip_while(|l| !l.trim_start().starts_with("------"))
        .skip(1)
        .filter_map(|l| {
            let row: Option<Vec<f64>> = l.split_whitespace().map(|v| v.parse().ok()).collect();
            row.filter(|r| r.len() >= 9)
        })
        .collect()
}

/// A per-test directory under `work/` (XFoil filenames are limited to 64
/// characters, so paths stay relative), removed when dropped.
struct Scratch(PathBuf);
impl Scratch {
    fn new(tag: &str) -> Scratch {
        let dir = Path::new("work").join(format!("{tag}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        Scratch(dir)
    }
    fn path(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn polars_match_mit_699_including_split_ncrit_wake_lag_and_transition_indices() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/xfoil699");
    let scratch = Scratch::new("t699");
    let _quiet = foiley::quiet_scope();
    for name in ["n0012", "n0012_240", "n23012", "n2412_split", "n2412_wake", "n2412_lag"] {
        let output = scratch.path(&format!("{name}.pol"));
        let script = std::fs::read_to_string(root.join(format!("{name}.in"))).unwrap().replace("@POL", &output);
        let mut s = Session::new(Input::from_string(&script));
        s.run(None);
        let expected = rows(&std::fs::read_to_string(root.join(format!("{name}.pol"))).unwrap());
        let actual = rows(&std::fs::read_to_string(&output).unwrap());
        assert_eq!(expected.len(), 5, "{name} reference");
        assert_eq!(actual.len(), expected.len(), "{name} convergence");
        // every column as printed: half a unit in its last digit
        let tolerances = [5e-4, 5e-5, 5e-6, 5e-6, 5e-5, 5e-5, 5e-5, 5e-5, 5e-5];
        for (r, a) in expected.iter().zip(&actual) {
            for (column, (&x, &y)) in r.iter().zip(a).enumerate() {
                assert!((x - y).abs() <= tolerances[column], "{name} alpha={} column {column}: reference {x}, Foiley {y}", r[0]);
            }
        }
    }
}

#[test]
fn vpar_controls_follow_xfoil_699() {
    let mut s = Session::new(Input::from_string("NACA 2412\nOPER\nVPAR\nNT 5\nNB 11\nWAKE 2.5\nLAG 5.6 1 1.1\n\n\nQUIT\n"));
    s.run(None);
    assert_eq!((s.xf.ncrit(1), s.xf.ncrit(2)), (5.0, 11.0));
    assert_eq!(s.xf.waklen, 2.5);
    assert_eq!(s.xf.bl.par.dlcon, 1.1);
    // N sets both surfaces; NT then changes only the top
    s.input = Input::from_string("OPER\nVPAR\nN 7\nNT 3\n\n\nQUIT\n");
    s.run(None);
    assert_eq!((s.xf.ncrit(1), s.xf.ncrit(2)), (3.0, 7.0));
    // values XFoil accepts are not range-checked
    s.input = Input::from_string("OPER\nVPAR\nN 0.05\nWAKE 0.5\n\n\nQUIT\n");
    s.run(None);
    assert_eq!((s.xf.ncrit(1), s.xf.ncrit(2)), (0.05, 0.05));
    assert_eq!(s.xf.waklen, 0.5);
}

#[test]
fn wake_length_invalidates_the_wake_and_bl_state() {
    let mut xf = XFoil::new();
    (xf.lwake, xf.lblini, xf.lipan, xf.lvconv, xf.lwdij) = (true, true, true, true, true);
    xf.set_wake_length(3.0);
    assert!(!xf.lwake && !xf.lblini && !xf.lipan && !xf.lvconv && !xf.lwdij);
}

#[test]
fn wake_count_follows_each_versions_rule_and_array_cap() {
    for (version, length, panels, expected) in [
        (SolverVersion::Xfoil697, 1.0, 160, 22),
        (SolverVersion::Xfoil699, 1.0, 160, 23),
        (SolverVersion::Xfoil699, 1.0, 240, 30),
        (SolverVersion::Xfoil699, 2.5, 160, 33),
        (SolverVersion::Xfoil699, 20.0, 160, foiley::xfoil::IWX),
    ] {
        let mut xf = XFoil::new();
        xf.solver_version = version;
        xf.npan = panels;
        xf.naca(12);
        xf.set_wake_length(length);
        xf.specal();
        xf.xywake();
        assert_eq!(xf.nw, expected);
        assert!(xf.x[xf.n + xf.nw].is_finite());
    }
}

#[test]
fn naca_seed_count_follows_the_solver_version() {
    for (version, nb) in [(SolverVersion::Xfoil697, 239), (SolverVersion::Xfoil699, 245)] {
        let mut xf = XFoil::new();
        xf.solver_version = version;
        xf.naca(12);
        assert_eq!(xf.nb, nb);
    }
}

#[test]
fn repeated_endpoint_input_collapses_without_underflow() {
    let mut n = 7;
    let mut x = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
    let mut y = [0.0, 10.0, 11.0, 12.0, 20.0, 30.0, 31.0, 32.0];
    foiley::gdes_geom::fixdup(&mut n, &mut x, &mut y);
    assert_eq!(&x[1..=n], &[0.0, 0.5, 1.0]);
    assert_eq!(&y[1..=n], &[12.0, 20.0, 32.0]);
    n = 3;
    x[1..=3].fill(0.0);
    foiley::gdes_geom::fixdup(&mut n, &mut x, &mut y);
    assert_eq!(n, 1);
}

#[test]
fn opposite_corner_normals_get_a_finite_fallback() {
    let x = [0.0, 1.0, 0.3, 0.3, 1.0];
    let y = [0.0; 5];
    let s = [0.0, 0.0, 0.7, 0.7, 1.4];
    let (mut nx, mut ny) = ([0.0; 5], [0.0; 5]);
    foiley::panel::ncalc(&x, &y, &s, 4, &mut nx, &mut ny);
    assert_eq!((nx[2], ny[2]), (-1.0, 0.0));
    assert_eq!((nx[3], ny[3]), (-1.0, 0.0));
}

#[test]
fn polar_parser_keeps_699_side_fields_and_legacy_ncrit() {
    use foiley::polar::{polread, Polar, JNC, JTI, JTN, JTP};
    let source = include_str!("data/xfoil699/n2412_split.pol");
    let mut polar = Polar::default();
    assert!(!polread(source, &mut polar));
    assert_eq!(polar.na(), 5);
    assert_eq!(polar.acrit, 5.0);
    assert_eq!(polar.acrit_bot, Some(11.0));
    assert_eq!(polar.cpolsd[0][JNC][2], 11.0);
    assert_eq!(polar.cpolsd[0][JTN][1], 0.5435);
    assert_eq!(polar.cpolsd[0][JTI][2], 145.0519);
    let mut legacy = Polar::default();
    assert!(!polread(&source.replace("5.000 11.000", "5.000"), &mut legacy));
    assert_eq!(legacy.cpolsd[0][JNC][2], 5.0);
    // trip columns next to transition columns: Xtr must not match Xtrip
    let head = &source[..source.find("   alpha").unwrap()];
    let both = format!(
        "{head}   alpha    CL        CD       CDp       CM    Top_Xtrip Bot_Xtrip  Top_Xtr  Bot_Xtr  Top_Itr  Bot_Itr\n  ------ -------- --------- --------- -------- --------- --------- -------- -------- -------- --------\n   0.000   0.2265   0.00587   0.00053  -0.0499    0.2000    0.3000   0.1900   0.2900  30.1734 145.0519\n"
    );
    let mut tripped = Polar::default();
    assert!(!polread(&both, &mut tripped));
    assert_eq!(tripped.na(), 1);
    assert_eq!((tripped.cpolsd[0][JTP][1], tripped.cpolsd[0][JTP][2]), (0.2, 0.3));
    assert_eq!((tripped.cpolsd[0][JTN][1], tripped.cpolsd[0][JTN][2]), (0.19, 0.29));
    assert_eq!((tripped.cpolsd[0][JTI][1], tripped.cpolsd[0][JTI][2]), (30.1734, 145.0519));
}

#[test]
fn defaults_file_round_trips_both_ncrit_values() {
    let scratch = Scratch::new("def699");
    let path = scratch.path("defaults.def");
    let mut source = Session::new(Input::from_string(""));
    source.xf.acrit = 5.0;
    source.xf.acrit_bot = Some(11.0);
    let text = foiley::opextra::wrtdef(&source);
    std::fs::write(&path, &text).unwrap();
    let mut restored = Session::new(Input::from_string(""));
    foiley::opextra::getdef(&mut restored, &path, false);
    assert_eq!(restored.xf.ncrit(1), 5.0);
    assert_eq!(restored.xf.ncrit(2), 11.0);
    // Single-Ncrit files must not consume the next trip-location record.
    let old = text
        .lines()
        .map(|line| if line.contains("REtype") { " 1 1.0 7.0 | REtype Re/10^6 Ncrit".to_string() } else { line.to_string() })
        .collect::<Vec<_>>()
        .join("\n");
    std::fs::write(&path, old).unwrap();
    foiley::opextra::getdef(&mut restored, &path, false);
    assert_eq!(restored.xf.ncrit(2), 7.0);
    assert_eq!(restored.xf.xstrip[1], 1.0);
}

#[test]
fn exports_use_the_699_layouts() {
    let scratch = Scratch::new("e699");
    let base = scratch.path("e");
    let script = format!(
        "NACA 2412\nOPER\nVISC 1e6\nITER 100\nALFA 4\nDUMP {base}.dump\nDMP {base}.dmp\nVPLO\nUE\nDUMP {base}.x\nXPLO\nDUMP {base}.s\n\n\nGDES\nCAMB\nWRTT {base}.thick\nWRTC {base}.camb\n\n\nQUIT\n"
    );
    let _quiet = foiley::quiet_scope();
    let mut s = Session::new(Input::from_string(&script));
    s.run(None);
    assert!(s.xf.lvconv);
    let data = |suffix: &str| std::fs::read_to_string(format!("{base}.{suffix}")).unwrap();
    let columns = |text: &str| -> Vec<usize> { text.lines().filter(|l| !l.starts_with('#') && !l.trim().is_empty()).map(|l| l.split_whitespace().count()).collect() };
    // DUMP: H*, P, m and K on the surface rows only, as BLDUMP writes them
    let dump = columns(&data("dump"));
    assert_eq!(dump.iter().filter(|&&c| c == 12).count(), s.xf.n);
    assert_eq!(dump.iter().filter(|&&c| c == 8).count(), s.xf.nw);
    let first = data("dump").lines().nth(1).unwrap().trim_end_matches('\r').to_string();
    assert_eq!(first.len(), 1 + 4 * 9 + 3 * 10 + 10 + 10 + 3 * 9, "FORMAT(1X, 4F9.5, 3F10.6, F10.4, F10.4, 3F9.5)");
    // DMP: 24 fixed-format columns, a blank line after the top side and after the bottom surface
    let dmp = data("dmp");
    assert!(columns(&dmp).iter().all(|&c| c == 24));
    assert_eq!(dmp.lines().filter(|l| l.trim().is_empty()).count(), 2);
    // XPLO: arc length runs on from the top trailing edge into the wake
    let top_s: Vec<f64> = data("s")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .take_while(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().next().unwrap().parse().unwrap())
        .collect();
    assert!(top_s.len() > s.xf.iblte[1]);
    // (the first wake node sits at the trailing edge: equal s there)
    assert!(top_s.windows(2).all(|w| w[1] >= w[0]), "top-side s must not fall back at the wake");
    assert_ne!(data("x"), data("s"));
    // WRTT/WRTC: FORMAT(1X, 3G15.7)
    let thickness = data("thick");
    assert!(thickness.starts_with("Thickness:"));
    assert!(thickness.lines().skip(1).all(|l| l.len() == 1 + 3 * 15));
    assert!(data("camb").lines().skip(1).all(|l| l.len() == 1 + 2 * 15));
}
