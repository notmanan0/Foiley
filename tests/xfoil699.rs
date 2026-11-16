use foiley::cli::{Input, Session};
use foiley::xfoil::{SolverVersion, XFoil};
use std::path::Path;

fn rows(text: &str) -> Vec<Vec<f64>> {
    text.lines().skip_while(|l| !l.trim_start().starts_with("------")).skip(1)
        .filter_map(|l| {
            let row: Option<Vec<f64>> = l.split_whitespace().map(|v| v.parse().ok()).collect();
            row.filter(|r| r.len() >= 9)
        }).collect()
}

#[test]
fn polars_match_mit_699_including_split_ncrit_wake_lag_and_transition_indices() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/xfoil699");
    // XFOIL-compatible filenames are limited to 64 characters. Keep this
    // argument relative even when the sandbox's temp path is very long.
    let scratch = Path::new("work").join(format!("t699-{}", std::process::id()));
    std::fs::create_dir_all(&scratch).unwrap();
    let _quiet = foiley::quiet_scope();
    for name in ["n0012", "n0012_240", "n23012", "n2412_split", "n2412_wake", "n2412_lag"] {
        let output = scratch.join(format!("{name}.pol"));
        let script = std::fs::read_to_string(root.join(format!("{name}.in"))).unwrap().replace("@POL", &output.to_string_lossy());
        // Input scripts use a filename argument, so an absolute path with spaces
        // is kept as a single COMARG. Remove the previous file to avoid prompts.
        if output.exists() { std::fs::remove_file(&output).unwrap(); }
        let mut s = Session::new(Input::from_string(&script));
        s.run(None);
        let expected = rows(&std::fs::read_to_string(root.join(format!("{name}.pol"))).unwrap());
        let actual = rows(&std::fs::read_to_string(&output).unwrap());
        assert_eq!(expected.len(), 5, "{name} reference");
        assert_eq!(actual.len(), expected.len(), "{name} convergence");
        let tolerances = [5e-4, 5e-5, 5e-6, 5e-6, 5e-5, 5e-4, 5e-4, 2e-2, 2e-2];
        for (r, a) in expected.iter().zip(&actual) {
            for (column, (&x, &y)) in r.iter().zip(a).enumerate() {
                assert!((x-y).abs() <= tolerances[column], "{name} alpha={} column {column}: reference {x}, Foiley {y}", r[0]);
            }
        }
    }
}

#[test]
fn commands_preserve_independent_ncrit_and_invalidate_wake_state() {
    let mut s = Session::new(Input::from_string("NACA 2412\nOPER\nVPAR\nNT 5\nNB 11\nWAKE 2.5\nLAG 5.6 1 1.1\n\n\nQUIT\n"));
    s.run(None);
    assert_eq!(s.xf.ncrit(1), 5.0);
    assert_eq!(s.xf.ncrit(2), 11.0);
    assert_eq!(s.xf.waklen, 2.5);
    assert_eq!(s.xf.bl.par.dlcon, 1.1);
    s.xf.lwake = true;
    s.xf.lblini = true;
    s.xf.lipan = true;
    s.xf.lvconv = true;
    assert!(s.xf.set_wake_length(3.0).is_ok());
    assert!(!s.xf.lwake && !s.xf.lblini && !s.xf.lipan && !s.xf.lvconv);
    assert!(s.xf.set_wake_length(f64::NAN).is_err());
    s.input = Input::from_string("OPER\nVPAR\nN 7\n\n\nQUIT\n");
    s.run(None);
    assert_eq!(s.xf.ncrit(1), 7.0);
    assert_eq!(s.xf.ncrit(2), 7.0);
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
        xf.set_wake_length(length).unwrap();
        xf.specal();
        xf.xywake();
        assert_eq!(xf.nw, expected);
        assert!(xf.x[xf.n+xf.nw].is_finite());
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
    let (mut nx, mut ny) = ([0.0;5], [0.0;5]);
    foiley::panel::ncalc(&x, &y, &s, 4, &mut nx, &mut ny);
    assert_eq!((nx[2],ny[2]), (-1.0,0.0));
    assert_eq!((nx[3],ny[3]), (-1.0,0.0));
}

#[test]
fn polar_parser_keeps_699_side_fields_and_legacy_ncrit() {
    use foiley::polar::{polread, Polar, JNC, JTI, JTN};
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
}

#[test]
fn defaults_file_round_trips_both_ncrit_values() {
    let mut source = Session::new(Input::from_string(""));
    source.xf.acrit = 5.0;
    source.xf.acrit_bot = Some(11.0);
    let text = foiley::opextra::wrtdef(&source);
    let path = Path::new("work").join(format!("defaults699-{}.def", std::process::id()));
    std::fs::create_dir_all("work").unwrap();
    std::fs::write(&path, &text).unwrap();
    let mut restored = Session::new(Input::from_string(""));
    foiley::opextra::getdef(&mut restored, path.to_str().unwrap(), false);
    assert_eq!(restored.xf.ncrit(1), 5.0);
    assert_eq!(restored.xf.ncrit(2), 11.0);
    // Single-Ncrit files must not consume the next trip-location record.
    let old = text.lines().map(|line| if line.contains("REtype") {
        " 1 1.0 7.0 | REtype Re/10^6 Ncrit".to_string()
    } else { line.to_string() }).collect::<Vec<_>>().join("\n");
    std::fs::write(&path, old).unwrap();
    foiley::opextra::getdef(&mut restored, path.to_str().unwrap(), false);
    assert_eq!(restored.xf.ncrit(2), 7.0);
    assert_eq!(restored.xf.xstrip[1], 1.0);
}

#[test]
fn exports_have_consistent_columns_and_arc_length_toggle() {
    let base = format!("work/e699-{}", std::process::id());
    std::fs::create_dir_all("work").unwrap();
    let script = format!("NACA 2412\nOPER\nVISC 1e6\nITER 100\nALFA 4\nDUMP {base}.dump\nDMP {base}.dmp\nVPLO\nUE\nDUMP {base}.x\nXPLO\nDUMP {base}.s\n\n\nGDES\nCAMB\nWRTT {base}.thick\n\n\nQUIT\n");
    let _quiet = foiley::quiet_scope();
    let mut s = Session::new(Input::from_string(&script));
    s.run(None);
    assert!(s.xf.lvconv);
    for (suffix, columns) in [("dump",12), ("dmp",24)] {
        let data = std::fs::read_to_string(format!("{base}.{suffix}")).unwrap();
        assert!(data.lines().filter(|l| !l.starts_with('#') && !l.trim().is_empty()).all(|l| l.split_whitespace().count() == columns), "{suffix}");
    }
    let x = std::fs::read_to_string(format!("{base}.x")).unwrap();
    let arc = std::fs::read_to_string(format!("{base}.s")).unwrap();
    assert!(arc.contains("#    s"));
    assert_ne!(x, arc);
    let thickness = std::fs::read_to_string(format!("{base}.thick")).unwrap();
    assert!(thickness.starts_with("Thickness:"));
    assert!(thickness.lines().skip(1).all(|l| l.split_whitespace().count() == 3));
}

#[test]
fn legacy_version_can_be_selected_before_naca_generation() {
    let mut s = Session::new(Input::from_string("VERS 6.97\nNACA 0012\nQUIT\n"));
    s.run(None);
    assert_eq!(s.xf.solver_version, SolverVersion::Xfoil697);
    assert_eq!(s.xf.nb, 239);
}
