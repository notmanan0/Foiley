use foiley_core::Progress;
use foiley_foil::{analyze, AnalysisSettings, Foil, FoilPolar, Sweep};

#[test]
fn settings_round_trip_and_keep_their_legacy_json() {
    let old: AnalysisSettings = serde_json::from_str(r#"{"ncrit":7.0}"#).unwrap();
    assert_eq!(old.ncrit_bot, None);
    assert_eq!((old.wake_length, old.wake_dissipation), (1.0, 0.9));
    // the 6.99 controls stay out of the JSON at their defaults, so settings
    // (and the identities hashed from them) read and write as before
    let json = serde_json::to_string(&AnalysisSettings::default()).unwrap();
    assert!(!json.contains("ncrit_bot") && !json.contains("wake_length") && !json.contains("wake_dissipation"), "{json}");
    let new = AnalysisSettings { ncrit_bot: Some(11.0), wake_length: 2.5, wake_dissipation: 1.1, ..Default::default() };
    assert_eq!(new, serde_json::from_str(&serde_json::to_string(&new).unwrap()).unwrap());
    assert!(new.validate().is_ok());
    assert!(AnalysisSettings { ncrit_bot: Some(f64::NAN), ..new.clone() }.validate().is_err());
    assert!(AnalysisSettings { wake_length: 0.0, ..new.clone() }.validate().is_err());
    assert!(AnalysisSettings { wake_dissipation: f64::NAN, ..new }.validate().is_err());
}

#[test]
fn a_bottom_ncrit_equal_to_the_top_is_no_split() {
    let same = AnalysisSettings { ncrit_bot: Some(9.0), ..Default::default() };
    assert_eq!(same.bottom_ncrit(), None);
    assert_eq!(same.polar_name(), AnalysisSettings::default().polar_name());
    let split = AnalysisSettings { ncrit_bot: Some(11.0), ..Default::default() };
    assert_eq!(split.bottom_ncrit(), Some(11.0));
    assert!(split.polar_name().ends_with("_N9.0_NB11.0"));
    let polar = |st: &AnalysisSettings| FoilPolar { foil: "f".into(), ncrit: 9.0, ncrit_bot: st.bottom_ncrit(), settings: Some(st.clone()), ..Default::default() };
    assert_eq!(polar(&same).identity(), polar(&AnalysisSettings::default()).identity());
    assert_ne!(polar(&split).identity(), polar(&AnalysisSettings::default()).identity());
}

#[test]
fn imports_old_and_new_ncrit_headers() {
    let source = include_str!("../../../tests/data/xfoil699/n2412_split.pol");
    let new = FoilPolar::from_xfoil_text(source, "fallback").unwrap();
    assert_eq!(new.ncrit, 5.0);
    assert_eq!(new.ncrit_bot, Some(11.0));
    // 6.99 CDp is integrated: friction is unknown without a CDf column
    assert!(new.points.iter().all(|p| p.cdf == 0.0));
    let again = FoilPolar::from_xfoil_text(&new.to_xfoil_text(), "fallback").unwrap();
    assert_eq!(again.ncrit_bot, Some(11.0));
    assert_eq!(again.identity(), new.identity());
    let old = FoilPolar::from_xfoil_text(&source.replace("5.000 11.000", "5.000"), "fallback").unwrap();
    assert_eq!(old.ncrit_bot, None);
    assert_ne!(old.identity(), new.identity());
    let equal = FoilPolar::from_xfoil_text(&source.replace("5.000 11.000", "5.000  5.000"), "fallback").unwrap();
    assert_eq!(equal.ncrit_bot, None);
}

#[test]
fn only_the_version_line_marks_a_699_polar() {
    // a 6.97 polar whose name and rows contain "6.99": CDf = CD - CDp still
    let source = include_str!("../../../tests/data/xfoil699/n2412_split.pol")
        .replace("Version 6.99", "Version 6.97")
        .replace("NACA 2412", "NACA 2412 rev 6.99")
        .replace("   8.000 ", "  16.990 ");
    let old = FoilPolar::from_xfoil_text(&source, "fallback").unwrap();
    assert_eq!(old.foil, "NACA 2412 rev 6.99");
    assert_eq!(old.points.last().unwrap().alpha, 16.99);
    assert!(old.points.iter().all(|p| p.cdf > 0.0 && (p.cdf - (p.cd - p.cdp)).abs() < 1e-12));
}

/// The library path (`analyze`, as Studio and the CLI run it) against the
/// MIT 6.99 split-Ncrit polar.  `Foil::naca` keeps 120 seed points per side
/// (6.99 script NACA uses 123), so the repaneled nodes differ slightly: the
/// coefficients still agree to their printed digits, the fractional
/// transition panel index to about 2e-3.
#[test]
fn library_analysis_matches_the_699_split_ncrit_polar() {
    let reference: Vec<Vec<f64>> = include_str!("../../../tests/data/xfoil699/n2412_split.pol")
        .lines()
        .skip_while(|l| !l.trim_start().starts_with("------"))
        .skip(1)
        .map(|l| l.split_whitespace().map(|v| v.parse().unwrap()).collect())
        .collect();
    let st = AnalysisSettings { ncrit: 5.0, ncrit_bot: Some(11.0), parallel: false, from_zero: false, ..Default::default() };
    let (pol, _) = analyze(&Foil::naca("2412").unwrap(), &st, &Sweep::Alpha { from: 0.0, to: 8.0, step: 2.0 }, &Progress::none()).unwrap();
    assert_eq!(pol.points.len(), reference.len());
    for (p, r) in pol.points.iter().zip(&reference) {
        for (name, got, want, tol) in [("CL", p.cl, r[1], 1e-4), ("CD", p.cd, r[2], 1e-5), ("CDp", p.cdp, r[3], 1e-5), ("CM", p.cm, r[4], 1e-4), ("Xtr top", p.xtr_top, r[5], 1e-4), ("Xtr bot", p.xtr_bot, r[6], 1e-4), ("Itr top", p.itr_top, r[7], 3e-3), ("Itr bot", p.itr_bot, r[8], 3e-3)] {
            assert!((got - want).abs() <= tol, "alpha {}: {name} {got} vs 6.99 {want}", r[0]);
        }
    }
}
