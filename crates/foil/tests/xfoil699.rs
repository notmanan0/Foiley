use foiley_foil::analysis::SolverVersion;
use foiley_foil::{AnalysisSettings, Foil, FoilPolar};

#[test]
fn old_settings_reproduce_697_and_new_settings_round_trip_699() {
    let old: AnalysisSettings = serde_json::from_str(r#"{"ncrit":7.0}"#).unwrap();
    assert_eq!(old.solver_version, SolverVersion::Xfoil697);
    assert_eq!(old.ncrit_bot, None);
    assert_eq!(old.wake_length, 1.0);
    assert_eq!(old.wake_dissipation, 0.9);
    let new = AnalysisSettings { ncrit_bot: Some(11.0), wake_length: 2.5, ..Default::default() };
    assert_eq!(new.solver_version, SolverVersion::Xfoil699);
    assert_eq!(new, serde_json::from_str(&serde_json::to_string(&new).unwrap()).unwrap());
    assert!(new.validate().is_ok());
    assert!(AnalysisSettings { ncrit_bot: Some(f64::NAN), ..new.clone() }.validate().is_err());
    assert!(AnalysisSettings { wake_length: 0.0, ..new }.validate().is_err());
}

#[test]
fn imports_old_and_new_ncrit_headers() {
    let source = include_str!("../../../tests/data/xfoil699/n2412_split.pol");
    let new = FoilPolar::from_xfoil_text(source, "fallback").unwrap();
    assert_eq!(new.ncrit, 5.0);
    assert_eq!(new.ncrit_bot, Some(11.0));
    let output = new.to_xfoil_text();
    let again = FoilPolar::from_xfoil_text(&output, "fallback").unwrap();
    assert_eq!(again.ncrit_bot, Some(11.0));
    let old = FoilPolar::from_xfoil_text(&source.replace("5.000 11.000", "5.000"), "fallback").unwrap();
    assert_eq!(old.ncrit_bot, None);
    assert_ne!(old.identity(), new.identity());
}

#[test]
fn normal_offset_changes_thickness_and_keeps_a_valid_foil() {
    let foil = Foil::naca("0012").unwrap();
    let offset = foil.with_normal_offset(0.001).unwrap();
    offset.validate().unwrap();
    assert!(offset.props().thickness > foil.props().thickness + 0.001);
    assert!(foil.with_normal_offset(f64::NAN).is_err());
}
