//! F017: 2D points report the engine's Newton iteration count and how they
//! were recovered.
use foiley_core::Progress;
use foiley_foil::analysis::analyze_point;
use foiley_foil::{analyze, AnalysisSettings, Foil, Recovery, Sweep};

#[test]
fn f017_iterations_nonzero() {
    let f = Foil::naca("0012").unwrap();
    let st = AnalysisSettings { re: 1e6, ..Default::default() };
    // single point
    let d = analyze_point(&f, &st, 3.0, true, &Progress::none()).unwrap().expect("converged");
    assert!(d.point.iterations > 0, "single point reports 0 iterations");
    assert!(d.point.iterations as usize <= st.max_iter * 2);
    // sweep path (parallel chunks)
    let (pol, _) = analyze(&f, &st, &Sweep::Alpha { from: -4.0, to: 10.0, step: 0.5 }, &Progress::none()).unwrap();
    assert!(!pol.points.is_empty());
    for p in &pol.points {
        assert!(p.iterations > 0, "alpha {} reports 0 iterations", p.alpha);
        assert_eq!(p.recovered, p.recovery != Recovery::None);
    }
}

#[test]
fn f017_inviscid_iterations_zero() {
    let f = Foil::naca("0012").unwrap();
    let st = AnalysisSettings { viscous: false, ..Default::default() };
    let (pol, _) = analyze(&f, &st, &Sweep::Alpha { from: 0.0, to: 4.0, step: 2.0 }, &Progress::none()).unwrap();
    assert!(pol.points.iter().all(|p| p.iterations == 0 && p.recovery == Recovery::None));
}

#[test]
fn f017_recovered_points_say_how() {
    // the hard sweep of the unit tests needs recovery near stall
    let f = Foil::naca("4412").unwrap();
    let st = AnalysisSettings { re: 2e5, ..Default::default() };
    let (pol, _) = analyze(&f, &st, &Sweep::Alpha { from: -6.0, to: 18.0, step: 0.5 }, &Progress::none()).unwrap();
    for p in pol.points.iter().filter(|p| p.recovered) {
        assert_ne!(p.recovery, Recovery::None, "alpha {}", p.alpha);
        if let Recovery::March { substeps, attempts } = p.recovery {
            assert!(substeps >= 2 && attempts >= 1);
        }
    }
}
