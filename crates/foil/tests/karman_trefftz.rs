//! XFoil's inviscid panel solution against the exact potential flow about
//! Karman-Trefftz sections (conformal map of the circle).

use foiley_core::Progress;
use foiley_foil::analysis::analyze_point;
use foiley_foil::gen::kt::KarmanTrefftz;
use foiley_foil::{AnalysisSettings, Foil, Recipe};

#[test]
fn inviscid_lift_and_pressure_match_the_exact_solution() {
    let st = AnalysisSettings { viscous: false, ..Default::default() };
    for kt in [KarmanTrefftz::default(), KarmanTrefftz { x_center: 0.10, y_center: 0.0, te_angle: 15.0 }, KarmanTrefftz { x_center: 0.05, y_center: 0.04, te_angle: 20.0 }] {
        let foil = Foil::from_recipe(&Recipe::KarmanTrefftz(kt), None).unwrap();
        for alpha in [0.0, 4.0] {
            let exact = kt.exact(alpha).unwrap();
            let p = analyze_point(&foil, &st, alpha, true, &Progress::none()).unwrap().expect("converged").point;
            let err = (p.cl - exact.cl).abs() / exact.cl.abs().max(0.1);
            println!("{kt:?} alpha {alpha}: XFoil CL {:.5}, exact {:.5} ({:.2}%), Cm {:.5} vs {:.5}", p.cl, exact.cl, 100.0 * err, p.cm, exact.cm);
            assert!(err < 0.005, "{kt:?} alpha {alpha}: XFoil {} exact {}", p.cl, exact.cl);
            assert!((p.cm - exact.cm).abs() < 1e-3, "{kt:?} alpha {alpha}: Cm XFoil {} exact {}", p.cm, exact.cm);
        }
    }
}
