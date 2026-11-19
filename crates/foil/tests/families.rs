//! One section of every generated family through a viscous XFoil polar.

use foiley_core::Progress;
use foiley_foil::analysis::analyze_point;
use foiley_foil::gen::{cst::Cst, kt::KarmanTrefftz, parsec::Parsec, plate::Plate};
use foiley_foil::{AnalysisSettings, Foil, Recipe};

#[test]
fn every_family_converges_in_a_viscous_analysis() {
    let mut foils: Vec<Foil> = ["23112", "43012", "2412 closed TE", "0012-64", "16-212", "63-206", "64(2)-415", "65-415 a=0.5", "64A010", "63A210"].iter().map(|c| Foil::naca(c).unwrap()).collect();
    for r in [Recipe::Cst(Cst::default()), Recipe::Parsec(Parsec::default()), Recipe::KarmanTrefftz(KarmanTrefftz::default())] {
        foils.push(Foil::from_recipe(&r, None).unwrap());
    }
    let st = AnalysisSettings { re: 1.0e6, ..Default::default() };
    for f in &foils {
        assert!(f.recipe.is_some(), "{}", f.name);
        for alpha in [0.0, 4.0] {
            let p = analyze_point(f, &st, alpha, true, &Progress::none()).unwrap().unwrap_or_else(|| panic!("{} did not converge at {alpha}", f.name)).point;
            assert!(p.cd > 0.002 && p.cd < 0.03, "{} alpha {alpha}: CD {}", f.name, p.cd);
            println!("{:<22} alpha {alpha}: CL {:.4} CD {:.5}", f.name, p.cl, p.cd);
        }
    }
    // a thin cambered plate at low Reynolds number
    let plate = Foil::from_recipe(&Recipe::Plate(Plate::default()), None).unwrap();
    let low = AnalysisSettings { re: 1.0e5, ..Default::default() };
    let p = analyze_point(&plate, &low, 3.0, true, &Progress::none()).unwrap().expect("plate converged").point;
    println!("{:<22} alpha 3: CL {:.4} CD {:.5}", plate.name, p.cl, p.cd);
    assert!(p.cl > 0.4, "plate CL {}", p.cl);
}
