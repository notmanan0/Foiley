//! The NACA equations (`NacaSection::surfaces_at` / `half_thickness`)
//! against the ordinates printed by NASA Langley's original NACA4 and NACA6
//! programs (tests/data/naca, 5 decimals).  This checks the analytic
//! construction; the codes `Foil::naca` hands to XFoil's own routine
//! (4-digit, 210xx-250xx, such as 23012 here) are checked bit for bit
//! against that routine in `gen::naca::tests`.

use foiley_foil::gen::naca;
use std::path::Path;

fn table(name: &str) -> Vec<Vec<f64>> {
    let p = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/naca").join(format!("{name}.txt"));
    std::fs::read_to_string(&p)
        .unwrap_or_else(|e| panic!("{}: {e}", p.display()))
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(|v| v.parse().unwrap()).collect())
        .collect()
}

/// Largest difference between our ordinates and NASA's within the nose
/// (x/c < 0.0125, where NASA's own spline of the mapping and of the
/// logarithmic mean-line slope differ most) and aft of it.
fn compare(code: &str) -> (f64, f64) {
    let sec = naca::parse(code).unwrap();
    let rows = table(code);
    let x: Vec<f64> = rows.iter().map(|r| r[0]).collect();
    let (mut nose, mut rest): (f64, f64) = (0.0, 0.0);
    let mut note = |x: f64, e: f64| {
        assert!(e.is_finite(), "NACA {code}: non-finite ordinate at x/c {x}");
        if x < 0.0125 {
            nose = nose.max(e);
        } else {
            rest = rest.max(e);
        }
    };
    if rows[0].len() == 2 {
        let yt = sec.half_thickness(&x).unwrap();
        for (r, y) in rows.iter().zip(&yt) {
            note(r[0], (r[1] - y).abs());
        }
    } else {
        let (up, lo) = sec.surfaces_at(&x).unwrap();
        for (i, r) in rows.iter().enumerate() {
            for (a, b) in [(r[1], up[i].0), (r[2], up[i].1), (r[3], lo[i].0), (r[4], lo[i].1)] {
                note(r[0], (a - b).abs());
            }
        }
    }
    (nose, rest)
}

#[test]
fn naca4_program_sections() {
    // 5-digit (standard and reflexed), 4-digit modified and 16-series.
    // NASA's 4-digit-modified runs were given the leading-edge radius
    // rounded to 3 figures (0.0110 for 0010-64; Abbott's formula gives
    // 0.011019), which accounts for up to ~2e-5 near the nose.
    for code in ["23012", "23112", "0010-64", "16-009"] {
        let (nose, rest) = compare(code);
        println!("NACA {code}: max difference {nose:.1e} at the nose, {rest:.1e} aft");
        assert!(nose <= 2e-5 && rest <= 2e-5, "NACA {code}: {nose:.1e} {rest:.1e}");
    }
}

#[test]
fn naca6_program_sections() {
    // 6- and 6A-series thickness, a=1.0 and modified a=0.8 mean lines
    for code in ["63-006", "63A010", "63-206", "63A409"] {
        let (nose, rest) = compare(code);
        println!("NACA {code}: max difference {nose:.1e} at the nose, {rest:.1e} aft");
        assert!(nose <= 2.5e-4 && rest <= 2e-5, "NACA {code}: {nose:.1e} {rest:.1e}");
    }
}
