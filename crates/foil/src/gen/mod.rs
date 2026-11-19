//! Airfoil generators.
//!
//! * [`naca`]: the NACA 4-digit, 5-digit (standard, reflexed, any design
//!   CL), 4-digit modified, 16-, 6- and 6A-series families.
//! * [`cst`]: Kulfan class-shape transformation, generated or fitted to
//!   any foil.
//! * [`parsec`]: Sobieczky's 11-parameter PARSEC sections.
//! * [`kt`]: Karman-Trefftz sections, with their exact inviscid solution.
//! * [`plate`]: cambered plates of constant thickness with round edges.
//!
//! Each generator returns a [`Foil`] whose [`Recipe`] records how it was
//! made, so an editor can reopen the parameters, and rejects parameters
//! whose contour would cross itself.

pub mod cst;
pub mod kt;
pub mod naca;
mod naca6_data;
pub mod parsec;
pub mod plate;
mod spline;

use crate::{Foil, FoilError};
use serde::{Deserialize, Serialize};

/// Upper and lower surface points, each from the leading to the trailing edge.
pub type Surfaces = (Vec<(f64, f64)>, Vec<(f64, f64)>);

/// Points per surface, leading to trailing edge (as `Foil::naca` has
/// always used: 239 points in all).
pub const NSIDE: usize = 120;

/// How a foil's shape was generated.  Stored with the foil, and cleared by
/// any operation that changes its coordinates (including repaneling, since
/// regenerating would not reproduce the new points).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "family", rename_all = "snake_case")]
pub enum Recipe {
    /// A NACA designation, as accepted by [`naca::parse`].
    Naca { designation: String },
    Cst(cst::Cst),
    Parsec(parsec::Parsec),
    KarmanTrefftz(kt::KarmanTrefftz),
    Plate(plate::Plate),
}

impl Recipe {
    /// Short description of the family ("NACA", "CST order 8", ...).
    pub fn family_label(&self) -> String {
        match self {
            Recipe::Naca { .. } => "NACA".into(),
            Recipe::Cst(c) => format!("CST order {}", c.order()),
            Recipe::Parsec(_) => "PARSEC".into(),
            Recipe::KarmanTrefftz(_) => "Karman-Trefftz".into(),
            Recipe::Plate(_) => "cambered plate".into(),
        }
    }

    /// Generate the foil, named `name` (the NACA designation's own name if
    /// `None`).
    pub fn generate(&self, name: Option<&str>) -> Result<Foil, FoilError> {
        // every family's generator attaches its own recipe
        let mut f = match self {
            Recipe::Naca { designation } => naca::parse(designation)?.generate()?,
            Recipe::Cst(c) => c.generate(name.unwrap_or("CST"))?,
            Recipe::Parsec(p) => p.generate(name.unwrap_or("PARSEC"))?,
            Recipe::KarmanTrefftz(k) => k.generate(name.unwrap_or("Karman-Trefftz"))?,
            Recipe::Plate(p) => p.generate(name.unwrap_or("Plate"))?,
        };
        if let Some(n) = name {
            f.name = n.to_string();
        }
        Ok(f)
    }

    /// Check the parameters without generating (for editors loading a saved
    /// recipe).
    pub fn validate(&self) -> Result<(), FoilError> {
        match self {
            Recipe::Naca { designation } => naca::parse(designation).map(|_| ()),
            Recipe::Cst(c) => c.validate(),
            Recipe::Parsec(p) => p.validate(),
            Recipe::KarmanTrefftz(k) => k.validate(),
            Recipe::Plate(p) => p.validate(),
        }
    }
}

/// Chord stations bunched at the leading edge only (half cosine), leading
/// (0) to trailing edge (1): about 1e-4 c at the nose and 0.013 c at the
/// trailing edge, close to the NACA tabulations.  The 6-series and a-type
/// mean lines have a logarithmic slope singularity at the trailing edge;
/// bunched stations there would hand the solver a hooked, near-cusped
/// trailing edge that its sharp-edge Kutta treatment cannot resolve.
pub(crate) fn nose_spacing(nside: usize) -> Vec<f64> {
    (0..nside).map(|i| if i + 1 == nside { 1.0 } else { 1.0 - (0.5 * std::f64::consts::PI * i as f64 / (nside - 1) as f64).cos() }).collect()
}

/// Cosine chord stations, leading (0) to trailing edge (1).
pub(crate) fn cosine_spacing(nside: usize) -> Vec<f64> {
    (0..nside).map(|i| 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / (nside - 1) as f64).cos())).collect()
}

/// A foil from upper and lower surfaces, each given from the leading edge
/// (shared first point) to the trailing edge.
pub(crate) fn from_surfaces(name: impl Into<String>, upper: &[(f64, f64)], lower: &[(f64, f64)]) -> Result<Foil, FoilError> {
    // A closed trailing edge must be one point: two points 1e-12 apart in
    // the wrong order give XFoil an inverted trailing-edge panel.  A real
    // inversion (upper below lower) is an error.
    let (mut upper, mut lower) = (upper.to_vec(), lower.to_vec());
    if let (Some(u), Some(l)) = (upper.last().copied(), lower.last().copied()) {
        if (u.0 - l.0).hypot(u.1 - l.1) >= 1e-9 && u.1 < l.1 {
            return Err(FoilError::Geometry(format!("the upper trailing edge ({:.6}) is below the lower one ({:.6})", u.1, l.1)));
        }
        if (u.0 - l.0).hypot(u.1 - l.1) < 1e-9 {
            let m = (0.5 * (u.0 + l.0), 0.5 * (u.1 + l.1));
            *upper.last_mut().unwrap() = m;
            *lower.last_mut().unwrap() = m;
        }
    }
    let mut x = Vec::with_capacity(upper.len() + lower.len());
    let mut y = Vec::with_capacity(upper.len() + lower.len());
    for &(a, b) in upper.iter().rev() {
        x.push(a);
        y.push(b);
    }
    for &(a, b) in &lower[1..] {
        x.push(a);
        y.push(b);
    }
    if x.iter().chain(&y).any(|v| !v.is_finite()) {
        return Err(FoilError::Geometry("the generator produced non-finite coordinates".into()));
    }
    Foil::new(name, x, y)
}

/// Reject a contour that crosses itself (any two non-adjacent segments
/// intersecting), naming the generator in the error.
pub(crate) fn check_contour(f: &Foil, what: &str) -> Result<(), FoilError> {
    let n = f.n();
    let p = |i: usize| (f.x[i], f.y[i]);
    let closed = (f.x[0] - f.x[n - 1]).hypot(f.y[0] - f.y[n - 1]) < 1e-12;
    let cross = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0);
    // segment boxes, to skip most pairs cheaply
    let boxes: Vec<[f64; 4]> = (0..n - 1).map(|i| [f.x[i].min(f.x[i + 1]), f.x[i].max(f.x[i + 1]), f.y[i].min(f.y[i + 1]), f.y[i].max(f.y[i + 1])]).collect();
    for i in 0..n - 1 {
        for j in i + 2..n - 1 {
            if closed && i == 0 && j == n - 2 {
                continue; // both end at the closed trailing edge
            }
            let (bi, bj) = (&boxes[i], &boxes[j]);
            if bi[1] < bj[0] || bj[1] < bi[0] || bi[3] < bj[2] || bj[3] < bi[2] {
                continue;
            }
            let (a, b, c, d) = (p(i), p(i + 1), p(j), p(j + 1));
            let (d1, d2, d3, d4) = (cross(a, b, c), cross(a, b, d), cross(c, d, a), cross(c, d, b));
            if d1 * d2 < 0.0 && d3 * d4 < 0.0 {
                return Err(FoilError::Geometry(format!("{what}: the contour crosses itself near x/c {:.3}", 0.25 * (a.0 + b.0 + c.0 + d.0))));
            }
        }
    }
    Ok(())
}

/// Whether the upper surface stays above the lower one (sampled at the
/// stations both surfaces share).
pub(crate) fn surfaces_cross(upper: &[(f64, f64)], lower: &[(f64, f64)]) -> bool {
    upper.iter().zip(lower).skip(1).take(upper.len().saturating_sub(2)).any(|(u, l)| u.1 < l.1 - 1e-12)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recipes_round_trip_through_json() {
        let recipes = vec![
            Recipe::Naca { designation: "64A010".into() },
            Recipe::Cst(cst::Cst::default()),
            Recipe::Parsec(parsec::Parsec::default()),
            Recipe::KarmanTrefftz(kt::KarmanTrefftz::default()),
            Recipe::Plate(plate::Plate::default()),
        ];
        for r in recipes {
            let j = serde_json::to_string(&r).unwrap();
            let back: Recipe = serde_json::from_str(&j).unwrap();
            assert_eq!(back, r, "{j}");
            let f = r.generate(Some("x")).unwrap();
            assert_eq!(f.recipe.as_ref(), Some(&r), "{j}");
            assert_eq!(f.name, "x");
            f.validate().unwrap();
            r.validate().unwrap();
        }
        // the family generators attach the recipe themselves
        assert!(cst::Cst::default().generate("c").unwrap().recipe.is_some());
        assert!(parsec::Parsec::default().generate("p").unwrap().recipe.is_some());
        assert!(kt::KarmanTrefftz::default().generate("k").unwrap().recipe.is_some());
        assert!(plate::Plate::default().generate("p").unwrap().recipe.is_some());
        assert!(naca::parse("64A010").unwrap().generate().unwrap().recipe.is_some());
    }

    #[test]
    fn crossing_contours_are_found() {
        let f = Foil::new("loop", vec![1.0, 0.5, 0.0, 0.5, 1.0, 0.4, 0.3, 0.6, 0.9, 1.0], vec![0.0, 0.1, 0.0, -0.1, -0.05, 0.2, -0.2, 0.02, 0.03, 0.0]).unwrap();
        assert!(check_contour(&f, "test").is_err());
        let g = naca::parse("2412").unwrap().generate().unwrap();
        check_contour(&g, "NACA 2412").unwrap();
    }

    #[test]
    fn spacings_span_the_chord() {
        for v in [cosine_spacing(NSIDE), nose_spacing(NSIDE)] {
            assert_eq!(v.len(), NSIDE);
            assert_eq!(v[0], 0.0);
            assert_eq!(*v.last().unwrap(), 1.0);
            assert!(v.windows(2).all(|w| w[1] > w[0]));
        }
    }
}
