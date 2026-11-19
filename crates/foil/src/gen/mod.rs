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
//! made, so an editor can reopen the parameters.

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
/// any operation that changes the shape (repaneling keeps it).
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
        let (mut f, recipe) = match self {
            Recipe::Naca { designation } => {
                let s = naca::parse(designation)?;
                (s.generate()?, Recipe::Naca { designation: s.designation })
            }
            Recipe::Cst(c) => (c.generate(name.unwrap_or("CST"))?, self.clone()),
            Recipe::Parsec(p) => (p.generate(name.unwrap_or("PARSEC"))?, self.clone()),
            Recipe::KarmanTrefftz(k) => (k.generate(name.unwrap_or("Karman-Trefftz"))?, self.clone()),
            Recipe::Plate(p) => (p.generate(name.unwrap_or("Plate"))?, self.clone()),
        };
        if let Some(n) = name {
            f.name = n.to_string();
        }
        f.recipe = Some(recipe);
        Ok(f)
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
    // the wrong order give XFoil an inverted trailing-edge panel.
    let (mut upper, mut lower) = (upper.to_vec(), lower.to_vec());
    if let (Some(u), Some(l)) = (upper.last().copied(), lower.last().copied()) {
        if (u.0 - l.0).hypot(u.1 - l.1) < 1e-9 || u.1 < l.1 {
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
        }
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
