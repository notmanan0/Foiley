//! Kulfan's class-shape transformation (CST).
//!
//! Each surface is
//!
//! ```text
//! y(x) = C(x) S(x) + x dz,   C(x) = x^N1 (1-x)^N2,
//! S(x) = sum_i A_i K_i x^i (1-x)^(n-i),   K_i = n! / (i! (n-i)!)
//! ```
//!
//! with N1 = 0.5 and N2 = 1 for a round-nosed airfoil (B. Kulfan, "Universal
//! parametric geometry representation method", J. Aircraft 45(1), 2008).
//! A0 sets the leading-edge radius (r = A0^2 / 2), A_n the trailing-edge
//! wedge angle, and dz the trailing-edge offset of that surface.

use super::{check_contour, check_thickness, cosine_spacing, from_surfaces, surfaces_cross, Recipe, Surfaces, NSIDE};
use crate::{Foil, FoilError};
use foiley_core::linalg::{lstsq, DenseMatrix};
use serde::{Deserialize, Serialize};

fn half() -> f64 {
    0.5
}
fn one() -> f64 {
    1.0
}

/// CST weights of both surfaces (same order on both).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Cst {
    /// Upper-surface weights A_0..A_n (positive for a normal section).
    pub upper: Vec<f64>,
    /// Lower-surface weights (negative for a normal section).
    pub lower: Vec<f64>,
    /// Class-function exponents.
    #[serde(default = "half")]
    pub n1: f64,
    #[serde(default = "one")]
    pub n2: f64,
    /// Trailing-edge y of the upper and lower surface (chord fractions).
    #[serde(default)]
    pub te_upper: f64,
    #[serde(default)]
    pub te_lower: f64,
}

impl Default for Cst {
    /// A 12% symmetric section close to NACA 0012 (order 5).
    fn default() -> Self {
        let a = vec![0.1703, 0.1603, 0.1458, 0.1473, 0.1371, 0.1460];
        Cst { lower: a.iter().map(|v| -v).collect(), upper: a, n1: 0.5, n2: 1.0, te_upper: 0.00063, te_lower: -0.00063 }
    }
}

/// How well a CST fit follows the foil it was fitted to.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct CstFit {
    /// RMS and largest vertical distance between the foil's points and the
    /// CST surfaces, chord fractions (of the normalized foil).
    pub rms: f64,
    pub max: f64,
}

fn binomial(n: usize, k: usize) -> f64 {
    (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
}

/// The class function times each Bernstein term at x: the columns of the
/// fit and the effect of each weight.
pub fn basis(order: usize, n1: f64, n2: f64, x: f64) -> Vec<f64> {
    let c = if x <= 0.0 || x >= 1.0 { 0.0 } else { x.powf(n1) * (1.0 - x).powf(n2) };
    (0..=order).map(|i| c * binomial(order, i) * x.powi(i as i32) * (1.0 - x).powi((order - i) as i32)).collect()
}

/// One CST surface at x.
pub fn surface(w: &[f64], n1: f64, n2: f64, dz: f64, x: f64) -> f64 {
    if w.is_empty() {
        return x * dz;
    }
    basis(w.len() - 1, n1, n2, x).iter().zip(w).map(|(b, a)| b * a).sum::<f64>() + x * dz
}

impl Cst {
    /// Polynomial order n (weights per surface minus one).
    pub fn order(&self) -> usize {
        self.upper.len().max(1) - 1
    }

    pub fn validate(&self) -> Result<(), FoilError> {
        if self.upper.is_empty() || self.upper.len() != self.lower.len() {
            return Err(FoilError::Geometry("CST: both surfaces need the same number of weights (at least one)".into()));
        }
        if self.upper.len() > 31 {
            return Err(FoilError::Geometry("CST: at most 30th order".into()));
        }
        // N1 = 0 or N2 = 0 would leave the class function nonzero at the
        // nose or the tail, so the surfaces would not meet there
        if !(self.n1 > 0.0 && self.n1 <= 1.0 && self.n2 > 0.0 && self.n2 <= 2.0) {
            return Err(FoilError::Geometry("CST: N1 must be in (0, 1] and N2 in (0, 2]".into()));
        }
        if self.upper.iter().chain(&self.lower).chain([&self.te_upper, &self.te_lower, &self.n1, &self.n2]).any(|v| !v.is_finite()) {
            return Err(FoilError::Geometry("CST: non-finite parameter".into()));
        }
        if self.te_upper < self.te_lower {
            return Err(FoilError::Geometry(format!("CST: the upper trailing edge ({}) is below the lower one ({})", self.te_upper, self.te_lower)));
        }
        Ok(())
    }

    /// Surface points at the chord stations `x` (leading to trailing edge).
    pub fn surfaces_at(&self, x: &[f64]) -> Surfaces {
        let up = x.iter().map(|&x| (x, surface(&self.upper, self.n1, self.n2, self.te_upper, x))).collect();
        let lo = x.iter().map(|&x| (x, surface(&self.lower, self.n1, self.n2, self.te_lower, x))).collect();
        (up, lo)
    }

    pub fn generate(&self, name: &str) -> Result<Foil, FoilError> {
        self.validate()?;
        let (up, lo) = self.surfaces_at(&cosine_spacing(NSIDE));
        if surfaces_cross(&up, &lo) {
            return Err(FoilError::Geometry("CST: the upper surface dips below the lower one".into()));
        }
        check_thickness("CST", |x| surface(&self.upper, self.n1, self.n2, self.te_upper, x) - surface(&self.lower, self.n1, self.n2, self.te_lower, x))?;
        let mut f = from_surfaces(name, &up, &lo)?;
        check_contour(&f, "CST")?;
        f.note = format!("CST order {}", self.order());
        f.recipe = Some(Recipe::Cst(self.clone()));
        Ok(f)
    }

    /// Least-squares CST weights of `order` for `foil` (normalized first:
    /// leading edge at the origin, unit chord, chord line horizontal).
    pub fn fit(foil: &Foil, order: usize) -> Result<(Cst, CstFit), FoilError> {
        if order > 30 {
            return Err(FoilError::Geometry("CST: at most 30th order".into()));
        }
        foil.validate()?;
        let f = foil.normalized();
        let n = f.n();
        let ile = f.le_index();
        // upper: TE..LE (indices 0..=ile), lower: LE..TE
        let upper: Vec<(f64, f64)> = (0..=ile).rev().map(|i| (f.x[i], f.y[i])).collect();
        let lower: Vec<(f64, f64)> = (ile..n).map(|i| (f.x[i], f.y[i])).collect();
        let (n1, n2) = (0.5, 1.0);
        let fit_side = |pts: &[(f64, f64)], dz: f64| -> Result<Vec<f64>, FoilError> {
            let rows: Vec<&(f64, f64)> = pts.iter().filter(|p| p.0 > 1e-9 && p.0 < 1.0 - 1e-9).collect();
            if rows.len() < order + 1 {
                return Err(FoilError::Geometry(format!("CST: {} points on a surface cannot fit order {order}", rows.len())));
            }
            let mut a = DenseMatrix::zeros(rows.len(), order + 1);
            let mut b = DenseMatrix::zeros(rows.len(), 1);
            for (r, p) in rows.iter().enumerate() {
                for (j, v) in basis(order, n1, n2, p.0).into_iter().enumerate() {
                    a.set(r, j, v);
                }
                b.set(r, 0, p.1 - p.0 * dz);
            }
            let w = lstsq(&a, &b).map_err(|e| FoilError::Geometry(format!("CST fit: {e:?}")))?;
            Ok(w.col(0).to_vec())
        };
        let te_upper = upper.last().map(|p| p.1).unwrap_or(0.0);
        let te_lower = lower.last().map(|p| p.1).unwrap_or(0.0);
        let cst = Cst { upper: fit_side(&upper, te_upper)?, lower: fit_side(&lower, te_lower)?, n1, n2, te_upper, te_lower };
        // a high-order fit can oscillate until its surfaces cross: report
        // that instead of returning weights `generate` would reject
        cst.generate("fit").map_err(|e| FoilError::Geometry(format!("the order-{order} CST fit is not a valid foil ({e}); try a lower order")))?;
        let mut sq = 0.0;
        let mut max: f64 = 0.0;
        let mut cnt = 0usize;
        for (pts, w, dz) in [(&upper, &cst.upper, te_upper), (&lower, &cst.lower, te_lower)] {
            for p in pts.iter() {
                let e = (surface(w, n1, n2, dz, p.0.clamp(0.0, 1.0)) - p.1).abs();
                sq += e * e;
                max = max.max(e);
                cnt += 1;
            }
        }
        Ok((cst, CstFit { rms: (sq / cnt.max(1) as f64).sqrt(), max }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn leading_edge_radius_is_a0_squared_over_two() {
        let c = Cst::default();
        // y ~ A0 sqrt(x) near the nose -> r = A0^2 / 2
        let x = 1e-8;
        let y = surface(&c.upper, 0.5, 1.0, 0.0, x);
        let r = y * y / (2.0 * x);
        assert!((r - c.upper[0] * c.upper[0] / 2.0).abs() < 1e-6, "{r}");
    }

    #[test]
    fn fitting_naca_sections_is_accurate() {
        for (code, order, tol) in [("2412", 8, 3e-4), ("0012", 6, 3e-4), ("4415", 10, 3e-4)] {
            let f = Foil::naca(code).unwrap();
            let (c, fit) = Cst::fit(&f, order).unwrap();
            assert!(fit.max < tol && fit.rms < tol / 4.0, "{code}: {fit:?}");
            assert_eq!(c.order(), order);
            let g = c.generate("fit").unwrap();
            let (p, q) = (f.props(), g.props());
            assert!((p.thickness - q.thickness).abs() < 5e-4, "{code} {p:?} {q:?}");
            assert!((p.camber - q.camber).abs() < 5e-4, "{code}");
        }
    }

    #[test]
    fn a_fit_of_a_cst_foil_recovers_its_weights() {
        let c = Cst::default();
        let f = c.generate("x").unwrap();
        let (back, fit) = Cst::fit(&f, c.order()).unwrap();
        assert!(fit.max < 1e-9, "{fit:?}");
        for (a, b) in back.upper.iter().zip(&c.upper).chain(back.lower.iter().zip(&c.lower)) {
            assert!((a - b).abs() < 1e-6, "{a} {b}");
        }
    }

    #[test]
    fn crossing_surfaces_are_rejected() {
        let mut c = Cst::default();
        c.lower = c.upper.iter().map(|v| v + 0.05).collect();
        assert!(c.generate("x").is_err());
    }

    #[test]
    fn unsupported_class_exponents_and_inverted_edges_are_rejected() {
        for (n1, n2) in [(0.5, 0.0), (0.0, 1.0)] {
            assert!(Cst { n1, n2, ..Cst::default() }.generate("x").is_err(), "{n1} {n2}");
        }
        assert!(Cst { te_upper: -1e-8, te_lower: 0.0, ..Cst::default() }.generate("x").is_err());
        assert!(Cst { upper: vec![], lower: vec![], ..Cst::default() }.generate("x").is_err());
        assert!(Cst { upper: vec![0.1, 0.2], lower: vec![-0.1], ..Cst::default() }.validate().is_err());
        // order 0 is a valid (if crude) section
        assert!(Cst { upper: vec![0.15], lower: vec![-0.15], ..Cst::default() }.generate("x").is_ok());
    }

    #[test]
    fn a_fit_that_would_cross_is_refused_and_lower_orders_regenerate() {
        let f = Foil::naca("63A210").unwrap();
        assert!(Cst::fit(&f, 30).is_err(), "the order-30 fit of 63A210 oscillates through itself");
        for order in [4, 8, 12] {
            let (c, _) = Cst::fit(&f, order).unwrap_or_else(|e| panic!("order {order}: {e}"));
            c.generate("x").unwrap_or_else(|e| panic!("order {order}: {e}"));
        }
    }

    #[test]
    fn coincident_and_between_station_crossings_are_rejected() {
        assert!(Cst { upper: vec![0.0], lower: vec![0.0], te_upper: 0.0, te_lower: 0.0, ..Cst::default() }.generate("x").is_err());
        let up = vec![0.249975, -0.250025, 0.249975];
        let c = Cst { lower: up.iter().map(|v| -v).collect(), upper: up, te_upper: 0.0, te_lower: 0.0, ..Cst::default() };
        assert!(c.generate("x").is_err());
    }
}
