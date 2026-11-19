//! PARSEC sections (H. Sobieczky, "Parametric airfoils and wings", Notes on
//! Numerical Fluid Mechanics 68, 1998).
//!
//! Each surface is `z(x) = sum_{n=1..6} a_n x^(n - 1/2)`; its six
//! coefficients follow from six conditions: the leading-edge radius, the
//! crest position, height and curvature, and the trailing-edge height and
//! slope.

use super::{check_contour, check_thickness, cosine_spacing, from_surfaces, surfaces_cross, Recipe, Surfaces, NSIDE};
use crate::{Foil, FoilError};
use foiley_core::linalg::solve_dense;
use serde::{Deserialize, Serialize};

/// The 11 PARSEC parameters (lengths are chord fractions, angles degrees).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Parsec {
    /// Leading-edge radius.
    pub r_le: f64,
    /// Upper crest: position, height and curvature z'' (negative).
    pub x_up: f64,
    pub z_up: f64,
    pub zxx_up: f64,
    /// Lower crest: position, height (negative) and curvature (positive).
    pub x_lo: f64,
    pub z_lo: f64,
    pub zxx_lo: f64,
    /// Height of the trailing-edge midpoint and trailing-edge thickness.
    pub z_te: f64,
    pub dz_te: f64,
    /// Direction of the trailing-edge bisector, degrees below the chord
    /// line (positive for an aft-cambered section).
    pub alpha_te: f64,
    /// Trailing-edge wedge angle, degrees.
    pub beta_te: f64,
}

impl Default for Parsec {
    /// Close to NACA 0012.
    fn default() -> Self {
        Parsec { r_le: 0.0155, x_up: 0.30, z_up: 0.060, zxx_up: -0.45, x_lo: 0.30, z_lo: -0.060, zxx_lo: 0.45, z_te: 0.0, dz_te: 0.00252, alpha_te: 0.0, beta_te: 16.0 }
    }
}

/// Coefficients a_1..a_6 of one surface.
fn coefficients(sign: f64, r_le: f64, xc: f64, zc: f64, zxx: f64, z_te: f64, slope_te: f64) -> Result<[f64; 6], FoilError> {
    // exponent of a_{k+1}: k + 1/2
    let e: Vec<f64> = (0..6).map(|k| k as f64 + 0.5).collect();
    let a = vec![
        vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        vec![1.0; 6],
        e.clone(),
        e.iter().map(|&e| xc.powf(e)).collect(),
        e.iter().map(|&e| e * xc.powf(e - 1.0)).collect(),
        e.iter().map(|&e| e * (e - 1.0) * xc.powf(e - 2.0)).collect(),
    ];
    let mut b = vec![0.0; 6];
    b[0] = sign * (2.0 * r_le).sqrt();
    b[1] = z_te;
    b[2] = slope_te;
    b[3] = zc;
    b[4] = 0.0;
    b[5] = zxx;
    let x = solve_dense(&a, &b).map_err(|e| FoilError::Geometry(format!("PARSEC: the conditions are inconsistent ({e:?})")))?;
    let mut out = [0.0; 6];
    out.copy_from_slice(&x);
    if out.iter().any(|v| !v.is_finite()) {
        return Err(FoilError::Geometry("PARSEC: the conditions are inconsistent".into()));
    }
    Ok(out)
}

fn eval(a: &[f64; 6], x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    a.iter().enumerate().map(|(k, c)| c * x.powf(k as f64 + 0.5)).sum()
}

impl Parsec {
    pub fn validate(&self) -> Result<(), FoilError> {
        let v = [self.r_le, self.x_up, self.z_up, self.zxx_up, self.x_lo, self.z_lo, self.zxx_lo, self.z_te, self.dz_te, self.alpha_te, self.beta_te];
        if v.iter().any(|v| !v.is_finite()) {
            return Err(FoilError::Geometry("PARSEC: non-finite parameter".into()));
        }
        if self.r_le <= 0.0 || self.r_le > 0.2 {
            return Err(FoilError::Geometry("PARSEC: the leading-edge radius must be in (0, 0.2]".into()));
        }
        if !(0.01..=0.99).contains(&self.x_up) || !(0.01..=0.99).contains(&self.x_lo) {
            return Err(FoilError::Geometry("PARSEC: crest positions must be in 0.01..0.99".into()));
        }
        if self.z_up <= self.z_lo {
            return Err(FoilError::Geometry("PARSEC: the upper crest must be above the lower crest".into()));
        }
        if self.zxx_up >= 0.0 || self.zxx_lo <= 0.0 {
            return Err(FoilError::Geometry("PARSEC: a crest needs curvature z'' < 0 on the upper surface and > 0 on the lower (otherwise it is a trough)".into()));
        }
        if self.dz_te < 0.0 || self.beta_te < 0.0 || self.beta_te >= 90.0 || self.alpha_te.abs() >= 45.0 {
            return Err(FoilError::Geometry("PARSEC: need TE thickness >= 0, wedge angle 0..90 deg, TE direction within +-45 deg".into()));
        }
        Ok(())
    }

    /// Coefficients of the upper and lower surface.
    pub fn coefficients(&self) -> Result<([f64; 6], [f64; 6]), FoilError> {
        self.validate()?;
        let (al, be) = (self.alpha_te.to_radians(), self.beta_te.to_radians());
        let up = coefficients(1.0, self.r_le, self.x_up, self.z_up, self.zxx_up, self.z_te + 0.5 * self.dz_te, -(al + 0.5 * be).tan())?;
        let lo = coefficients(-1.0, self.r_le, self.x_lo, self.z_lo, self.zxx_lo, self.z_te - 0.5 * self.dz_te, -(al - 0.5 * be).tan())?;
        Ok((up, lo))
    }

    /// Surface points at the chord stations `x` (leading to trailing edge).
    pub fn surfaces_at(&self, x: &[f64]) -> Result<Surfaces, FoilError> {
        let (u, l) = self.coefficients()?;
        Ok((x.iter().map(|&x| (x, eval(&u, x))).collect(), x.iter().map(|&x| (x, eval(&l, x))).collect()))
    }

    pub fn generate(&self, name: &str) -> Result<Foil, FoilError> {
        let (up, lo) = self.surfaces_at(&cosine_spacing(NSIDE))?;
        if surfaces_cross(&up, &lo) {
            return Err(FoilError::Geometry("PARSEC: the upper surface dips below the lower one".into()));
        }
        // the crests must be the extreme points of their surfaces: the
        // conditions only make them stationary, and far-off crests can let
        // the polynomial run away elsewhere
        let (cu, cl) = self.coefficients()?;
        for k in 0..=2000 {
            let x = k as f64 / 2000.0;
            let (zu, zl) = (eval(&cu, x), eval(&cl, x));
            if zu > self.z_up + 1e-9 || zl < self.z_lo - 1e-9 {
                return Err(FoilError::Geometry(format!(
                    "PARSEC: the {} crest is not the extreme of its surface (it reaches {:.4} at x/c {x:.3})",
                    if zu > self.z_up + 1e-9 { "upper" } else { "lower" },
                    if zu > self.z_up + 1e-9 { zu } else { zl }
                )));
            }
        }
        check_thickness("PARSEC", |x| eval(&cu, x) - eval(&cl, x))?;
        let mut f = from_surfaces(name, &up, &lo)?;
        check_contour(&f, "PARSEC")?;
        f.note = "PARSEC".into();
        f.recipe = Some(Recipe::Parsec(*self));
        Ok(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn d(a: &[f64; 6], x: f64, order: i32) -> f64 {
        a.iter()
            .enumerate()
            .map(|(k, c)| {
                let e = k as f64 + 0.5;
                match order {
                    1 => c * e * x.powf(e - 1.0),
                    _ => c * e * (e - 1.0) * x.powf(e - 2.0),
                }
            })
            .sum()
    }

    #[test]
    fn the_surfaces_meet_every_condition() {
        let p = Parsec { r_le: 0.012, x_up: 0.35, z_up: 0.075, zxx_up: -0.6, x_lo: 0.25, z_lo: -0.04, zxx_lo: 0.5, z_te: -0.01, dz_te: 0.002, alpha_te: 6.0, beta_te: 10.0 };
        let (u, l) = p.coefficients().unwrap();
        assert!((eval(&u, p.x_up) - p.z_up).abs() < 1e-12);
        assert!(d(&u, p.x_up, 1).abs() < 1e-12);
        assert!((d(&u, p.x_up, 2) - p.zxx_up).abs() < 1e-10);
        assert!((eval(&l, p.x_lo) - p.z_lo).abs() < 1e-12);
        assert!((d(&l, p.x_lo, 2) - p.zxx_lo).abs() < 1e-10);
        assert!((eval(&u, 1.0) - eval(&l, 1.0) - p.dz_te).abs() < 1e-12);
        assert!((0.5 * (eval(&u, 1.0) + eval(&l, 1.0)) - p.z_te).abs() < 1e-12);
        let wedge = (d(&l, 1.0, 1).atan() - d(&u, 1.0, 1).atan()).to_degrees();
        assert!((wedge - p.beta_te).abs() < 1e-9, "{wedge}");
        assert!((u[0] * u[0] / 2.0 - p.r_le).abs() < 1e-12);
        let f = p.generate("p").unwrap();
        assert!((f.props().le_radius - p.r_le).abs() < 2e-3, "{}", f.props().le_radius);
    }

    #[test]
    fn the_default_is_close_to_naca_0012() {
        let f = Parsec::default().generate("p").unwrap();
        let p = f.props();
        assert!((p.thickness - 0.12).abs() < 2e-3, "{p:?}");
        assert!(p.camber.abs() < 1e-6, "{p:?}");
    }

    #[test]
    fn bad_parameters_are_rejected() {
        let p = Parsec { z_up: -0.1, ..Default::default() };
        assert!(p.generate("p").is_err());
        let p = Parsec { r_le: 0.0, ..Default::default() };
        assert!(p.generate("p").is_err());
        // a "crest" with positive curvature is a trough
        assert!(Parsec { zxx_up: 0.1, ..Default::default() }.generate("p").is_err());
        assert!(Parsec { zxx_lo: -0.1, ..Default::default() }.generate("p").is_err());
        // crests at x/c 0.99 let the surfaces run away to ~1e5 chords
        let p = Parsec { x_up: 0.99, x_lo: 0.99, ..Default::default() };
        assert!(p.generate("p").is_err());
    }
}
