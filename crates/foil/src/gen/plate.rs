//! Cambered plates: a constant-thickness plate bent to a circular-arc or
//! parabolic camber line, with a semicircular leading edge and a square-cut
//! or semicircular trailing edge.  Typical of model and micro-air-vehicle
//! wings at very low Reynolds number.
//!
//! XFoil's viscous solution needs the square trailing edge (a blunt edge cut
//! vertically at x/c = 1): its boundary layer cannot follow a round cap
//! that closes to a point, and no viscous point converges.

use super::{check_contour, cosine_spacing, Recipe, NSIDE};
use crate::{Foil, FoilError};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PlateCamber {
    /// Circular arc, maximum camber at mid-chord.
    CircularArc,
    /// Two parabolas meeting at the point of maximum camber (the NACA
    /// 4-digit mean line).
    Parabolic,
}

/// Shape of the plate's trailing edge.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PlateEdge {
    /// Cut vertically at x/c = 1: a blunt trailing edge about as tall as the
    /// plate is thick (what XFoil can analyse).
    #[default]
    Square,
    /// Semicircular, like the leading edge (geometry only: XFoil's viscous
    /// solution does not converge on it).
    Round,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Plate {
    /// Plate thickness / chord.
    pub thickness: f64,
    /// Maximum camber / chord.
    pub camber: f64,
    /// Position of maximum camber (x/c); parabolic camber only.
    pub x_camber: f64,
    pub kind: PlateCamber,
    #[serde(default)]
    pub trailing_edge: PlateEdge,
}

impl Default for Plate {
    fn default() -> Self {
        Plate { thickness: 0.02, camber: 0.04, x_camber: 0.5, kind: PlateCamber::CircularArc, trailing_edge: PlateEdge::Square }
    }
}

impl Plate {
    pub fn validate(&self) -> Result<(), FoilError> {
        if ![self.thickness, self.camber, self.x_camber].iter().all(|v| v.is_finite()) {
            return Err(FoilError::Geometry("plate: non-finite parameter".into()));
        }
        if !(0.002..=0.15).contains(&self.thickness) {
            return Err(FoilError::Geometry("plate: thickness must be 0.2%..15% of the chord".into()));
        }
        if self.camber.abs() > 0.15 {
            return Err(FoilError::Geometry("plate: camber must be within 15% of the chord".into()));
        }
        if self.kind == PlateCamber::Parabolic && !(0.1..=0.9).contains(&self.x_camber) {
            return Err(FoilError::Geometry("plate: the position of maximum camber must be 0.1..0.9".into()));
        }
        // the inner side of the plate folds over where its half-thickness
        // exceeds the radius of curvature of the camber line
        let rmin = self.min_camber_radius();
        if 0.5 * self.thickness >= rmin {
            return Err(FoilError::Geometry(format!(
                "plate: half the thickness ({:.4}) must be less than the smallest radius of curvature of the camber line ({rmin:.4}); reduce the thickness or the camber",
                0.5 * self.thickness
            )));
        }
        Ok(())
    }

    /// Smallest radius of curvature of the camber line.
    fn min_camber_radius(&self) -> f64 {
        let h = self.camber.abs();
        if h == 0.0 {
            return f64::INFINITY;
        }
        match self.kind {
            PlateCamber::CircularArc => 0.25 / (2.0 * h) + 0.5 * h,
            // y'' = -2h/p^2 ahead of the crest, -2h/(1-p)^2 behind it; the
            // radius (1 + y'^2)^1.5 / |y''| is smallest at the crest (y' = 0)
            PlateCamber::Parabolic => {
                let p = self.x_camber;
                p.min(1.0 - p).powi(2) / (2.0 * h)
            }
        }
    }

    /// Camber line ordinate and slope at x.
    pub fn camber_line(&self, x: f64) -> (f64, f64) {
        let h = self.camber;
        if h == 0.0 {
            return (0.0, 0.0);
        }
        match self.kind {
            PlateCamber::CircularArc => {
                // with the curvature k = 1/rho = 2h / (0.25 + h^2):
                // y = h - k d^2 / (1 + sqrt(1 - k^2 d^2)), y' = -k d / sqrt(1 - k^2 d^2)
                // (no rho^2, so a tiny camber neither overflows nor cancels)
                let k = 2.0 * h / (0.25 + h * h);
                let d = x - 0.5;
                let w = (1.0 - k * k * d * d).max(0.0).sqrt();
                (h - k * d * d / (1.0 + w), -k * d / w)
            }
            PlateCamber::Parabolic => {
                let p = self.x_camber;
                if x < p {
                    (h / (p * p) * (2.0 * p * x - x * x), 2.0 * h / (p * p) * (p - x))
                } else {
                    let q = 1.0 - p;
                    (h / (q * q) * ((1.0 - 2.0 * p) + 2.0 * p * x - x * x), 2.0 * h / (q * q) * (p - x))
                }
            }
        }
    }

    pub fn generate(&self, name: &str) -> Result<Foil, FoilError> {
        self.validate()?;
        let r = 0.5 * self.thickness;
        let round_te = self.trailing_edge == PlateEdge::Round;
        let normal = |x: f64| {
            let (y, s) = self.camber_line(x);
            let l = (1.0 + s * s).sqrt();
            ((x, y), (-s / l, 1.0 / l), (1.0 / l, s / l))
        };
        // Camber stations of the last upper and lower points.  A round edge
        // ends both at the centre of its circle; a square edge is cut
        // vertically at x/c = 1, so each surface ends where its offset point
        // reaches x = 1 (solve x -/+ r n_x(x) = 1 by fixed point).
        let end = |side: f64| {
            if round_te {
                return 1.0 - r;
            }
            let mut x = 1.0;
            for _ in 0..60 {
                x = 1.0 - side * r * normal(x).1 .0;
            }
            x
        };
        let (xu_end, xl_end) = (end(1.0), end(-1.0));
        let stations = |x_end: f64| -> Vec<f64> { cosine_spacing(NSIDE).iter().map(|s| r + (x_end - r) * s).collect() };
        let (xs_up, xs_lo) = (stations(xu_end), stations(xl_end));
        let ncap = 24;
        let mut pts: Vec<(f64, f64)> = Vec::new();
        let (c1, n1, t1) = normal(1.0 - r);
        if round_te {
            // trailing-edge cap, from its tip round to the upper surface
            for k in 0..ncap {
                let a = 0.5 * PI * k as f64 / ncap as f64;
                pts.push((c1.0 + r * (a.cos() * t1.0 + a.sin() * n1.0), c1.1 + r * (a.cos() * t1.1 + a.sin() * n1.1)));
            }
        }
        // upper surface, trailing to leading edge
        for &x in xs_up.iter().rev() {
            let (c, n, _) = normal(x);
            pts.push((c.0 + r * n.0, c.1 + r * n.1));
        }
        // leading-edge cap: angle measured from the upper normal round the nose
        let (c0, n0, t0) = normal(r);
        for k in 1..2 * ncap {
            let b = PI * k as f64 / (2 * ncap) as f64;
            pts.push((c0.0 + r * (b.cos() * n0.0 - b.sin() * t0.0), c0.1 + r * (b.cos() * n0.1 - b.sin() * t0.1)));
        }
        // lower surface, leading to trailing edge
        for &x in &xs_lo {
            let (c, n, _) = normal(x);
            pts.push((c.0 - r * n.0, c.1 - r * n.1));
        }
        if round_te {
            // trailing-edge cap from the lower surface back to the tip
            for k in 1..ncap {
                let a = 0.5 * PI * k as f64 / ncap as f64;
                pts.push((c1.0 + r * (a.sin() * t1.0 - a.cos() * n1.0), c1.1 + r * (a.sin() * t1.1 - a.cos() * n1.1)));
            }
            pts.push(pts[0]);
        } else {
            // the cut is exactly at x/c = 1
            let last = pts.len() - 1;
            pts[0].0 = 1.0;
            pts[last].0 = 1.0;
        }
        pts.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-12 && (a.1 - b.1).abs() < 1e-12);
        if pts.iter().any(|p| !p.0.is_finite() || !p.1.is_finite()) {
            return Err(FoilError::Geometry("plate: non-finite coordinates".into()));
        }
        let (x, y): (Vec<f64>, Vec<f64>) = pts.into_iter().unzip();
        let mut f = Foil::new(name, x, y)?;
        check_contour(&f, "plate")?;
        f.recipe = Some(Recipe::Plate(*self));
        f.note = format!(
            "{} plate, t/c {:.1}%, camber {:.1}%, {} trailing edge",
            match self.kind {
                PlateCamber::CircularArc => "circular-arc",
                PlateCamber::Parabolic => "parabolic",
            },
            self.thickness * 100.0,
            self.camber * 100.0,
            if round_te { "round" } else { "square" }
        );
        Ok(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn folding_plates_are_rejected() {
        let p = Plate { thickness: 0.15, camber: 0.15, x_camber: 0.1, kind: PlateCamber::Parabolic, trailing_edge: PlateEdge::Square };
        assert!(p.generate("p").is_err());
        // the extremes Studio allows that do not fold still generate
        for kind in [PlateCamber::CircularArc, PlateCamber::Parabolic] {
            for te in [PlateEdge::Square, PlateEdge::Round] {
                for t in [0.002, 0.05, 0.15] {
                    for c in [-0.15, 0.0, 0.08, 0.15] {
                        for xc in [0.1, 0.5, 0.9] {
                            let p = Plate { thickness: t, camber: c, x_camber: xc, kind, trailing_edge: te };
                            if p.validate().is_ok() {
                                p.generate("p").unwrap_or_else(|e| panic!("{p:?}: {e}"));
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn steep_and_tiny_cambers_stay_well_defined() {
        // a steep aft-loaded plate near the fold limit: the square cut is
        // vertical at x/c = 1, no point lies aft of it, and the reported
        // maximum thickness is no larger than the cut itself (before, the
        // slanted cut reported 33% at x/c 1.021)
        let p = Plate { thickness: 0.066, camber: 0.15, x_camber: 0.9, kind: PlateCamber::Parabolic, trailing_edge: PlateEdge::Square };
        let f = p.generate("p").unwrap();
        let pr = f.props();
        assert!(pr.thickness.is_finite() && pr.thickness <= pr.te_gap + 1e-3 && pr.x_thickness < 1.01, "{pr:?}");
        assert!(f.x.iter().all(|&x| x <= 1.0 + 1e-12));
        // a vanishing circular-arc camber neither overflows nor cancels
        for h in [1e-160, -1e-160, 1e-9] {
            let f = Plate { camber: h, ..Plate::default() }.generate("p").unwrap();
            assert!(f.x.iter().chain(&f.y).all(|v| v.is_finite()), "{h}");
        }
        let (y, s) = Plate { camber: 0.06, ..Plate::default() }.camber_line(0.0);
        assert!(y.abs() < 1e-15 && s > 0.2, "{y} {s}");
    }

    #[test]
    fn camber_lines_peak_where_asked() {
        let p = Plate { camber: 0.06, ..Default::default() };
        let (y, s) = p.camber_line(0.5);
        assert!((y - 0.06).abs() < 1e-12 && s.abs() < 1e-12);
        assert!(p.camber_line(0.0).0.abs() < 1e-12 && p.camber_line(1.0).0.abs() < 1e-12);
        let q = Plate { camber: 0.05, x_camber: 0.3, kind: PlateCamber::Parabolic, ..Default::default() };
        let (y, s) = q.camber_line(0.3);
        assert!((y - 0.05).abs() < 1e-12 && s.abs() < 1e-12);
    }

    #[test]
    fn the_plate_has_constant_thickness_and_round_edges() {
        for (kind, trailing_edge) in [(PlateCamber::CircularArc, PlateEdge::Round), (PlateCamber::Parabolic, PlateEdge::Round), (PlateCamber::CircularArc, PlateEdge::Square)] {
            let p = Plate { thickness: 0.03, camber: 0.05, x_camber: 0.4, kind, trailing_edge };
            let f = p.generate("plate").unwrap();
            let pr = f.props();
            assert!((pr.thickness - 0.03).abs() < 1e-3, "{kind:?} {pr:?}");
            assert!((pr.camber - 0.05).abs() < 3e-3, "{kind:?} {pr:?}");
            assert!((pr.le_radius - 0.015).abs() < 3e-3, "{kind:?} {pr:?}");
            let n = f.n();
            let gap = (f.x[0] - f.x[n - 1]).hypot(f.y[0] - f.y[n - 1]);
            match trailing_edge {
                PlateEdge::Round => assert_eq!(gap, 0.0),
                PlateEdge::Square => {
                    // cut vertically at x/c = 1: at least the plate thickness
                    assert_eq!((f.x[0], f.x[n - 1]), (1.0, 1.0));
                    assert!((0.03 - 1e-12..0.04).contains(&gap), "{gap}");
                }
            }
            assert!(f.x.iter().all(|&x| x <= 1.0 + 1e-12), "{kind:?} {trailing_edge:?}");
            let xmin = f.x.iter().cloned().fold(f64::MAX, f64::min);
            let xmax = f.x.iter().cloned().fold(f64::MIN, f64::max);
            assert!(xmin.abs() < 0.01 && (xmax - 1.0).abs() < 0.01, "{xmin} {xmax}");
        }
    }
}
