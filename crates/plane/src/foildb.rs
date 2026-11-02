//! Foils and 2D polars available to the 3D analyses, with fast sampled
//! shapes for meshing and polar interpolation for viscous drag.

use foiley_foil::{Foil, FoilPoint, FoilPolar};
use std::collections::HashMap;
use std::sync::Arc;

/// Sampled foil shape: upper / lower / camber y at cosine-spaced x/c.
#[derive(Clone, Debug)]
pub struct FoilShape {
    pub xs: Vec<f64>,
    pub upper: Vec<f64>,
    pub lower: Vec<f64>,
    pub hinge: Option<(f64, f64)>,
    /// cross-section area / c^2
    pub area: f64,
}

impl FoilShape {
    pub fn from_foil(f: &Foil) -> FoilShape {
        let g = f.normalized();
        let (xs, c, t) = g.camber_thickness(120);
        let upper: Vec<f64> = c.iter().zip(&t).map(|(c, t)| c + 0.5 * t).collect();
        let lower: Vec<f64> = c.iter().zip(&t).map(|(c, t)| c - 0.5 * t).collect();
        let mut area = 0.0;
        for i in 1..xs.len() {
            area += 0.5 * ((upper[i] - lower[i]) + (upper[i - 1] - lower[i - 1])) * (xs[i] - xs[i - 1]);
        }
        FoilShape { xs, upper, lower, hinge: f.hinge, area }
    }
    /// A thin flat plate (used when a foil is missing and the user allows it).
    pub fn flat_plate() -> FoilShape {
        let xs = foiley_core::interp::distribution(foiley_core::Spacing::Cosine, 40);
        let z = vec![0.0; xs.len()];
        FoilShape { xs, upper: z.clone(), lower: z, hinge: None, area: 0.0 }
    }
    pub fn upper_at(&self, x: f64) -> f64 {
        foiley_core::interp::lerp_table(&self.xs, &self.upper, x)
    }
    pub fn lower_at(&self, x: f64) -> f64 {
        foiley_core::interp::lerp_table(&self.xs, &self.lower, x)
    }
    pub fn camber_at(&self, x: f64) -> f64 {
        0.5 * (self.upper_at(x) + self.lower_at(x))
    }
}

/// The set of foils and 2D polars known to an analysis.
#[derive(Clone, Debug, Default)]
pub struct FoilDb {
    foils: HashMap<String, Arc<FoilShape>>,
    polars: HashMap<String, Vec<FoilPolar>>,
}

impl FoilDb {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn add_foil(&mut self, f: &Foil) {
        self.foils.insert(f.name.clone(), Arc::new(FoilShape::from_foil(f)));
    }
    /// Add a polar (only type-1 viscous polars are used for interpolation).
    pub fn add_polar(&mut self, p: &FoilPolar) {
        if p.points.len() < 2 {
            return;
        }
        let v = self.polars.entry(p.foil.clone()).or_default();
        v.retain(|q| q.name != p.name);
        v.push(p.clone());
        v.sort_by(|a, b| a.re.total_cmp(&b.re));
    }
    pub fn shape(&self, name: &str) -> Option<Arc<FoilShape>> {
        self.foils.get(name).cloned()
    }
    pub fn has_foil(&self, name: &str) -> bool {
        self.foils.contains_key(name)
    }
    pub fn polars_of(&self, foil: &str) -> &[FoilPolar] {
        self.polars.get(foil).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Usable polars for viscous interpolation (type 1, viscous).
    fn t1(&self, foil: &str) -> Vec<&FoilPolar> {
        self.polars_of(foil).iter().filter(|p| p.viscous && p.polar_type == 1 && p.points.len() >= 2).collect()
    }

    /// Re range covered by the polars of a foil.
    pub fn re_range(&self, foil: &str) -> Option<(f64, f64)> {
        let v = self.t1(foil);
        Some((v.first()?.re, v.last()?.re))
    }

    /// Interpolate a 2D point at (Re, CL) for one foil: linear in Re
    /// between the bracketing polars.  Returns (point, out_of_envelope).
    pub fn at_re_cl(&self, foil: &str, re: f64, cl: f64) -> Option<(FoilPoint, bool)> {
        self.interp_re(foil, re, |p| p.at_cl(cl))
    }

    /// Interpolate at (Re, alpha).
    pub fn at_re_alpha(&self, foil: &str, re: f64, alpha: f64) -> Option<(FoilPoint, bool)> {
        self.interp_re(foil, re, |p| p.at_alpha(alpha))
    }

    fn interp_re(&self, foil: &str, re: f64, f: impl Fn(&FoilPolar) -> Option<(FoilPoint, bool)>) -> Option<(FoilPoint, bool)> {
        let v = self.t1(foil);
        if v.is_empty() {
            return None;
        }
        if re <= v[0].re {
            let (p, o) = f(v[0])?;
            return Some((p, o || re < v[0].re * 0.98));
        }
        let last = v[v.len() - 1];
        if re >= last.re {
            let (p, o) = f(last)?;
            return Some((p, o || re > last.re * 1.02));
        }
        let j = v.partition_point(|p| p.re < re);
        let (a, b) = (v[j - 1], v[j]);
        let (pa, oa) = f(a)?;
        let (pb, ob) = f(b)?;
        let t = (re - a.re) / (b.re - a.re);
        Some((foiley_foil::polar::lerp_point(&pa, &pb, t), oa || ob))
    }

    /// Linearised lift: (zero-lift angle deg, slope per degree) at a Re.
    pub fn linear_lift(&self, foil: &str, re: f64) -> Option<(f64, f64)> {
        let v = self.t1(foil);
        if v.is_empty() {
            return None;
        }
        let pick = |p: &FoilPolar| Some((p.alpha0()?, p.cl_alpha()?));
        if re <= v[0].re {
            return pick(v[0]);
        }
        if re >= v[v.len() - 1].re {
            return pick(v[v.len() - 1]);
        }
        let j = v.partition_point(|p| p.re < re);
        let (a0a, sa) = pick(v[j - 1])?;
        let (a0b, sb) = pick(v[j])?;
        let t = (re - v[j - 1].re) / (v[j].re - v[j - 1].re);
        Some((a0a + t * (a0b - a0a), sa + t * (sb - sa)))
    }
}
