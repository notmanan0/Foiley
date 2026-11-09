//! Foils and 2D polars available to the 3D analyses, with fast sampled
//! shapes for meshing and polar interpolation for viscous drag.

use foiley_core::hash::{Hash32, Hasher};
use foiley_foil::{Foil, FoilPoint, FoilPolar};
use serde::{Deserialize, Serialize};
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

/// Where a viscous interpolation landed relative to the 2D polar data.
/// Stored per strip as a `u8` (`StripResult.envelope`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Envelope {
    #[default]
    Inside = 0,
    /// Reynolds number outside the polars' range (clamped to the nearest)
    ClampedRe = 1,
    /// CL (or angle of attack, when interpolating by angle) outside a polar
    ClampedCl = 2,
    ClampedBoth = 3,
    /// no usable polar for the foil
    Missing = 4,
}

impl Envelope {
    pub const KINDS: [Envelope; 5] = [Envelope::Inside, Envelope::ClampedRe, Envelope::ClampedCl, Envelope::ClampedBoth, Envelope::Missing];
    pub fn from_flags(re: bool, cl: bool) -> Envelope {
        match (re, cl) {
            (false, false) => Envelope::Inside,
            (true, false) => Envelope::ClampedRe,
            (false, true) => Envelope::ClampedCl,
            (true, true) => Envelope::ClampedBoth,
        }
    }
    pub fn from_u8(v: u8) -> Envelope {
        Envelope::KINDS.get(v as usize).copied().unwrap_or(Envelope::Missing)
    }
    pub fn re_clamped(self) -> bool {
        matches!(self, Envelope::ClampedRe | Envelope::ClampedBoth)
    }
    pub fn cl_clamped(self) -> bool {
        matches!(self, Envelope::ClampedCl | Envelope::ClampedBoth)
    }
    pub fn is_clamped(self) -> bool {
        self.re_clamped() || self.cl_clamped()
    }
    /// The worse of two (Missing dominates, clamp flags add up).
    pub fn combine(self, o: Envelope) -> Envelope {
        if self == Envelope::Missing || o == Envelope::Missing {
            Envelope::Missing
        } else {
            Envelope::from_flags(self.re_clamped() || o.re_clamped(), self.cl_clamped() || o.cl_clamped())
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Envelope::Inside => "inside",
            Envelope::ClampedRe => "clamped in Re",
            Envelope::ClampedCl => "clamped in CL",
            Envelope::ClampedBoth => "clamped in Re and CL",
            Envelope::Missing => "no polar",
        }
    }
}

/// Hash of a foil's normalized coordinates (unit chord, LE at the origin,
/// de-rotated), quantized to 1e-9 chord.  The name is not part of it.
pub fn foil_geom_hash(f: &Foil) -> Hash32 {
    f.geom_hash()
}

/// Quantize a family parameter for comparison and hashing (1e-6).
fn qv(x: f64) -> i64 {
    let r = (x * 1e6).round();
    if r == 0.0 || !r.is_finite() {
        0
    } else {
        r as i64
    }
}

/// Everything that makes two 2D polars of one foil belong to different
/// families: only the Reynolds number varies inside a family.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct FamilyKey {
    /// hash of the foil's normalized coordinates (`foil_geom_hash`)
    pub foil_geom_hash: Hash32,
    pub polar_type: u8,
    pub mach: f64,
    pub ncrit: f64,
    pub xtr_top: f64,
    pub xtr_bot: f64,
    /// trailing-edge flap deflection of the polar's geometry, degrees (+ down)
    pub flap_deg: f64,
    /// moment reference x/c of the polar's Cm
    pub mom_ref_x: f64,
}

impl FamilyKey {
    /// Moment reference used when a polar does not record one.
    pub const DEFAULT_MOM_REF_X: f64 = 0.25;
    pub fn of(p: &FoilPolar, foil_geom_hash: Hash32) -> FamilyKey {
        let mom_ref_x = p.settings.as_ref().map(|s| s.cm_ref.0).unwrap_or(FamilyKey::DEFAULT_MOM_REF_X);
        FamilyKey { foil_geom_hash, polar_type: p.polar_type, mach: p.mach, ncrit: p.ncrit, xtr_top: p.xtr_top, xtr_bot: p.xtr_bot, flap_deg: p.flap_deg, mom_ref_x }
    }
    fn q(&self) -> (Hash32, u8, [i64; 6]) {
        (self.foil_geom_hash, self.polar_type, [qv(self.mach), qv(self.ncrit), qv(self.xtr_top), qv(self.xtr_bot), qv(self.flap_deg), qv(self.mom_ref_x)])
    }
    /// The same family apart from the flap angle.
    pub fn same_except_flap(&self, o: &FamilyKey) -> bool {
        let (a, b) = (self.q(), o.q());
        a.0 == b.0 && a.1 == b.1 && a.2[0..4] == b.2[0..4] && a.2[5] == b.2[5]
    }
    pub fn describe(&self) -> String {
        format!(
            "type {} M={} Ncrit={} xtr=({}, {}){}{} geometry {}",
            self.polar_type,
            self.mach,
            self.ncrit,
            self.xtr_top,
            self.xtr_bot,
            if qv(self.flap_deg) != 0 { format!(" flap {:+} deg", self.flap_deg) } else { String::new() },
            if qv(self.mom_ref_x) != qv(FamilyKey::DEFAULT_MOM_REF_X) { format!(" Cm about x/c={}", self.mom_ref_x) } else { String::new() },
            self.foil_geom_hash.short()
        )
    }
    /// Stable hash of the key (for manifests).
    pub fn hash(&self) -> Hash32 {
        let mut h = Hasher::new();
        h.tag("foiley.family.v1").bytes(&self.foil_geom_hash.0).u64(self.polar_type as u64);
        for v in [self.mach, self.ncrit, self.xtr_top, self.xtr_bot, self.flap_deg, self.mom_ref_x] {
            h.f64(v);
        }
        h.finish()
    }
}

impl PartialEq for FamilyKey {
    fn eq(&self, o: &Self) -> bool {
        self.q() == o.q()
    }
}
impl Eq for FamilyKey {}
impl std::hash::Hash for FamilyKey {
    fn hash<H: std::hash::Hasher>(&self, h: &mut H) {
        self.q().hash(h)
    }
}

/// Index of a family in `FoilDb::families()`.  Ids are stable while the
/// database is not modified.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct FamilyId(pub u32);

/// The polars of one foil that differ only in Reynolds number.
#[derive(Clone, Debug)]
pub struct Family {
    pub foil: String,
    pub key: FamilyKey,
    /// sorted by Re
    pub polars: Vec<FoilPolar>,
}

impl Family {
    /// Usable for viscous interpolation: type-1 viscous polars.
    pub fn usable(&self) -> bool {
        self.key.polar_type == 1 && !self.polars.is_empty()
    }
    pub fn re_range(&self) -> Option<(f64, f64)> {
        Some((self.polars.first()?.re, self.polars.last()?.re))
    }
}

/// Which family of each foil a 3D analysis uses.  `None` fields match any
/// value; the selection must leave exactly one family per foil.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct FamilySelector {
    pub mach: Option<f64>,
    pub ncrit: Option<f64>,
    pub xtr_top: Option<f64>,
    pub xtr_bot: Option<f64>,
}

impl FamilySelector {
    pub fn matches(&self, k: &FamilyKey) -> bool {
        let m = |s: Option<f64>, v: f64| s.is_none_or(|s| qv(s) == qv(v));
        m(self.mach, k.mach) && m(self.ncrit, k.ncrit) && m(self.xtr_top, k.xtr_top) && m(self.xtr_bot, k.xtr_bot)
    }
    pub fn describe(&self) -> String {
        let f = |n: &str, v: Option<f64>| v.map(|v| format!("{n}={v}"));
        let v: Vec<String> = [f("M", self.mach), f("Ncrit", self.ncrit), f("xtr_top", self.xtr_top), f("xtr_bot", self.xtr_bot)].into_iter().flatten().collect();
        if v.is_empty() {
            "any family".into()
        } else {
            v.join(" ")
        }
    }
}

/// The set of foils and 2D polars known to an analysis.  Polars are grouped
/// into families (`FamilyKey`); interpolation works inside one family.
#[derive(Clone, Debug, Default)]
pub struct FoilDb {
    foils: HashMap<String, (Arc<FoilShape>, Hash32)>,
    families: Vec<Family>,
    index: HashMap<(String, FamilyKey), FamilyId>,
    /// polars computed on another geometry than the loaded foil's (their
    /// `foil_hash` does not match): kept out of every family
    stale: Vec<FoilPolar>,
}

impl FoilDb {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn add_foil(&mut self, f: &Foil) {
        let h = foil_geom_hash(f);
        let old = self.foils.insert(f.name.clone(), (Arc::new(FoilShape::from_foil(f)), h));
        if old.is_none_or(|o| o.1 != h) && self.families.iter().any(|fam| fam.foil == f.name) {
            // polars of this foil are re-keyed to the new geometry
            self.rebuild();
        }
    }
    fn rebuild(&mut self) {
        let mut polars: Vec<FoilPolar> = std::mem::take(&mut self.families).into_iter().flat_map(|f| f.polars).collect();
        polars.append(&mut self.stale);
        self.index.clear();
        for p in &polars {
            self.insert(p);
        }
    }
    /// Hash of the loaded geometry of a foil.
    pub fn foil_hash(&self, name: &str) -> Option<Hash32> {
        self.foils.get(name).map(|f| f.1)
    }
    /// Add a polar. A polar with the same content identity replaces the old
    /// one. Rounded display names may be shared by distinct Re stations.
    /// Polars with fewer than two
    /// points and inviscid polars are not kept.
    pub fn add_polar(&mut self, p: &FoilPolar) {
        if p.points.len() < 2 || !p.viscous {
            return;
        }
        let identity = p.identity();
        for fam in self.families.iter_mut().filter(|f| f.foil == p.foil) {
            fam.polars.retain(|q| q.re != p.re || q.identity() != identity);
        }
        self.stale.retain(|q| q.foil != p.foil || q.re != p.re || q.identity() != identity);
        self.insert(p);
    }
    /// Polars left out because their geometry hash does not match the
    /// loaded foil (the foil was edited after they were computed).
    pub fn stale_polars(&self) -> &[FoilPolar] {
        &self.stale
    }
    fn insert(&mut self, p: &FoilPolar) {
        let current = self.foil_hash(&p.foil);
        if let (Some(c), Some(h)) = (current, p.foil_hash) {
            if c != h {
                self.stale.push(p.clone());
                return;
            }
        }
        let key = FamilyKey::of(p, current.or(p.foil_hash).unwrap_or(Hash32::ZERO));
        let families = &mut self.families;
        let id = *self.index.entry((p.foil.clone(), key)).or_insert_with(|| {
            families.push(Family { foil: p.foil.clone(), key, polars: vec![] });
            FamilyId(families.len() as u32 - 1)
        });
        let v = &mut self.families[id.0 as usize].polars;
        v.push(p.clone());
        v.sort_by(|a, b| a.re.total_cmp(&b.re));
    }
    pub fn shape(&self, name: &str) -> Option<Arc<FoilShape>> {
        self.foils.get(name).map(|f| f.0.clone())
    }
    pub fn has_foil(&self, name: &str) -> bool {
        self.foils.contains_key(name)
    }
    pub fn families(&self) -> &[Family] {
        &self.families
    }
    pub fn family(&self, id: FamilyId) -> &Family {
        &self.families[id.0 as usize]
    }
    /// Usable (type-1 viscous) families of a foil, in a fixed order (by key,
    /// not by insertion).
    pub fn families_of(&self, foil: &str) -> Vec<FamilyId> {
        let mut v: Vec<FamilyId> = (0..self.families.len()).map(|i| FamilyId(i as u32)).filter(|&i| self.family(i).foil == foil && self.family(i).usable()).collect();
        v.sort_by_key(|&i| {
            let k = self.family(i).key;
            (qv(k.flap_deg), qv(k.mach), qv(k.ncrit), qv(k.xtr_top), qv(k.xtr_bot), qv(k.mom_ref_x), k.foil_geom_hash)
        });
        v
    }
    /// Whether a foil has any usable polar.
    pub fn has_polars(&self, foil: &str) -> bool {
        self.families.iter().any(|f| f.foil == foil && f.usable())
    }
    /// Resolve the undeflected family of `foil` that `sel` selects.
    /// `Ok(None)`: the foil has no usable polar at all.  Errors on no match
    /// or several matches, listing the candidates.
    pub fn resolve(&self, foil: &str, sel: &FamilySelector) -> Result<Option<FamilyId>, String> {
        let all = self.families_of(foil);
        if all.is_empty() {
            return Ok(None);
        }
        let list = |v: &[FamilyId]| v.iter().map(|&i| format!("  {}", self.family(i).key.describe())).collect::<Vec<_>>().join("\n");
        let hit: Vec<FamilyId> = all.iter().copied().filter(|&i| qv(self.family(i).key.flap_deg) == 0 && sel.matches(&self.family(i).key)).collect();
        match hit.len() {
            1 => Ok(Some(hit[0])),
            0 => Err(format!("no 2D polar family of foil '{foil}' matches the selection ({}); candidates:\n{}", sel.describe(), list(&all))),
            _ => Err(format!("{} 2D polar families of foil '{foil}' match the selection ({}): set the family's Mach, Ncrit or trips to choose one; candidates:\n{}", hit.len(), sel.describe(), list(&hit))),
        }
    }

    /// The families that differ from `base` only in flap angle (`base`
    /// included), sorted by flap angle: what a deflected strip brackets.
    pub fn flap_families(&self, base: FamilyId) -> Vec<(f64, FamilyId)> {
        let k = self.family(base).key;
        let foil = &self.family(base).foil;
        let mut v: Vec<(f64, FamilyId)> = self.families_of(foil).into_iter().filter(|&i| self.family(i).key.same_except_flap(&k)).map(|i| (self.family(i).key.flap_deg, i)).collect();
        v.sort_by(|a, b| a.0.total_cmp(&b.0));
        v
    }

    /// Re range covered by a family.
    pub fn re_range(&self, id: FamilyId) -> Option<(f64, f64)> {
        self.family(id).re_range()
    }

    /// Interpolate a 2D point at (Re, CL) in one family: linear in Re
    /// between the bracketing polars.  Returns the point and where it fell
    /// relative to the data (never `Missing`: that is `None`).
    pub fn at_re_cl(&self, id: FamilyId, re: f64, cl: f64) -> Option<(FoilPoint, Envelope)> {
        self.interp_re(id, re, |p| p.at_cl(cl))
    }

    /// Interpolate at (Re, alpha).
    pub fn at_re_alpha(&self, id: FamilyId, re: f64, alpha: f64) -> Option<(FoilPoint, Envelope)> {
        self.interp_re(id, re, |p| p.at_alpha(alpha))
    }

    fn interp_re(&self, id: FamilyId, re: f64, f: impl Fn(&FoilPolar) -> Option<(FoilPoint, bool)>) -> Option<(FoilPoint, Envelope)> {
        let v = &self.family(id).polars;
        if v.is_empty() {
            return None;
        }
        if re <= v[0].re {
            let (p, o) = f(&v[0])?;
            return Some((p, Envelope::from_flags(re < v[0].re * 0.98, o)));
        }
        let last = &v[v.len() - 1];
        if re >= last.re {
            let (p, o) = f(last)?;
            return Some((p, Envelope::from_flags(re > last.re * 1.02, o)));
        }
        let j = v.partition_point(|p| p.re < re);
        let (a, b) = (&v[j - 1], &v[j]);
        let (pa, oa) = f(a)?;
        let (pb, ob) = f(b)?;
        let t = (re - a.re) / (b.re - a.re);
        Some((foiley_foil::polar::lerp_point(&pa, &pb, t), Envelope::from_flags(false, oa || ob)))
    }

    /// Linearised lift: (zero-lift angle deg, slope per degree) at a Re.
    pub fn linear_lift(&self, id: FamilyId, re: f64) -> Option<(f64, f64)> {
        let v = &self.family(id).polars;
        if v.is_empty() {
            return None;
        }
        let pick = |p: &FoilPolar| Some((p.alpha0()?, p.cl_alpha()?));
        if re <= v[0].re {
            return pick(&v[0]);
        }
        if re >= v[v.len() - 1].re {
            return pick(&v[v.len() - 1]);
        }
        let j = v.partition_point(|p| p.re < re);
        let (a0a, sa) = pick(&v[j - 1])?;
        let (a0b, sb) = pick(&v[j])?;
        let t = (re - v[j - 1].re) / (v[j].re - v[j - 1].re);
        Some((a0a + t * (a0b - a0a), sa + t * (sb - sa)))
    }
}
