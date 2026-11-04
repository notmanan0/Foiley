//! Airfoil geometry.

use crate::FoilError;
use foiley::cli::{Input, Session};
use foiley::geom::{geopar, naca4, naca5};
use foiley::spline::{scalc, segspl, seval, sinvrt};
use foiley::xfoil::{XFoil, IBX, IQX};
use serde::{Deserialize, Serialize};

/// An airfoil: coordinates ordered from the trailing edge along the upper
/// surface to the leading edge and back along the lower surface
/// (counter-clockwise, "Selig" order), chord normally 1.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Foil {
    pub name: String,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Trailing-edge flap hinge (x/c, y/t in 0..1), if this foil is used as a
    /// control surface in a 3D wing.  Deflection is applied by the 3D mesher.
    #[serde(default)]
    pub hinge: Option<(f64, f64)>,
    /// Free-form description / provenance ("NACA 2412", "flap 10 deg"...).
    #[serde(default)]
    pub note: String,
}

/// Geometric properties of a foil.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct FoilProps {
    pub points: usize,
    pub chord: f64,
    /// max thickness / chord
    pub thickness: f64,
    pub x_thickness: f64,
    /// max camber / chord
    pub camber: f64,
    pub x_camber: f64,
    /// leading-edge radius / chord
    pub le_radius: f64,
    /// trailing-edge gap / chord
    pub te_gap: f64,
    /// trailing-edge included angle, degrees
    pub te_angle: f64,
    /// cross-section area / chord^2
    pub area: f64,
}

/// 1-based copy (index 0 unused) as the XFoil engine expects.
fn one_based(v: &[f64]) -> Vec<f64> {
    let mut o = Vec::with_capacity(v.len() + 1);
    o.push(0.0);
    o.extend_from_slice(v);
    o
}

fn quiet_engine() -> XFoil {
    foiley::set_quiet(true);
    let mut xf = XFoil::new();
    xf.verbose = false;
    xf.par_inner = false;
    xf
}

impl Foil {
    pub fn new(name: impl Into<String>, x: Vec<f64>, y: Vec<f64>) -> Result<Foil, FoilError> {
        let f = Foil { name: name.into(), x, y, hinge: None, note: String::new() };
        f.validate()?;
        Ok(f)
    }

    pub fn n(&self) -> usize {
        self.x.len()
    }

    /// Check that the geometry is usable by the solvers.
    pub fn validate(&self) -> Result<(), FoilError> {
        let n = self.x.len();
        if n != self.y.len() {
            return Err(FoilError::Geometry("x and y have different lengths".into()));
        }
        if n < 10 {
            return Err(FoilError::Geometry(format!("only {n} points; at least 10 are needed")));
        }
        if n > IBX - 2 {
            return Err(FoilError::Geometry(format!("{n} points; the maximum is {}", IBX - 2)));
        }
        if self.x.iter().chain(&self.y).any(|v| !v.is_finite()) {
            return Err(FoilError::Geometry("coordinates contain NaN or infinite values".into()));
        }
        let chord = self.x.iter().cloned().fold(f64::MIN, f64::max) - self.x.iter().cloned().fold(f64::MAX, f64::min);
        if chord <= 1e-9 {
            return Err(FoilError::Geometry("zero chord".into()));
        }
        Ok(())
    }

    /// NACA 4- or 5-digit section ("2412", "0012", "23012").
    pub fn naca(code: &str) -> Result<Foil, FoilError> {
        let digits: String = code.chars().filter(|c| c.is_ascii_digit()).collect();
        let ides: i32 = digits.parse().map_err(|_| FoilError::Parse(format!("'{code}' is not a NACA 4/5-digit designation")))?;
        let _q = foiley::quiet_scope();
        let nside = IQX / 3;
        let r = match digits.len() {
            1..=4 => Some(naca4(ides, nside)),
            5 => naca5(ides, nside),
            _ => None,
        };
        let (xb, yb, nb, _name) = r.ok_or_else(|| FoilError::Parse(format!("NACA {digits} is not a supported designation (4-digit, or 5-digit 210xx-250xx)")))?;
        let mut f = Foil::new(format!("NACA {digits}"), xb[1..=nb].to_vec(), yb[1..=nb].to_vec())?;
        f.note = format!("NACA {digits} generated");
        Ok(f)
    }

    /// Parse a coordinate file (Selig or Lednicer .dat; blanks, tabs or commas).
    pub fn from_dat(text: &str, fallback_name: &str) -> Result<Foil, FoilError> {
        let mut name = String::new();
        let mut pts: Vec<(f64, f64)> = Vec::new();
        for line in text.lines() {
            let l = line.trim();
            if l.is_empty() || l.starts_with('#') {
                continue;
            }
            let nums: Vec<f64> = l.split(|c: char| c.is_whitespace() || c == ',' || c == ';').filter(|t| !t.is_empty()).map_while(|t| t.parse::<f64>().ok()).collect();
            if nums.len() >= 2 {
                pts.push((nums[0], nums[1]));
            } else if pts.is_empty() && name.is_empty() {
                name = l.to_string();
            }
        }
        if name.is_empty() {
            name = fallback_name.to_string();
        }
        // Lednicer format: first numeric line holds the point counts (> 1).
        if let Some(&(a, b)) = pts.first() {
            if a > 1.5 && b > 1.5 && a.fract() == 0.0 && b.fract() == 0.0 {
                let (nu, nl) = (a as usize, b as usize);
                let body = &pts[1..];
                if body.len() >= nu + nl {
                    let upper = &body[..nu];
                    let lower = &body[nu..nu + nl];
                    let mut out: Vec<(f64, f64)> = upper.iter().rev().cloned().collect();
                    let skip = if lower.first().map(|p| (p.0 - upper[0].0).abs() < 1e-9 && (p.1 - upper[0].1).abs() < 1e-9).unwrap_or(false) { 1 } else { 0 };
                    out.extend(lower[skip..].iter().cloned());
                    pts = out;
                } else {
                    pts.remove(0);
                }
            }
        }
        // drop consecutive duplicates
        pts.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-12 && (a.1 - b.1).abs() < 1e-12);
        let (mut x, mut y): (Vec<f64>, Vec<f64>) = pts.into_iter().unzip();
        // make counter-clockwise
        let n = x.len();
        let mut area = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            area += 0.5 * (y[i] + y[j]) * (x[i] - x[j]);
        }
        if area < 0.0 {
            x.reverse();
            y.reverse();
        }
        let mut f = Foil::new(name.trim().to_string(), x, y)?;
        f.note = "loaded from file".into();
        Ok(f)
    }

    pub fn load(path: &std::path::Path) -> Result<Foil, FoilError> {
        let text = std::fs::read_to_string(path).map_err(|e| FoilError::Io(format!("{}: {e}", path.display())))?;
        let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
        Foil::from_dat(&text, &stem)
    }

    /// Selig-format .dat text.
    pub fn to_dat(&self) -> String {
        let mut s = format!("{}\n", self.name);
        for (x, y) in self.x.iter().zip(&self.y) {
            s.push_str(&format!("{:10.6} {:10.6}\n", x, y));
        }
        s
    }

    pub fn save(&self, path: &std::path::Path) -> Result<(), FoilError> {
        std::fs::write(path, self.to_dat()).map_err(|e| FoilError::Io(format!("{}: {e}", path.display())))
    }

    /// Build an XFoil engine instance with this foil as the buffer airfoil,
    /// paneled with `npan` panels (XFoil PANGEN).
    pub fn engine(&self, pan: &Paneling) -> Result<XFoil, FoilError> {
        self.validate()?;
        let mut xf = quiet_engine();
        let n = self.n();
        xf.set_buffer(&one_based(&self.x), &one_based(&self.y), n, &self.name);
        xf.npan = pan.panels.clamp(20, IQX - 6);
        xf.cvpar = pan.le_bunching;
        xf.cterat = pan.te_bunching;
        xf.ctrrat = pan.refined_ratio;
        xf.pangen(false);
        if xf.n < 10 {
            return Err(FoilError::Geometry("paneling failed".into()));
        }
        Ok(xf)
    }

    /// Build an XFoil engine with the panel nodes `input` asks for: PANGEN
    /// ([`Foil::engine`]) or the coordinates themselves (XFoil LOAD + PCOP).
    pub fn engine_nodes(&self, input: &NodeInput) -> Result<XFoil, FoilError> {
        if let NodeInput::Pangen(p) = input {
            return self.engine(p);
        }
        self.validate()?;
        let n = self.n();
        if n > IQX - 5 {
            return Err(FoilError::Geometry(format!("'{}' has {n} points: preserving them as panel nodes allows at most {}; use PANGEN paneling", self.name, IQX - 5)));
        }
        let mut xf = quiet_engine();
        xf.set_buffer(&one_based(&self.x), &one_based(&self.y), n, &self.name);
        xf.abcopy(false);
        if xf.n < 10 {
            return Err(FoilError::Geometry(format!("'{}': fewer than 10 distinct panel nodes", self.name)));
        }
        Ok(xf)
    }

    /// Repanel with XFoil's curvature-based distribution.
    pub fn repaneled(&self, pan: &Paneling) -> Result<Foil, FoilError> {
        let xf = self.engine(pan)?;
        let mut f = self.clone();
        f.x = xf.x[1..=xf.n].to_vec();
        f.y = xf.y[1..=xf.n].to_vec();
        f.note = format!("repaneled ({} panels)", xf.n);
        Ok(f)
    }

    /// Geometric properties (via XFoil GEOPAR on a spline of the coordinates).
    pub fn props(&self) -> FoilProps {
        let n = self.n();
        if n < 3 {
            return FoilProps::default();
        }
        let x = one_based(&self.x);
        let y = one_based(&self.y);
        let mut s = vec![0.0; n + 1];
        let mut xp = vec![0.0; n + 1];
        let mut yp = vec![0.0; n + 1];
        scalc(&x, &y, &mut s, n);
        segspl(&x, &mut xp, &s, n);
        segspl(&y, &mut yp, &s, n);
        let _q = foiley::quiet_scope();
        let g = geopar(&x, &xp, &y, &yp, &s, n, false);
        let gap = ((self.x[0] - self.x[n - 1]).powi(2) + (self.y[0] - self.y[n - 1]).powi(2)).sqrt();
        let ch = if g.chord > 0.0 { g.chord } else { 1.0 };
        FoilProps {
            points: n,
            chord: g.chord,
            thickness: g.thick / ch,
            x_thickness: g.xthick,
            camber: g.cambr / ch,
            x_camber: g.xcambr,
            le_radius: g.radle / ch,
            te_gap: gap / ch,
            te_angle: g.angte,
            area: g.area / (ch * ch),
        }
    }

    /// Index of the leading-edge point (minimum x).
    pub fn le_index(&self) -> usize {
        let mut k = 0;
        for i in 1..self.n() {
            if self.x[i] < self.x[k] {
                k = i;
            }
        }
        k
    }

    /// y of the upper and lower surfaces at chord position `xc` (0..1, in
    /// the foil's own x range), by spline inversion.
    pub fn surface_y(&self, xc: f64) -> (f64, f64) {
        let _q = foiley::quiet_scope();
        let n = self.n();
        let x = one_based(&self.x);
        let y = one_based(&self.y);
        let mut s = vec![0.0; n + 1];
        let mut xp = vec![0.0; n + 1];
        let mut yp = vec![0.0; n + 1];
        scalc(&x, &y, &mut s, n);
        segspl(&x, &mut xp, &s, n);
        segspl(&y, &mut yp, &s, n);
        let ile = self.le_index() + 1;
        let xmin = self.x[ile - 1];
        let xmax = self.x[0].max(self.x[n - 1]);
        let xi = xmin + xc.clamp(0.0, 1.0) * (xmax - xmin);
        let mut su = s[1] + (s[ile] - s[1]) * (1.0 - xc);
        let mut sl = s[ile] + (s[n] - s[ile]) * xc;
        sinvrt(&mut su, xi, &x, &xp, &s, n);
        sinvrt(&mut sl, xi, &x, &xp, &s, n);
        let su = su.clamp(s[1], s[ile]);
        let sl = sl.clamp(s[ile], s[n]);
        (seval(su, &y, &yp, &s, n), seval(sl, &y, &yp, &s, n))
    }

    /// Camber and thickness lines at `m+1` cosine-spaced stations:
    /// returns (x, camber, thickness).
    pub fn camber_thickness(&self, m: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let xs = foiley_core::interp::distribution(foiley_core::Spacing::Cosine, m);
        let mut c = Vec::with_capacity(xs.len());
        let mut t = Vec::with_capacity(xs.len());
        for &xc in &xs {
            let (yu, yl) = self.surface_y(xc);
            c.push(0.5 * (yu + yl));
            t.push(yu - yl);
        }
        (xs, c, t)
    }

    /// Mid-camber line y at chord fraction `xc` (used by the thin-surface 3D solvers).
    pub fn camber_at(&self, xc: f64) -> f64 {
        let (u, l) = self.surface_y(xc);
        0.5 * (u + l)
    }

    /// Normalise: de-rotate, unit chord, LE at origin.
    pub fn normalized(&self) -> Foil {
        let ile = self.le_index();
        let (xl, yl) = (self.x[ile], self.y[ile]);
        let n = self.n();
        let xt = 0.5 * (self.x[0] + self.x[n - 1]);
        let yt = 0.5 * (self.y[0] + self.y[n - 1]);
        let ang = (yt - yl).atan2(xt - xl);
        let ch = ((xt - xl).powi(2) + (yt - yl).powi(2)).sqrt().max(1e-12);
        let (s, c) = (-ang).sin_cos();
        let mut f = self.clone();
        for i in 0..n {
            let dx = self.x[i] - xl;
            let dy = self.y[i] - yl;
            f.x[i] = (c * dx - s * dy) / ch;
            f.y[i] = (s * dx + c * dy) / ch;
        }
        f.note = "normalized".into();
        f
    }

    /// Linear blend between two foils (`t` = 0 gives `a`, 1 gives `b`),
    /// done on camber and thickness so the result stays well formed.
    pub fn interpolate(a: &Foil, b: &Foil, t: f64, name: &str) -> Result<Foil, FoilError> {
        let m = 100;
        let (xs, ca, ta) = a.normalized().camber_thickness(m);
        let (_, cb, tb) = b.normalized().camber_thickness(m);
        let mut x = Vec::new();
        let mut y = Vec::new();
        for i in (0..=m).rev() {
            let c = ca[i] + t * (cb[i] - ca[i]);
            let th = ta[i] + t * (tb[i] - ta[i]);
            x.push(xs[i]);
            y.push(c + 0.5 * th);
        }
        for i in 1..=m {
            let c = ca[i] + t * (cb[i] - ca[i]);
            let th = ta[i] + t * (tb[i] - ta[i]);
            x.push(xs[i]);
            y.push(c - 0.5 * th);
        }
        let mut f = Foil::new(name, x, y)?;
        f.note = format!("{:.0}% {} / {:.0}% {}", (1.0 - t) * 100.0, a.name, t * 100.0, b.name);
        Ok(f)
    }

    /// Run XFoil GDES commands on this foil's coordinates and return the result.
    fn gdes(&self, cmds: &str, note: String) -> Result<Foil, FoilError> {
        self.validate()?;
        let script = format!("GDES\n{cmds}\n\nQUIT\n");
        let mut s = Session::new(Input::from_string(&script));
        foiley::set_quiet(true);
        s.xf.verbose = false;
        let n = self.n();
        s.xf.set_buffer(&one_based(&self.x), &one_based(&self.y), n, &self.name);
        s.xf.abcopy(false);
        s.run(None);
        let nb = s.xf.nb;
        let mut f = self.clone();
        f.x = s.xf.xb[1..=nb].to_vec();
        f.y = s.xf.yb[1..=nb].to_vec();
        f.note = note;
        f.validate()?;
        Ok(f)
    }

    /// Deflect a trailing-edge flap: hinge at `x_hinge` (x/c) and relative
    /// height `y_rel` (0 = lower surface, 1 = upper), `deg` positive down.
    pub fn with_flap(&self, x_hinge: f64, y_rel: f64, deg: f64) -> Result<Foil, FoilError> {
        if deg == 0.0 {
            return Ok(self.clone());
        }
        if !(0.05..0.98).contains(&x_hinge) {
            return Err(FoilError::Geometry("flap hinge x/c must be between 0.05 and 0.98".into()));
        }
        let mut f = self.gdes(&format!("FLAP {x_hinge} 999\n{y_rel}\n{deg}"), format!("flap {deg:+.1} deg at x/c={x_hinge}"))?;
        f.hinge = self.hinge;
        Ok(f)
    }
    /// Set the trailing-edge gap (fraction of chord), blended over `blend` (x/c).
    pub fn with_te_gap(&self, gap: f64, blend: f64) -> Result<Foil, FoilError> {
        self.gdes(&format!("TGAP {gap} {blend}"), format!("TE gap {gap}"))
    }
    /// Scale the leading-edge radius by `factor`, blended to `blend` (x/c).
    pub fn with_le_radius(&self, factor: f64, blend: f64) -> Result<Foil, FoilError> {
        self.gdes(&format!("LERA {factor} {blend}"), format!("LE radius x{factor}"))
    }
    /// Scale thickness and camber.
    pub fn scaled_tc(&self, tfac: f64, cfac: f64) -> Result<Foil, FoilError> {
        self.gdes(&format!("TFAC {tfac} {cfac}"), format!("thickness x{tfac}, camber x{cfac}"))
    }
    /// Set absolute max thickness and camber (chord fractions); `None` keeps.
    pub fn with_tc(&self, t: Option<f64>, c: Option<f64>) -> Result<Foil, FoilError> {
        let p = self.props();
        let tf = t.map(|v| if p.thickness > 0.0 { v / p.thickness } else { 1.0 }).unwrap_or(1.0);
        let cf = c.map(|v| if p.camber.abs() > 1e-9 { v / p.camber } else { 1.0 }).unwrap_or(1.0);
        if c.is_some() && p.camber.abs() <= 1e-9 && c.unwrap().abs() > 1e-9 {
            return Err(FoilError::Geometry("cannot set the camber of a symmetric foil by scaling; use interpolation with a cambered foil".into()));
        }
        self.scaled_tc(tf, cf)
    }
    /// Move the points of max thickness / max camber to new chord positions.
    pub fn with_high_points(&self, x_thick: f64, x_camber: f64) -> Result<Foil, FoilError> {
        self.gdes(&format!("HIGH {x_thick} {x_camber}"), format!("max t at {x_thick}, max c at {x_camber}"))
    }
    /// Rotate so that the chord line is horizontal.
    pub fn derotated(&self) -> Result<Foil, FoilError> {
        self.gdes("DERO", "derotated".into())
    }
}

/// XFoil PANGEN paneling parameters.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Paneling {
    /// number of panel nodes (XFoil default 160, max 354)
    pub panels: usize,
    /// panel bunching by curvature (LE), 0..~2 (default 1)
    pub le_bunching: f64,
    /// TE/LE panel density ratio (default 0.15)
    pub te_bunching: f64,
    /// refined-area/LE panel density ratio (default 0.2)
    pub refined_ratio: f64,
}

impl Default for Paneling {
    fn default() -> Self {
        Paneling { panels: 160, le_bunching: 1.0, te_bunching: 0.15, refined_ratio: 0.2 }
    }
}

/// Where the solver's panel nodes come from (`AnalysisSettings::paneling`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NodeInput {
    /// XFoil PANE: PANGEN repanels the input coordinates by curvature
    Pangen(Paneling),
    /// XFoil LOAD + PCOP: the input coordinates are the panel nodes
    /// (repeated consecutive points are dropped, as PCOP does)
    Preserve,
}

impl Default for NodeInput {
    fn default() -> Self {
        NodeInput::Pangen(Paneling::default())
    }
}

impl NodeInput {
    /// PANGEN parameters, if any.
    pub fn pangen(&self) -> Option<&Paneling> {
        match self {
            NodeInput::Pangen(p) => Some(p),
            NodeInput::Preserve => None,
        }
    }
}

/// Serialized as the PANGEN parameters (as `Paneling` was before) or as the
/// string "preserve".
#[derive(Serialize, Deserialize)]
#[serde(untagged)]
enum NodeInputRepr {
    Pangen(Paneling),
    Mode(String),
}

impl Serialize for NodeInput {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        match self {
            NodeInput::Pangen(p) => NodeInputRepr::Pangen(*p),
            NodeInput::Preserve => NodeInputRepr::Mode("preserve".into()),
        }
        .serialize(s)
    }
}

impl<'de> Deserialize<'de> for NodeInput {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        match NodeInputRepr::deserialize(d)? {
            NodeInputRepr::Pangen(p) => Ok(NodeInput::Pangen(p)),
            NodeInputRepr::Mode(m) if m.eq_ignore_ascii_case("preserve") => Ok(NodeInput::Preserve),
            NodeInputRepr::Mode(m) => Err(serde::de::Error::custom(format!("unknown paneling '{m}' (expected PANGEN parameters or \"preserve\")"))),
        }
    }
}

/// FNV-1a (64 bit) over the IEEE bits of the node coordinates, in order
/// x1, y1, x2, y2, ...: identical nodes give the same hash on every platform.
pub fn node_hash(x: &[f64], y: &[f64]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for (a, b) in x.iter().zip(y) {
        for v in [a.to_bits(), b.to_bits()] {
            for byte in v.to_le_bytes() {
                h ^= byte as u64;
                h = h.wrapping_mul(0x100000001b3);
            }
        }
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn naca_props() {
        let f = Foil::naca("2412").unwrap();
        let p = f.props();
        assert!((p.thickness - 0.12).abs() < 2e-3, "{p:?}");
        assert!((p.camber - 0.02).abs() < 1e-3, "{p:?}");
        assert!((p.x_camber - 0.4).abs() < 0.03, "{p:?}");
        let (u, l) = f.surface_y(0.3);
        assert!(u > 0.07 && l < -0.03, "{u} {l}");
    }
    #[test]
    fn dat_roundtrip_and_ops() {
        let f = Foil::naca("0012").unwrap();
        let g = Foil::from_dat(&f.to_dat(), "x").unwrap();
        assert_eq!(g.n(), f.n());
        let t = f.scaled_tc(1.5, 1.0).unwrap();
        assert!((t.props().thickness - 0.18).abs() < 3e-3);
        let fl = f.with_flap(0.75, 0.5, 10.0).unwrap();
        // flap down -> TE moves down
        assert!(fl.y[0] < -0.03, "{}", fl.y[0]);
        let r = f.repaneled(&Paneling { panels: 120, ..Default::default() }).unwrap();
        assert_eq!(r.n(), 120);
        let m = Foil::interpolate(&f, &Foil::naca("4412").unwrap(), 0.5, "mix").unwrap();
        assert!((m.props().camber - 0.02).abs() < 2e-3);
    }
    #[test]
    fn lednicer() {
        let txt = "TEST\n3. 3.\n\n0 0\n0.5 0.1\n1 0\n\n0 0\n0.5 -0.1\n1 0\n";
        let r = Foil::from_dat(txt, "x");
        // too few points -> error but parsed without panic
        assert!(r.is_err());
    }
}
