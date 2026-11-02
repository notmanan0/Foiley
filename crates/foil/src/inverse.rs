//! Inverse design (XFoil QDES mixed-inverse and MDES full-inverse) as plain
//! functions: get the surface-speed distribution, modify it, and solve for
//! the airfoil that produces it.  Runs the XFoil engine's own routines.

use crate::foil::{Foil, Paneling};
use crate::FoilError;
use foiley::cli::{Input, Session};
use serde::{Deserialize, Serialize};
use std::panic::{catch_unwind, AssertUnwindSafe};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum InverseMode {
    /// QDES: modify a segment of the surface; the rest is kept
    #[default]
    Mixed,
    /// MDES: conformal-mapping full inverse over the whole surface
    Full,
}

impl InverseMode {
    pub fn label(self) -> &'static str {
        match self {
            InverseMode::Mixed => "Mixed inverse (QDES): change one segment",
            InverseMode::Full => "Full inverse (MDES): whole surface",
        }
    }
}

/// Inviscid surface speed Q/Vinf at the design points.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SpeedDist {
    /// chordwise position of each point (x/c)
    pub x: Vec<f64>,
    /// normalised arc length (0 upper TE .. 1 lower TE)
    pub s: Vec<f64>,
    /// speed (positive on the upper surface in XFoil's convention)
    pub q: Vec<f64>,
    /// arc position of the leading edge
    pub s_le: f64,
    pub alpha: f64,
    pub cl: f64,
}

impl SpeedDist {
    pub fn is_upper(&self, i: usize) -> bool {
        self.s[i] < self.s_le
    }
}

/// A change to the target speed distribution.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Modification {
    pub x_from: f64,
    pub x_to: f64,
    pub upper: bool,
    pub kind: ModKind,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum ModKind {
    /// multiply the speed by a factor (blended at the segment ends)
    Scale(f64),
    /// add a smooth bump of `amplitude` (Q/Vinf) centred in the segment
    Bump(f64),
    /// smooth the speed in the segment (moving average)
    Smooth,
}

fn session(foil: &Foil, pan: &Paneling, script: &str) -> Result<Session, FoilError> {
    let xf = foil.engine(pan)?;
    let mut s = Session::new(Input::from_string(script));
    s.xf = xf;
    s.xf.verbose = false;
    Ok(s)
}

fn run(s: &mut Session) -> Result<(), FoilError> {
    let _q = foiley::quiet_scope();
    catch_unwind(AssertUnwindSafe(|| s.run(None))).map_err(|_| FoilError::Geometry("the inverse solver failed for this target (try a gentler modification)".into()))
}

fn enter(foil: &Foil, alpha: f64, mode: InverseMode, pan: &Paneling) -> Result<Session, FoilError> {
    let menu = match mode {
        InverseMode::Mixed => "QDES",
        InverseMode::Full => "MDES",
    };
    let mut s = session(foil, pan, &format!("OPER\nALFA {alpha}\n\n{menu}\n\nQUIT\n"))?;
    run(&mut s)?;
    if s.qdes_state.nsp < 10 {
        return Err(FoilError::Geometry("could not initialise the speed distribution".into()));
    }
    Ok(s)
}

fn dist_of(s: &Session, alpha: f64) -> SpeedDist {
    let q = &s.qdes_state;
    let n = q.nsp;
    SpeedDist { x: q.xspoc[1..=n].to_vec(), s: q.sspec[1..=n].to_vec(), q: q.qspec[1][1..=n].to_vec(), s_le: q.ssple, alpha, cl: s.xf.cl }
}

/// Current inviscid speed distribution of `foil` at `alpha` (degrees).
pub fn speed_distribution(foil: &Foil, alpha: f64, mode: InverseMode, pan: &Paneling) -> Result<SpeedDist, FoilError> {
    let s = enter(foil, alpha, mode, pan)?;
    Ok(dist_of(&s, alpha))
}

/// Apply a modification; returns the new target and the modified index range.
pub fn modify(d: &SpeedDist, m: &Modification) -> (Vec<f64>, Option<(usize, usize)>) {
    let (lo, hi) = (m.x_from.min(m.x_to), m.x_from.max(m.x_to));
    let idx: Vec<usize> = (0..d.q.len()).filter(|&i| d.is_upper(i) == m.upper && d.x[i] >= lo && d.x[i] <= hi).collect();
    let mut q = d.q.clone();
    if idx.len() < 2 {
        return (q, None);
    }
    let (i1, i2) = (*idx.first().unwrap(), *idx.last().unwrap());
    let w = |x: f64| {
        // smooth window: 0 at the ends, 1 in the middle
        let t = ((x - lo) / (hi - lo).max(1e-12)).clamp(0.0, 1.0);
        (std::f64::consts::PI * t).sin().powi(2)
    };
    match m.kind {
        ModKind::Scale(f) => {
            for &i in &idx {
                q[i] *= 1.0 + (f - 1.0) * w(d.x[i]);
            }
        }
        ModKind::Bump(a) => {
            for &i in &idx {
                q[i] += a * w(d.x[i]) * q[i].signum();
            }
        }
        ModKind::Smooth => {
            let src = q.clone();
            for _ in 0..3 {
                let cur = q.clone();
                for &i in &idx[1..idx.len() - 1] {
                    q[i] = 0.25 * cur[i - 1] + 0.5 * cur[i] + 0.25 * cur[i + 1];
                }
            }
            let _ = src;
        }
    }
    (q, Some((i1.min(i2) + 1, i1.max(i2) + 1)))
}

/// Solve for the airfoil producing `target` (same points as the
/// distribution).  Mixed mode only changes the segment `seg` (1-based
/// indices from [`modify`]).
pub fn solve(foil: &Foil, alpha: f64, mode: InverseMode, pan: &Paneling, target: &[f64], seg: Option<(usize, usize)>, iterations: usize) -> Result<(Foil, SpeedDist), FoilError> {
    let mut s = enter(foil, alpha, mode, pan)?;
    let n = s.qdes_state.nsp;
    if target.len() != n {
        return Err(FoilError::Geometry(format!("target has {} points, expected {n}", target.len())));
    }
    if target.iter().any(|v| !v.is_finite()) {
        return Err(FoilError::Geometry("target speed contains NaN".into()));
    }
    for i in 1..=n {
        s.qdes_state.qspec[1][i] = target[i - 1];
    }
    let script = match mode {
        InverseMode::Mixed => {
            let (i1, i2) = seg.ok_or_else(|| FoilError::Geometry("mixed inverse needs a modified segment".into()))?;
            format!("QDES\nMARK {i1} {i2}\nEXEC {}\n\nQUIT\n", iterations.max(1))
        }
        InverseMode::Full => "MDES\nEXEC\n\nQUIT\n".to_string(),
    };
    s.input = Input::from_string(&script);
    run(&mut s)?;
    let xf = &s.xf;
    let nn = xf.n;
    let x: Vec<f64> = xf.x[1..=nn].to_vec();
    let y: Vec<f64> = xf.y[1..=nn].to_vec();
    if x.iter().chain(&y).any(|v| !v.is_finite()) {
        return Err(FoilError::Geometry("inverse design diverged".into()));
    }
    let mut f = Foil::new(format!("{} inv", foil.name), x, y)?;
    f.note = format!("{} at alpha {alpha:.2}", mode.label());
    let d = dist_of(&s, alpha);
    Ok((f, d))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mixed_inverse_changes_segment() {
        let f = Foil::naca("0012").unwrap();
        let pan = Paneling::default();
        let d = speed_distribution(&f, 2.0, InverseMode::Mixed, &pan).unwrap();
        assert!(d.q.len() > 100);
        let m = Modification { x_from: 0.3, x_to: 0.8, upper: true, kind: ModKind::Scale(0.95) };
        let (t, seg) = modify(&d, &m);
        let (g, _) = solve(&f, 2.0, InverseMode::Mixed, &pan, &t, seg, 10).unwrap();
        // slower upper surface in the segment -> foil changes there
        let (yu0, _) = f.surface_y(0.55);
        let (yu1, _) = g.surface_y(0.55);
        assert!((yu1 - yu0).abs() > 1e-4, "{yu0} -> {yu1}");
        let (_, yl0) = f.surface_y(0.5);
        let (_, yl1) = g.surface_y(0.5);
        assert!((yl1 - yl0).abs() < 2e-3, "lower surface should be nearly unchanged");
    }
    #[test]
    fn full_inverse_runs() {
        let f = Foil::naca("2412").unwrap();
        let pan = Paneling::default();
        let d = speed_distribution(&f, 3.0, InverseMode::Full, &pan).unwrap();
        let m = Modification { x_from: 0.2, x_to: 0.7, upper: false, kind: ModKind::Bump(0.03) };
        let (t, _) = modify(&d, &m);
        let (g, _) = solve(&f, 3.0, InverseMode::Full, &pan, &t, None, 1).unwrap();
        assert!(g.n() > 50);
    }
}
