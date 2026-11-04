//! Triangulation of STEP faces.
//!
//! Interim (F022a): the previous behaviour on top of the new topology —
//! B-spline surfaces on a parameter grid over their whole domain, planar
//! faces as a fan of their outer loop; other faces are reported, not
//! triangulated.  F022c replaces this with trimmed, watertight tessellation.

use super::trim::{Model, StepOptions};
use super::{clean, orient, Surface, Tri};
use std::collections::BTreeMap;

pub(crate) fn legacy(text: &str, density: f64) -> Result<(Vec<Tri>, Vec<String>), String> {
    let model = Model::from_step(text, &StepOptions::default())?;
    let mut tris = Vec::new();
    let mut notes = model.notes.clone();
    let mut n_surf = 0;
    let mut skipped: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mut done_surfaces = std::collections::HashSet::new();
    for (fi, f) in model.faces.iter().enumerate() {
        match &f.surface {
            Surface::BSpline(s) => {
                if !done_surfaces.insert(f.surface_id) {
                    continue;
                }
                let nu = s.ctrl.len();
                let nv = s.ctrl.first().map(|r| r.len()).unwrap_or(0);
                let seg = |n: usize| (((n.max(2) - 1) as f64 * density).round() as usize).max(2);
                tris.extend(s.tessellate(seg(nu), seg(nv)));
                n_surf += 1;
            }
            Surface::Plane(_) => {
                let Ok(d) = model.face_domain(fi) else { continue };
                let Some(outer) = d.loops.iter().max_by(|a, b| a.signed_area().abs().total_cmp(&b.signed_area().abs())) else { continue };
                let poly: Vec<_> = outer.ids.iter().map(|&i| model.vertex(&d, i)).collect();
                for k in 1..poly.len().saturating_sub(1) {
                    tris.push([poly[0], poly[k], poly[k + 1]]);
                }
                n_surf += 1;
            }
            s => *skipped.entry(s.kind()).or_default() += 1,
        }
    }
    for (k, n) in skipped {
        notes.push(format!("{n} face(s) on {k} surfaces are not tessellated and were skipped"));
    }
    let mut t = clean(tris);
    if t.is_empty() {
        return Err("no supported surfaces found (B-spline/NURBS surfaces and planar faces are supported)".into());
    }
    orient(&mut t);
    notes.insert(0, format!("{n_surf} surface(s), unit scale {} m, {} triangles", model.length_scale, t.len()));
    Ok((t, notes))
}
