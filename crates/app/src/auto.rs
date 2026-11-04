//! Automatic 2D polar meshes for 3D viscous analyses.

use foiley_foil::batch::BatchCase;
use foiley_foil::{AnalysisSettings, Foil, Sweep};
use foiley_plane::analysis::{PolarKind, G};
use foiley_plane::{FoilDb, Plane, PlanePolarDef};

/// Reynolds range spanned by the wings of `plane` for this analysis.
pub fn reynolds_range(plane: &Plane, db: &FoilDb, def: &PlanePolarDef) -> (f64, f64) {
    let nu = def.fluid.nu();
    let chords: Vec<f64> = plane.wings.iter().flat_map(|w| w.sections.iter().map(|s| s.chord)).collect();
    let cmin = chords.iter().cloned().fold(f64::MAX, f64::min).max(1e-3);
    let cmax = chords.iter().cloned().fold(0.0, f64::max).max(cmin);
    let (vmin, vmax) = match def.kind {
        PolarKind::FixedAlpha => {
            let v = def.range.values();
            (v.iter().cloned().fold(f64::MAX, f64::min), v.iter().cloned().fold(0.0, f64::max))
        }
        PolarKind::FixedLift | PolarKind::Glide | PolarKind::Stability => {
            let m = def.inertia.map(|i| i.mass).unwrap_or_else(|| plane.inertia(db).mass).max(1e-3);
            let s = plane.reference(true).0.max(1e-6);
            let v = |cl: f64| (2.0 * m * G / (def.fluid.rho * s * cl)).sqrt();
            (v(1.4), v(0.1))
        }
        PolarKind::Control => {
            let r = def.ranges.iter().find(|r| r.target == foiley_plane::analysis::ControlTarget::Velocity);
            match r {
                Some(r) => (r.min.min(r.max), r.min.max(r.max)),
                None => (def.velocity, def.velocity),
            }
        }
        _ => (def.velocity, def.velocity),
    };
    let lo = (0.6 * cmin * vmin.max(0.1) / nu).max(2e4);
    let hi = (1.5 * cmax * vmax.max(0.1) / nu).max(lo * 2.0);
    (lo, hi)
}

/// Foils of the plane that have no type-1 viscous polars yet.
pub fn foils_needing_polars<'a>(plane: &Plane, foils: &'a [Foil], db: &FoilDb) -> Vec<&'a Foil> {
    let names = plane.foil_names();
    foils.iter().filter(|f| names.contains(&f.name) && !db.has_polars(&f.name)).collect()
}

/// Cases and sweep for an automatic polar mesh.
pub fn polar_mesh_cases(plane: &Plane, db: &FoilDb, def: &PlanePolarDef, count: usize) -> (Vec<BatchCase>, AnalysisSettings, Sweep) {
    let (lo, hi) = reynolds_range(plane, db, def);
    let res = foiley_foil::batch::re_list_geometric(lo, hi, count.max(2));
    let mach = 0.0;
    let cases = res.into_iter().map(|re| BatchCase { re, mach, ncrit: 9.0 }).collect();
    (cases, AnalysisSettings::default(), Sweep::Alpha { from: -10.0, to: 18.0, step: 0.5 })
}
