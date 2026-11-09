//! Project model and file exchange.

use foiley_foil::{Foil, FoilPolar};
use foiley_plane::{FoilDb, Plane, PlanePolar};
use serde::{Deserialize, Serialize};
use std::path::Path;

/// Project file format.  2: polar provenance (settings snapshot, geometry
/// hash, run manifest) and explicit families; 1 files migrate on load.
pub const FORMAT_VERSION: u32 = 2;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Project {
    #[serde(default)]
    pub format: u32,
    #[serde(default)]
    pub name: String,
    #[serde(default)]
    pub foils: Vec<Foil>,
    #[serde(default)]
    pub foil_polars: Vec<FoilPolar>,
    #[serde(default)]
    pub planes: Vec<Plane>,
    #[serde(default)]
    pub plane_polars: Vec<PlanePolar>,
}

impl Project {
    pub fn new(name: &str) -> Self {
        Project { format: FORMAT_VERSION, name: name.into(), ..Default::default() }
    }

    /// A starter project: two NACA foils and the default glider.
    pub fn starter() -> Self {
        let mut p = Project::new("Untitled");
        for c in ["2412", "0009"] {
            if let Ok(f) = Foil::naca(c) {
                p.foils.push(f);
            }
        }
        p.planes.push(Plane::default_plane("NACA 2412", "NACA 0009"));
        p
    }

    pub fn load(path: &Path) -> Result<Project, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let mut p: Project = serde_json::from_str(&text).map_err(|e| format!("{}: not a Foiley project ({e})", path.display()))?;
        if p.format > FORMAT_VERSION {
            return Err(format!("{} was written by a newer Foiley (format {})", path.display(), p.format));
        }
        if p.format < 2 {
            p.migrate_v1();
        }
        p.fill_polar_hashes();
        p.format = FORMAT_VERSION;
        Ok(p)
    }

    /// Format 1 -> 2: every 2D polar without a settings snapshot or manifest
    /// becomes a legacy member of the family its stored Mach, Ncrit and trip
    /// fields name; its geometry hash is that of the foil it is stored with.
    fn migrate_v1(&mut self) {
        for q in &mut self.foil_polars {
            if q.settings.is_none() && q.manifest.is_none() {
                q.legacy = true;
            }
        }
    }

    /// Give every 2D polar without a geometry hash the hash of its foil
    /// (legacy and imported polars), so later geometry edits are detected.
    fn fill_polar_hashes(&mut self) {
        if self.foil_polars.iter().all(|q| q.foil_hash.is_some()) { return; }
        let hashes: Vec<(String, foiley_core::hash::Hash32)> = self.foils.iter().map(|f| (f.name.clone(), f.geom_hash())).collect();
        for q in &mut self.foil_polars {
            if q.foil_hash.is_none() {
                q.foil_hash = hashes.iter().find(|h| h.0 == q.foil).map(|h| h.1);
            }
        }
    }

    /// Drop the 2D polars that no longer describe their foil: computed on
    /// other coordinates (hash mismatch), or deflected about a hinge that
    /// has moved.  Returns how many were dropped.
    pub fn drop_stale_polars(&mut self, foil: &str, old_hinge: Option<(f64, f64)>) -> usize {
        let Some(f) = self.foil(foil) else { return 0 };
        let (h, hinge) = (f.geom_hash(), f.hinge);
        let n = self.foil_polars.len();
        self.foil_polars.retain(|q| q.foil != foil || (q.foil_hash.is_none_or(|qh| qh == h) && (q.flap_deg == 0.0 || hinge == old_hinge)));
        n - self.foil_polars.len()
    }

    /// Save atomically (write to a temp file, then rename).  Always writes
    /// the current format.
    pub fn save(&self, path: &Path) -> Result<(), String> {
        let text = if self.format == FORMAT_VERSION {
            serde_json::to_string(self)
        } else {
            serde_json::to_string(&Project { format: FORMAT_VERSION, ..self.clone() })
        }
        .map_err(|e| e.to_string())?;
        let tmp = path.with_extension("foiley.tmp");
        std::fs::write(&tmp, text).map_err(|e| format!("{}: {e}", tmp.display()))?;
        std::fs::rename(&tmp, path).map_err(|e| format!("{}: {e}", path.display()))
    }

    pub fn foil(&self, name: &str) -> Option<&Foil> {
        self.foils.iter().find(|f| f.name == name)
    }

    /// Unique name derived from `base`.
    pub fn unique_foil_name(&self, base: &str) -> String {
        if self.foil(base).is_none() {
            return base.to_string();
        }
        (2..).map(|i| format!("{base} ({i})")).find(|n| self.foil(n).is_none()).unwrap()
    }

    /// Add or replace a foil (by name).  Replacing a foil's coordinates
    /// drops the polars computed on the old ones (geometry hash mismatch),
    /// and moving its hinge drops its deflected polars; a renamed note
    /// keeps them.
    pub fn upsert_foil(&mut self, f: Foil) {
        if let Some(i) = self.foils.iter().position(|g| g.name == f.name) {
            let old_hinge = self.foils[i].hinge;
            let name = f.name.clone();
            self.fill_polar_hashes();
            self.foils[i] = f;
            self.drop_stale_polars(&name, old_hinge);
        } else {
            self.foils.push(f);
        }
    }

    pub fn remove_foil(&mut self, name: &str) {
        self.foils.retain(|f| f.name != name);
        self.foil_polars.retain(|p| p.foil != name);
    }

    /// Add or merge a 2D polar (same geometry and analysis settings). A
    /// polar without a geometry hash gets its foil's.
    pub fn add_foil_polar(&mut self, p: FoilPolar) -> (String, String) {
        // a deflected copy from an automatic batch goes under its base foil
        let mut p = crate::auto::adopt(p, |base| self.foil(base).map(|f| f.geom_hash()));
        if p.foil_hash.is_none() {
            p.foil_hash = self.foil(&p.foil).map(|f| f.geom_hash());
        }
        self.fill_polar_hashes();
        let key = p.key();
        let identity = p.identity();
        if let Some(q) = self.foil_polars.iter_mut().find(|q| q.foil == p.foil && q.re == p.re && q.identity() == identity) {
            for pt in p.points {
                q.insert(pt);
            }
            q.failed = p.failed;
        } else {
            self.foil_polars.push(p);
        }
        key
    }

    pub fn add_plane_polar(&mut self, p: PlanePolar) {
        if let Some(q) = self.plane_polars.iter_mut().find(|q| q.plane == p.plane && q.name == p.name) {
            *q = p;
        } else {
            self.plane_polars.push(p);
        }
    }

    /// Add a fuselage from an STL or STEP file to plane `plane`.
    /// `scale` converts file units to metres (None = automatic: STEP files
    /// carry their units, STL files are assumed to be in metres).
    pub fn import_fuselage(&mut self, path: &Path, plane: usize, scale: Option<f64>, target_triangles: usize) -> Result<String, String> {
        let ext = path.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default();
        let (mut tris, mut notes) = match ext.as_str() {
            "stl" => (foiley_plane::cad::read_stl(&std::fs::read(path).map_err(|e| e.to_string())?)?, vec![]),
            "step" | "stp" => foiley_plane::cad::read_step_target(&std::fs::read_to_string(path).map_err(|e| e.to_string())?, target_triangles)?,
            _ => return Err("fuselage files must be .stl, .step or .stp".into()),
        };
        if let Some(s) = scale {
            for t in &mut tris {
                for p in t.iter_mut() {
                    *p = *p * s;
                }
            }
            notes.push(format!("scaled by {s}"));
        }
        let n = tris.len();
        let name = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_else(|| "Fuselage".into());
        let pl = self.planes.get_mut(plane).ok_or("no such plane")?;
        let f = foiley_plane::Fuse::from_triangles(&name, tris, 0.0);
        let len = f.length();
        pl.fuses.clear();
        pl.fuses.push(f);
        let mut msg = format!("fuselage '{name}': {n} triangles, length {len:.3} m");
        if n > 6000 {
            msg.push_str(" (large: panel analyses will be slow; consider a coarser export)");
        }
        for nt in notes {
            msg.push_str(&format!("; {nt}"));
        }
        Ok(msg)
    }

    /// Foil database (shapes + polars) for the 3D analyses.
    pub fn foil_db(&self) -> FoilDb {
        let mut db = FoilDb::new();
        for f in &self.foils {
            db.add_foil(f);
        }
        for p in &self.foil_polars {
            db.add_polar(p);
        }
        db
    }

    /// Import any supported file: .dat (foil), .xml (plane), .txt/.pol
    /// (XFoil polar), .foiley/.json (project, merged).  Returns a summary.
    pub fn import(&mut self, path: &Path) -> Result<String, String> {
        let ext = path.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default();
        if matches!(ext.as_str(), "xfl" | "fl5" | "plr") {
                let imp = crate::legacy_io::read_any(path)?;
                let (nf, np, npl) = (imp.foils.len(), imp.polars.len(), imp.planes.len());
                let mut skipped = 0;
                for f in imp.foils {
                    // keep names: planes and polars reference foils by name
                    if self.foil(&f.name).is_none() {
                        self.foils.push(f);
                    } else {
                        skipped += 1;
                    }
                }
                for q in imp.polars {
                    self.add_foil_polar(q);
                }
                for mut pl in imp.planes {
                    if self.planes.iter().any(|p| p.name == pl.name) {
                        pl.name = format!("{} (imported)", pl.name);
                    }
                    self.planes.push(pl);
                }
                let mut msg = format!("{nf} foils, {np} polars, {npl} planes from {}", path.file_name().map(|f| f.to_string_lossy().to_string()).unwrap_or_default());
                if skipped > 0 {
                    msg.push_str(&format!(" ({skipped} foils already present were kept)"));
                }
                for w in imp.warnings.iter().take(5) {
                    msg.push_str(&format!("; {w}"));
                }
            return Ok(msg);
        }
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
        match ext.as_str() {
            "xml" => {
                let mut pl = foiley_plane::xml::read_plane(&text).map_err(|e| e.to_string())?;
                // load foils referenced as files next to the xml
                let dir = path.parent().unwrap_or(Path::new("."));
                let mut loaded = Vec::new();
                for name in pl.foil_names() {
                    if self.foil(&name).is_none() {
                        let f = dir.join(format!("{name}.dat"));
                        if let Ok(foil) = Foil::load(&f) {
                            let mut foil = foil;
                            foil.name = name.clone();
                            self.foils.push(foil);
                            loaded.push(name);
                        }
                    }
                }
                // NACA designations can be generated
                for name in pl.foil_names() {
                    if self.foil(&name).is_none() {
                        let l = name.to_lowercase().replace([' ', '-', '_'], "");
                        if let Some(Ok(mut f)) = l.strip_prefix("naca").map(Foil::naca) {
                            f.name = name.clone();
                            self.foils.push(f);
                            loaded.push(format!("{name} (generated)"));
                        }
                    }
                }
                let missing: Vec<String> = pl.foil_names().into_iter().filter(|n| self.foil(n).is_none()).collect();
                if self.planes.iter().any(|p| p.name == pl.name) {
                    pl.name = format!("{} (imported)", pl.name);
                }
                let mut msg = format!("plane '{}' with {} wing(s)", pl.name, pl.wings.len());
                if !loaded.is_empty() {
                    msg.push_str(&format!(", loaded foils {}", loaded.join(", ")));
                }
                if !missing.is_empty() {
                    msg.push_str(&format!("; foils still missing: {} (load their .dat files)", missing.join(", ")));
                }
                self.planes.push(pl);
                Ok(msg)
            }
            "txt" | "pol" if text.contains("alpha") && text.contains("CL") => {
                let p = FoilPolar::from_xfoil_text(&text, &stem).ok_or("no polar data found")?;
                let msg = format!("polar {} for foil '{}' ({} points)", p.name, p.foil, p.points.len());
                self.add_foil_polar(p);
                Ok(msg)
            }
            "foiley" | "json" => {
                let other = Project::load(path)?;
                let msg = format!("{} foils, {} polars, {} planes, {} plane polars", other.foils.len(), other.foil_polars.len(), other.planes.len(), other.plane_polars.len());
                let mut foil_names = std::collections::HashMap::new();
                for f in other.foils {
                    let mut f = f;
                    let renamed = self.unique_foil_name(&f.name);
                    foil_names.insert(f.name.clone(), renamed.clone());
                    f.name = renamed;
                    self.foils.push(f);
                }
                for mut p in other.foil_polars {
                    if let Some(name) = foil_names.get(&p.foil) { p.foil = name.clone(); }
                    self.add_foil_polar(p);
                }
                let mut plane_names = std::collections::HashMap::new();
                for mut plane in other.planes {
                    for section in plane.wings.iter_mut().flat_map(|w| &mut w.sections) {
                        for foil in [&mut section.foil_left, &mut section.foil_right] {
                            if let Some(name) = foil_names.get(foil) { *foil = name.clone(); }
                        }
                    }
                    let base = plane.name.clone();
                    if self.planes.iter().any(|p| p.name == base) {
                        plane.name = (2..).map(|i| format!("{base} ({i})")).find(|name| self.planes.iter().all(|p| p.name != *name)).unwrap();
                    }
                    plane_names.insert(base, plane.name.clone());
                    self.planes.push(plane);
                }
                for mut p in other.plane_polars {
                    if let Some(name) = plane_names.get(&p.plane) { p.plane = name.clone(); }
                    self.add_plane_polar(p);
                }
                Ok(msg)
            }
            _ => {
                let mut f = Foil::load(path).map_err(|e| e.to_string())?;
                f.name = self.unique_foil_name(&f.name);
                let msg = format!("foil '{}' ({} points)", f.name, f.n());
                self.foils.push(f);
                Ok(msg)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn save_load_roundtrip() {
        let p = Project::starter();
        let dir = std::env::temp_dir().join(format!("foiley_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let f = dir.join("p.foiley");
        p.save(&f).unwrap();
        let q = Project::load(&f).unwrap();
        assert_eq!(q.foils.len(), 2);
        assert_eq!(q.planes.len(), 1);
        std::fs::remove_dir_all(&dir).ok();
    }
}
