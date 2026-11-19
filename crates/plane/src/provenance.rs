//! Run manifests for 3D polars (F016): plane and mesh hashes and the 2D
//! polar families a run used.  Computed once per run.

use crate::analysis::{mesh_options, Ctx, Method};
use crate::geom::Plane;
use crate::mesh::{build, Mesh};
use foiley_core::hash::{Hash32, Hasher};
use foiley_core::provenance::{hash_json, FamilyRef, RunManifest};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};

/// Hash of a plane's geometry (its full serialized definition).
pub fn plane_hash(p: &Plane) -> Hash32 {
    hash_json("foiley.plane.v1", p)
}

/// Hash of a mesh: panel nodes, kinds and owners, and the strips' foils,
/// surfaces and flap angles.
pub fn mesh_hash(m: &Mesh) -> Hash32 {
    mesh_hash_as(m, &BTreeMap::new())
}

/// [`mesh_hash`] with the strips' foils under other names: `foils` maps a
/// current foil name to the name to hash, as [`RenamedLabels::foils`] does.
pub fn mesh_hash_as(m: &Mesh, foils: &BTreeMap<String, String>) -> Hash32 {
    fn name<'a>(foils: &'a BTreeMap<String, String>, foil: &'a str) -> &'a str {
        foils.get(foil).map_or(foil, String::as_str)
    }
    let mut h = Hasher::new();
    h.tag("foiley.mesh.v1").u64(m.panels.len() as u64).u64(m.thick as u64);
    for p in &m.panels {
        for n in &p.nodes {
            h.f64(n.x).f64(n.y).f64(n.z);
        }
        h.u64(p.kind as u64).u64(p.wing as u64).u64(p.strip.map(|s| s as u64 + 1).unwrap_or(0)).u64(p.flap.map(|s| s as u64 + 1).unwrap_or(0));
    }
    h.u64(m.strips.len() as u64);
    for s in &m.strips {
        h.u64(s.wing as u64).u64(s.left as u64).u64(s.surf as u64).str(name(foils, &s.foil_a)).str(name(foils, &s.foil_b)).f64(s.tau).f64(s.flap_deg);
    }
    h.finish()
}

/// Labels of a plane renamed after a 3D run, mapped back to the names the
/// run's manifest hashed.  A rename of the plane, its description, a wing, a
/// body, a point mass or a foil changes no geometry, so matching a polar to
/// its plane puts these names back before comparing hashes.
///
/// Wing, body and mass names are kept by index; foils by name, because a
/// section refers to its foil by name.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct RenamedLabels {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plane: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    /// wing index -> name
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub wings: BTreeMap<usize, String>,
    /// body (`Plane::fuses`) index -> name
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub bodies: BTreeMap<usize, String>,
    /// point mass index -> tag
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub masses: BTreeMap<usize, String>,
    /// current foil name -> name at run time
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub foils: BTreeMap<String, String>,
}

impl RenamedLabels {
    pub fn is_empty(&self) -> bool {
        *self == RenamedLabels::default()
    }

    /// Record a rename of the plane itself.
    pub fn follow_plane_rename(&mut self, before: &str, after: &str) {
        track(&mut self.plane, before, after);
    }

    /// Record the labels an edit of the plane from `before` to `after`
    /// renamed.  Geometry is not tracked: a geometry edit still makes the
    /// restored plane differ from the one the manifest hashed.
    pub fn follow_plane_edit(&mut self, before: &Plane, after: &Plane) {
        track(&mut self.plane, &before.name, &after.name);
        track(&mut self.description, &before.description, &after.description);
        // Labels are kept by index, which an added or removed item shifts.
        // The geometry then differs anyway, so such an edit is not tracked.
        if before.wings.len() == after.wings.len() {
            for (i, (b, a)) in before.wings.iter().zip(&after.wings).enumerate() {
                track_at(&mut self.wings, i, &b.name, &a.name);
            }
        }
        if before.fuses.len() == after.fuses.len() {
            for (i, (b, a)) in before.fuses.iter().zip(&after.fuses).enumerate() {
                track_at(&mut self.bodies, i, &b.name, &a.name);
            }
        }
        if before.point_masses.len() == after.point_masses.len() {
            for (i, (b, a)) in before.point_masses.iter().zip(&after.point_masses).enumerate() {
                track_at(&mut self.masses, i, &b.tag, &a.tag);
            }
        }
    }

    /// Record foil renames made at the same time (old name -> new name):
    /// every reference to a renamed foil now uses its new name.
    pub fn follow_foil_renames(&mut self, renames: &HashMap<String, String>) {
        let old = std::mem::take(&mut self.foils);
        let mut foils = BTreeMap::new();
        for (from, to) in renames {
            foils.insert(to.clone(), old.get(from).unwrap_or(from).clone());
        }
        for (current, original) in old {
            if !renames.contains_key(&current) {
                foils.entry(current).or_insert(original);
            }
        }
        foils.retain(|current, original| current != original);
        self.foils = foils;
    }

    /// `plane` with these labels put back, or None when they do not fit it:
    /// a renamed wing, body or mass is gone, or two foils the plane uses now
    /// map to one name (a section was switched to a different foil).
    pub fn restore(&self, plane: &Plane) -> Option<Plane> {
        let mut p = plane.clone();
        if let Some(name) = &self.plane {
            p.name.clone_from(name);
        }
        if let Some(description) = &self.description {
            p.description.clone_from(description);
        }
        for (&i, name) in &self.wings {
            p.wings.get_mut(i)?.name.clone_from(name);
        }
        for (&i, name) in &self.bodies {
            p.fuses.get_mut(i)?.name.clone_from(name);
        }
        for (&i, tag) in &self.masses {
            p.point_masses.get_mut(i)?.tag.clone_from(tag);
        }
        if !self.foils.is_empty() {
            // run-time name -> the one current name that maps to it
            let mut source: HashMap<String, String> = HashMap::new();
            for section in p.wings.iter_mut().flat_map(|w| &mut w.sections) {
                for foil in [&mut section.foil_left, &mut section.foil_right] {
                    let original = self.foils.get(foil.as_str()).unwrap_or(foil).clone();
                    if *source.entry(original.clone()).or_insert_with(|| foil.clone()) != *foil {
                        return None;
                    }
                    *foil = original;
                }
            }
        }
        Some(p)
    }
}

/// Keep the first label of a slot renamed from `before` to `after`; forget it
/// when it is renamed back.
fn track(slot: &mut Option<String>, before: &str, after: &str) {
    if before != after {
        let original = slot.take().unwrap_or_else(|| before.to_owned());
        if original != after {
            *slot = Some(original);
        }
    }
}

fn track_at(slots: &mut BTreeMap<usize, String>, i: usize, before: &str, after: &str) {
    let mut slot = slots.remove(&i);
    track(&mut slot, before, after);
    if let Some(original) = slot {
        slots.insert(i, original);
    }
}

/// The manifest of a 3D run: settings, plane, base-configuration mesh (not
/// for the lifting line) and the resolved polar families.
pub fn run_manifest(ctx: &Ctx) -> RunManifest {
    let mut m = RunManifest::new(hash_json("foiley.def3d.v1", ctx.def), plane_hash(ctx.plane));
    if ctx.def.method != Method::Llt {
        if let Ok(mesh) = build(ctx.plane, ctx.db, &mesh_options(ctx.def, &HashMap::new(), &HashMap::new())) {
            m.mesh_hash = Some(mesh_hash(&mesh));
        }
    }
    for (foil, id) in &ctx.foil_family {
        let Some(id) = *id else { continue };
        let set = ctx.flap_sets.get(id.0 as usize).filter(|s| !s.is_empty()).cloned().unwrap_or_else(|| vec![(0.0, id)]);
        for (flap, fid) in set {
            let k = ctx.db.family(fid).key;
            m.families.push(FamilyRef { foil: foil.clone(), family_id: fid.0, key_hash: k.hash(), foil_geom_hash: k.foil_geom_hash, flap_deg: flap });
        }
    }
    m
}
