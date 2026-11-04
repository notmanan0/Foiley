//! Run manifests for 3D polars (F016): plane and mesh hashes and the 2D
//! polar families a run used.  Computed once per run.

use crate::analysis::{mesh_options, Ctx, Method};
use crate::geom::Plane;
use crate::mesh::{build, Mesh};
use foiley_core::hash::{Hash32, Hasher};
use foiley_core::provenance::{hash_json, FamilyRef, RunManifest};
use std::collections::HashMap;

/// Hash of a plane's geometry (its full serialized definition).
pub fn plane_hash(p: &Plane) -> Hash32 {
    hash_json("foiley.plane.v1", p)
}

/// Hash of a mesh: panel nodes, kinds and owners, and the strips' foils,
/// surfaces and flap angles.
pub fn mesh_hash(m: &Mesh) -> Hash32 {
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
        h.u64(s.wing as u64).u64(s.left as u64).u64(s.surf as u64).str(&s.foil_a).str(&s.foil_b).f64(s.tau).f64(s.flap_deg);
    }
    h.finish()
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
