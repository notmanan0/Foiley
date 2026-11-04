//! F016: provenance and project format v2.
use foiley_app::project::FORMAT_VERSION;
use foiley_app::Project;
use foiley_core::provenance::hash_json;
use foiley_core::Progress;
use foiley_foil::{analyze, AnalysisSettings, Foil, FoilPolar, Sweep};
use foiley_plane::analysis::Range;
use foiley_plane::provenance::{mesh_hash, plane_hash};
use foiley_plane::{run_polar, Method, Plane, PlanePolarDef};
use std::path::{Path, PathBuf};

fn projects_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/projects")
}

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("foiley_f016_{}_{name}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn json(p: &Project) -> serde_json::Value {
    serde_json::to_value(p).unwrap()
}

fn polar_of(f: &Foil, re: f64) -> FoilPolar {
    let st = AnalysisSettings { re, ..Default::default() };
    analyze(f, &st, &Sweep::Alpha { from: -2.0, to: 6.0, step: 2.0 }, &Progress::none()).unwrap().0
}

/// A small project with computed 2D polars and one 3D polar.
fn computed_project() -> Project {
    let mut p = Project::starter();
    let foils = p.foils.clone();
    for f in &foils {
        for re in [2e5, 6e5] {
            p.add_foil_polar(polar_of(f, re));
        }
    }
    let db = p.foil_db();
    let def = PlanePolarDef { method: Method::Vlm2, velocity: 15.0, range: Range { from: 2.0, to: 4.0, step: 2.0 }, ..Default::default() };
    let pp = run_polar(&p.planes[0], &db, &def, &Progress::none(), &|_| {}).unwrap();
    p.add_plane_polar(pp);
    p
}

#[test]
fn project_roundtrip() {
    let dir = projects_dir();
    let inv = std::fs::read_to_string(dir.join("INVENTORY.txt")).expect("tests/data/projects/INVENTORY.txt");
    let files: Vec<&str> = inv.lines().map(str::trim).filter(|l| !l.is_empty() && !l.starts_with('#')).map(|l| l.split_whitespace().next().unwrap()).collect();
    assert!(files.len() >= 6, "INVENTORY.txt lists {} files (need at least 6)", files.len());
    let out = tmp("roundtrip");
    let mut v1 = 0;
    for f in &files {
        let path = dir.join(f);
        assert!(path.is_file(), "INVENTORY.txt names a missing file: {f}");
        let raw: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
        let was = raw.get("format").and_then(|v| v.as_u64()).unwrap_or(0);
        v1 += (was < 2) as usize;
        let a = Project::load(&path).unwrap_or_else(|e| panic!("{f}: {e}"));
        assert_eq!(a.format, FORMAT_VERSION);
        if was < 2 {
            assert!(a.foil_polars.iter().all(|q| q.legacy && q.foil_hash.is_some()), "{f}: v1 polars not migrated");
        }
        let saved = out.join("x.foiley");
        a.save(&saved).unwrap();
        let text = std::fs::read_to_string(&saved).unwrap();
        assert!(text.contains(&format!("\"format\":{FORMAT_VERSION}")), "{f}: saved in another format");
        let b = Project::load(&saved).unwrap();
        assert_eq!(json(&a), json(&b), "{f}: load-save-load changed the project");
        assert_eq!((a.foils.len(), a.foil_polars.len(), a.planes.len(), a.plane_polars.len()), (b.foils.len(), b.foil_polars.len(), b.planes.len(), b.plane_polars.len()));
    }
    assert!(v1 >= 2, "the inventory needs at least two format-1 files (has {v1})");
    std::fs::remove_dir_all(&out).ok();
}

#[test]
fn f016_manifest_present() {
    let f = Foil::naca("2412").unwrap();
    let q = polar_of(&f, 3e5);
    let m = q.manifest.as_ref().expect("2D polar manifest");
    assert!(!m.git_commit.is_empty() && !m.foiley_version.is_empty() && !m.build_profile.is_empty() && !m.started_at.is_empty());
    assert!(m.thread_config.rayon_threads >= 1);
    assert_eq!(m.geometry_hash, f.geom_hash());
    assert_eq!(q.foil_hash, Some(f.geom_hash()));
    let st = q.settings.as_ref().expect("settings snapshot");
    assert_eq!(st.re, 3e5);
    assert_eq!(m.settings, hash_json("foiley.settings2d.v1", st));

    let p = computed_project();
    let pp = &p.plane_polars[0];
    let m = pp.manifest.as_ref().expect("3D polar manifest");
    assert_eq!(m.geometry_hash, plane_hash(&p.planes[0]));
    assert_eq!(m.settings, hash_json("foiley.def3d.v1", pp.def.as_ref().unwrap()));
    assert!(m.mesh_hash.is_some(), "VLM run without a mesh hash");
    // one family per foil of the plane, with the foils' geometry
    let names = p.planes[0].foil_names();
    assert_eq!(m.families.len(), names.len(), "{:?}", m.families);
    for fr in &m.families {
        assert_eq!(Some(fr.foil_geom_hash), p.foil(&fr.foil).map(|f| f.geom_hash()));
    }
    // the lifting line has no mesh
    let db = p.foil_db();
    let llt = run_polar(&p.planes[0], &db, &PlanePolarDef { method: Method::Llt, range: Range { from: 2.0, to: 2.0, step: 1.0 }, ..Default::default() }, &Progress::none(), &|_| {}).unwrap();
    assert!(llt.manifest.as_ref().is_some_and(|m| m.mesh_hash.is_none() && !m.families.is_empty()));
}

#[test]
fn g03_hashes_stable_across_save_load() {
    let p = computed_project();
    let dir = tmp("g03");
    let path = dir.join("p.foiley");
    p.save(&path).unwrap();
    let q = Project::load(&path).unwrap();
    for (a, b) in p.foil_polars.iter().zip(&q.foil_polars) {
        assert_eq!(a.foil_hash, b.foil_hash);
        assert_eq!(a.manifest, b.manifest);
        assert_eq!(b.foil_hash, q.foil(&b.foil).map(|f| f.geom_hash()), "stored hash differs from the reloaded geometry");
    }
    let (a, b) = (&p.plane_polars[0], &q.plane_polars[0]);
    assert_eq!(a.manifest, b.manifest);
    let m = b.manifest.as_ref().unwrap();
    let pl: &Plane = &q.planes[0];
    assert_eq!(m.geometry_hash, plane_hash(pl), "plane hash changed across save/load");
    let mesh = foiley_plane::mesh::build(pl, &q.foil_db(), &foiley_plane::analysis::mesh_options(b.def.as_ref().unwrap(), &Default::default(), &Default::default())).unwrap();
    assert_eq!(m.mesh_hash, Some(mesh_hash(&mesh)), "mesh hash changed across save/load");
    // family hashes recomputed from the reloaded database agree
    let db = q.foil_db();
    for fr in &m.families {
        let id = db.resolve(&fr.foil, &Default::default()).unwrap().unwrap();
        assert_eq!(db.family(id).key.hash(), fr.key_hash);
    }
    // a v1 file: hashes computed on load are the same after a v2 round trip
    let v1 = Project::load(&projects_dir().join("starter_v1.foiley")).unwrap();
    v1.save(&path).unwrap();
    let v2 = Project::load(&path).unwrap();
    let h = |p: &Project| p.foil_polars.iter().map(|q| (q.foil_hash, q.legacy)).collect::<Vec<_>>();
    assert_eq!(h(&v1), h(&v2));
    assert!(h(&v2).iter().all(|(h, l)| h.is_some() && *l));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn f016_geometry_edit_invalidates_family() {
    let mut p = Project::new("edit");
    let mut f = Foil::naca("0012").unwrap();
    f.hinge = Some((0.75, 0.5));
    p.foils.push(f.clone());
    for re in [2e5, 6e5] {
        p.add_foil_polar(polar_of(&f, re));
    }
    // a deflected family on the same foil
    let mut d = polar_of(&f.with_flap(0.75, 0.5, 10.0).unwrap(), 6e5);
    d.name.push_str(" flap 10");
    d.foil = f.name.clone();
    d.flap_deg = 10.0;
    d.foil_hash = Some(f.geom_hash());
    p.add_foil_polar(d);
    let fam = |p: &Project| p.foil_db().families_of(&f.name).len();
    assert_eq!(fam(&p), 2);

    // an edit outside upsert_foil: the families go stale on hash mismatch
    let mut q = p.clone();
    q.foils[0].y[10] += 0.002;
    let db = q.foil_db();
    assert_eq!(db.families_of(&f.name).len(), 0, "edited geometry still has families");
    assert_eq!(db.stale_polars().len(), 3);

    // a note-only change keeps everything
    let mut q = p.clone();
    q.upsert_foil(Foil { note: "renamed note".into(), ..f.clone() });
    assert_eq!(q.foil_polars.len(), 3);
    // a moved hinge drops the deflected polar only
    let mut q = p.clone();
    q.upsert_foil(Foil { hinge: Some((0.7, 0.5)), ..f.clone() });
    assert_eq!(q.foil_polars.len(), 2);
    assert!(q.foil_polars.iter().all(|x| x.flap_deg == 0.0));
    // new coordinates drop every polar of the foil
    let mut q = p.clone();
    let g = Foil { name: f.name.clone(), hinge: f.hinge, ..Foil::naca("2412").unwrap() };
    q.upsert_foil(g);
    assert!(q.foil_polars.is_empty());
    assert_eq!(fam(&q), 0);
}
