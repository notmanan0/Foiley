//! A cosmetic rename (of the plane, its description, a wing, a body, a point
//! mass or a foil) keeps a 3D polar matched to its plane, so its results are
//! still shown; a geometry edit does not.
use foiley_app::Project;
use foiley_core::Progress;
use foiley_foil::Foil;
use foiley_plane::analysis::{mesh_options, Range};
use foiley_plane::geom::Fuse;
use foiley_plane::mesh::build;
use foiley_plane::{run_polar, Method, Plane, PlanePolarDef};
use std::path::PathBuf;

/// The starter project, with a body, and one inviscid VLM polar of its plane.
fn project() -> Project {
    let mut p = Project::starter();
    p.planes[0].fuses.push(Fuse::default_pod());
    for wing in &mut p.planes[0].wings {
        for section in &mut wing.sections {
            section.nx = 3;
            section.ny = 3;
        }
    }
    let def = PlanePolarDef { method: Method::Vlm2, viscous: false, range: Range { from: 2.0, to: 2.0, step: 1.0 }, ..Default::default() };
    let polar = run_polar(&p.planes[0], &p.foil_db(), &def, &Progress::none(), &|_| {}).unwrap();
    p.add_plane_polar(polar);
    assert!(matches(&p, 0), "a fresh polar must match its plane");
    p
}

/// Whether 3D polar `i` matches its plane: plane hash and mesh hash.
fn matches(p: &Project, i: usize) -> bool {
    let polar = &p.plane_polars[i];
    let Some(plane) = p.planes.iter().find(|pl| pl.name == polar.plane) else { return false };
    if !polar.matches_plane_geometry(plane) {
        return false;
    }
    let options = mesh_options(polar.def.as_ref().unwrap(), &Default::default(), &Default::default());
    let mesh = build(plane, &p.foil_db(), &options).unwrap();
    polar.manifest.as_ref().unwrap().mesh_hash == Some(polar.mesh_hash(&mesh))
}

/// Edit plane 0 the way the plane editor does.
fn edit(p: &mut Project, change: impl FnOnce(&mut Plane)) {
    let before = p.planes[0].clone();
    let mut after = before.clone();
    change(&mut after);
    for polar in p.plane_polars.iter_mut().filter(|polar| polar.plane == before.name) {
        polar.follow_plane_edit(&before, &after);
    }
    p.planes[0] = after;
}

fn reloaded(p: &Project) -> Project {
    serde_json::from_str(&serde_json::to_string(p).unwrap()).unwrap()
}

fn tmp(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("foiley_renames_{}_{name}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn a_renamed_foil_keeps_the_results() {
    let mut p = project();
    p.rename_foil("NACA 2412", "Root");
    assert!(p.planes[0].foil_names().contains(&"Root".to_string()));
    assert!(p.foil("Root").is_some());
    assert!(matches(&p, 0), "a foil rename hid the results");
    assert!(matches(&reloaded(&p), 0), "the rename record did not survive saving");
    // renamed back: nothing left to record
    p.rename_foil("Root", "NACA 2412");
    assert!(p.plane_polars[0].renamed.is_empty());
    assert!(matches(&p, 0));
}

#[test]
fn foils_swapped_by_renames_match_but_swapped_sections_do_not() {
    let mut p = project();
    // three renames trade the two foils' names
    p.rename_foil("NACA 2412", "spare");
    p.rename_foil("NACA 0009", "NACA 2412");
    p.rename_foil("spare", "NACA 0009");
    assert!(matches(&p, 0));

    // a section switched to the other foil is a geometry change
    let mut q = project();
    edit(&mut q, |plane| plane.wings[0].sections[0].foil_left = "NACA 0009".into());
    assert!(!matches(&q, 0));
}

#[test]
fn a_new_foil_under_a_renamed_foils_old_name_is_not_the_old_foil() {
    let mut p = project();
    p.rename_foil("NACA 2412", "Root");
    let mut other = Foil::naca("4415").unwrap();
    other.name = "NACA 2412".into();
    p.foils.push(other);
    edit(&mut p, |plane| plane.wings[0].sections[0].foil_left = "NACA 2412".into());
    assert!(!matches(&p, 0), "a different foil passed for the renamed one");
}

#[test]
fn renamed_labels_keep_the_results_and_geometry_edits_do_not() {
    let mut p = project();
    let original = p.planes[0].clone();
    edit(&mut p, |plane| {
        plane.name = "Glider".into();
        plane.description = "for the club".into();
        plane.wings[0].name = "Main wing".into();
        plane.fuses[0].name = "Pod".into();
        plane.point_masses[0].tag = "battery".into();
    });
    assert_eq!(p.plane_polars[0].plane, "Glider");
    assert!(matches(&p, 0), "a rename hid the results");
    assert!(matches(&reloaded(&p), 0));

    edit(&mut p, |plane| plane.wings[0].sections[0].twist += 2.0);
    assert!(!matches(&p, 0), "a geometry edit kept the results");
    edit(&mut p, |plane| plane.wings[0].sections[0].twist -= 2.0);
    assert!(matches(&p, 0), "undoing the edit by hand must match again");

    // reordering wings keeps every name but not the geometry
    edit(&mut p, |plane| plane.wings.swap(0, 1));
    assert!(!matches(&p, 0));
    edit(&mut p, |plane| plane.wings.swap(0, 1));
    assert!(matches(&p, 0));

    // every label renamed back: the record is empty again
    edit(&mut p, |plane| *plane = original);
    assert!(p.plane_polars[0].renamed.is_empty(), "{:?}", p.plane_polars[0].renamed);
    assert!(matches(&p, 0));
}

#[test]
fn imported_polars_match_their_renamed_plane_and_foils() {
    let dir = tmp("import");
    let path = dir.join("incoming.foiley");
    project().save(&path).unwrap();
    // the plane and both foil names are taken: the import renames all three
    let mut p = Project::starter();
    p.import(&path).unwrap();
    assert_eq!(p.planes[1].name, "Plane (2)");
    assert!(p.planes[1].foil_names().iter().all(|name| name.ends_with("(2)")), "{:?}", p.planes[1].foil_names());
    let i = p.plane_polars.iter().position(|polar| polar.plane == "Plane (2)").unwrap();
    assert!(matches(&p, i), "the imported polar lost its results");
    std::fs::remove_dir_all(dir).ok();
}
