//! F004b: the fuselage tab edits every body (UI state only, no rendering).
#[path = "../src/body_list.rs"]
mod body_list;

use body_list::{apply, selected, selected_mut, BodyAction};
use foiley_app::CadImportResult;
use foiley_core::Vec3;
use foiley_plane::geom::Plane;

/// A closed octahedron, in millimetres.
fn octahedron_mm() -> Vec<[Vec3; 3]> {
    let v = [Vec3::new(1000.0, 0.0, 0.0), Vec3::new(-1000.0, 0.0, 0.0), Vec3::new(0.0, 100.0, 0.0), Vec3::new(0.0, -100.0, 0.0), Vec3::new(0.0, 0.0, 100.0), Vec3::new(0.0, 0.0, -100.0)];
    let mut t = Vec::new();
    for &x in &[0, 1] {
        for &y in &[2, 3] {
            for &z in &[4, 5] {
                t.push([v[x], v[y], v[z]]);
            }
        }
    }
    t
}

#[test]
fn f004b_gui_state_edits_second_body() {
    let mut p = Plane::default_plane("NACA 0012", "NACA 0012");
    let mut sel = 0usize;
    assert_eq!(selected(&p, &mut sel), None);
    apply(&mut p, &mut sel, BodyAction::AddPod);
    apply(&mut p, &mut sel, BodyAction::AddPod);
    assert_eq!(p.fuses.len(), 2);
    assert_eq!(sel, 1, "a new body is selected");
    assert_eq!((p.fuses[0].name.as_str(), p.fuses[1].name.as_str()), ("Fuselage", "Fuselage 2"));
    assert_eq!(body_list::label(&p, 1), "2. Fuselage 2");
    assert!(p.fuses[1].position.y - p.fuses[0].position.y >= body_list::NEW_BODY_SPACING - 1e-12, "a new pod is placed beside the last one");

    // edit the second body through the selection; the first one is untouched
    let first = p.fuses[0].clone();
    let f = selected_mut(&mut p, &mut sel).unwrap();
    f.position.x = 0.4;
    f.frames[2].half_width = 0.08;
    f.mass = 1.5;
    assert_eq!(p.fuses[0], first);
    assert_eq!((p.fuses[1].position.x, p.fuses[1].frames[2].half_width, p.fuses[1].mass), (0.4, 0.08, 1.5));

    // per-body CAD import (the CAD job's result, scaled to metres by the
    // reader) goes into the selected body only, keeping its mass
    let target = selected(&p, &mut sel).unwrap();
    let metres: Vec<[Vec3; 3]> = octahedron_mm().into_iter().map(|t| t.map(|v| v * 1e-3)).collect();
    let import = |body: usize, triangles: Vec<[Vec3; 3]>| CadImportResult { body, name: "pod_cad".into(), triangles, notes: vec![] };
    import(target, metres).apply(&mut p).unwrap();
    assert_eq!(p.fuses[0], first);
    assert_eq!(p.fuses[1].name, "pod_cad");
    assert_eq!(p.fuses[1].triangles.len(), 8);
    assert_eq!(p.fuses[1].mass, 1.5);
    assert!((p.fuses[1].length() - 2.0).abs() < 1e-12, "mm scaled to m: {}", p.fuses[1].length());
    assert!(import(5, octahedron_mm()).apply(&mut p).unwrap_err().contains("body 6 does not exist"));

    // duplicate, select and remove
    apply(&mut p, &mut sel, BodyAction::Select(0));
    apply(&mut p, &mut sel, BodyAction::Duplicate);
    assert_eq!(p.fuses.len(), 3);
    assert_eq!((sel, p.fuses[2].name.as_str()), (2, "Fuselage copy"));
    apply(&mut p, &mut sel, BodyAction::Select(1));
    apply(&mut p, &mut sel, BodyAction::Remove);
    assert_eq!(p.fuses.iter().map(|f| f.name.as_str()).collect::<Vec<_>>(), ["Fuselage", "Fuselage copy"]);
    assert_eq!(sel, 1);
    apply(&mut p, &mut sel, BodyAction::Select(9));
    assert_eq!(sel, 1, "selection clamps to the list");
    apply(&mut p, &mut sel, BodyAction::Remove);
    apply(&mut p, &mut sel, BodyAction::Remove);
    assert!(p.fuses.is_empty());
    assert_eq!(selected(&p, &mut sel), None);
}
