//! Generated foils in projects: recipes survive saving and loading, files
//! without recipes still load, and planes naming any NACA family get their
//! foils generated on import.
use foiley_app::Project;
use foiley_foil::gen::cst::Cst;
use foiley_foil::{Foil, Recipe};
use foiley_plane::Plane;

fn temp(name: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("foiley_families_{}_{name}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn recipes_survive_a_project_round_trip() {
    let mut p = Project::new("families");
    p.foils.push(Foil::naca("64A010").unwrap());
    p.foils.push(Foil::naca("2412").unwrap());
    p.foils.push(Foil::from_recipe(&Recipe::Cst(Cst::default()), Some("wing root")).unwrap());
    let path = temp("roundtrip").join("p.foiley");
    p.save(&path).unwrap();
    let q = Project::load(&path).unwrap();
    for (a, b) in p.foils.iter().zip(&q.foils) {
        assert_eq!(a.recipe, b.recipe, "{}", a.name);
        assert!(a.recipe.is_some(), "{}", a.name);
        assert_eq!(a.geom_hash(), b.geom_hash(), "{}", a.name);
    }
    // the recipe regenerates the stored coordinates exactly
    for f in &q.foils {
        let g = Foil::from_recipe(f.recipe.as_ref().unwrap(), Some(&f.name)).unwrap();
        assert_eq!((&g.x, &g.y), (&f.x, &f.y), "{}", f.name);
    }
}

#[test]
fn foils_without_a_recipe_still_load_and_edits_clear_it() {
    let f = Foil::naca("0012").unwrap();
    let mut v = serde_json::to_value(&f).unwrap();
    v.as_object_mut().unwrap().remove("recipe");
    let g: Foil = serde_json::from_value(v).unwrap();
    assert_eq!(g.recipe, None);
    assert_eq!((g.x, g.y), (f.x.clone(), f.y.clone()));
    // a shape change clears the recipe; the unchanged foil keeps it
    assert!(f.with_te_gap(0.002, 0.8).unwrap().recipe.is_none());
    assert!(f.with_flap(0.75, 0.5, 0.0).unwrap().recipe.is_some());
}

#[test]
fn imported_planes_generate_every_naca_family_by_name() {
    let d = temp("import");
    let plane = Plane::default_plane("NACA 63-206", "NACA 0012-64");
    let xml = d.join("glider.xml");
    std::fs::write(&xml, foiley_plane::xml::write_plane(&plane)).unwrap();
    let mut p = Project::new("import");
    let msg = p.import(&xml).unwrap();
    for name in ["NACA 63-206", "NACA 0012-64"] {
        let f = p.foil(name).unwrap_or_else(|| panic!("{name} not generated: {msg}"));
        assert!(matches!(&f.recipe, Some(Recipe::Naca { .. })), "{name}");
    }
}
