//! F005: stability derivatives from the full (viscous + extra drag) loads.
mod common;
use common::*;
use foiley_plane::analysis::*;
use foiley_plane::geom::Plane;

#[test]
#[ignore = "F005a: added parasite drag leaves Xu unchanged"]
fn f005_fullload_xu_moves_with_drag() {
    let f = naca0012();
    let f9 = foiley_foil::Foil::naca("0009").unwrap();
    let mut db = db_with(&f);
    db.add_foil(&f9);
    let glider = Plane::default_plane(&f.name, &f9.name);
    let sd = PlanePolarDef { viscous: false, derivatives: true, range: Range { from: 4.0, to: 4.0, step: 1.0 }, ..Default::default() };
    let clean = run1(&glider, &db, &sd);
    let draggy = run1(&glider, &db, &PlanePolarDef { extra_drag: vec![ExtraDrag { name: "audit drag".into(), area: 1.0, cd: 0.1 }], ..sd });
    assert!(draggy.cd > clean.cd);
    let (a, b) = (clean.stability.unwrap().dim.xu, draggy.stability.unwrap().dim.xu);
    assert!((a - b).abs() > 1e-6, "Xu {a} vs {b}");
}
