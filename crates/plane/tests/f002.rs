//! F002: explicit polar families.
mod common;
use common::*;
use foiley_plane::analysis::*;
use foiley_plane::{FamilySelector, PlaneError};

fn ncrit(n: f64) -> FamilySelector {
    FamilySelector { ncrit: Some(n), ..Default::default() }
}

/// Two families (Ncrit 9 and 5) at the same Re: the selected family decides,
/// never the order in which the polars were added.
#[test]
fn f002_polar_order_irrelevant() {
    let f = naca0012();
    let db = db_with(&f);
    let pa = synth_polar(&f.name, "N9", 1e6, 9.0, 0.01);
    let pb = synth_polar(&f.name, "N5", 1e6, 5.0, 0.03);
    let mut ab = db.clone();
    ab.add_polar(&pa);
    ab.add_polar(&pb);
    let mut ba = db.clone();
    ba.add_polar(&pb);
    ba.add_polar(&pa);
    assert_eq!(ab.families_of(&f.name).len(), 2, "Ncrit 9 and 5 are two families");
    for (n, want) in [(9.0, 0.01), (5.0, 0.03)] {
        let cd = |db: &foiley_plane::FoilDb| db.at_re_cl(db.resolve(&f.name, &ncrit(n)).unwrap().unwrap(), 1e6, 0.5).unwrap().0.cd;
        assert_eq!(cd(&ab), cd(&ba), "CD depends on polar insertion order");
        assert_eq!(cd(&ab), want, "Ncrit {n} selects the wrong family");
    }
    // and through a 3D analysis
    let p = review_rect(&f);
    let def = PlanePolarDef { viscous: true, free_wake: 0, family: ncrit(9.0), ..review_def() };
    let (x, y) = (run1(&p, &ab, &def), run1(&p, &ba, &def));
    assert_eq!(x.cdv, y.cdv);
    assert!((x.cdv - 0.01).abs() < 2e-3, "profile drag {} is not from the Ncrit 9 family", x.cdv);
}

/// Several families match the selection: the run fails and names them.
#[test]
fn g07_ambiguous_family_errors() {
    let f = naca0012();
    let mut db = db_with(&f);
    db.add_polar(&synth_polar(&f.name, "N9", 1e6, 9.0, 0.01));
    db.add_polar(&synth_polar(&f.name, "N5", 1e6, 5.0, 0.03));
    let p = review_rect(&f);
    let def = PlanePolarDef { viscous: true, free_wake: 0, ..review_def() };
    let run = |d: &PlanePolarDef| run_polar(&p, &db, d, &foiley_core::Progress::none(), &|_| {});
    match run(&def) {
        Err(PlaneError::Input(m)) => {
            assert!(m.contains("2 2D polar families") && m.contains("Ncrit=9") && m.contains("Ncrit=5"), "message does not list the candidates: {m}");
        }
        other => panic!("ambiguous families must fail, got {:?}", other.map(|p| p.points.len())),
    }
    // no match also fails, listing what exists
    match run(&PlanePolarDef { family: ncrit(7.0), ..def.clone() }) {
        Err(PlaneError::Input(m)) => assert!(m.contains("no 2D polar family") && m.contains("Ncrit=9"), "{m}"),
        other => panic!("no matching family must fail, got {:?}", other.map(|p| p.points.len())),
    }
    // the LLT resolves through the same selection
    assert!(matches!(run(&PlanePolarDef { method: Method::Llt, ..def.clone() }), Err(PlaneError::Input(_))));
    // a unique selection runs; an inviscid run does not need polars at all
    assert!(run(&PlanePolarDef { family: ncrit(5.0), ..def.clone() }).is_ok());
    assert!(run(&PlanePolarDef { viscous: false, ..def }).is_ok());
}

/// Families are keyed by geometry too: replacing a foil's coordinates
/// re-keys its polars.
#[test]
fn f002_family_key_has_geometry_hash() {
    let f = naca0012();
    let mut db = db_with(&f);
    db.add_polar(&synth_polar(&f.name, "N9", 1e6, 9.0, 0.01));
    let id = db.resolve(&f.name, &FamilySelector::default()).unwrap().unwrap();
    let h0 = db.family(id).key.foil_geom_hash;
    assert_eq!(h0, foiley_plane::foildb::foil_geom_hash(&f));
    let mut g = foiley_foil::Foil::naca("2412").unwrap();
    g.name = f.name.clone();
    db.add_foil(&g);
    let id = db.resolve(&f.name, &FamilySelector::default()).unwrap().unwrap();
    assert_ne!(db.family(id).key.foil_geom_hash, h0);
}
