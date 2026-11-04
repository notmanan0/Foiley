//! F010: AVL export fidelity.
mod common;
use common::*;

#[test]
#[ignore = "F010a: asymmetric wings are exported with YDUPLICATE"]
fn f010_asymmetric_export_two_surfaces() {
    let f = naca0012();
    let f9 = foiley_foil::Foil::naca("0009").unwrap();
    let mut db = db_with(&f);
    db.add_foil(&f9);
    let mut asym = review_rect(&f);
    asym.wings[0].symmetric = false;
    for sec in &mut asym.wings[0].sections {
        sec.foil_left = f9.name.clone();
    }
    let avl = foiley_plane::export::avl(&asym, &db, 0.0);
    assert!(avl.contains("NACA\n0012"));
    assert!(avl.contains("NACA\n0009"), "left-side foil lost");
    assert!(!avl.contains("YDUPLICATE"), "asymmetric wing exported as mirrored");
}
