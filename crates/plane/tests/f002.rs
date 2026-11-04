//! F002: explicit polar families.
mod common;
use common::*;

#[test]
#[ignore = "F002: polar lookup depends on insertion order"]
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
    let x = ab.at_re_cl(&f.name, 1e6, 0.5).unwrap().0.cd;
    let y = ba.at_re_cl(&f.name, 1e6, 0.5).unwrap().0.cd;
    assert_eq!(x, y, "CD depends on polar insertion order");
}
