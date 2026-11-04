//! F004: multi-body meshing.
mod common;
use common::*;
use foiley_plane::geom::{Fuse, Plane};
use foiley_plane::mesh::{build, MeshOptions, PanelKind};

#[test]
#[ignore = "F004a: only the first fuselage is meshed"]
fn f004_two_pods_double_panels_mirror_loads() {
    let f = naca0012();
    let db = db_with(&f);
    let mut pods = review_rect(&f);
    pods.fuses.push(Fuse::default_pod());
    let opts = MeshOptions { thick: true, include_fuse: true, ..Default::default() };
    let count = |p: &Plane| build(p, &db, &opts).unwrap().panels.iter().filter(|p| p.kind == PanelKind::Fuse).count();
    let one = count(&pods);
    let mut second = Fuse::default_pod();
    second.position.y = 3.0;
    pods.fuses.push(second);
    assert_eq!(count(&pods), 2 * one);
}
