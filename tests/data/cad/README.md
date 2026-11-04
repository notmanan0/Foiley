Example CAD bodies from flow5's `API_examples` (PlaneRun5 STL, PlaneRun6 STEP),
copyright André Deperrois, distributed under the GNU GPL v3 (see flow5's
LICENSE). Used by `crates/plane/tests/validation.rs`.

The other `.stp` files are hand-made STEP (AP214) B-rep fixtures for the
CAD importer (`crates/plane/tests/f022.rs`), written by `make_fixtures.py`
(part of Foiley, same licence as Foiley; run
`python tests/data/cad/make_fixtures.py` to regenerate them). Each has a
closed-form area and volume:

| file | part |
|---|---|
| `box_hole.stp` | 2 x 1 x 0.5 m box with a 0.2 m radius through hole (file in mm) |
| `cyl_seam.stp`, `cyl_pcurve.stp` | closed cylinder along x with a tilted seam, without and with PCURVEs |
| `frustum_rev.stp` | frustum as a surface of revolution of a line, planar caps |
| `cyl_partial.stp` | open cylindrical patch crossing the u = 0 seam |
| `cyl_annulus.stp` | cylindrical band bounded by two circles, no seam edge |
| `cone_apex.stp` | cone with its apex (singular point), planar base |
| `sphere.stp` | sphere bounded by one meridian seam |
| `torus.stp` | torus with u and v seams |
| `ellipse_prism.stp` | linear extrusion of an ellipse, planar caps |
| `bspline_patch.stp` | bicubic B-spline patch S = (u, v, u^2 v) split along a B-spline diagonal curve, no PCURVEs |
| `bspline_plate_hole.stp` | bilinear B-spline plate with a rational B-spline circle hole |
| `loft_fuselage.stp` | closed rational B-spline loft: cone, cylinder, cone, poles at nose and tail |
| `cyl_fuselage.stp` | cylinder along x with hemispherical nose and tail, each hemisphere closing through its pole |
| `nacelle.stp` | annular nacelle: an ellipse revolved at a radius, one doubly periodic face |
| `multi_shell.stp` | two solids in one file: the box with a hole and a sphere beside it |
