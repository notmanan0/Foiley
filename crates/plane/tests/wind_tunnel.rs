//! NACA Report 530, Table VI (p. 444), full-scale rectangular wing.
//! https://ntrs.nasa.gov/citations/19930091603
//! Table VI gives tunnel-corrected finite-span alpha and CL and reports Cm
//! about the measured aerodynamic center. Figure 4 gives test 42-2's speed.

use foiley_core::{Progress, Vec3};
use foiley_foil::Foil;
use foiley_plane::geom::Inertia;
use foiley_plane::{run_polar, FoilDb, Method, Plane, PlanePolarDef, PolarKind, Wing, WingKind, WingSection};

// (finite-wing alpha in degrees, CL, Cm about the reported aerodynamic center)
const MEASUREMENTS: [(f64, f64, f64); 5] = [
    (-2.6, -0.1, -0.005),
    (0.3, 0.1, -0.006),
    (3.1, 0.3, -0.006),
    (5.8, 0.5, -0.008),
    (8.4, 0.7, -0.008),
];

#[test]
fn naca_23012_full_scale_wing() {
    let foot = 0.3048;
    let chord = 6.0 * foot;
    let half_span = 18.0 * foot;
    let speed = 85.5 * foot;
    let foil = Foil::naca("23012").unwrap();
    let mut db = FoilDb::new();
    db.add_foil(&foil);
    let mut root = WingSection::new(0.0, chord, 0.0, &foil.name);
    root.nx = 12;
    root.ny = 18;
    let tip = WingSection::new(half_span, chord, 0.0, &foil.name);
    let wing = Wing::new("NACA tunnel model", WingKind::Main, vec![root, tip]);
    let plane = Plane { name: foil.name, description: String::new(), wings: vec![wing], fuses: vec![], point_masses: vec![] };

    // The report gives x forward and y upward from quarter chord (p. 438).
    // Foiley's origin is the root leading edge, with x aft and z upward.
    let reference = Vec3::new((0.25 - 0.0191) * chord, 0.0, 0.0887 * chord);
    let mut def = PlanePolarDef {
        kind: PolarKind::List,
        viscous: false,
        list: MEASUREMENTS.iter().map(|&(alpha, _, _)| (alpha, 0.0, speed)).collect(),
        // A positive placeholder mass selects this reference in a prescribed-speed run.
        inertia: Some(Inertia { mass: 1.0, cog: reference, ..Default::default() }),
        ..Default::default()
    };
    // Per-method absolute regression limits for this mesh and attached-flow range.
    // The experimental Cm precision is +/-0.003; the solvers have larger errors.
    for (method, cl_limit, cm_limit) in [
        (Method::Vlm2, 0.020, 0.013),
        (Method::Panel, 0.035, 0.009),
        (Method::TriLinear, 0.065, 0.020),
    ] {
        def.method = method;
        let polar = run_polar(&plane, &db, &def, &Progress::none(), &|_| {}).unwrap();
        assert_eq!(polar.points.len(), MEASUREMENTS.len(), "{method:?}: {:?}", polar.failed);
        for (point, (alpha, cl, cm)) in polar.points.iter().zip(MEASUREMENTS) {
            assert_eq!(point.alpha, alpha);
            println!("{method:?}, {alpha:5.1} deg: CL {:.4} (tunnel {cl:.3}), Cm {:.4} (tunnel {cm:.3})", point.cl, point.cm);
            assert!((point.cl - cl).abs() < cl_limit, "{method:?}, alpha {alpha}: CL {} vs {cl}", point.cl);
            assert!(point.cm < 0.0 && (point.cm - cm).abs() < cm_limit, "{method:?}, alpha {alpha}: Cm {} vs {cm}", point.cm);
        }
    }
}
