//! Regressions for the lifting-surface failure that kernel-only checks miss.
use foiley_core::{Progress, Spacing, Vec3};
use foiley_foil::Foil;
use foiley_plane::analysis::*;
use foiley_plane::geom::*;
use foiley_plane::{mesh, FoilDb};

fn wing(foil: &str, nx: usize, ny: usize) -> (Plane, FoilDb) {
    let mut db = FoilDb::new();
    let f = Foil::naca(foil).unwrap();
    db.add_foil(&f);
    let mut a = WingSection::new(0.0, 1.0, 0.0, &f.name);
    a.nx = nx;
    a.ny = ny;
    a.ydist = Spacing::Uniform;
    let b = WingSection::new(4.0, 1.0, 0.0, &f.name);
    (
        Plane {
            name: "linear test".into(),
            description: String::new(),
            wings: vec![Wing::new("wing", WingKind::Main, vec![a, b])],
            fuses: vec![],
            point_masses: vec![],
        },
        db,
    )
}

fn definition(method: Method) -> PlanePolarDef {
    PlanePolarDef { method, viscous: false, keep_cp: true, reference: RefDims::Planform, range: Range { from: -5.0, to: 5.0, step: 5.0 }, ..Default::default() }
}

/// Symmetric polar at -5, 0 and 5 degrees with per-mesh symmetry checks.
/// Returns CL, |pressure CL / wake CL - 1| and the RMS TE Cp jump at 5 degrees.
fn refine_step(nx: usize, ny: usize) -> (f64, f64, f64) {
    let (p, db) = wing("0012", nx, ny);
    let d = definition(Method::TriLinear);
    let mesh = mesh::build(&p, &db, &mesh_options(&d, &Default::default(), &Default::default())).unwrap();
    let r = run_polar(&p, &db, &d, &Progress::none(), &|_| {}).unwrap();
    assert!(r.failed.is_empty(), "{:?}", r.failed);
    assert_eq!(r.points.len(), 3);
    assert_eq!(r.unknowns, 3 * mesh.panels.len());
    let (neg, zero, pos) = (&r.points[0], &r.points[1], &r.points[2]);
    println!(
        "{nx}x{ny} ({} unknowns): CL={} CDi={} Cm={} zero={} yaw={} cp range={:?}",
        r.unknowns,
        pos.cl,
        pos.cdi,
        pos.cm,
        zero.cl,
        pos.cn.unwrap(),
        (pos.cp.iter().copied().fold(f32::INFINITY, f32::min), pos.cp.iter().copied().fold(f32::NEG_INFINITY, f32::max))
    );
    assert!(zero.cl.abs() < 1e-5);
    assert!((pos.cl + neg.cl).abs() < 2e-5);
    assert!(pos.cn.unwrap().abs() < 1e-6 && pos.cy.unwrap().abs() < 1e-6 && pos.cr.unwrap().abs() < 1e-6);
    assert!(pos.cdi > 0.0 && (pos.cdi - neg.cdi).abs() < 1e-6);
    assert!(pos.cp.iter().all(|c| c.is_finite()));
    assert!((0.39..0.46).contains(&pos.cl), "CL={}", pos.cl);
    let cl_pressure: f64 = mesh.panels.iter().zip(&pos.cp).map(|(p, c)| -(*c as f64) * p.area * p.normal.dot(wind_normal(5.0)) / 8.0).sum();
    let te_error = (mesh.strips.iter().map(|s| (pos.cp[s.te_up] - pos.cp[s.te_low]).powi(2) as f64).sum::<f64>() / mesh.strips.len() as f64).sqrt();
    println!("pressure CL={cl_pressure} TE Cp jump rms={te_error}");
    (pos.cl, (cl_pressure / pos.cl - 1.0).abs(), te_error)
}

#[test]
fn lifting_wing_refines_and_preserves_symmetry() {
    let [a, b, c] = [(6, 8), (10, 12), (14, 16)].map(|(nx, ny)| refine_step(nx, ny));
    // Lift converges (differences shrink about fivefold); pressure lift and
    // the TE traces approach the wake values.
    assert!((c.0 - b.0).abs() < 0.3 * (b.0 - a.0).abs() && (c.0 / b.0 - 1.0).abs() < 0.002, "{a:?} {b:?} {c:?}");
    assert!(c.1 < 0.05 && c.1 < 0.55 * a.1, "{a:?} {c:?}");
    assert!(c.2 < 0.02 && c.2 < 0.45 * a.2, "{a:?} {c:?}");
}

/// The 24x24 mesh (14,100 unknowns, a 1.6 GB dense matrix) and the 40x40
/// quad reference are too heavy for the default suite. Run with
/// `cargo test --release -p foiley-plane --test trilinear -- --ignored`.
#[test]
#[ignore]
fn lifting_wing_fine_mesh_convergence() {
    let [coarse, medium, fine] = [(6, 8), (14, 16), (24, 24)].map(|(nx, ny)| refine_step(nx, ny));
    assert!((fine.0 / medium.0 - 1.0).abs() < 0.005, "{medium:?} {fine:?}");
    assert!(fine.1 < 0.035 && fine.1 < 0.5 * coarse.1, "{coarse:?} {fine:?}");
    assert!(fine.2 < 0.01 && fine.2 < 0.25 * coarse.2, "{coarse:?} {fine:?}");
    let (p, db) = wing("0012", 40, 40);
    let q = run_polar(&p, &db, &definition(Method::Panel), &Progress::none(), &|_| {}).unwrap();
    println!("fine quad CL={}", q.points[2].cl);
    // Different trial/test spaces retain a ~3.6% offset in this benchmark.
    // Check independent thin-limit theory, pressure convergence, and TE traces
    // as well, rather than treating the uniform method as an exact reference.
    assert!((fine.0 / q.points[2].cl - 1.0).abs() < 0.04, "linear {} quad {}", fine.0, q.points[2].cl);
}

#[test]
fn camber_ground_and_rotational_unit_solutions() {
    let (mut p, db) = wing("2412", 8, 10);
    p.wings[0].sections[1].offset = 0.6;
    let d = definition(Method::TriLinear);
    let m = mesh::build(&p, &db, &mesh_options(&d, &Default::default(), &Default::default())).unwrap();
    let s = make_system(&d, &m, Vec3::new(0.25, 0.0, 0.0), 1.0, &Progress::none()).unwrap();
    let zero = s.loads(Vec3::X, Vec3::ZERO);
    let pos = s.loads(wind_dir(5.0, 0.0), Vec3::ZERO);
    assert!(zero.force.z > 0.0 && pos.force.z > zero.force.z);
    let rot = s.loads(Vec3::X, Vec3::new(0.01, -0.02, 0.01));
    assert!(rot.force.is_finite() && rot.moment.is_finite());
    let roll = s.loads(Vec3::X, Vec3::X * 0.01);
    assert!(roll.moment.x < zero.moment.x);
    let mut gd = d.clone();
    gd.ground_height = Some(0.4);
    let ground = make_system(&gd, &m, Vec3::new(0.25, 0.0, 0.0), 1.0, &Progress::none()).unwrap().loads(wind_dir(5.0, 0.0), Vec3::ZERO);
    println!("camber free lift {} drag {} ground lift {} drag {}", pos.force.z, pos.induced_force.x, ground.force.z, ground.induced_force.x);
    assert!(ground.force.z > pos.force.z);
    assert!(ground.induced_force.dot(wind_dir(5.0, 0.0)) < pos.induced_force.dot(wind_dir(5.0, 0.0)));
}

#[test]
fn thin_limit_matches_lifting_surface_theory() {
    let (p, mut db) = wing("0012", 14, 20);
    let f = Foil::naca("0012").unwrap().scaled_tc(0.01, 1.0).unwrap();
    db.add_foil(&f);
    let d = definition(Method::TriLinear);
    let r = run_polar(&p, &db, &d, &Progress::none(), &|_| {}).unwrap();
    let slope = r.points[2].cl / 5f64.to_radians().sin();
    println!("thin wing lift slope {slope}");
    assert!((4.5..4.75).contains(&slope));
}

#[test]
fn fuselage_carry_through_and_flapped_wing() {
    let (mut p, mut db) = wing("2412", 6, 6);
    let mut f = Foil::naca("2412").unwrap();
    f.hinge = Some((0.7, 0.5));
    db.add_foil(&f);
    p.fuses.push(Fuse::default_pod());
    p.wings[0].sections[1].chord = 0.6;
    p.wings[0].sections[1].offset = 0.5;
    p.wings[0].sections[0].dihedral = 4.0;
    let mut d = definition(Method::TriLinear);
    d.include_fuse = true;
    d.range = Range { from: 2.0, to: 2.0, step: 1.0 };
    let base = run_polar(&p, &db, &d, &Progress::none(), &|_| {}).unwrap();
    d.flaps = vec![((0, 0), 5.0), ((0, 1), 5.0)];
    let flap = run_polar(&p, &db, &d, &Progress::none(), &|_| {}).unwrap();
    assert!(base.failed.is_empty() && flap.failed.is_empty());
    let (a, b) = (&base.points[0], &flap.points[0]);
    println!("fuselage CL={} CDi={} flap CL={} CDi={} Cm={} -> {}", a.cl, a.cdi, b.cl, b.cdi, a.cm, b.cm);
    assert!(b.cl > a.cl && b.cm < a.cm);
    for r in [a, b] {
        assert!(r.cdi > 0.0 && r.cn.unwrap().abs() < 1e-5 && r.cy.unwrap().abs() < 1e-5);
        assert!(r.cp.iter().all(|c| c.is_finite()));
    }
}
