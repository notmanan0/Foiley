//! F010b: .xfl / .fl5 import of hand-written minimal files (written by
//! tests/data/interchange/make_minimal.py), and the conversion report.
use foiley_app::legacy_io::{import_into, read_any, ConversionReport};
use foiley_app::Project;
use foiley_plane::analysis::{Method, PolarKind, RefDims};
use foiley_plane::WingKind;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/interchange").join(name)
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-12
}

#[test]
fn f010b_minimal_xfl_and_fl5_load() {
    // --- xflr5 v6 ---
    let x = read_any(&data("minimal_v6.xfl")).unwrap();
    assert_eq!(x.foils.len(), 1);
    assert_eq!((x.foils[0].name.as_str(), x.foils[0].hinge), ("M1", Some((0.75, 0.5))));
    assert_eq!((x.polars.len(), x.polars[0].points.len()), (1, 3));
    assert_eq!(x.planes.len(), 1);
    let p = &x.planes[0];
    assert_eq!(p.name, "Mini");
    assert_eq!(p.wings.iter().map(|w| (w.name.as_str(), w.kind)).collect::<Vec<_>>(), [("Main", WingKind::Main), ("Stab", WingKind::Elevator)]);
    let m = &p.wings[0];
    assert!(close(m.ry, 1.0) && close(m.sections[1].y, 1.0) && close(m.sections[1].chord, 0.15) && close(m.sections[0].dihedral, 3.0) && close(m.sections[1].twist, -2.0));
    assert!(close(p.wings[1].position.x, 0.8) && close(p.wings[1].ry, -1.5));
    assert_eq!(p.point_masses.len(), 1);
    // the plane polar definition (the one for a plane not in the file is dropped)
    assert_eq!(x.plane_polars.len(), 1);
    let pp = &x.plane_polars[0];
    assert_eq!((pp.plane.as_str(), pp.name.as_str()), ("Mini", "T2 trimmed"));
    assert!(pp.points.is_empty());
    let d = pp.def.as_ref().unwrap();
    assert_eq!((d.kind, d.method), (PolarKind::FixedLift, Method::Vlm2));
    assert!(close(d.velocity, 12.0) && close(d.alpha, 3.0) && d.viscous && d.include_fuse);
    assert_eq!(d.reference, RefDims::Custom { area: 0.42, span: 2.1, chord: 0.2 });
    let i = d.inertia.unwrap();
    assert!(close(i.mass, 1.1) && close(i.cog.x, 0.08) && close(i.cog.z, 0.01) && close(i.ixx, 0.05) && close(i.ixz, 0.001));
    assert_eq!(d.ground_height, Some(0.5));
    assert!(close(d.fluid.rho, 1.2) && close(d.fluid.mu, 1.2 * 1.6e-5), "kinematic viscosity -> dynamic");
    assert_eq!((d.wake_panels, d.wake_length, d.wake_ratio), (5, 50.0, 1.1));
    assert_eq!(d.extra_drag.len(), 1);
    assert!(close(d.extra_drag[0].area, 0.01) && close(d.extra_drag[0].cd, 0.5));

    // --- flow5 v7.54 ---
    let f = read_any(&data("minimal_v7.fl5")).unwrap();
    assert_eq!(f.foils.iter().map(|f| f.name.as_str()).collect::<Vec<_>>(), ["N1"]);
    assert_eq!(f.planes.len(), 1);
    let p = &f.planes[0];
    assert_eq!(p.name, "Flow mini");
    assert_eq!(p.wings.iter().map(|w| (w.kind, w.two_sided)).collect::<Vec<_>>(), [(WingKind::Main, true), (WingKind::Fin, false)]);
    assert!(close(p.wings[0].ry, 1.5) && close(p.wings[0].sections[1].y, 1.2));
    assert_eq!(f.plane_polars.len(), 1);
    let d = f.plane_polars[0].def.as_ref().unwrap();
    assert_eq!((f.plane_polars[0].plane.as_str(), d.name.as_str()), ("Flow mini", "T6 flap sweep"));
    assert_eq!((d.kind, d.method), (PolarKind::Control, Method::TriPanel));
    assert!(close(d.velocity, 15.0) && close(d.alpha, 2.0) && !d.include_fuse && !d.viscous_from_cl);
    assert_eq!(d.flaps, vec![((0, 1), 5.0)]);
    assert_eq!(d.reference, RefDims::Planform);
    assert_eq!(d.ground_height, Some(0.4));
    assert!(close(d.inertia.unwrap().mass, 1.3));
    assert_eq!(d.extra_drag.iter().map(|e| e.name.as_str()).collect::<Vec<_>>(), ["wheels"]);

    // --- into a project, twice: the second plane is renamed and its polar follows ---
    let mut proj = Project::new("t");
    for _ in 0..2 {
        import_into(&mut proj, &data("minimal_v6.xfl")).unwrap();
    }
    assert_eq!(proj.planes.iter().map(|p| p.name.as_str()).collect::<Vec<_>>(), ["Mini", "Mini (imported)"]);
    let mut owners: Vec<&str> = proj.plane_polars.iter().map(|q| q.plane.as_str()).collect();
    owners.sort();
    assert_eq!(owners, ["Mini", "Mini (imported)"]);
    assert!(proj.plane_polars.iter().all(|q| q.def.is_some()));
    assert_eq!(proj.foils.len(), 1, "the second import keeps the project's foil");

    // --- AVL, with its AFILE foil next to it ---
    let a = read_any(&data("minimal.avl")).unwrap();
    let w = &a.planes[0].wings;
    assert_eq!(w.iter().map(|w| (w.name.as_str(), w.two_sided)).collect::<Vec<_>>(), [("Wing", true), ("Stab", true)]);
    assert!(close(w[0].sections[0].twist, 1.0) && close(w[0].sections[1].twist, -1.0), "ANGLE adds to every section");
    assert!(close(w[1].frames(false)[0].le().x, 0.9), "TRANSLATE");
    let mini = a.foils.iter().find(|f| f.name == "mini_foil").expect("AFILE foil loaded from the file's folder");
    assert_eq!(mini.hinge, Some((0.75, 0.5)), "hinged for its CONTROL");
    assert!(a.foils.iter().any(|f| f.name == "NACA 0009"));
}

fn has(entries: &[foiley_app::legacy_io::ReportEntry], item: &str, reason: &str) -> bool {
    entries.iter().any(|e| e.item.contains(item) && e.reason.contains(reason))
}

fn well_formed(r: &ConversionReport) {
    for e in r.approximated.iter().chain(&r.skipped) {
        assert!(!e.item.trim().is_empty() && !e.reason.trim().is_empty(), "entry without item or reason: {e:?}");
    }
    let t = r.text();
    assert!(t.contains(&format!("mapped: {}, approximated: {}, skipped: {}", r.mapped.len(), r.approximated.len(), r.skipped.len())), "{t}");
    for e in &r.skipped {
        assert!(t.contains(&e.reason), "{t}");
    }
}

#[test]
fn g45_conversion_report_lists_skips() {
    let x = read_any(&data("minimal_v6.xfl")).unwrap().report;
    well_formed(&x);
    assert!(x.mapped.iter().any(|m| m.contains("plane 'Mini'")) && x.mapped.iter().any(|m| m.contains("'Mini' / 'T2 trimmed' -> T2 VLM2")), "{}", x.text());
    assert!(has(&x.skipped, "'Ghost' / 'orphan'", "its plane was not imported"), "{}", x.text());
    assert!(has(&x.skipped, "tilted geometry", "freestream"), "{}", x.text());
    assert!(has(&x.skipped, "2 control gain(s)", "T7 gains"), "{}", x.text());
    assert!(has(&x.skipped, "1 computed point(s)", "not imported"), "{}", x.text());
    assert!(has(&x.skipped, "foil operating points", "not imported"), "{}", x.text());
    assert!(has(&x.approximated, "foil 'M1'", "flap deflection"), "{}", x.text());
    assert!(has(&x.approximated, "sweep range", "default range"), "{}", x.text());

    let f = read_any(&data("minimal_v7.fl5")).unwrap().report;
    well_formed(&f);
    assert!(has(&f.skipped, "'Nope' / 'orphan'", "its plane was not imported"), "{}", f.text());
    assert!(has(&f.skipped, "bank angle 5 deg", "no bank angle"), "{}", f.text());
    assert!(has(&f.skipped, "1 control range(s)", "T6 ranges"), "{}", f.text());
    assert!(has(&f.skipped, "2 computed point(s)", "not imported"), "{}", f.text());
    assert!(has(&f.skipped, "plane operating points", "not imported"), "{}", f.text());
    assert!(has(&f.approximated, "thin-surface panels", "thick"), "{}", f.text());
    assert!(has(&f.approximated, "computed on the fly", "2D polars"), "{}", f.text());
    assert!(f.mapped.iter().any(|m| m.contains("preset flap angle")), "{}", f.text());

    let a = read_any(&data("minimal.avl")).unwrap().report;
    well_formed(&a);
    for (item, reason) in [("BODY 'Fuse pod'", "slender bodies"), ("NOWAKE", "no Foiley equivalent"), ("CLAF", "2D polars"), ("CDp", "extra drag"), ("Sref", "reference dimensions")] {
        assert!(has(&a.skipped, item, reason), "{item}: {}", a.text());
    }
    assert!(has(&a.approximated, "mini_foil", "hinge set at x/c = 0.75"), "{}", a.text());
    assert!(a.mapped.iter().any(|m| m.contains("CONTROL 'aileron' -> Wing flap 1 (left side, gain -1)")), "{}", a.text());

    // import_into reports name clashes too
    let mut proj = Project::new("t");
    import_into(&mut proj, &data("minimal_v7.fl5")).unwrap();
    let again = import_into(&mut proj, &data("minimal_v7.fl5")).unwrap();
    assert_eq!(again.source, "minimal_v7.fl5");
    assert!(has(&again.approximated, "plane 'Flow mini'", "imported as 'Flow mini (imported)'"), "{}", again.text());
    assert!(has(&again.approximated, "foil 'N1'", "project's one is kept"), "{}", again.text());
}
