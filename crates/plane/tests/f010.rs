//! F010: AVL export fidelity and importer.
mod common;
use common::*;
use foiley_core::Vec3;
use foiley_foil::Foil;
use foiley_plane::analysis::ControlTarget;
use foiley_plane::export::{avl, avl_with_gains, import_avl};
use foiley_plane::geom::{Plane, Wing, WingKind, WingSection};
use foiley_plane::FoilDb;

#[test]
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
    assert_eq!(avl.matches("\nSURFACE\n").count(), 2, "one surface per side");
    // the symmetric original keeps its single mirrored surface
    let sym = foiley_plane::export::avl(&review_rect(&f), &db, 0.0);
    assert_eq!((sym.matches("\nSURFACE\n").count(), sym.matches("YDUPLICATE\n0.0").count()), (1, 1));
}

/// Every kind of wing the exporter distinguishes.
fn zoo() -> (Plane, FoilDb) {
    let names = ["2412", "0012", "0009", "4415"];
    let mut db = FoilDb::new();
    for n in names {
        db.add_foil(&Foil::naca(n).unwrap());
    }
    let mut p = Plane::default_plane("NACA 2412", "NACA 0009");
    // main wing: three sections, sweep, taper, dihedral break and washout, mirrored
    let main = &mut p.wings[0];
    let mut mid = WingSection::new(0.8, 0.26, 0.03, "NACA 2412");
    mid.dihedral = 4.0;
    mid.twist = -0.7;
    main.sections.insert(1, mid);
    main.sections[0].dihedral = 1.5;
    main.sections[2].twist = -2.5;
    main.sections[2].foil_right = "NACA 4415".into();
    main.sections[2].foil_left = "NACA 4415".into();
    main.ry = 1.25;
    // asymmetric canard: centred, other foil on the left
    let mut canard = Wing::new("canard", WingKind::Other, vec![WingSection::new(0.0, 0.12, 0.0, "NACA 0012"), WingSection::new(0.35, 0.08, 0.02, "NACA 0012")]);
    canard.position = Vec3::new(-0.45, 0.0, 0.04);
    canard.symmetric = false;
    canard.sections[0].dihedral = 3.0;
    canard.sections[1].twist = 1.0;
    for s in &mut canard.sections {
        s.foil_left = "NACA 0009".into();
    }
    p.wings.push(canard);
    // offset two-sided wing (mirrored about its own y = 0.6 plane), with a gap at the root
    let mut pod = Wing::new("pod wing", WingKind::Other, vec![WingSection::new(0.05, 0.1, 0.0, "NACA 0012"), WingSection::new(0.25, 0.07, 0.015, "NACA 0012")]);
    pod.position = Vec3::new(0.3, 0.6, -0.05);
    pod.sections[0].dihedral = -6.0;
    pod.sections[1].twist = 2.0;
    pod.ry = -0.5;
    p.wings.push(pod);
    (p, db)
}

fn sides(w: &Wing) -> Vec<bool> {
    if w.two_sided {
        vec![false, true]
    } else {
        vec![false]
    }
}

#[test]
fn f010_avl_roundtrip_1e9() {
    let (p, db) = zoo();
    let text = avl(&p, &db, 0.0);
    let imp = import_avl(&text, &[]).unwrap();
    assert_eq!(imp.plane.wings.len(), p.wings.len(), "{}", imp.report.text());
    for (a, b) in p.wings.iter().zip(&imp.plane.wings) {
        assert_eq!(a.name, b.name);
        assert_eq!(a.kind, b.kind, "{}", a.name);
        assert_eq!(a.two_sided, b.two_sided, "{}", a.name);
        assert_eq!(a.sections.len(), b.sections.len(), "{}", a.name);
        for left in sides(a) {
            let (fa, fb) = (a.frames(left), b.frames(left));
            for k in 0..a.sections.len() {
                let what = format!("{} section {} {}", a.name, k + 1, if left { "left" } else { "right" });
                let d = (fa[k].le() - fb[k].le()).norm();
                assert!(d < 1e-9, "{what}: LE {:?} vs {:?} ({d:e})", fa[k].le(), fb[k].le());
                assert!((a.sections[k].chord - b.sections[k].chord).abs() < 1e-9, "{what}: chord");
                let (ia, ib) = (a.sections[k].twist + a.ry, b.sections[k].twist + b.ry);
                assert!((ia - ib).abs() < 1e-9, "{what}: incidence {ia} vs {ib}");
                assert_eq!(a.foil(k, left), b.foil(k, left), "{what}: foil");
            }
        }
        assert_eq!(a.sections.iter().map(|s| s.ny).collect::<Vec<_>>(), b.sections.iter().map(|s| s.ny).collect::<Vec<_>>(), "{}: spanwise panels", a.name);
        assert_eq!(a.sections[0].nx, b.sections[0].nx, "{}: chordwise panels", a.name);
    }
    // a second trip changes nothing
    let mut db2 = FoilDb::new();
    for f in &imp.foils {
        db2.add_foil(f);
    }
    let again = import_avl(&avl(&imp.plane, &db2, 0.0), &[]).unwrap();
    for (a, b) in imp.plane.wings.iter().zip(&again.plane.wings) {
        for left in sides(a) {
            for (x, y) in a.frames(left).iter().zip(b.frames(left)) {
                assert!((x.le() - y.le()).norm() < 1e-12 && (x.te() - y.te()).norm() < 1e-12, "{}", a.name);
            }
        }
    }
    // only the header items that belong to the analysis definition are reported
    assert!(imp.report.approximated.is_empty(), "{}", imp.report.text());
    assert!(imp.report.skipped.iter().all(|e| e.item.starts_with("Sref") || e.item.starts_with("Xref")), "{}", imp.report.text());
}

fn hinged(code: &str, name: &str, xh: f64) -> Foil {
    let mut f = Foil::naca(code).unwrap();
    f.name = name.into();
    f.hinge = Some((xh, 0.5));
    f
}

/// Flapped main wing (two flaps per side), elevator, rudder.
fn flapped() -> (Plane, FoilDb, Vec<Foil>) {
    let foils = vec![hinged("2412", "wing flapped", 0.75), hinged("0009", "tail hinged", 0.7)];
    let mut db = FoilDb::new();
    for f in &foils {
        db.add_foil(f);
    }
    let mut p = Plane::default_plane("wing flapped", "tail hinged");
    let main = &mut p.wings[0];
    main.split_section(0);
    (p, db, foils)
}

fn controls(avl: &str) -> Vec<(String, f64, f64, f64)> {
    let l: Vec<&str> = avl.lines().collect();
    l.iter()
        .enumerate()
        .filter(|(_, x)| **x == "CONTROL")
        .map(|(i, _)| {
            let t: Vec<&str> = l[i + 2].split_whitespace().collect();
            (t[0].to_string(), t[1].parse().unwrap(), t[2].parse().unwrap(), t[6].parse().unwrap())
        })
        .collect()
}

#[test]
fn f010_control_signs_and_names() {
    let (p, db, foils) = flapped();
    let flap = |index| ControlTarget::Flap { wing: 0, index };
    // flap numbering: left side tip to root (0, 1), then right side root to tip (2, 3);
    // the outer pair is mixed as ailerons, the inner pair as flaps
    let gains = vec![(flap(0), -1.0), (flap(3), 1.0), (flap(1), 0.5), (flap(2), 0.5)];
    let text = avl_with_gains(&p, &db, 0.0, &gains);
    assert_eq!(text.matches("YDUPLICATE").count(), 2, "main wing and elevator stay mirrored:\n{text}");
    let c = controls(&text);
    let mut names: Vec<&str> = c.iter().map(|x| x.0.as_str()).collect();
    names.sort();
    names.dedup();
    assert_eq!(names, ["Elevator_elevator1", "Fin_rudder1", "Main_wing_aileron2", "Main_wing_flap1"], "{text}");
    for (name, gain, xh, sgn) in &c {
        assert_eq!(c.iter().filter(|x| &x.0 == name).count(), 2, "{name} spans exactly two sections");
        let (g, s, h) = match name.as_str() {
            "Main_wing_aileron2" => (1.0, -1.0, 0.75),
            "Main_wing_flap1" => (0.5, 1.0, 0.75),
            "Elevator_elevator1" => (1.0, 1.0, 0.7),
            _ => (1.0, 1.0, 0.7),
        };
        assert_eq!((*gain, *sgn, *xh), (g, s, h), "{name}");
    }
    // the importer reads the mixing back
    let imp = import_avl(&text, &foils).unwrap();
    let mut mix: Vec<(String, ControlTarget, f64)> = imp.mixing.clone();
    mix.sort_by(|a, b| format!("{:?}", a.1).cmp(&format!("{:?}", b.1)));
    let main: Vec<_> = mix.iter().filter(|m| matches!(m.1, ControlTarget::Flap { wing: 0, .. })).map(|m| (m.1, m.2)).collect();
    assert_eq!(main, vec![(flap(0), -1.0), (flap(1), 0.5), (flap(2), 0.5), (flap(3), 1.0)]);
    assert!(imp.report.approximated.is_empty(), "{}", imp.report.text());

    // unequal left/right gains are not a mirror image: explicit halves, one control per flap
    let lopsided = vec![(flap(0), -2.0), (flap(3), 1.0)];
    let text = avl_with_gains(&p, &db, 0.0, &lopsided);
    assert_eq!(text.matches("YDUPLICATE").count(), 1, "only the elevator stays mirrored");
    let c = controls(&text);
    let mut names: Vec<&str> = c.iter().map(|x| x.0.as_str()).collect();
    names.dedup();
    for n in ["Main_wing_flap1_R", "Main_wing_flap2_R", "Main_wing_flap1_L", "Main_wing_flap2_L"] {
        assert!(names.contains(&n), "{n} missing: {names:?}");
    }
    assert!(c.iter().all(|x| x.3 == 1.0), "no SgnDup without YDUPLICATE");
    let tip_l = c.iter().find(|x| x.0 == "Main_wing_flap2_L").unwrap();
    assert_eq!(tip_l.1, -2.0);
    let imp = import_avl(&text, &foils).unwrap();
    assert!(imp.mixing.iter().any(|m| m.1 == flap(0) && m.2 == -2.0) && imp.mixing.iter().any(|m| m.1 == flap(3) && m.2 == 1.0), "{:?}", imp.mixing);
}
