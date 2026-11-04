//! F021: coverage-aware automatic 2D polars.
use foiley_app::auto::{self, deflected_name, AutoPlan};
use foiley_app::Project;
use foiley_core::Progress;
use foiley_foil::batch::{re_list_geometric, run_batch};
use foiley_foil::{Foil, FoilPoint, FoilPolar};
use foiley_plane::analysis::{ControlRange, ControlTarget, PolarKind, Range};
use foiley_plane::foildb::Envelope;
use foiley_plane::{run_polar, FamilySelector, Method, PlanePolarDef, WingKind};

fn synth(foil: &str, re: f64, mach: f64, ncrit: f64) -> FoilPolar {
    FoilPolar {
        foil: foil.into(),
        name: format!("synth Re{re:.0} M{mach} N{ncrit}"),
        re,
        mach,
        ncrit,
        xtr_top: 1.0,
        xtr_bot: 1.0,
        polar_type: 1,
        viscous: true,
        points: (-20..=20).map(|a| FoilPoint { alpha: a as f64, cl: 0.1 * a as f64, cd: 0.01, re, ..Default::default() }).collect(),
        ..Default::default()
    }
}

/// Run every batch of a plan and file the polars in the project.
fn execute(p: &mut Project, plan: &AutoPlan) {
    for b in &plan.batches {
        for r in run_batch(&b.foils, &b.cases, &b.settings, &b.sweep, &Progress::none(), &|_| {}) {
            if let Ok(q) = r.polar {
                p.add_foil_polar(q);
            }
        }
    }
}

/// The starter glider with a hinged elevator foil.
fn flapped_project() -> Project {
    let mut p = Project::starter();
    let mut ef = Foil::naca("0009").unwrap();
    ef.name = "elevator 0009".into();
    ef.hinge = Some((0.7, 0.5));
    for w in p.planes[0].wings.iter_mut().filter(|w| w.kind == WingKind::Elevator) {
        for s in &mut w.sections {
            s.foil_left = ef.name.clone();
            s.foil_right = ef.name.clone();
        }
    }
    p.foils.push(ef);
    p
}

#[test]
fn g08_auto_covers_t8_envelope() {
    let mut p = Project::starter();
    let plane = p.planes[0].clone();
    let def = PlanePolarDef { kind: PolarKind::List, method: Method::Vlm2, list: vec![(0.0, 0.0, 8.0), (4.0, 0.0, 18.0), (6.0, 0.0, 30.0)], ..Default::default() };
    let plan = auto::plan(&plane, &p.foils, &p.foil_db(), &def, 5).unwrap();
    let cov = plan.coverage.clone().unwrap();
    assert_eq!(cov.states, 3);
    // Re per strip over the min / max chords and the slowest / fastest state
    let chords: Vec<f64> = plane.wings.iter().flat_map(|w| w.sections.iter().map(|s| s.chord)).collect();
    let nu = def.fluid.nu();
    let (cmin, cmax) = (chords.iter().cloned().fold(f64::MAX, f64::min), chords.iter().cloned().fold(0.0, f64::max));
    let lo = cov.foils.iter().map(|c| c.re.0).fold(f64::MAX, f64::min);
    let hi = cov.foils.iter().map(|c| c.re.1).fold(0.0, f64::max);
    assert!(lo <= 8.0 * cmin / nu && hi >= 30.0 * cmax / nu * 0.95, "Re envelope {lo:.0}..{hi:.0} misses 8 m/s x {cmin} m .. 30 m/s x {cmax} m");
    assert!(!plan.is_empty());
    execute(&mut p, &plan);
    // the 3D run then never clamps in Re and finds a polar for every strip
    let o = run_polar(&plane, &p.foil_db(), &def, &Progress::none(), &|_| {}).unwrap();
    assert_eq!(o.points.len(), 3, "{:?}", o.failed);
    for pt in &o.points {
        for s in &pt.strips {
            let e = Envelope::from_u8(s.envelope);
            assert!(!e.re_clamped() && e != Envelope::Missing, "V={} strip y={:.3} Re={:.0}: {}", pt.velocity, s.span_pos, s.re, e.label());
        }
    }
    // and nothing is left to generate
    let again = auto::plan(&plane, &p.foils, &p.foil_db(), &def, 5).unwrap();
    assert!(again.is_empty(), "second plan still wants {}", again.describe());
}

#[test]
fn f021_generates_deflected_families() {
    let mut p = flapped_project();
    // undeflected families that already cover every Re station
    for f in p.foils.clone() {
        for re in re_list_geometric(1e4, 1e7, 40) {
            p.add_foil_polar(synth(&f.name, re, 0.0, 9.0));
        }
    }
    let plane = p.planes[0].clone();
    let iw = plane.wings.iter().position(|w| w.kind == WingKind::Elevator).unwrap();
    let ranges = (0..2).map(|index| ControlRange { target: ControlTarget::Flap { wing: iw, index }, min: -10.0, max: 10.0 }).collect();
    let def = PlanePolarDef { kind: PolarKind::Control, method: Method::Vlm2, velocity: 15.0, alpha: 3.0, range: Range { from: 0.0, to: 1.0, step: 0.5 }, ranges, ..Default::default() };
    let plan = auto::plan(&plane, &p.foils, &p.foil_db(), &def, 2).unwrap();
    // only the deflected elevator geometries, at the two nonzero angles
    let mut names: Vec<String> = plan.batches.iter().flat_map(|b| b.foils.iter().map(|f| f.name.clone())).collect();
    names.sort();
    let mut want = vec![deflected_name("elevator 0009", -10.0), deflected_name("elevator 0009", 10.0)];
    want.sort();
    assert_eq!(names, want);
    let e = p.foil("elevator 0009").unwrap().clone();
    for b in &plan.batches {
        assert!(b.cases.iter().all(|c| c.mach == 0.0 && c.ncrit == 9.0));
        for f in &b.foils {
            assert_ne!(f.geom_hash(), e.geom_hash(), "deflected copy has the undeflected geometry");
        }
    }
    execute(&mut p, &plan);
    let db = p.foil_db();
    let base = db.resolve("elevator 0009", &FamilySelector::default()).unwrap().unwrap();
    let flaps: Vec<f64> = db.flap_families(base).iter().map(|x| x.0).collect();
    assert_eq!(flaps, vec![-10.0, 0.0, 10.0]);
    assert!(p.foil_polars.iter().filter(|q| q.flap_deg != 0.0).all(|q| q.foil == "elevator 0009" && q.foil_hash == Some(e.geom_hash())));
    // the control polar's deflected strips now find their family
    let o = run_polar(&plane, &db, &def, &Progress::none(), &|_| {}).unwrap();
    assert_eq!(o.points.len(), 3, "{:?}", o.failed);
    assert!(o.points.iter().all(|pt| pt.strips.iter().all(|s| !s.flap_clamped)));
    assert!(auto::plan(&plane, &p.foils, &db, &def, 2).unwrap().is_empty());
}

#[test]
fn f021_uses_selected_family_mach_ncrit() {
    let mut p = Project::starter();
    for f in p.foils.clone() {
        // M 0.1 / Ncrit 5 covers only low Re; Ncrit 9 everything
        for re in [3e4, 1e5] {
            p.add_foil_polar(synth(&f.name, re, 0.1, 5.0));
        }
        for re in re_list_geometric(1e4, 1e7, 40) {
            p.add_foil_polar(synth(&f.name, re, 0.0, 9.0));
        }
    }
    let plane = p.planes[0].clone();
    let db = p.foil_db();
    let base = PlanePolarDef { method: Method::Vlm2, velocity: 20.0, range: Range { from: 0.0, to: 6.0, step: 3.0 }, ..Default::default() };
    // two families and no selection: ambiguous
    assert!(auto::plan(&plane, &p.foils, &db, &base, 4).is_err());
    // the selected family's Mach and Ncrit, missing stations only
    let def = PlanePolarDef { family: FamilySelector { ncrit: Some(5.0), ..Default::default() }, ..base.clone() };
    let plan = auto::plan(&plane, &p.foils, &db, &def, 4).unwrap();
    assert!(!plan.is_empty());
    for b in &plan.batches {
        assert_eq!((b.settings.mach, b.settings.ncrit, b.settings.xtr_top, b.settings.xtr_bot), (0.1, 5.0, 1.0, 1.0));
        assert!(b.cases.iter().all(|c| c.mach == 0.1 && c.ncrit == 5.0 && c.re > 1.2e5), "{:?}", b.cases);
    }
    // the Ncrit 9 family covers everything already
    let n9 = PlanePolarDef { family: FamilySelector { ncrit: Some(9.0), ..Default::default() }, ..base.clone() };
    assert!(auto::plan(&plane, &p.foils, &db, &n9, 4).unwrap().is_empty());
    // a selected family that does not exist yet is generated with its values
    let n7 = PlanePolarDef { family: FamilySelector { ncrit: Some(7.0), ..Default::default() }, ..base };
    let plan = auto::plan(&plane, &p.foils, &db, &n7, 4).unwrap();
    assert_eq!(plan.polar_count(), 2 * 4);
    assert!(plan.batches.iter().all(|b| b.cases.iter().all(|c| c.ncrit == 7.0 && c.mach == 0.0)));
}
