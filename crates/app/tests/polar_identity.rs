//! Regression: rounded polar labels are presentation, never storage identity.
use foiley_app::Project;
use foiley_foil::{AnalysisSettings, Foil, FoilPoint, FoilPolar};

fn polar(foil: &Foil, settings: AnalysisSettings) -> FoilPolar {
    let mut p = FoilPolar {
        name: settings.polar_name(), foil: foil.name.clone(), foil_hash: Some(foil.geom_hash()),
        polar_type: settings.polar_type.number(), re: settings.re, mach: settings.mach,
        ncrit: settings.ncrit, xtr_top: settings.xtr_top, xtr_bot: settings.xtr_bot, viscous: settings.viscous,
        settings: Some(settings), ..Default::default()
    };
    for alpha in [0.0, 1.0] {
        p.insert(FoilPoint { alpha, re: p.re, cl: alpha * 0.1, cd: p.re / 1e8, ..Default::default() });
    }
    p
}

#[test]
fn near_re_runs_survive_merge_save_database_and_export() {
    let foil = Foil::naca("0012").unwrap();
    let a = polar(&foil, AnalysisSettings { re: 300_000.0, ..Default::default() });
    let b = polar(&foil, AnalysisSettings { re: 300_100.0, ..Default::default() });
    assert_eq!(a.name, b.name, "exercise the rounded-name collision");
    assert_ne!(a.identity(), b.identity());
    assert_ne!(a.key(), b.key());
    assert_ne!(a.display_name(), b.display_name());
    assert_ne!(a.export_stem(), b.export_stem());
    let mut project = Project::new("exact Re");
    project.foils.push(foil);
    project.add_foil_polar(a.clone());
    project.add_foil_polar(b.clone());
    assert_eq!(project.foil_polars.len(), 2);
    for p in &project.foil_polars {
        assert!(p.points.iter().all(|pt| pt.re == p.re));
    }
    let dir = std::env::temp_dir().join(format!("foiley_polar_identity_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("near_re.foiley");
    project.save(&path).unwrap();
    let mut loaded = Project::load(&path).unwrap();
    assert_eq!(loaded.foil_polars.len(), 2);
    assert_eq!(loaded.foil_polars[0].key(), a.key());
    assert_eq!(loaded.foil_polars[1].key(), b.key());
    let db = loaded.foil_db();
    assert_eq!(db.families().len(), 1);
    assert_eq!(db.families()[0].polars.len(), 2, "3D database must retain both Re stations");
    for p in &loaded.foil_polars {
        std::fs::write(dir.join(format!("{}.csv", p.export_stem())), p.to_csv()).unwrap();
    }
    assert_eq!(std::fs::read_dir(&dir).unwrap().filter_map(Result::ok).filter(|e| e.path().extension().is_some_and(|e| e == "csv")).count(), 2);
    let mut extension = a.clone();
    extension.name = "renamed display label".into();
    extension.settings.as_mut().unwrap().keep_details = true;
    extension.points = vec![FoilPoint { alpha: 2.0, re: a.re, ..Default::default() }];
    loaded.add_foil_polar(extension);
    assert_eq!(loaded.foil_polars.len(), 2, "same settings extend the existing sweep");
    assert_eq!(loaded.foil_polars[0].points.len(), 3);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn geometry_and_solver_settings_distinguish_same_label() {
    let foil = Foil::naca("0012").unwrap();
    let a = polar(&foil, AnalysisSettings::default());
    let mut p = Project::new("settings");
    p.foils.push(foil.clone());
    p.add_foil_polar(a.clone());
    let mut settings = a.clone();
    settings.settings.as_mut().unwrap().cm_ref = (0.5, 0.0);
    assert_eq!(a.name, settings.name);
    assert_ne!(a.identity(), settings.identity());
    p.add_foil_polar(settings);
    let mut geometry = a.clone();
    geometry.foil_hash = Some(Foil::naca("2412").unwrap().geom_hash());
    p.add_foil_polar(geometry);
    assert_eq!(p.foil_polars.len(), 3);
}

#[test]
fn fractional_reynolds_numbers_remain_exact_in_identity_and_csv() {
    let foil = Foil::naca("0012").unwrap();
    let a = polar(&foil, AnalysisSettings { re: 300_000.125, ..Default::default() });
    let b = polar(&foil, AnalysisSettings { re: 300_000.25, ..Default::default() });
    assert_eq!(a.name, b.name);
    assert_ne!(a.identity(), b.identity());
    assert!(a.to_csv().contains(",300000.125,"));
    assert!(b.to_csv().contains(",300000.25,"));
}

#[test]
fn legacy_named_polars_merge_only_matching_exact_parameters() {
    let foil = Foil::naca("0012").unwrap();
    let mut a = polar(&foil, AnalysisSettings { re: 300_000.0, ..Default::default() });
    a.settings = None;
    a.foil_hash = None;
    let mut p = Project::new("legacy");
    p.foils.push(foil);
    p.add_foil_polar(a.clone());
    let mut extension = a.clone();
    extension.points = vec![FoilPoint { alpha: 2.0, re: a.re, ..Default::default() }];
    p.add_foil_polar(extension);
    assert_eq!(p.foil_polars.len(), 1);
    assert_eq!(p.foil_polars[0].points.len(), 3);
    let mut near = a.clone();
    near.re = 300_100.0;
    p.add_foil_polar(near);
    let mut separately_named = a;
    separately_named.name = "wind tunnel import".into();
    p.add_foil_polar(separately_named);
    assert_eq!(p.foil_polars.len(), 3);
}

#[test]
fn repeated_project_import_renames_foil_and_plane_references_consistently() {
    let foil = Foil::naca("0012").unwrap();
    let mut source = Project::new("source");
    source.foils.push(foil.clone());
    source.add_foil_polar(polar(&foil, AnalysisSettings::default()));
    let plane = foiley_plane::Plane::default_plane(&foil.name, &foil.name);
    source.plane_polars.push(foiley_plane::PlanePolar { plane: plane.name.clone(), name: "analysis".into(), ..Default::default() });
    source.planes.push(plane);
    let path = std::env::temp_dir().join(format!("foiley_import_identity_{}.foiley", std::process::id()));
    source.save(&path).unwrap();
    let mut target = source.clone();
    target.import(&path).unwrap();
    target.import(&path).unwrap();
    assert_eq!((target.foils.len(), target.foil_polars.len(), target.planes.len(), target.plane_polars.len()), (3, 3, 3, 3));
    let keys: std::collections::HashSet<_> = target.foil_polars.iter().map(FoilPolar::key).collect();
    assert_eq!(keys.len(), 3);
    for i in 0..3 {
        assert_eq!(target.foil_polars[i].foil, target.foils[i].name);
        assert_eq!(target.planes[i].foil_names(), vec![target.foils[i].name.clone()]);
        assert_eq!(target.plane_polars[i].plane, target.planes[i].name);
    }
    std::fs::remove_file(path).unwrap();
}
