#![cfg(feature = "snap")]

use foiley_gui::snap::{self, Options, Preset, SnapError};
use serde_json::Value;
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn starter_views_show_the_aircraft_inside_the_canvas() {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let out = std::env::temp_dir().join(format!("foiley-views-{}-{stamp}", std::process::id()));
    let files = match snap::run(Options {
        preset: Some(Preset::Views),
        out: Some(out.clone()),
        grid: Some(false),
        ..Options::default()
    }) {
        Ok(files) => files,
        Err(SnapError::Adapter(error)) if !cfg!(windows) => {
            eprintln!("Skipping snap views: {error}");
            return;
        }
        Err(error) => panic!("{error}"),
    };
    assert_eq!(files.len(), 10);
    let manifest: Value =
        serde_json::from_slice(&std::fs::read(files.last().unwrap()).unwrap()).unwrap();
    assert_eq!(manifest["viewport_renderer"], "cpu");
    let captures = manifest["captures"].as_array().unwrap();
    assert_eq!(captures.len(), 8);
    let expected = [
        (135.0, 0.52_f64.to_degrees()),
        (0.0, 90.0),
        (180.0, 0.0),
        (90.0, 0.0),
    ];
    for (i, capture) in captures.iter().enumerate() {
        let file = capture["file"].as_str().unwrap();
        let viewport = &capture["viewport"];
        assert_eq!(viewport["animating"], false, "{file}");
        assert_eq!(viewport["perspective"], i < 4, "{file}");
        assert!(
            (viewport["yaw"].as_f64().unwrap().to_degrees() - expected[i % 4].0).abs() < 1e-6,
            "{file}"
        );
        assert!(
            (viewport["pitch"].as_f64().unwrap().to_degrees() - expected[i % 4].1).abs() < 1e-6,
            "{file}"
        );
        let rect: Vec<_> = viewport["rect"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap().round() as u32)
            .collect();
        let image = image::open(out.join(file)).unwrap().to_rgb8();
        let (width, height) = image.dimensions();
        assert_eq!(
            (width, height),
            (rect[2] - rect[0], rect[3] - rect[1]),
            "{file}"
        );
        let mut bounds = [width, height, 0, 0];
        let (mut sx, mut sy, mut count) = (0.0, 0.0, 0);
        for (x, y, pixel) in image.enumerate_pixels() {
            // Exclude viewport chrome and count the blue surfaces on the slate backdrop.
            if y < 38 || y + 24 > height || (x < 97 && y + 110 > height) {
                continue;
            }
            let [r, g, b] = pixel.0;
            if b.saturating_sub(r) > 33 && b.saturating_sub(g) > 13 && r > 24 {
                bounds = [
                    bounds[0].min(x),
                    bounds[1].min(y),
                    bounds[2].max(x),
                    bounds[3].max(y),
                ];
                sx += f64::from(x);
                sy += f64::from(y);
                count += 1;
            }
        }
        assert!(count > 70, "no aircraft surface in {file}");
        assert!(
            bounds[0] > 2 && bounds[1] > 2 && bounds[2] + 3 < width && bounds[3] + 3 < height,
            "{file}: aircraft reaches canvas edge"
        );
        let (cx, cy) = (
            sx / f64::from(count) / f64::from(width),
            sy / f64::from(count) / f64::from(height),
        );
        assert!(
            (0.25..=0.75).contains(&cx) && (0.25..=0.75).contains(&cy),
            "{file}: aircraft centroid ({cx:.2}, {cy:.2})"
        );
    }
    for file in files {
        std::fs::remove_file(file).unwrap();
    }
    std::fs::remove_dir(out).unwrap();
}

#[test]
fn a_preset_that_needs_a_wing_reports_a_plane_without_one() {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("foiley-wingless-{}-{stamp}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut plane = foiley_plane::Plane::default_plane("NACA 0009", "NACA 0009");
    plane.wings.clear();
    plane.fuses.push(foiley_plane::geom::Fuse::default_pod());
    let mut project = foiley_app::Project::new("body only");
    project.planes.push(plane);
    let file = dir.join("body_only.foiley");
    project.save(&file).unwrap();
    let result = snap::run(Options {
        launch: foiley_gui::Launch {
            files: vec![file],
            ..Default::default()
        },
        preset: Some(Preset::Narrow),
        out: Some(dir.join("out")),
        ..Options::default()
    });
    std::fs::remove_dir_all(&dir).ok();
    match result {
        Err(SnapError::Adapter(error)) if !cfg!(windows) => {
            eprintln!("Skipping wingless snap: {error}");
        }
        // reported as an input error (exit status 1), never a panic
        Err(SnapError::Input(reason)) => assert!(reason.contains("wing"), "{reason}"),
        other => panic!("expected an input error, got {other:?}"),
    }
}
