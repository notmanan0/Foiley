use std::path::Path;
use std::process::{Command, Output};

use foiley_plane::analysis::G;
use foiley_plane::PlanePolar;

fn analyze(project: &str, elevator: &str, args: &[&str]) -> Output {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data/projects")
        .join(project);
    Command::new(env!("CARGO_BIN_EXE_foiley"))
        .args(["--threads", "2", "--quiet", "plane", "analyze"])
        .arg(fixture)
        .args(["--type", "t9", "--elevator", elevator, "-f", "json"])
        .args(args)
        .output()
        .unwrap()
}

fn polar(output: Output) -> PlanePolar {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}

#[test]
fn elevator_trim_balances_pitch_and_glide_force_at_requested_speeds() {
    let p = polar(analyze(
        "flapped.foiley",
        "1,-10,10",
        &["--inviscid", "--range", "8:16:4"],
    ));
    assert!(p.failed.is_empty(), "{:?}", p.failed);
    assert_eq!(
        p.points.iter().map(|o| o.velocity).collect::<Vec<_>>(),
        [8.0, 12.0, 16.0]
    );
    for o in &p.points {
        assert!((-10.0..=10.0).contains(&o.elevator_deg.unwrap()));
        assert!(o.cm.abs() <= 1e-7, "Cm = {}", o.cm);
        assert!((o.lift.hypot(o.drag) / (o.mass * G) - 1.0).abs() <= 1e-6);
        assert!(o.balance_residual.unwrap() <= 1e-6);
        assert!(o.vz() > 0.0 && o.ld() > 0.0);
    }
    assert_ne!(p.points[0].elevator_deg, p.points[2].elevator_deg);
    assert!(p
        .to_csv()
        .lines()
        .next()
        .unwrap()
        .ends_with(",elevator_deg"));
    assert!(p.trimmed_glide_extrema().is_some());
}

#[test]
fn trim_reports_elevator_travel_and_2d_coverage_failures() {
    let limited = polar(analyze(
        "flapped.foiley",
        "1,8,10",
        &["--inviscid", "--range", "12:12:1"],
    ));
    assert!(limited.points.is_empty());
    assert_eq!(limited.failed.len(), 1);
    assert!(
        limited.failed[0].contains("V=12.000 m/s") && limited.failed[0].contains("elevator travel")
    );

    let uncovered = polar(analyze(
        "flapped.foiley",
        "1,-10,10",
        &["--auto-polars", "0", "--range", "12:12:1"],
    ));
    assert!(uncovered.points.is_empty());
    assert_eq!(uncovered.failed.len(), 1);
    assert!(uncovered.failed[0].contains("outside the 2D polars"));
}

#[test]
fn trim_requires_a_hinged_elevator_and_positive_speeds() {
    let missing = analyze(
        "starter.foiley",
        "1,-10,10",
        &["--inviscid", "--range", "12:12:1"],
    );
    assert!(!missing.status.success());
    assert!(String::from_utf8_lossy(&missing.stderr).contains("has no hinged flap"));

    let bad_speed = analyze(
        "flapped.foiley",
        "1,-10,10",
        &["--inviscid", "--range", "-4:10:1"],
    );
    assert!(!bad_speed.status.success());
    assert!(String::from_utf8_lossy(&bad_speed.stderr)
        .contains("airspeeds must be finite and positive"));
}
