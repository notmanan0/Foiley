//! Regression tests for multi-point wing-section optimisation, output-path
//! safety, hard drag caps, and JSON script validation.

use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

fn temp_dir(tag: &str) -> PathBuf {
    static N: AtomicU64 = AtomicU64::new(0);
    let d = std::env::temp_dir().join(format!("foiley_opt_test_{}_{}_{}", tag, std::process::id(), N.fetch_add(1, Ordering::Relaxed)));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn foiley() -> Command {
    Command::new(env!("CARGO_BIN_EXE_foiley"))
}

#[test]
fn save_does_not_delete_tmp_input_or_pareto_or_leave_device_junk() {
    let dir = temp_dir("save_guards");
    // Build a minimal valid .foiley project stored at `study.foiley.tmp`
    let in_tmp = dir.join("study.foiley.tmp");
    let proj_json = r#"{"format":3,"name":"t","foils":[],"foil_polars":[],"planes":[],"plane_polars":[]}"#;
    std::fs::write(&in_tmp, proj_json).unwrap();
    let dest = dir.join("study.foiley");
    let out_dat = dir.join("best.dat");

    // 1. --project study.foiley.tmp --save study.foiley must preserve the input file!
    let r = foiley()
        .args([
            "foil",
            "optimize",
            "naca2412",
            "--project",
            in_tmp.to_str().unwrap(),
            "--op",
            "Cruise:1.8e5:0.50:min",
            "--shape",
            "cst",
            "--cst-order",
            "2",
            "--population",
            "2",
            "--iterations",
            "0",
            "--save",
            dest.to_str().unwrap(),
            "-o",
            out_dat.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(r.status.success(), "stderr: {}", String::from_utf8_lossy(&r.stderr));
    assert!(in_tmp.exists(), "input study.foiley.tmp was deleted!");
    assert_eq!(std::fs::read_to_string(&in_tmp).unwrap(), proj_json);
    assert!(dest.exists() && out_dat.exists());

    // 2. --pareto study2.foiley.tmp --save study2.foiley must keep the CSV intact!
    let dest2 = dir.join("study2.foiley");
    let pareto_tmp = dir.join("study2.foiley.tmp");
    let out_dat2 = dir.join("best2.dat");
    let r2 = foiley()
        .args([
            "foil",
            "optimize",
            "naca2412",
            "--op",
            "Cruise:1.8e5:0.50:min",
            "--shape",
            "cst",
            "--cst-order",
            "2",
            "--population",
            "2",
            "--iterations",
            "0",
            "--pareto",
            pareto_tmp.to_str().unwrap(),
            "--save",
            dest2.to_str().unwrap(),
            "-o",
            out_dat2.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(r2.status.success(), "stderr: {}", String::from_utf8_lossy(&r2.stderr));
    assert!(pareto_tmp.exists(), "pareto CSV at study2.foiley.tmp was lost!");
    assert!(std::fs::read_to_string(&pareto_tmp).unwrap().starts_with("index,violation,spar_tc"));

    // 3. Invalid --out parent dir must fail up front without creating --save
    let bad_out = dir.join("missing_dir").join("best.dat");
    let untouched_save = dir.join("should_not_exist.foiley");
    let r3 = foiley()
        .args([
            "foil",
            "optimize",
            "naca2412",
            "--op",
            "Cruise:1.8e5:0.50:min",
            "--population",
            "2",
            "--iterations",
            "0",
            "--save",
            untouched_save.to_str().unwrap(),
            "-o",
            bad_out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!r3.status.success());
    assert!(!untouched_save.exists(), "--save must not be written when --out fails preflight");

    // 4. --save NUL on Windows must fail up front and leave no NUL.foiley.tmp
    #[cfg(windows)]
    {
        let r4 = foiley()
            .current_dir(&dir)
            .args(["foil", "optimize", "naca2412", "--op", "Cruise:1.8e5:0.50:min", "--population", "2", "--iterations", "0", "--save", "NUL"])
            .output()
            .unwrap();
        assert!(!r4.status.success());
        assert!(!dir.join("NUL.foiley.tmp").exists());
    }

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn impossible_cap_and_bad_args_and_script_guards_fail_cleanly() {
    let dir = temp_dir("opt_errors");
    let out_dat = dir.join("impossible.dat");
    let out_pareto = dir.join("impossible.csv");
    let out_save = dir.join("impossible.foiley");

    // Impossible drag cap must exit non-zero and produce NO output files
    let r = foiley()
        .args([
            "foil",
            "optimize",
            "naca2412",
            "--op",
            "Cap:1.8e5:0.50:<=1e-7",
            "--shape",
            "cst",
            "--cst-order",
            "2",
            "--population",
            "3",
            "--iterations",
            "1",
            "-o",
            out_dat.to_str().unwrap(),
            "--pareto",
            out_pareto.to_str().unwrap(),
            "--save",
            out_save.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!r.status.success(), "impossible cap should fail");
    assert!(!out_dat.exists() && !out_pareto.exists() && !out_save.exists());

    // Extra trailing --op field is rejected
    let r2 = foiley()
        .args(["foil", "optimize", "naca2412", "--op", "cruise:175000:0.50:min:1:extra"])
        .output()
        .unwrap();
    assert!(!r2.status.success());

    // Script guards: copy_plane duplicate name, out-of-range wing, and unmatched foil all fail
    let base_proj = dir.join("base.foiley");
    let p = foiley_app::Project::starter();
    let existing_plane = p.planes[0].name.clone();
    p.save(&base_proj).unwrap();

    let dup_plane_step = format!(r#"{{"copy_plane":{{"name":"{existing_plane}"}}}}"#);
    for (tag, step_json) in [
        ("dup_plane", dup_plane_step.as_str()),
        ("bad_wing", r#"{"copy_plane":{"name":"Glider 2","wing":99}}"#),
        ("no_match", r#"{"copy_plane":{"name":"Glider 2","wing":0,"replace_foil":["NoSuchFoil","NACA 0009"]}}"#),
        ("dup_foil", r#"{"optimize":{"foil":"NACA 2412","name":"NACA 2412","settings":{"population":2,"iterations":0}}}"#),
        ("empty_3d", r#"{"optimize":{"foil":"NACA 2412","from_3d":{},"settings":{"population":2,"iterations":0}}}"#),
    ] {
        let sc_path = dir.join(format!("{tag}.json"));
        let sc = format!(r#"{{"project":"{}","steps":[{step_json}]}}"#, base_proj.to_string_lossy().replace('\\', "\\\\"));
        std::fs::write(&sc_path, sc).unwrap();
        let out = foiley().args(["run", sc_path.to_str().unwrap()]).output().unwrap();
        assert!(!out.status.success(), "script '{tag}' unexpectedly succeeded: {}", String::from_utf8_lossy(&out.stderr));
    }

    let _ = std::fs::remove_dir_all(&dir);
}
