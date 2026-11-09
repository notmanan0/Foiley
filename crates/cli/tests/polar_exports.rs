use foiley_app::Project;
use serde_json::json;
use std::{fs, process::Command};

#[test]
fn nearby_reynolds_numbers_survive_project_and_csv_export() {
    let root = std::env::temp_dir().join(format!("foiley-cli-polars-{}", std::process::id()));
    fs::create_dir(&root).unwrap();
    for mode in ["script", "script-batch", "batch"] {
        let dir = root.join(mode);
        fs::create_dir(&dir).unwrap();
        let mut cli = Command::new(env!("CARGO_BIN_EXE_foiley"));
        cli.current_dir(&dir).args(["--threads", "2"]);
        if mode == "batch" {
            cli.args(["foil", "batch", "naca0012", "--re", "300000,300100", "--alpha", "1", "--out-dir", "csv", "--save", "out.foiley"]);
        } else {
            let sweep = json!({"AlphaList": [1.0]});
            let mut steps = vec![json!({"naca": "0012"})];
            if mode == "script-batch" {
                steps.push(json!({"batch": {"foils": ["NACA 0012"], "re": [300000, 300100], "sweep": sweep}}));
            } else {
                for re in [300000, 300100] {
                    steps.push(json!({"foil_polar": {"foil": "NACA 0012", "settings": {"re": re}, "sweep": sweep}}));
                }
            }
            fs::write(dir.join("study.json"), json!({"steps": steps, "output_dir": "csv", "save": "out.foiley"}).to_string()).unwrap();
            cli.args(["run", "study.json"]);
        }
        let result = cli.output().unwrap();
        assert!(result.status.success(), "{mode}: {}", String::from_utf8_lossy(&result.stderr));
        let project = Project::load(&dir.join("out.foiley")).unwrap();
        let mut project_re = Vec::new();
        for polar in project.foil_polars {
            assert_eq!(polar.points.len(), 1, "{mode}");
            assert_eq!(polar.points[0].re, polar.re, "{mode}");
            project_re.push(polar.re);
        }
        project_re.sort_by(f64::total_cmp);
        assert_eq!(project_re, [300000.0, 300100.0], "{mode}");

        let mut csv_re = Vec::new();
        for entry in fs::read_dir(dir.join("csv")).unwrap() {
            let csv = fs::read_to_string(entry.unwrap().path()).unwrap();
            let mut rows = csv.lines();
            let col = rows.next().unwrap().split(',').position(|s| s == "Re").unwrap();
            csv_re.push(rows.next().unwrap().split(',').nth(col).unwrap().parse::<u32>().unwrap());
            assert!(rows.next().is_none(), "{mode}");
        }
        csv_re.sort_unstable();
        assert_eq!(csv_re, [300000, 300100], "{mode}");
    }
    fs::remove_dir_all(root).unwrap();
}
