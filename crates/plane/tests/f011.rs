//! F011: a stability failure is reported on the point, never dropped.
mod common;
use common::*;
use foiley_plane::analysis::*;
use foiley_plane::geom::Inertia;

/// Mass but no moments of inertia: `stability::analyse` refuses it.
fn no_inertia() -> Inertia {
    Inertia { ixx: 0.0, iyy: 0.0, izz: 0.0, ..review_inertia() }
}

#[test]
fn f011_stability_failure_reported() {
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    let base = PlanePolarDef { free_wake: 0, derivatives: true, inertia: Some(no_inertia()), ..review_def() };
    // T1 (solve_point) and T6 (run_control)
    let t6 = PlanePolarDef { kind: PolarKind::Control, range: Range { from: 0.0, to: 0.0, step: 1.0 }, ..base.clone() };
    for d in [base.clone(), t6] {
        let o = run1(&p, &db, &d);
        assert!(o.stability.is_none());
        let e = o.stability_error.as_deref().unwrap_or_else(|| panic!("{:?}: stability silently missing", d.kind));
        assert!(e.contains("moments of inertia"), "{e}");
        assert!(o.warnings.iter().any(|w| w.contains("stability")), "{:?}", o.warnings);
    }
    // no mass at all
    let massless = PlanePolarDef { inertia: Some(Inertia { mass: 0.0, ..no_inertia() }), ..base.clone() };
    let o = run1(&p, &db, &massless);
    assert!(o.stability_error.as_deref().unwrap().contains("no mass"));
    // a good inertia gives stability and no error
    let ok = run1(&p, &db, &PlanePolarDef { inertia: Some(review_inertia()), ..base });
    assert!(ok.stability.is_some() && ok.stability_error.is_none());
}

/// Body of `fn name(` in `src`, by brace matching.
fn fn_body<'a>(src: &'a str, name: &str) -> &'a str {
    let start = src.find(&format!("fn {name}(")).unwrap_or_else(|| panic!("fn {name} not found"));
    let open = start + src[start..].find('{').unwrap();
    let mut depth = 0;
    for (i, c) in src[open..].char_indices() {
        match c {
            '{' => depth += 1,
            '}' => {
                depth -= 1;
                if depth == 0 {
                    return &src[start..open + i + 1];
                }
            }
            _ => {}
        }
    }
    panic!("unbalanced braces in fn {name}");
}

/// Calls whose `Result`/`Option` carries a solver outcome.
const SOLVER_CALLS: [&str; 7] = ["loads(", "solve", "analyse(", "trim(", "build(", "make_system(", "add_control_derivatives("];

/// Silent result-swallowing in the 3D drivers.  Later items may only lower it.
const MAX_SILENT: usize = 0;

#[test]
fn no_silent_ok() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let analysis = std::fs::read_to_string(root.join("src/analysis.rs")).unwrap();
    let llt = std::fs::read_to_string(root.join("src/llt.rs")).unwrap();
    let mut hits = Vec::new();
    for (file, src, names) in [("analysis.rs", &analysis, &["solve_point", "run_polar", "run_control", "run_stability"][..]), ("llt.rs", &llt, &["run"][..])] {
        for name in names {
            for (k, line) in fn_body(src, name).lines().enumerate() {
                let solver = SOLVER_CALLS.iter().any(|c| line.contains(c));
                let silent = line.contains(".ok()") || line.contains("if let Ok(") || (solver && line.contains("unwrap_or"));
                if silent {
                    hits.push(format!("{file} fn {name} +{k}: {}", line.trim()));
                }
            }
        }
    }
    assert!(hits.len() <= MAX_SILENT, "{} silent result(s) in the 3D drivers:\n{}", hits.len(), hits.join("\n"));
}
