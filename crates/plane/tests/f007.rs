//! F007 / F019: one wake-aware evaluator for every driver.
mod common;
use common::*;
use foiley_plane::analysis::*;

#[test]
fn f007_t2_lift_equals_weight() {
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    let t2 = run1(&p, &db, &PlanePolarDef { kind: PolarKind::FixedLift, ..review_def() });
    let r = t2.lift / (review_inertia().mass * G);
    assert!((r - 1.0).abs() < 1e-4, "lift/weight = {r}");
    assert!(!t2.wake.is_empty(), "relaxed T2 keeps no wake lines");
}

#[test]
fn f007_t3_relaxed_differs_from_rigid() {
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    let t3 = run1(&p, &db, &PlanePolarDef { kind: PolarKind::Glide, ..review_def() });
    let rigid = run1(&p, &db, &PlanePolarDef { kind: PolarKind::Glide, free_wake: 0, ..review_def() });
    assert!((t3.cl - rigid.cl).abs() > 1e-6, "delta_CL = {}", t3.cl - rigid.cl);
    assert!(!t3.wake.is_empty(), "relaxed T3 keeps no wake lines");
    assert!(t3.balance_residual.unwrap() < GLIDE_TOL);
}

#[test]
fn f019_t1_t6_equivalent() {
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    let t1 = run1(&p, &db, &review_def());
    let t6 = run1(&p, &db, &PlanePolarDef { kind: PolarKind::Control, range: Range { from: 0.0, to: 0.0, step: 1.0 }, ..review_def() });
    assert!((t1.cl - t6.cl).abs() < 1e-9, "CL T1 {} vs T6 {}", t1.cl, t6.cl);
    assert!((t1.cm - t6.cm).abs() < 1e-9, "Cm T1 {} vs T6 {}", t1.cm, t6.cm);
    assert!((t1.cd - t6.cd).abs() < 1e-9, "CD T1 {} vs T6 {}", t1.cd, t6.cd);
    assert_eq!(t1.wake.len(), t6.wake.len());
    // and with derivatives (frozen relaxed wake) the stability results agree
    let d = PlanePolarDef { derivatives: true, ..review_def() };
    let s1 = run1(&p, &db, &d).stability.unwrap();
    let s6 = run1(&p, &db, &PlanePolarDef { kind: PolarKind::Control, range: Range { from: 0.0, to: 0.0, step: 1.0 }, ..d }).stability.unwrap();
    assert!((s1.dim.zw - s6.dim.zw).abs() < 1e-9 * s1.dim.zw.abs(), "Zw T1 {} vs T6 {}", s1.dim.zw, s6.dim.zw);
}

#[test]
fn f007_wake_report_residuals() {
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    // converged: both residuals below the tolerance, in fewer steps than the cap
    let ok = run1(&p, &db, &PlanePolarDef { free_wake: 100, wake_tol: 1e-3, ..review_def() });
    let w = ok.diagnostics.as_ref().and_then(|d| d.wake.clone()).expect("relaxed point carries its wake report");
    assert!(w.relaxed && w.converged, "{w:?}");
    assert!(w.iterations >= 1 && w.iterations < 100, "{w:?}");
    assert!(w.displacement < 1e-3 && w.circulation < 1e-3, "{w:?}");
    assert!(w.displacement > 0.0 && w.circulation > 0.0, "{w:?}");
    assert!(!ok.warnings.iter().any(|s| s.contains("relaxed wake")), "{:?}", ok.warnings);
    // capped: not converged, residuals reported, a warning on the point
    let capped = PlanePolarDef { free_wake: 1, wake_tol: 1e-12, ..review_def() };
    let o = run1(&p, &db, &capped);
    let w1 = o.diagnostics.as_ref().and_then(|d| d.wake.clone()).unwrap();
    assert!(!w1.converged && w1.iterations == 1, "{w1:?}");
    assert!(w1.displacement > w.displacement && w1.circulation > w.circulation, "{w1:?} vs {w:?}");
    assert!(o.warnings.iter().any(|s| s.contains("relaxed wake not converged")), "{:?}", o.warnings);
    // strict mode fails the point
    let pol = run_polar1(&p, &db, &PlanePolarDef { strict_wake: true, ..capped });
    assert!(pol.points.is_empty() && pol.failed.len() == 1, "{:?}", pol.failed);
    assert!(pol.failed[0].contains("strict wake"), "{}", pol.failed[0]);
    // segments and relaxation are honoured
    let fine = run1(&p, &db, &PlanePolarDef { free_wake: 100, wake_nseg: 40, wake_relax: 0.7, ..review_def() });
    assert_eq!(fine.wake[0].len(), 41);
    assert_eq!(ok.wake[0].len(), 25);
    assert!((fine.cl - ok.cl).abs() > 0.0);
    // a rigid wake has no report
    let rigid = run1(&p, &db, &PlanePolarDef { free_wake: 0, ..review_def() });
    assert!(rigid.diagnostics.is_none() && rigid.wake.is_empty());
}

#[test]
fn f007_relaxed_derivatives() {
    // frozen-wake derivatives (default) and re-relaxed ones are close, and the
    // frozen wake differs from the rigid one
    let f = naca0012();
    let (p, db) = (review_rect(&f), db_with(&f));
    let d = PlanePolarDef { derivatives: true, free_wake: 100, ..review_def() };
    let frozen = run1(&p, &db, &d).stability.unwrap();
    let relaxed = run1(&p, &db, &PlanePolarDef { relaxed_derivatives: true, ..d.clone() }).stability.unwrap();
    let rigid = run1(&p, &db, &PlanePolarDef { free_wake: 0, ..d }).stability.unwrap();
    let rel = |a: f64, b: f64| (a - b).abs() / b.abs();
    eprintln!("Zw frozen {} relaxed {} rigid {}; Mq {} {} {}", frozen.dim.zw, relaxed.dim.zw, rigid.dim.zw, frozen.dim.mq, relaxed.dim.mq, rigid.dim.mq);
    assert!(rel(frozen.dim.zw, relaxed.dim.zw) < 0.02, "Zw frozen {} relaxed {}", frozen.dim.zw, relaxed.dim.zw);
    assert!(rel(frozen.dim.mq, relaxed.dim.mq) < 0.05, "Mq frozen {} relaxed {}", frozen.dim.mq, relaxed.dim.mq);
    assert!(frozen.dim.zw != rigid.dim.zw);
}

/// The block starting at `head` in `src`, by brace matching.
fn block<'a>(src: &'a str, head: &str) -> &'a str {
    let start = src.find(head).unwrap_or_else(|| panic!("{head} not found"));
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
    panic!("unbalanced braces after {head}");
}

/// No driver calls a system's loads (or the free-wake solver) directly:
/// everything goes through `Evaluator::loads`.
#[test]
fn no_direct_sys_loads() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let read = |f: &str| std::fs::read_to_string(root.join(f)).unwrap();
    let analysis = read("analysis.rs");
    // the System dispatch itself is the one place allowed to call loads
    let system_impl = block(&analysis, "impl<'a> System<'a> {");
    let analysis = analysis.replace(system_impl, "");
    let mut hits = Vec::new();
    for (file, src) in [("analysis.rs", analysis), ("stability.rs", read("stability.rs"))] {
        for (k, line) in src.lines().enumerate() {
            let code = line.split("//").next().unwrap();
            let direct_loads = code.match_indices(".loads(").any(|(i, _)| !code[..i].ends_with("ev"));
            if direct_loads || code.contains("solve_free_wake(") || code.contains(".solve(") {
                hits.push(format!("{file}:{}: {}", k + 1, line.trim()));
            }
        }
    }
    assert!(hits.is_empty(), "direct solver calls outside the evaluator:\n{}", hits.join("\n"));
    // and the evaluator really is used by every driver
    let analysis = read("analysis.rs");
    for name in ["fn solve_point(", "fn run_control(", "fn run_stability(", "fn evaluate("] {
        let body = block(&analysis, name);
        assert!(body.contains("ev.loads(") || body.contains("evaluate(") || body.contains("Evaluator::new("), "{name} bypasses the evaluator");
    }
    let stab = read("stability.rs");
    for name in ["fn cm_at(", "fn trim(", "fn forces("] {
        assert!(block(&stab, name).contains("ev.loads"), "{name} bypasses the evaluator");
    }
}
