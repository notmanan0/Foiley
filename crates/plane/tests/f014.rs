//! F014: above 2,000 unknowns the solvers factor with the cancellable
//! blocked LU, and get the same answer as the single faer call.
mod common;
use common::*;
use foiley_core::linalg::{set_blocked_lu_min, BLOCKED_LU_MIN};
use foiley_core::Progress;
use foiley_plane::analysis::*;

#[test]
fn f014_blocked_matches_faer_on_a_wing() {
    let f = naca0012();
    let db = db_with(&f);
    let mut p = review_rect(&f);
    p.wings[0].sections[0].nx = 10;
    p.wings[0].sections[0].ny = 110; // 2 x 10 x 110 = 2,200 VLM panels
    let d = PlanePolarDef { viscous: false, range: Range { from: 2.0, to: 6.0, step: 4.0 }, ..Default::default() };
    // the only test in this binary, so the process-wide switch is safe here
    set_blocked_lu_min(usize::MAX);
    let mono = run_polar1(&p, &db, &d);
    set_blocked_lu_min(BLOCKED_LU_MIN);
    let blocked = run_polar1(&p, &db, &d);
    assert_eq!(mono.unknowns, 2200);
    for (a, b) in mono.points.iter().zip(&blocked.points) {
        assert!((a.cl / b.cl - 1.0).abs() < 1e-10, "CL {} vs {}", a.cl, b.cl);
        assert!((a.cm / b.cm - 1.0).abs() < 1e-10, "Cm {} vs {}", a.cm, b.cm);
    }
    // and cancelling mid-factorization ends the polar as Cancelled
    let prog = Progress::new();
    let p2 = prog.clone();
    let (pc, dc, dbc) = (p.clone(), d.clone(), db.clone());
    let h = std::thread::spawn(move || run_polar(&pc, &dbc, &dc, &p2, &|_| {}));
    while !prog.stage_text().contains("Factorising") {
        std::thread::yield_now();
    }
    prog.cancel();
    assert!(matches!(h.join().unwrap(), Err(foiley_plane::PlaneError::Cancelled)));
}
