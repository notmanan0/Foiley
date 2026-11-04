//! F018: node-preserving 2D input (XFoil LOAD + PCOP instead of PANGEN).
use foiley_core::Progress;
use foiley_foil::analysis::analyze_point;
use foiley_foil::foil::{node_hash, NodeInput, Paneling};
use foiley_foil::{analyze, AnalysisSettings, Foil, Sweep};

fn e387() -> Foil {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/foils/e387.dat");
    Foil::load(&p).unwrap()
}

fn settings(paneling: NodeInput, viscous: bool) -> AnalysisSettings {
    AnalysisSettings { re: 2e5, viscous, paneling, keep_details: true, parallel: false, ..Default::default() }
}

#[test]
fn f018_preserve_keeps_nodes() {
    let f = e387();
    for viscous in [false, true] {
        let st = settings(NodeInput::Preserve, viscous);
        let d = analyze_point(&f, &st, 2.0, true, &Progress::none()).unwrap().expect("converged");
        assert_eq!(d.x.len(), f.x.len(), "every input point is a panel node");
        assert!(d.x.iter().zip(&f.x).all(|(a, b)| a.to_bits() == b.to_bits()), "x nodes moved");
        assert!(d.y.iter().zip(&f.y).all(|(a, b)| a.to_bits() == b.to_bits()), "y nodes moved");
        // the sweep reports the solved node count and hash
        let (pol, det) = analyze(&f, &st, &Sweep::Alpha { from: 0.0, to: 2.0, step: 1.0 }, &Progress::none()).unwrap();
        assert_eq!(pol.points.len(), 3);
        assert_eq!(pol.nodes, f.x.len());
        assert_eq!(pol.node_hash, node_hash(&f.x, &f.y));
        assert!(det.iter().all(|d| d.x == f.x && d.y == f.y));
    }
    // PANGEN repanels: other nodes, other count
    let st = settings(NodeInput::Pangen(Paneling::default()), false);
    let (pol, det) = analyze(&f, &st, &Sweep::Alpha { from: 0.0, to: 0.0, step: 1.0 }, &Progress::none()).unwrap();
    assert_eq!(pol.nodes, 160);
    assert_eq!(det[0].x.len(), 160);
    assert_ne!(pol.node_hash, node_hash(&f.x, &f.y));
    // repeated points are dropped like XFoil's PCOP; more points than the engine holds is an error
    let mut dup = f.clone();
    dup.x.insert(10, dup.x[10]);
    dup.y.insert(10, dup.y[10]);
    let (pol, _) = analyze(&dup, &settings(NodeInput::Preserve, false), &Sweep::Alpha { from: 1.0, to: 1.0, step: 1.0 }, &Progress::none()).unwrap();
    assert_eq!(pol.nodes, f.x.len());
    assert_eq!(pol.node_hash, node_hash(&f.x, &f.y));
    let n = 500;
    let big = Foil::new("big", (0..n).map(|k| 0.5 + 0.5 * (std::f64::consts::TAU * k as f64 / n as f64).cos()).collect(), (0..n).map(|k| 0.06 * (std::f64::consts::TAU * k as f64 / n as f64).sin()).collect());
    let e = analyze(&big.unwrap(), &settings(NodeInput::Preserve, false), &Sweep::Alpha { from: 1.0, to: 1.0, step: 1.0 }, &Progress::none()).unwrap_err().to_string();
    assert!(e.contains("at most") && e.contains("PANGEN"), "{e}");
}

#[test]
fn g09_node_hash_reported_and_stable() {
    let f = e387();
    let pre = settings(NodeInput::Preserve, true);
    let pan = settings(NodeInput::Pangen(Paneling::default()), true);
    let sweep = Sweep::Alpha { from: 1.0, to: 3.0, step: 1.0 };
    let run = |st: &AnalysisSettings| analyze(&f, st, &sweep, &Progress::none()).unwrap().0;
    let (a, b) = (run(&pre), run(&pre));
    assert_eq!((a.nodes, a.node_hash), (b.nodes, b.node_hash));
    // the hash of the tracked e387.dat coordinates, the same on every platform
    assert_eq!(a.node_hash, 0xaea0aee61e2f6591, "e387 preserved-node hash");
    let (c, d) = (run(&pan), run(&pan));
    assert_eq!((c.nodes, c.node_hash), (d.nodes, d.node_hash), "PANGEN nodes are repeatable");
    assert_ne!(a.node_hash, c.node_hash);
    // parallel chunks solve on the same nodes
    let par = AnalysisSettings { parallel: true, ..pan.clone() };
    let e = analyze(&f, &par, &Sweep::Alpha { from: -4.0, to: 8.0, step: 0.5 }, &Progress::none()).unwrap().0;
    assert_eq!((e.nodes, e.node_hash), (c.nodes, c.node_hash));
    // settings round-trip through JSON, and old files (PANGEN parameters only) still read
    for st in [&pre, &pan] {
        let j = serde_json::to_string(st).unwrap();
        let back: AnalysisSettings = serde_json::from_str(&j).unwrap();
        assert_eq!(&back, st);
    }
    let old: AnalysisSettings = serde_json::from_str(r#"{"paneling": {"panels": 200, "le_bunching": 1.0, "te_bunching": 0.15, "refined_ratio": 0.2}}"#).unwrap();
    assert_eq!(old.paneling, NodeInput::Pangen(Paneling { panels: 200, ..Default::default() }));
    let keep: AnalysisSettings = serde_json::from_str(r#"{"paneling": "preserve"}"#).unwrap();
    assert_eq!(keep.paneling, NodeInput::Preserve);
}
