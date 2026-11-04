//! F004: multi-body meshing.
mod common;
use common::*;
use foiley_core::Vec3;
use foiley_plane::analysis::{Method, PlaneOpp, PlanePolarDef, Range};
use foiley_plane::geom::{Fuse, FuseFrame, Plane, Wing, WingKind, WingSection};
use foiley_plane::mesh::{build, Mesh, MeshOptions, PanelKind};

fn thick() -> MeshOptions {
    MeshOptions { thick: true, include_fuse: true, ..Default::default() }
}

/// Inviscid rigid-wake panel method at alpha = 4 deg, bodies included.
fn panel_def() -> PlanePolarDef {
    let mut d = review_def();
    d.method = Method::Panel;
    d.free_wake = 0;
    d.include_fuse = true;
    d.alpha = 4.0;
    d.range = Range { from: 4.0, to: 4.0, step: 1.0 };
    d
}

fn body_panels(m: &Mesh, body: usize) -> usize {
    m.panels.iter().filter(|p| p.kind == PanelKind::Fuse && p.body_id == Some(body)).count()
}

fn pod_at(y: f64, z: f64) -> Fuse {
    let mut f = Fuse::default_pod();
    f.position.y = y;
    f.position.z = z;
    f
}

fn lateral(o: &PlaneOpp) -> (f64, f64, f64) {
    (o.cy.unwrap(), o.cr.unwrap(), o.cn.unwrap())
}

#[test]
fn f004_two_pods_double_panels_mirror_loads() {
    let f = naca0012();
    let db = db_with(&f);
    let mut pods = review_rect(&f);
    pods.fuses.push(Fuse::default_pod());
    let opts = thick();
    let count = |p: &Plane| build(p, &db, &opts).unwrap().panels.iter().filter(|p| p.kind == PanelKind::Fuse).count();
    let one = count(&pods);
    let mut second = Fuse::default_pod();
    second.position.y = 3.0;
    pods.fuses.push(second);
    assert_eq!(count(&pods), 2 * one);
    let m = build(&pods, &db, &opts).unwrap();
    assert_eq!(body_panels(&m, 0), one);
    assert_eq!(body_panels(&m, 1), one);

    // a pod under the right wing and its mirror image under the left wing give mirrored loads
    let d = panel_def();
    let with = |pods: Vec<Fuse>| {
        let mut p = review_rect(&f);
        p.fuses = pods;
        run1(&p, &db, &d)
    };
    let (r, l) = (with(vec![pod_at(2.0, -0.25)]), with(vec![pod_at(-2.0, -0.25)]));
    let ((ry, rr, rn), (ly, lr, ln)) = (lateral(&r), lateral(&l));
    assert!(rr.abs() + rn.abs() > 1e-6, "a single offset pod must load the plane laterally: Cl {rr}, Cn {rn}");
    for (a, b, what) in [(ry, ly, "CY"), (rr, lr, "Cl"), (rn, ln, "Cn")] {
        assert!((a + b).abs() < 1e-9, "{what} not mirrored: {a} vs {b}");
    }
    for (a, b, what) in [(r.cl, l.cl, "CL"), (r.cdi, l.cdi, "CDi"), (r.cm, l.cm, "Cm")] {
        assert!((a - b).abs() < 1e-9, "{what} differs between mirror images: {a} vs {b}");
    }
    // both pods together: no lateral load
    let both = with(vec![pod_at(2.0, -0.25), pod_at(-2.0, -0.25)]);
    let (by, br, bn) = lateral(&both);
    assert!(by.abs() < 1e-9 && br.abs() < 1e-9 && bn.abs() < 1e-9, "symmetric pods: CY {by}, Cl {br}, Cn {bn}");
}

/// FNV-1a over the mesh data that existed before multi-body meshing. With `quantum` = None
/// every float goes in bit for bit; with Some(dq) floats are first rounded to multiples of dq,
/// so ulp-level differences between platform libms (glibc vs MSVC) do not change the hash.
fn mesh_hash(m: &Mesh, quantum: Option<f64>) -> u64 {
    let q = |x: f64| match quantum {
        None => x.to_bits(),
        Some(dq) => (x / dq).round() as i64 as u64,
    };
    let mut h: u64 = 0xcbf29ce484222325;
    let mut eat = |v: u64| {
        for b in v.to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
    };
    for p in &m.panels {
        for n in &p.nodes {
            eat(q(n.x));
            eat(q(n.y));
            eat(q(n.z));
        }
        eat(q(p.normal.x));
        eat(q(p.normal.y));
        eat(q(p.normal.z));
        eat(q(p.area));
        eat(p.wing as u64);
        eat(p.strip.map(|s| s as u64 + 1).unwrap_or(0));
        eat(p.kind as u64);
    }
    for s in &m.strips {
        eat(q(s.chord));
        eat(q(s.c4.y));
        eat(s.panels.len() as u64);
    }
    for &(a, b) in &m.carry {
        eat(a as u64);
        eat(b as u64);
    }
    // neighbour lists are built from a HashMap: compare them as sets
    for nb in &m.neighbours {
        let mut nb = nb.clone();
        nb.sort();
        eat(nb.len() as u64);
        for &j in &nb {
            eat(j as u64);
        }
    }
    eat(m.n_wing_panels as u64);
    h
}

/// The golden values were recorded with the single-body mesher of 10062ce in a Windows release
/// build. The bit-exact hashes only hold there: other libms (glibc on CI) and the debug profile
/// (no compile-time folding of libm calls) differ in the last ulp. Every build checks the hash
/// of the mesh rounded to 1e-9 m, recorded from that same bit-identical mesh.
const EXACT: bool = cfg!(all(windows, not(debug_assertions)));
const DQ: Option<f64> = Some(1e-9);

fn check_hash(m: &Mesh, exact: u64, coarse: u64) {
    if EXACT {
        assert_eq!(mesh_hash(m, None), exact, "bit-exact mesh hash");
    }
    assert_eq!(mesh_hash(m, DQ), coarse, "mesh hash at 1e-9 m");
}

#[test]
fn f004_single_body_bit_identical() {
    let f = naca0012();
    let db = db_with(&f);
    let mut p = Plane::default_plane(&f.name, &f.name);
    p.fuses.push(Fuse::default_pod());
    let quads = build(&p, &db, &thick()).unwrap();
    check_hash(&quads, 0x4e65f51a9895c45f, 0x231026c52ab3102c);
    assert_eq!(quads.carry, vec![(18, 19), (44, 45)]);
    assert!(quads.panels.iter().all(|q| (q.kind == PanelKind::Fuse) == (q.body_id == Some(0))));
    let tris = build(&p, &db, &MeshOptions { triangulate: true, ..thick() }).unwrap();
    check_hash(&tris, 0x16dd8f8fea190a43, 0x8970540f7cd3e87b);
    assert!(tris.panels.iter().all(|q| (q.kind == PanelKind::Fuse) == (q.body_id == Some(0))));
    let thin = build(&p, &db, &MeshOptions { include_fuse: true, fuse_in_thin: true, ..Default::default() }).unwrap();
    check_hash(&thin, 0xa37001f712cebe9e, 0x4b9022463490c5de);
    // the solve itself: 10062ce gave CL 0x3fd841166e9946e8, CDi 0x3f709591dfe0640e and
    // Cm 0x3fc8254bc689ecea (+-1 ulp from run to run: the neighbour lists came out of a
    // HashMap in random order; they are sorted now)
    let o = run1(&p, &db, &panel_def());
    let tol = if EXACT { 1e-13 } else { 1e-9 };
    for (v, bits) in [(o.cl, 0x3fd841166e9946e8u64), (o.cdi, 0x3f709591dfe0640e), (o.cm, 0x3fc8254bc689ecea)] {
        let g = f64::from_bits(bits);
        assert!((v - g).abs() <= tol * g.abs(), "{v} vs 10062ce {g}");
    }
    let again = run1(&p, &db, &panel_def());
    assert_eq!((o.cl.to_bits(), o.cdi.to_bits(), o.cm.to_bits()), (again.cl.to_bits(), again.cdi.to_bits(), again.cm.to_bits()), "repeatable to the bit");
}

fn boom(name: &str, y: f64) -> Fuse {
    let fr = |x: f64, r: f64| FuseFrame { x, z: 0.0, half_width: r, half_height: r };
    Fuse {
        name: name.into(),
        position: Vec3::new(-0.4, y, -0.2),
        frames: vec![fr(0.0, 0.005), fr(0.15, 0.05), fr(0.6, 0.06), fr(2.4, 0.04), fr(3.4, 0.005)],
        exponent: 2.0,
        nx: 24,
        nh: 6,
        mass: 0.2,
        nurbs: None,
        triangles: vec![],
    }
}

/// Central pod (trims the wing root) plus two tail booms below a raised tailplane
/// (kept out of the main-wing wake sheet).
fn twin_boom(f: &foiley_foil::Foil) -> Plane {
    let mut p = review_rect(f);
    p.fuses = vec![Fuse::default_pod(), boom("left boom", -1.5), boom("right boom", 1.5)];
    let mut root = WingSection::new(0.0, 0.5, 0.0, &f.name);
    root.nx = 3;
    root.ny = 6;
    let mut tail = Wing::new("tailplane", WingKind::Elevator, vec![root, WingSection::new(1.2, 0.5, 0.0, &f.name)]);
    tail.position = Vec3::new(2.5, 0.0, 0.3);
    p.wings.push(tail);
    p
}

#[test]
fn g23_twin_boom_fixture() {
    let f = naca0012();
    let db = db_with(&f);
    let p = twin_boom(&f);
    let m = build(&p, &db, &thick()).unwrap();
    let counts: Vec<usize> = (0..3).map(|b| body_panels(&m, b)).collect();
    assert!(counts.iter().all(|&c| c > 0), "every body is meshed: {counts:?}");
    assert_eq!(counts[1], counts[2], "the booms are mirror images");
    assert_eq!(m.panels.iter().filter(|q| q.kind == PanelKind::Fuse).count(), counts.iter().sum::<usize>());
    // the main wing root is trimmed on the pod (one carry-through pair); the
    // tailplane root is clear of every body
    assert_eq!(m.carry.len(), 1, "{:?}", m.carry);
    let (l, r) = m.carry[0];
    assert_eq!((m.strips[l].wing, m.strips[r].wing), (0, 0));
    // inboard leading-edge node of the right root strip sits on the pod side, not at y = 0
    let lead = m.strips[r].panels.iter().map(|&i| &m.panels[i]).find(|q| q.leading && q.kind == PanelKind::Top).unwrap();
    assert!(lead.nodes[0].y > 0.02 && lead.nodes[0].y < 0.06, "root trimmed at the pod side: {}", lead.nodes[0].y);

    let o = run1(&p, &db, &panel_def());
    let (cy, cr, cn) = lateral(&o);
    assert!(cy.abs() < 1e-9 && cr.abs() < 1e-9 && cn.abs() < 1e-9, "symmetric twin boom: CY {cy}, Cl {cr}, Cn {cn}");
    assert!(o.cl > 0.2 && o.cl.is_finite(), "CL {}", o.cl);

    // a boom moved into the pod is a clear error naming both bodies
    let mut clash = p.clone();
    clash.fuses[2].position = Vec3::new(-0.3, 0.06, -0.05);
    let e = build(&clash, &db, &thick()).unwrap_err().to_string();
    assert!(e.contains("interpenetrate") && e.contains("'Fuselage' (body 1)") && e.contains("'right boom' (body 3)"), "{e}");
    // a body entirely inside another one is an error too
    let mut nested = p.clone();
    let mut inner = Fuse::default_pod();
    inner.name = "inner".into();
    for fr in &mut inner.frames {
        fr.half_width *= 0.5;
        fr.half_height *= 0.5;
    }
    nested.fuses.push(inner);
    let e = build(&nested, &db, &thick()).unwrap_err().to_string();
    assert!(e.contains("'inner' (body 4) lies inside body 'Fuselage' (body 1)"), "{e}");
}
