//! Mesh checks for imported CAD bodies: manifoldness and open edges,
//! orientation consistency, signed volume, area against a high-density
//! reference tessellation, and self-intersection (BVH over triangles).

use super::tessellate::{tessellate, CadMesh};
use super::trim::{Model, StepOptions};
use super::step_parse::StepFile;
use foiley_core::{Progress, Vec3};
use rayon::prelude::*;
use std::collections::HashMap;

/// Result of [`check`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MeshReport {
    pub triangles: usize,
    pub vertices: usize,
    /// edges used by one triangle
    pub open_edges: usize,
    /// edges used by more than two triangles
    pub nonmanifold_edges: usize,
    /// edges shared by two triangles that traverse it in the same direction
    pub misoriented_edges: usize,
    /// connected parts (by shared edges)
    pub components: usize,
    /// triangles with (near) zero area
    pub degenerate: usize,
    pub area: f64,
    /// signed volume (positive = outward normals, for closed parts)
    pub volume: f64,
    /// pairs of triangles without a common vertex that intersect
    pub self_intersections: usize,
}

impl MeshReport {
    /// Closed two-manifold.
    pub fn watertight(&self) -> bool {
        self.open_edges == 0 && self.nonmanifold_edges == 0
    }
    /// Consistent winding with outward normals.
    pub fn oriented(&self) -> bool {
        self.misoriented_edges == 0 && self.volume > 0.0
    }
    pub fn summary(&self) -> String {
        format!(
            "{} triangles, {} vertices, {} part(s); open edges {}, non-manifold {}, misoriented {}, degenerate {}, self-intersections {}; area {:.6} m2, volume {:.6} m3",
            self.triangles, self.vertices, self.components, self.open_edges, self.nonmanifold_edges, self.misoriented_edges, self.degenerate, self.self_intersections, self.area, self.volume
        )
    }
}

/// Check an indexed mesh (see [`MeshReport`]).
pub fn check(vertices: &[Vec3], tris: &[[u32; 3]]) -> MeshReport {
    let mut r = MeshReport { triangles: tris.len(), vertices: vertices.len(), ..Default::default() };
    let p = |i: u32| vertices[i as usize];
    let (mut lo, mut hi) = (Vec3::new(f64::MAX, f64::MAX, f64::MAX), Vec3::new(f64::MIN, f64::MIN, f64::MIN));
    for v in vertices {
        lo = Vec3::new(lo.x.min(v.x), lo.y.min(v.y), lo.z.min(v.z));
        hi = Vec3::new(hi.x.max(v.x), hi.y.max(v.y), hi.z.max(v.z));
    }
    let diag = (hi - lo).norm().max(1e-300);
    // (forward uses, backward uses, triangles)
    let mut edges: HashMap<(u32, u32), (u32, u32, Vec<u32>)> = HashMap::new();
    for (i, t) in tris.iter().enumerate() {
        let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
        let n = (b - a).cross(c - a);
        r.area += 0.5 * n.norm();
        r.volume += a.dot(b.cross(c)) / 6.0;
        if n.norm() <= 1e-14 * diag * diag {
            r.degenerate += 1;
        }
        for k in 0..3 {
            let (x, y) = (t[k], t[(k + 1) % 3]);
            let e = edges.entry((x.min(y), x.max(y))).or_insert((0, 0, Vec::new()));
            if x < y {
                e.0 += 1;
            } else {
                e.1 += 1;
            }
            e.2.push(i as u32);
        }
    }
    for (f, b, _) in edges.values() {
        match f + b {
            1 => r.open_edges += 1,
            2 if *f == 2 || *b == 2 => r.misoriented_edges += 1,
            2 => {}
            _ => r.nonmanifold_edges += 1,
        }
    }
    // components (union-find over shared edges)
    let mut parent: Vec<u32> = (0..tris.len() as u32).collect();
    fn find(p: &mut [u32], mut i: u32) -> u32 {
        while p[i as usize] != i {
            p[i as usize] = p[p[i as usize] as usize];
            i = p[i as usize];
        }
        i
    }
    for (_, _, ts) in edges.values() {
        for w in ts.windows(2) {
            let (a, b) = (find(&mut parent, w[0]), find(&mut parent, w[1]));
            if a != b {
                parent[a.max(b) as usize] = a.min(b);
            }
        }
    }
    r.components = (0..tris.len() as u32).filter(|&i| find(&mut parent, i) == i).count();
    r.self_intersections = self_intersections(vertices, tris);
    r
}

/// Check a CAD import's mesh.
pub fn check_mesh(m: &CadMesh) -> MeshReport {
    check(&m.vertices, &m.triangles)
}

/// Area of the model's tessellation and of a reference tessellation with
/// tolerances `factor` times tighter: (area, reference area).
pub fn area_vs_reference(text: &str, opts: &StepOptions, factor: f64) -> Result<(f64, f64), String> {
    let file = StepFile::parse(text)?;
    let m = Model::from_file(&file, opts)?;
    let (mesh, _) = tessellate(&m, opts.strict, &Progress::none())?;
    let fine = StepOptions { chord_tol: m.tol.chord / factor, max_edge: m.tol.max_edge / factor.sqrt(), ..opts.clone() };
    let mr = Model::from_file(&file, &fine)?;
    let (reference, _) = tessellate(&mr, opts.strict, &Progress::none())?;
    Ok((mesh.area(), reference.area()))
}

// ---------------------------------------------------------------- BVH

#[derive(Clone, Copy)]
struct Aabb {
    lo: Vec3,
    hi: Vec3,
}

impl Aabb {
    fn of(t: &[Vec3; 3]) -> Aabb {
        let lo = Vec3::new(t[0].x.min(t[1].x).min(t[2].x), t[0].y.min(t[1].y).min(t[2].y), t[0].z.min(t[1].z).min(t[2].z));
        let hi = Vec3::new(t[0].x.max(t[1].x).max(t[2].x), t[0].y.max(t[1].y).max(t[2].y), t[0].z.max(t[1].z).max(t[2].z));
        Aabb { lo, hi }
    }
    fn join(a: Aabb, b: Aabb) -> Aabb {
        Aabb { lo: Vec3::new(a.lo.x.min(b.lo.x), a.lo.y.min(b.lo.y), a.lo.z.min(b.lo.z)), hi: Vec3::new(a.hi.x.max(b.hi.x), a.hi.y.max(b.hi.y), a.hi.z.max(b.hi.z)) }
    }
    fn overlaps(&self, o: &Aabb) -> bool {
        self.lo.x <= o.hi.x && o.lo.x <= self.hi.x && self.lo.y <= o.hi.y && o.lo.y <= self.hi.y && self.lo.z <= o.hi.z && o.lo.z <= self.hi.z
    }
}

enum Node {
    Leaf(Aabb, Vec<u32>),
    Inner(Aabb, Box<Node>, Box<Node>),
}

impl Node {
    fn bbox(&self) -> &Aabb {
        match self {
            Node::Leaf(b, _) | Node::Inner(b, _, _) => b,
        }
    }
    fn build(ids: Vec<u32>, boxes: &[Aabb], cent: &[Vec3]) -> Node {
        let b = ids.iter().skip(1).fold(boxes[ids[0] as usize], |a, &i| Aabb::join(a, boxes[i as usize]));
        if ids.len() <= 8 {
            return Node::Leaf(b, ids);
        }
        let d = b.hi - b.lo;
        let axis = if d.x >= d.y && d.x >= d.z { 0 } else if d.y >= d.z { 1 } else { 2 };
        let key = |i: &u32| {
            let c = cent[*i as usize];
            [c.x, c.y, c.z][axis]
        };
        let mut ids = ids;
        let mid = ids.len() / 2;
        ids.select_nth_unstable_by(mid, |a, b| key(a).total_cmp(&key(b)));
        let right = ids.split_off(mid);
        Node::Inner(b, Box::new(Node::build(ids, boxes, cent)), Box::new(Node::build(right, boxes, cent)))
    }
    fn query(&self, b: &Aabb, out: &mut Vec<u32>) {
        if !self.bbox().overlaps(b) {
            return;
        }
        match self {
            Node::Leaf(_, ids) => out.extend_from_slice(ids),
            Node::Inner(_, l, r) => {
                l.query(b, out);
                r.query(b, out);
            }
        }
    }
}

/// Does the segment p-q cross the interior of triangle t (strictly)?
fn segment_hits(p: Vec3, q: Vec3, t: &[Vec3; 3], eps: f64) -> bool {
    let d = q - p;
    let (e1, e2) = (t[1] - t[0], t[2] - t[0]);
    let h = d.cross(e2);
    let a = e1.dot(h);
    if a.abs() <= eps * d.norm() * e1.norm() * e2.norm() {
        return false; // parallel or coplanar: not counted
    }
    let f = 1.0 / a;
    let s = p - t[0];
    let u = f * s.dot(h);
    if u <= eps || u >= 1.0 - eps {
        return false;
    }
    let qv = s.cross(e1);
    let v = f * d.dot(qv);
    if v <= eps || u + v >= 1.0 - eps {
        return false;
    }
    let w = f * e2.dot(qv);
    w > eps && w < 1.0 - eps
}

fn tri_tri(a: &[Vec3; 3], b: &[Vec3; 3]) -> bool {
    let eps = 1e-9;
    (0..3).any(|k| segment_hits(a[k], a[(k + 1) % 3], b, eps)) || (0..3).any(|k| segment_hits(b[k], b[(k + 1) % 3], a, eps))
}

/// Number of intersecting triangle pairs that share no vertex.
pub fn self_intersections(vertices: &[Vec3], tris: &[[u32; 3]]) -> usize {
    if tris.len() < 2 {
        return 0;
    }
    let tv: Vec<[Vec3; 3]> = tris.iter().map(|t| [vertices[t[0] as usize], vertices[t[1] as usize], vertices[t[2] as usize]]).collect();
    let boxes: Vec<Aabb> = tv.iter().map(Aabb::of).collect();
    let cent: Vec<Vec3> = tv.iter().map(|t| (t[0] + t[1] + t[2]) / 3.0).collect();
    let root = Node::build((0..tris.len() as u32).collect(), &boxes, &cent);
    (0..tris.len())
        .into_par_iter()
        .map_init(Vec::new, |cand, i| {
            cand.clear();
            root.query(&boxes[i], cand);
            cand.iter()
                .filter(|&&j| (j as usize) > i)
                .filter(|&&j| {
                    let (a, b) = (&tris[i], &tris[j as usize]);
                    !a.iter().any(|x| b.contains(x)) && tri_tri(&tv[i], &tv[j as usize])
                })
                .count()
        })
        .sum()
}
