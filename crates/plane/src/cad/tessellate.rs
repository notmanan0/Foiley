//! Watertight tessellation of STEP faces.
//!
//! * Every edge is sampled once (in [`super::trim`]) and its samples are
//!   global vertices, so adjacent faces share their boundary vertices and
//!   the mesh is watertight by construction.
//! * Each face's trimmed parameter domain (outer loop, holes, seams, poles)
//!   is triangulated by constrained Delaunay ([`spade`]) in a parameter
//!   space scaled to the surface's mean metric, triangles outside the
//!   domain are dropped (parity flood fill across the constraints), then
//!   refined by inserting interior points until every triangle meets the
//!   chordal tolerance at its centroid and edge midpoints and the longest
//!   side is below the length limit.  Boundary (constraint) edges are never
//!   split, so the shared vertices stay shared.
//! * Faces are processed in parallel with Rayon.

use super::trim::{Domain, Model, StepOptions};
use super::{Surface, Tri};
use foiley_core::{Progress, Vec3};
use rayon::prelude::*;
use spade::handles::FixedVertexHandle;
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
use std::collections::{HashMap, VecDeque};

/// An indexed triangle mesh from a STEP import.
#[derive(Clone, Debug, Default)]
pub struct CadMesh {
    pub vertices: Vec<Vec3>,
    pub triangles: Vec<[u32; 3]>,
    /// source face (index into `Model::faces`) of each triangle
    pub face: Vec<u32>,
    /// shell (body) of each triangle
    pub shell: Vec<u32>,
}

impl CadMesh {
    pub fn to_tris(&self) -> Vec<Tri> {
        self.triangles.iter().map(|t| [self.vertices[t[0] as usize], self.vertices[t[1] as usize], self.vertices[t[2] as usize]]).collect()
    }
    pub fn area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let (a, b, c) = (self.vertices[t[0] as usize], self.vertices[t[1] as usize], self.vertices[t[2] as usize]);
                0.5 * (b - a).cross(c - a).norm()
            })
            .sum()
    }
    /// Signed enclosed volume (positive for outward-oriented closed shells).
    pub fn volume(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let (a, b, c) = (self.vertices[t[0] as usize], self.vertices[t[1] as usize], self.vertices[t[2] as usize]);
                a.dot(b.cross(c))
            })
            .sum::<f64>()
            / 6.0
    }
}

/// A tessellated face: parameter coordinates of its vertices and its
/// triangles (vertex numbering local to the face, `ids` maps to the model:
/// a global vertex, a `Domain::extra` vertex, or `u32::MAX` for a point
/// inserted inside the face).
#[derive(Clone, Debug, Default)]
pub struct FaceMesh {
    pub uv: Vec<[f64; 2]>,
    pub ids: Vec<u32>,
    pub pos: Vec<Vec3>,
    /// counter-clockwise in (u, v); flipped to the face normal on assembly
    pub tris: Vec<[u32; 3]>,
    /// local vertex pairs that are domain-boundary (constraint) edges
    pub boundary: Vec<[u32; 2]>,
    /// parameters that are free at the vertex (pole / apex)
    pub sing: Vec<[bool; 2]>,
    pub domain: Domain,
}

/// Result of [`import_step`].
pub struct StepImport {
    pub mesh: CadMesh,
    pub notes: Vec<String>,
    pub model: Model,
}

/// Parse, tessellate and assemble a STEP file.
pub fn import_step(text: &str, opts: &StepOptions, p: &Progress) -> Result<StepImport, String> {
    p.check().map_err(|e| e.to_string())?;
    let model = Model::from_step(text, opts)?;
    p.check().map_err(|e| e.to_string())?;
    let (mesh, notes) = tessellate(&model, opts.strict, p)?;
    Ok(StepImport { mesh, notes, model })
}

/// Tessellate every face of a model and assemble one indexed mesh.
pub fn tessellate(model: &Model, strict: bool, p: &Progress) -> Result<(CadMesh, Vec<String>), String> {
    let n = model.faces.len();
    let results: Vec<Result<FaceMesh, String>> = (0..n)
        .into_par_iter()
        .map(|fi| {
            p.check().map_err(|e| e.to_string())?;
            let r = face_mesh(model, fi);
            p.set((fi + 1) as f64 / n.max(1) as f64);
            r
        })
        .collect();
    p.check().map_err(|e| e.to_string())?;
    let mut notes = model.notes.clone();
    let nv = model.vertices.len();
    let mut map: Vec<u32> = vec![u32::MAX; nv];
    let mut mesh = CadMesh::default();
    let mut failed = 0;
    let mut seen: std::collections::HashSet<[u32; 3]> = std::collections::HashSet::new();
    let mut duplicates = 0usize;
    for (fi, r) in results.into_iter().enumerate() {
        let fm = match r {
            Ok(fm) => fm,
            Err(e) => {
                if strict {
                    return Err(e);
                }
                failed += 1;
                if failed <= 20 {
                    notes.push(format!("face #{} skipped: {e}", model.faces[fi].step_id));
                }
                continue;
            }
        };
        let mut extra_map: HashMap<u32, u32> = HashMap::new();
        let local: Vec<u32> = (0..fm.ids.len())
            .map(|k| {
                let id = fm.ids[k];
                let push = |m: &mut CadMesh, p: Vec3| {
                    m.vertices.push(p);
                    (m.vertices.len() - 1) as u32
                };
                if id == u32::MAX {
                    push(&mut mesh, fm.pos[k])
                } else if (id as usize) < nv {
                    if map[id as usize] == u32::MAX {
                        map[id as usize] = push(&mut mesh, model.vertices[id as usize]);
                    }
                    map[id as usize]
                } else {
                    *extra_map.entry(id).or_insert_with(|| push(&mut mesh, fm.pos[k]))
                }
            })
            .collect();
        let flip = !model.faces[fi].same_sense;
        for t in &fm.tris {
            let (a, b, c) = (local[t[0] as usize], local[t[1] as usize], local[t[2] as usize]);
            if a == b || b == c || a == c {
                continue; // collapsed at a pole or seam
            }
            // two faces that collapse to the same corner (untrimmed patches
            // meeting at a nose) can both emit the triangle on their shared
            // boundary; the second copy would make its edges non-manifold
            let mut key = [a, b, c];
            key.sort_unstable();
            if !seen.insert(key) {
                duplicates += 1;
                continue;
            }
            mesh.triangles.push(if flip { [a, c, b] } else { [a, b, c] });
            mesh.face.push(fi as u32);
            mesh.shell.push(model.faces[fi].shell as u32);
        }
    }
    if failed > 20 {
        notes.push(format!("... {} more faces skipped", failed - 20));
    }
    if duplicates > 0 {
        notes.push(format!("{duplicates} duplicate triangle(s) where faces meet at a collapsed corner dropped"));
    }
    if mesh.triangles.is_empty() {
        return Err("no face could be tessellated".into());
    }
    if model.bare_surfaces {
        // patches carry their own parametrisation: make their windings agree
        let flipped = super::orient_indexed(&mut mesh.triangles, &mesh.vertices);
        if flipped > 0 {
            notes.push(format!("{flipped} triangle(s) reoriented for a consistent outward winding"));
        }
    }
    notes.insert(0, format!("{} face(s) in {} body(ies), unit scale {} m, {} triangles, {} vertices", n - failed, model.n_shells, model.length_scale, mesh.triangles.len(), mesh.vertices.len()));
    Ok((mesh, notes))
}

/// Mean |dS/du| and |dS/dv| over the domain's bounding box.
fn metric(s: &Surface, lo: [f64; 2], hi: [f64; 2]) -> [f64; 2] {
    let mut acc = [0.0f64; 2];
    let mut n = 0.0;
    for i in 0..5 {
        for j in 0..5 {
            let u = lo[0] + (hi[0] - lo[0]) * (i as f64 + 0.5) / 5.0;
            let v = lo[1] + (hi[1] - lo[1]) * (j as f64 + 0.5) / 5.0;
            let (_, su, sv) = s.derivs(u, v);
            acc[0] += su.norm();
            acc[1] += sv.norm();
            n += 1.0;
        }
    }
    let m = [acc[0] / n, acc[1] / n];
    let big = m[0].max(m[1]).max(1e-300);
    [m[0].max(1e-6 * big), m[1].max(1e-6 * big)]
}

/// Tessellate one face (see the module notes).
pub fn face_mesh(model: &Model, fi: usize) -> Result<FaceMesh, String> {
    let face = &model.faces[fi];
    let s = &face.surface;
    let dom = model.face_domain(fi)?;
    let tol = model.tol;
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for l in &dom.loops {
        let (a, b) = l.bbox();
        for k in 0..2 {
            lo[k] = lo[k].min(a[k]);
            hi[k] = hi[k].max(b[k]);
        }
    }
    let sc = metric(s, lo, hi);
    let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> = ConstrainedDelaunayTriangulation::new();
    let mut fm = FaceMesh::default();
    // spade vertex index -> local vertex
    let mut of_handle: Vec<u32> = Vec::new();
    let add_vertex = |cdt: &mut ConstrainedDelaunayTriangulation<Point2<f64>>, fm: &mut FaceMesh, of_handle: &mut Vec<u32>, uv: [f64; 2], id: u32, pos: Vec3| -> Result<FixedVertexHandle, String> {
        // spade refuses non-zero coordinates below ~1e-43
        let snap = |x: f64| if x.abs() < 1e-30 { 0.0 } else { x };
        let h = cdt.insert(Point2::new(snap((uv[0] - lo[0]) * sc[0]), snap((uv[1] - lo[1]) * sc[1]))).map_err(|e| format!("face #{}: triangulation: {e:?}", face.step_id))?;
        if h.index() >= of_handle.len() {
            of_handle.resize(h.index() + 1, u32::MAX);
        }
        if of_handle[h.index()] == u32::MAX {
            of_handle[h.index()] = fm.uv.len() as u32;
            fm.sing.push(if id == u32::MAX { [false; 2] } else { s.singular(uv[0], uv[1]) });
            fm.uv.push(uv);
            fm.ids.push(id);
            fm.pos.push(pos);
        }
        Ok(h)
    };
    let mut constraints = Vec::new();
    for l in &dom.loops {
        let hs: Vec<FixedVertexHandle> = l.uv.iter().zip(&l.ids).map(|(uv, &id)| add_vertex(&mut cdt, &mut fm, &mut of_handle, *uv, id, model.vertex(&dom, id))).collect::<Result<_, _>>()?;
        for k in 0..hs.len() {
            let (a, b) = (hs[k], hs[(k + 1) % hs.len()]);
            if a != b {
                constraints.push((a, b));
            }
        }
    }
    for &(a, b) in &constraints {
        if !cdt.can_add_constraint(a, b) {
            return Err(format!("face #{}: its boundary crosses itself in parameter space", face.step_id));
        }
        cdt.add_constraint(a, b);
    }
    for &(a, b) in &constraints {
        fm.boundary.push([of_handle[a.index()], of_handle[b.index()]]);
    }
    // refinement: insert centroids / edge midpoints of inside triangles that
    // miss the chordal or length tolerance
    let max_pts = 2_000_000;
    for _pass in 0..40 {
        let inside = classify(&cdt);
        let mut cand: Vec<[f64; 2]> = Vec::new();
        for f in cdt.inner_faces() {
            if !inside[f.fix().index()] {
                continue;
            }
            // never refine a sliver (nearly collinear boundary points): its
            // centroid or midpoints would land on the boundary
            let ps = f.positions();
            let a2 = ((ps[1].x - ps[0].x) * (ps[2].y - ps[0].y) - (ps[2].x - ps[0].x) * (ps[1].y - ps[0].y)).abs();
            let l2 = (0..3).map(|k| (ps[(k + 1) % 3].x - ps[k].x).powi(2) + (ps[(k + 1) % 3].y - ps[k].y).powi(2)).fold(0.0, f64::max);
            if a2 <= 1e-9 * l2 {
                continue;
            }
            let vs = f.vertices();
            let l: [u32; 3] = [of_handle[vs[0].fix().index()], of_handle[vs[1].fix().index()], of_handle[vs[2].fix().index()]];
            // a triangle with two vertices at the same pole is dropped on
            // assembly and never refined; one spanning the two copies of a
            // seam (same ids, no pole) must still be refined until it no
            // longer does
            let ids = l.map(|i| fm.ids[i as usize]);
            let sg = l.map(|i| fm.sing.get(i as usize).copied().unwrap_or([false; 2]));
            if (0..3).any(|k| {
                let j = (k + 1) % 3;
                ids[k] != u32::MAX && ids[k] == ids[j] && (sg[k] != [false; 2] || sg[j] != [false; 2])
            }) {
                continue;
            }
            let dv = deviations(s, &fm, l);
            let q = l.map(|i| fm.uv[i as usize]);
            let edges = f.adjacent_edges();
            // spade: edge k runs from vertex k to vertex k + 1
            let mut worst: Option<(f64, [f64; 2])> = None;
            let mut consider = |score: f64, uv: [f64; 2]| {
                if score > 1.0 && worst.map(|w| score > w.0).unwrap_or(true) {
                    worst = Some((score, uv));
                }
            };
            for (k, e) in edges.iter().enumerate() {
                if e.as_undirected().is_constraint_edge() {
                    continue;
                }
                let (a, b) = (k, (k + 1) % 3);
                let mid = [0.5 * (q[a][0] + q[b][0]), 0.5 * (q[a][1] + q[b][1])];
                consider((dv.edges[k] / tol.chord).max(dv.lengths[k] / tol.max_edge), mid);
            }
            consider(dv.centroid / tol.chord, [(q[0][0] + q[1][0] + q[2][0]) / 3.0, (q[0][1] + q[1][1] + q[2][1]) / 3.0]);
            if let Some((_, uv)) = worst {
                cand.push(uv);
            }
        }
        if cand.is_empty() {
            break;
        }
        if fm.uv.len() + cand.len() > max_pts {
            return Err(format!("face #{}: refinement needs more than {max_pts} points (tolerance too tight?)", face.step_id));
        }
        for uv in cand {
            let pos = s.eval(uv[0], uv[1]);
            add_vertex(&mut cdt, &mut fm, &mut of_handle, uv, u32::MAX, pos)?;
        }
    }
    let inside = classify(&cdt);
    for f in cdt.inner_faces() {
        if !inside[f.fix().index()] {
            continue;
        }
        let vs = f.vertices();
        fm.tris.push([of_handle[vs[0].fix().index()], of_handle[vs[1].fix().index()], of_handle[vs[2].fix().index()]]);
    }
    fm.domain = dom;
    Ok(fm)
}

/// Chordal deviations of a face triangle (local vertices `t`).
#[derive(Clone, Copy, Debug)]
pub struct Deviations {
    /// distance from the surface at the parametric centroid to the
    /// triangle's plane
    pub centroid: f64,
    /// distance from the surface at each edge's parametric midpoint to the
    /// edge's line (edge k joins vertex k and k + 1)
    pub edges: [f64; 3],
    /// 3D side lengths
    pub lengths: [f64; 3],
}

/// Parameter coordinates of a triangle's vertices for evaluating its
/// interior: a vertex at a pole takes its free parameter from the others.
pub fn effective_uv(fm: &FaceMesh, t: [u32; 3]) -> [[f64; 2]; 3] {
    let mut q = t.map(|i| fm.uv[i as usize]);
    let sing = t.map(|i| fm.sing.get(i as usize).copied().unwrap_or([false; 2]));
    for d in 0..2 {
        for k in 0..3 {
            if !sing[k][d] {
                continue;
            }
            let others: Vec<f64> = (0..3).filter(|&j| j != k && !sing[j][d]).map(|j| q[j][d]).collect();
            if !others.is_empty() {
                q[k][d] = others.iter().sum::<f64>() / others.len() as f64;
            }
        }
    }
    q
}

/// See [`Deviations`].
pub fn deviations(s: &Surface, fm: &FaceMesh, t: [u32; 3]) -> Deviations {
    let p = t.map(|i| fm.pos[i as usize]);
    let q = effective_uv(fm, t);
    let n = (p[1] - p[0]).cross(p[2] - p[0]);
    let nn = n.norm();
    let c = [(q[0][0] + q[1][0] + q[2][0]) / 3.0, (q[0][1] + q[1][1] + q[2][1]) / 3.0];
    let sc = s.eval(c[0], c[1]);
    let centroid = if nn > 0.0 { ((sc - p[0]).dot(n) / nn).abs() } else { (sc - (p[0] + p[1] + p[2]) / 3.0).norm() };
    let sing = t.map(|i| fm.sing.get(i as usize).copied().unwrap_or([false; 2]));
    let mut edges = [0.0; 3];
    let mut lengths = [0.0; 3];
    for k in 0..3 {
        let (a, b) = (k, (k + 1) % 3);
        lengths[k] = (p[b] - p[a]).norm();
        let (mut qa, mut qb) = (fm.uv[t[a] as usize], fm.uv[t[b] as usize]);
        for d in 0..2 {
            if sing[a][d] && !sing[b][d] {
                qa[d] = qb[d];
            } else if sing[b][d] && !sing[a][d] {
                qb[d] = qa[d];
            }
        }
        let m = s.eval(0.5 * (qa[0] + qb[0]), 0.5 * (qa[1] + qb[1]));
        let ab = p[b] - p[a];
        let w = m - p[a];
        let l2 = ab.norm2();
        edges[k] = if l2 > 0.0 { (w - ab * (w.dot(ab) / l2)).norm() } else { w.norm() };
    }
    Deviations { centroid, edges, lengths }
}

/// Inside flags of the faces (by fixed face index): parity of the number
/// of constraint edges crossed from the outer face.
fn classify(cdt: &ConstrainedDelaunayTriangulation<Point2<f64>>) -> Vec<bool> {
    let mut parity: Vec<i8> = vec![-1; cdt.num_all_faces()];
    let mut queue = VecDeque::new();
    for f in cdt.inner_faces() {
        for e in f.adjacent_edges() {
            if e.rev().face().is_outer() {
                let i = f.fix().index();
                let p = e.as_undirected().is_constraint_edge() as i8;
                if parity[i] < 0 {
                    parity[i] = p;
                    queue.push_back(f.fix());
                }
            }
        }
    }
    while let Some(fh) = queue.pop_front() {
        let f = cdt.face(fh);
        let pf = parity[fh.index()];
        for e in f.adjacent_edges() {
            let Some(nb) = e.rev().face().as_inner() else { continue };
            let j = nb.fix().index();
            if parity[j] >= 0 {
                continue;
            }
            parity[j] = pf ^ (e.as_undirected().is_constraint_edge() as i8);
            queue.push_back(nb.fix());
        }
    }
    parity.into_iter().map(|p| p == 1).collect()
}
