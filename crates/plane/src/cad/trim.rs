//! B-rep topology of a STEP file (shells, faces, bounds, loops, edges,
//! vertices), edge sampling, and the trimmed parameter-space domain of
//! every face.
//!
//! Edges are sampled once, in 3D, and every sample is a global vertex, so
//! two faces sharing an edge share its vertices.  A face's loops are mapped
//! to its surface's (u, v) by its `PCURVE`s when present and otherwise by
//! inverting the 3D samples (Newton), with explicit seam handling on
//! periodic surfaces:
//! * values are unwrapped along each loop by continuity;
//! * an edge used twice by a face (a seam) is placed one period apart on
//!   its two uses, on the side that gives the loop its STEP orientation;
//! * singular points (cone apex, sphere poles) take their free parameter
//!   from their neighbours and are duplicated where the two neighbours
//!   disagree, so the pole becomes a segment of the domain boundary;
//! * a face bounded by two loops that each wind once around the period
//!   (a band without a seam edge) gets an artificial seam joining them, and
//!   a face bounded by one winding loop is closed through the surface pole.

use super::step_parse::{StepFile, Val};
use super::surfaces::{wrap, Curve, Frame, Surface};
use super::NurbsSurface;
use foiley_core::Vec3;
use std::collections::{BTreeMap, HashMap};

/// Import settings.
#[derive(Clone, Debug)]
pub struct StepOptions {
    /// chordal deviation tolerance in metres (0 = `rel_chord_tol` times
    /// the model's bounding-box diagonal)
    pub chord_tol: f64,
    /// chordal tolerance relative to the bounding-box diagonal
    pub rel_chord_tol: f64,
    /// largest turning angle between adjacent edge segments, degrees
    pub angle_tol_deg: f64,
    /// longest triangle / edge-segment side in metres (0 = `rel_max_edge`
    /// times the bounding-box diagonal)
    pub max_edge: f64,
    /// longest side relative to the bounding-box diagonal
    pub rel_max_edge: f64,
    /// fail on any unsupported entity instead of skipping it with a note
    pub strict: bool,
}

impl Default for StepOptions {
    fn default() -> Self {
        StepOptions { chord_tol: 0.0, rel_chord_tol: 1e-3, angle_tol_deg: 20.0, max_edge: 0.0, rel_max_edge: 0.1, strict: false }
    }
}

/// Resolved absolute tolerances.
#[derive(Clone, Copy, Debug)]
pub struct Tol {
    pub chord: f64,
    /// radians
    pub angle: f64,
    pub max_edge: f64,
}

fn seg_dev(p: Vec3, a: Vec3, b: Vec3) -> f64 {
    let ab = b - a;
    let s = ((p - a).dot(ab) / ab.norm2().max(1e-300)).clamp(0.0, 1.0);
    (a + ab * s - p).norm()
}

fn turn(a: Vec3, m: Vec3, b: Vec3) -> f64 {
    let (d1, d2) = (m - a, b - m);
    let (n1, n2) = (d1.norm(), d2.norm());
    if n1 < 1e-300 || n2 < 1e-300 {
        return 0.0;
    }
    (d1.dot(d2) / (n1 * n2)).clamp(-1.0, 1.0).acos()
}

/// Adaptive sampling of the curve `f` on [t0, t1] (both ends included),
/// first split uniformly into `n0` pieces, then bisected until each piece
/// meets the chordal, angle and length tolerances.  Symmetric under
/// reversal of the interval, so a curve shared by two patches is sampled
/// identically from either side.
pub fn sample_curve(f: &dyn Fn(f64) -> Vec3, t0: f64, t1: f64, n0: usize, tol: &Tol) -> Vec<f64> {
    fn rec(f: &dyn Fn(f64) -> Vec3, a: f64, pa: Vec3, b: f64, pb: Vec3, tol: &Tol, depth: usize, out: &mut Vec<f64>) {
        let m = 0.5 * (a + b);
        let pm = f(m);
        let ok = depth >= 18 || {
            let (q1, q3) = (f(a + 0.25 * (b - a)), f(a + 0.75 * (b - a)));
            seg_dev(pm, pa, pb) <= tol.chord && seg_dev(q1, pa, pb) <= tol.chord && seg_dev(q3, pa, pb) <= tol.chord && turn(pa, pm, pb) <= 0.5 * tol.angle && (pb - pa).norm() <= tol.max_edge
        };
        if ok {
            out.push(b);
        } else {
            rec(f, a, pa, m, pm, tol, depth + 1, out);
            rec(f, m, pm, b, pb, tol, depth + 1, out);
        }
    }
    let n0 = n0.max(1);
    let mut out = vec![t0];
    for k in 0..n0 {
        let a = t0 + (t1 - t0) * k as f64 / n0 as f64;
        let b = if k + 1 == n0 { t1 } else { t0 + (t1 - t0) * (k + 1) as f64 / n0 as f64 };
        rec(f, a, f(a), b, f(b), tol, 0, &mut out);
    }
    out
}

/// An edge: its curve, end vertices and samples (in the edge's direction,
/// from `v[0]` to `v[1]`).
#[derive(Clone, Debug)]
pub struct Edge {
    pub step_id: u64,
    pub curve: Curve,
    pub v: [u32; 2],
    /// (curve parameter, global vertex)
    pub samples: Vec<(f64, u32)>,
    /// parameter-space curves: (surface entity, 2D curve with z = 0)
    pub pcurves: Vec<(u64, Curve)>,
}

/// A face bound: oriented edges in loop order.
#[derive(Clone, Debug)]
pub struct Bound {
    pub edges: Vec<(usize, bool)>,
    /// FACE_OUTER_BOUND
    pub outer: bool,
    /// the bound's orientation flag
    pub orientation: bool,
}

#[derive(Clone, Debug)]
pub struct Face {
    pub step_id: u64,
    pub surface_id: u64,
    pub surface: Surface,
    /// face normal = surface normal (dS/du x dS/dv) when true
    pub same_sense: bool,
    pub bounds: Vec<Bound>,
    /// shell (body) index
    pub shell: usize,
    /// boundary of the whole parameter domain, for faces without edge loops
    pub natural: Option<ParamLoop>,
}

/// A closed polygon in a face's parameter space; `ids[i]` is the vertex of
/// `uv[i]` (global, or `>= Model::vertices.len()` for a vertex of
/// `Domain::extra`).  The closing segment (last to first) is implicit.
#[derive(Clone, Debug, Default)]
pub struct ParamLoop {
    pub uv: Vec<[f64; 2]>,
    pub ids: Vec<u32>,
}

impl ParamLoop {
    /// Signed area (positive counter-clockwise in (u, v)).
    pub fn signed_area(&self) -> f64 {
        let n = self.uv.len();
        (0..n)
            .map(|i| {
                let (a, b) = (self.uv[i], self.uv[(i + 1) % n]);
                a[0] * b[1] - b[0] * a[1]
            })
            .sum::<f64>()
            * 0.5
    }
    /// ([u_min, v_min], [u_max, v_max])
    pub fn bbox(&self) -> ([f64; 2], [f64; 2]) {
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for p in &self.uv {
            for k in 0..2 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        (lo, hi)
    }
}

/// The trimmed parameter domain of a face.
#[derive(Clone, Debug, Default)]
pub struct Domain {
    pub loops: Vec<ParamLoop>,
    /// face-local vertices (artificial seams, poles)
    pub extra: Vec<Vec3>,
    /// loop points placed from a PCURVE
    pub used_pcurves: usize,
    /// artificial seams inserted between winding loops
    pub artificial_seams: usize,
    /// loops closed through a surface pole
    pub pole_closures: usize,
}

impl Domain {
    /// Sum of the loops' signed areas.
    pub fn area(&self) -> f64 {
        self.loops.iter().map(|l| l.signed_area()).sum()
    }
}

/// A STEP B-rep model: vertices, edges and faces.
pub struct Model {
    pub vertices: Vec<Vec3>,
    pub edges: Vec<Edge>,
    pub faces: Vec<Face>,
    pub notes: Vec<String>,
    /// unsupported entity kind -> count of faces skipped because of it
    pub unsupported: BTreeMap<String, usize>,
    pub length_scale: f64,
    /// bounding-box diagonal of the model, metres
    pub diag: f64,
    pub tol: Tol,
    pub n_shells: usize,
    /// the file had no faces, only bare surfaces (GEOMETRIC_SET etc.)
    pub bare_surfaces: bool,
}

struct Builder<'a> {
    f: &'a StepFile,
    m: Model,
    vmap: HashMap<u64, u32>,
    emap: HashMap<u64, Result<usize, String>>,
    smap: HashMap<u64, Result<Surface, String>>,
}

/// The entity kind that makes `id` unsupported, for the notes.
fn kind_of(f: &StepFile, id: u64) -> String {
    let k = f.kind(id);
    if k.is_empty() {
        format!("missing entity #{id}")
    } else {
        k.to_string()
    }
}

/// Build a curve from a STEP curve entity (`length`: apply the length unit;
/// false for parameter-space curves).
pub(crate) fn build_curve(f: &StepFile, id: u64, length: bool) -> Result<Curve, String> {
    let e = f.get(id).ok_or_else(|| format!("missing entity #{id}"))?;
    let s = if length { f.length_scale } else { 1.0 };
    let pt = |r: Option<&Val>| -> Option<Vec3> { r.and_then(|v| v.rf()).and_then(|r| if length { f.point(r) } else { f.point_raw(r) }) };
    let frame = |r: Option<&Val>| -> Option<Frame> {
        let fr = f.axis2(r?.rf()?)?;
        if length {
            Some(fr)
        } else {
            // parameter-space placement: undo the unit applied to the origin
            Some(Frame { o: fr.o / f.length_scale, ..fr })
        }
    };
    let bad = || kind_of(f, id);
    if let Some(a) = e.part("LINE") {
        let p = pt(a.get(1)).ok_or_else(bad)?;
        let d = a.get(2).and_then(|v| v.rf()).and_then(|r| f.vector(r, length)).ok_or_else(bad)?;
        return Ok(Curve::Line { p, d });
    }
    if let Some(a) = e.part("CIRCLE") {
        let fr = frame(a.get(1)).ok_or_else(bad)?;
        let r = a.get(2).and_then(|v| v.num()).ok_or_else(bad)? * s;
        return Ok(Curve::Circle { f: fr, r });
    }
    if let Some(a) = e.part("ELLIPSE") {
        let fr = frame(a.get(1)).ok_or_else(bad)?;
        let r1 = a.get(2).and_then(|v| v.num()).ok_or_else(bad)? * s;
        let r2 = a.get(3).and_then(|v| v.num()).ok_or_else(bad)? * s;
        return Ok(Curve::Ellipse { f: fr, a: r1, b: r2 });
    }
    if let Some(a) = e.part("POLYLINE") {
        let pts: Option<Vec<Vec3>> = a.get(1).and_then(|v| v.list()).map(|l| l.iter().map(|v| pt(Some(v))).collect()).ok_or_else(bad)?;
        let pts = pts.ok_or_else(bad)?;
        if pts.len() < 2 {
            return Err(bad());
        }
        return Ok(Curve::Polyline(pts));
    }
    if let Some(a) = e.part("B_SPLINE_CURVE_WITH_KNOTS") {
        return bspline_curve(f, id, a, length).ok_or_else(bad);
    }
    if let Some(a) = e.part("TRIMMED_CURVE") {
        return build_curve(f, a.get(1).and_then(|v| v.rf()).ok_or_else(bad)?, length);
    }
    if let Some(a) = e.part("SURFACE_CURVE").or(e.part("SEAM_CURVE")).or(e.part("INTERSECTION_CURVE")) {
        return build_curve(f, a.get(1).and_then(|v| v.rf()).ok_or_else(bad)?, length);
    }
    Err(bad())
}

/// B_SPLINE_CURVE_WITH_KNOTS, simple or complex (rational) instance.
fn bspline_curve(f: &StepFile, id: u64, a: &[Val], length: bool) -> Option<Curve> {
    let e = f.get(id)?;
    // simple form: name, degree, ctrl, form, closed, self_int, mults, knots, spec
    // complex form: B_SPLINE_CURVE(degree, ctrl, form, closed, si) + B_SPLINE_CURVE_WITH_KNOTS(mults, knots, spec)
    let (deg, ctrl, mk) = if let Some(b) = e.part("B_SPLINE_CURVE") { (b.first()?.num()?, b.get(1)?, &a[..]) } else { (a.get(1)?.num()?, a.get(2)?, &a[6.min(a.len())..]) };
    let pts: Option<Vec<Vec3>> = ctrl.list()?.iter().map(|v| v.rf().and_then(|r| if length { f.point(r) } else { f.point_raw(r) })).collect();
    let pts = pts?;
    let (m, k) = (mk.first()?.list()?, mk.get(1)?.list()?);
    let mut knots = Vec::new();
    for (mm, kk) in m.iter().zip(k) {
        for _ in 0..mm.num()? as usize {
            knots.push(kk.num()?);
        }
    }
    let p = deg as usize;
    if pts.len() < 2 || knots.len() != pts.len() + p + 1 || p > super::MAX_DEGREE {
        return None;
    }
    let w = match e.part("RATIONAL_B_SPLINE_CURVE") {
        Some(r) => r.first()?.list()?.iter().map(|x| x.num().unwrap_or(1.0)).collect(),
        None => vec![1.0; pts.len()],
    };
    if w.len() != pts.len() {
        return None;
    }
    Some(Curve::BSpline(super::NurbsCurve { p, knots, ctrl: pts, w }))
}

/// Parameter-space curves of an edge geometry: (basis surface, 2D curve).
fn pcurves_of(f: &StepFile, geom: u64) -> Vec<(u64, Curve)> {
    let Some(e) = f.get(geom) else { return vec![] };
    let Some(a) = e.part("SURFACE_CURVE").or(e.part("SEAM_CURVE")).or(e.part("INTERSECTION_CURVE")) else { return vec![] };
    let mut out = Vec::new();
    for r in a.get(2).map(|v| v.refs()).unwrap_or_default() {
        let Some(pc) = f.part(r, "PCURVE") else { continue };
        let (Some(surf), Some(rep)) = (pc.get(1).and_then(|v| v.rf()), pc.get(2).and_then(|v| v.rf())) else { continue };
        let Some(items) = f.part(rep, "DEFINITIONAL_REPRESENTATION").and_then(|d| d.get(1)).map(|v| v.refs()) else { continue };
        if let Some(c) = items.first().and_then(|&c| build_curve(f, c, false).ok()) {
            out.push((surf, c));
        }
    }
    out
}

/// Build a surface from a STEP surface entity.
pub(crate) fn build_surface(f: &StepFile, id: u64) -> Result<Surface, String> {
    let e = f.get(id).ok_or_else(|| format!("missing entity #{id}"))?;
    let s = f.length_scale;
    let bad = || kind_of(f, id);
    let frame = |v: Option<&Val>| v.and_then(|v| v.rf()).and_then(|r| f.axis2(r));
    let num = |v: Option<&Val>| v.and_then(|v| v.num());
    if let Some(a) = e.part("PLANE") {
        return Ok(Surface::Plane(frame(a.get(1)).ok_or_else(bad)?));
    }
    if let Some(a) = e.part("CYLINDRICAL_SURFACE") {
        return Ok(Surface::Cylinder { f: frame(a.get(1)).ok_or_else(bad)?, r: num(a.get(2)).ok_or_else(bad)? * s });
    }
    if let Some(a) = e.part("CONICAL_SURFACE") {
        let ang = num(a.get(3)).ok_or_else(bad)? * f.angle_scale;
        return Ok(Surface::Cone { f: frame(a.get(1)).ok_or_else(bad)?, r: num(a.get(2)).ok_or_else(bad)? * s, tan_a: ang.tan() });
    }
    if let Some(a) = e.part("SPHERICAL_SURFACE") {
        return Ok(Surface::Sphere { f: frame(a.get(1)).ok_or_else(bad)?, r: num(a.get(2)).ok_or_else(bad)? * s });
    }
    if let Some(a) = e.part("TOROIDAL_SURFACE").or(e.part("DEGENERATE_TOROIDAL_SURFACE")) {
        return Ok(Surface::Torus { f: frame(a.get(1)).ok_or_else(bad)?, big: num(a.get(2)).ok_or_else(bad)? * s, small: num(a.get(3)).ok_or_else(bad)? * s });
    }
    if let Some(a) = e.part("SURFACE_OF_REVOLUTION") {
        let c = build_curve(f, a.get(1).and_then(|v| v.rf()).ok_or_else(bad)?, true)?;
        let (o, ax) = a.get(2).and_then(|v| v.rf()).and_then(|r| f.axis1(r)).ok_or_else(bad)?;
        return Ok(Surface::Revolution { curve: Box::new(c), o, a: ax });
    }
    if let Some(a) = e.part("SURFACE_OF_LINEAR_EXTRUSION") {
        let c = build_curve(f, a.get(1).and_then(|v| v.rf()).ok_or_else(bad)?, true)?;
        let d = a.get(2).and_then(|v| v.rf()).and_then(|r| f.vector(r, true)).ok_or_else(bad)?;
        return Ok(Surface::Extrusion { curve: Box::new(c), d });
    }
    if e.has("B_SPLINE_SURFACE_WITH_KNOTS") {
        return bspline_surface(f, id).map(Surface::bspline).map_err(|_| bad());
    }
    Err(bad())
}

/// B_SPLINE_SURFACE_WITH_KNOTS, simple or complex (rational) instance.
pub(crate) fn bspline_surface(f: &StepFile, id: u64) -> Result<NurbsSurface, String> {
    let e = f.get(id).ok_or("missing")?;
    let a = e.part("B_SPLINE_SURFACE_WITH_KNOTS").ok_or("not a B-spline surface")?;
    // simple form: name, pu, pv, ctrl, form, uc, vc, si, umult, vmult, uknots, vknots, spec
    // complex form: B_SPLINE_SURFACE(pu, pv, ctrl, ...) + B_SPLINE_SURFACE_WITH_KNOTS(umult, vmult, uknots, vknots, spec)
    let (pu, pv, net, mults) = if let Some(b) = e.part("B_SPLINE_SURFACE") {
        (b.first().and_then(|v| v.num()), b.get(1).and_then(|v| v.num()), b.get(2), &a[..])
    } else {
        (a.get(1).and_then(|v| v.num()), a.get(2).and_then(|v| v.num()), a.get(3), &a[8.min(a.len())..])
    };
    let (Some(pu), Some(pv), Some(net)) = (pu, pv, net) else { return Err("bad B-spline surface".into()) };
    let rows = net.list().ok_or("bad control net")?;
    let ctrl: Option<Vec<Vec<Vec3>>> = rows.iter().map(|r| r.list().and_then(|c| c.iter().map(|p| p.rf().and_then(|r| f.point(r))).collect())).collect();
    let ctrl = ctrl.ok_or("bad control net")?;
    let expand = |m: Option<&Val>, k: Option<&Val>| -> Option<Vec<f64>> {
        let m = m?.list()?;
        let k = k?.list()?;
        let mut out = Vec::new();
        for (mm, kk) in m.iter().zip(k) {
            for _ in 0..mm.num()? as usize {
                out.push(kk.num()?);
            }
        }
        Some(out)
    };
    let (Some(ku), Some(kv)) = (expand(mults.first(), mults.get(2)), expand(mults.get(1), mults.get(3))) else { return Err("bad knots".into()) };
    let nu = ctrl.len();
    let nv = ctrl.first().map(|r| r.len()).unwrap_or(0);
    if nu < 2 || nv < 2 || ctrl.iter().any(|r| r.len() != nv) || ku.len() != nu + pu as usize + 1 || kv.len() != nv + pv as usize + 1 {
        return Err("inconsistent knot vectors".into());
    }
    let w = match e.part("RATIONAL_B_SPLINE_SURFACE") {
        Some(r) => r.first().and_then(|l| l.list()).map(|rows| rows.iter().map(|row| row.list().map(|c| c.iter().map(|x| x.num().unwrap_or(1.0)).collect()).unwrap_or_default()).collect()).unwrap_or_else(|| vec![vec![1.0; nv]; nu]),
        None => vec![vec![1.0; nv]; nu],
    };
    Ok(NurbsSurface { pu: pu as usize, pv: pv as usize, ku, kv, ctrl, w })
}

impl<'a> Builder<'a> {
    fn vertex(&mut self, vp: u64) -> Option<u32> {
        if let Some(&i) = self.vmap.get(&vp) {
            return Some(i);
        }
        let p = self.f.part(vp, "VERTEX_POINT")?.get(1)?.rf().and_then(|r| self.f.point(r))?;
        let i = self.push_vertex(p);
        self.vmap.insert(vp, i);
        Some(i)
    }
    fn push_vertex(&mut self, p: Vec3) -> u32 {
        self.m.vertices.push(p);
        (self.m.vertices.len() - 1) as u32
    }

    /// EDGE_CURVE -> edge index (sampled), or the unsupported kind.
    fn edge(&mut self, id: u64) -> Result<usize, String> {
        if let Some(r) = self.emap.get(&id) {
            return r.clone();
        }
        let r = self.make_edge(id);
        self.emap.insert(id, r.clone());
        r
    }
    fn make_edge(&mut self, id: u64) -> Result<usize, String> {
        let f = self.f;
        let a = f.part(id, "EDGE_CURVE").ok_or_else(|| kind_of(f, id))?;
        let v0 = a.get(1).and_then(|v| v.rf()).and_then(|r| self.vertex(r)).ok_or("VERTEX_POINT")?;
        let v1 = a.get(2).and_then(|v| v.rf()).and_then(|r| self.vertex(r)).ok_or("VERTEX_POINT")?;
        let geom = a.get(3).and_then(|v| v.rf()).ok_or("EDGE_CURVE")?;
        let same = a.get(4).and_then(|v| v.boolean()).unwrap_or(true);
        let curve = build_curve(f, geom, true)?;
        let pcurves = pcurves_of(f, geom);
        let (p0, p1) = (self.m.vertices[v0 as usize], self.m.vertices[v1 as usize]);
        let closed = v0 == v1 || (p1 - p0).norm() <= 1e-9 * self.m.diag;
        let mut t0 = curve.param_of(p0);
        let mut t1 = curve.param_of(p1);
        if let Some(per) = curve.period() {
            if closed {
                t1 = if same { t0 + per } else { t0 - per };
            } else if same {
                let d = (t1 - t0).rem_euclid(per);
                t1 = t0 + if d <= 1e-12 { per } else { d };
            } else {
                let d = (t0 - t1).rem_euclid(per);
                t1 = t0 - if d <= 1e-12 { per } else { d };
            }
        } else if closed {
            if let Some((ra, rb)) = curve.range() {
                (t0, t1) = if same { (ra, rb) } else { (rb, ra) };
            }
        }
        let n0 = if closed { 4.max(curve.pieces()) } else { curve.pieces() };
        let c = &curve;
        let ts = if closed && curve.period().is_none() && curve.range().is_none() {
            vec![t0, t1] // degenerate (zero-length) edge
        } else {
            sample_curve(&|t| c.eval(t), t0, t1, n0, &self.m.tol)
        };
        let mut samples = Vec::with_capacity(ts.len());
        let last = ts.len() - 1;
        for (k, &t) in ts.iter().enumerate() {
            let vid = if k == 0 {
                v0
            } else if k == last {
                v1
            } else {
                let p = curve.eval(t);
                self.push_vertex(p)
            };
            samples.push((t, vid));
        }
        self.m.edges.push(Edge { step_id: id, curve, v: [v0, v1], samples, pcurves });
        Ok(self.m.edges.len() - 1)
    }

    fn surface(&mut self, id: u64) -> Result<Surface, String> {
        if let Some(s) = self.smap.get(&id) {
            return s.clone();
        }
        let s = build_surface(self.f, id);
        self.smap.insert(id, s.clone());
        s
    }

    fn skip(&mut self, what: String, face: u64) -> Result<(), String> {
        if self.m.notes.len() < 200 {
            self.m.notes.push(format!("face #{face}: {what} is not supported; face skipped"));
        }
        *self.m.unsupported.entry(what).or_default() += 1;
        Ok(())
    }

    fn face(&mut self, id: u64, shell: usize) -> Result<(), String> {
        let f = self.f;
        let Some(a) = f.part(id, "ADVANCED_FACE").or(f.part(id, "FACE_SURFACE")) else {
            return self.skip(kind_of(f, id), id);
        };
        let Some(sid) = a.get(2).and_then(|v| v.rf()) else { return self.skip("face without surface".into(), id) };
        let same_sense = a.get(3).and_then(|v| v.boolean()).unwrap_or(true);
        let surface = match self.surface(sid) {
            Ok(s) => s,
            Err(k) => return self.skip(k, id),
        };
        let mut bounds = Vec::new();
        for b in a.get(1).map(|v| v.refs()).unwrap_or_default() {
            let Some(be) = f.get(b) else { return self.skip(format!("missing bound #{b}"), id) };
            let outer = be.has("FACE_OUTER_BOUND");
            let Some(ba) = be.part("FACE_OUTER_BOUND").or(be.part("FACE_BOUND")) else { return self.skip(be.kind().to_string(), id) };
            let lp = ba.get(1).and_then(|v| v.rf()).unwrap_or(0);
            let orientation = ba.get(2).and_then(|v| v.boolean()).unwrap_or(true);
            let Some(le) = f.get(lp) else { return self.skip(format!("missing loop #{lp}"), id) };
            if le.has("VERTEX_LOOP") {
                continue; // a pole: the domain closes through it
            }
            let mut edges = Vec::new();
            if let Some(el) = le.part("EDGE_LOOP") {
                for oe in el.get(1).map(|v| v.refs()).unwrap_or_default() {
                    let Some(o) = f.part(oe, "ORIENTED_EDGE") else { return self.skip(kind_of(f, oe), id) };
                    let fwd = o.get(4).and_then(|v| v.boolean()).unwrap_or(true);
                    let Some(ec) = o.get(3).and_then(|v| v.rf()) else { return self.skip("ORIENTED_EDGE".into(), id) };
                    match self.edge(ec) {
                        Ok(ei) => edges.push((ei, fwd)),
                        Err(k) => return self.skip(k, id),
                    }
                }
            } else if let Some(pl) = le.part("POLY_LOOP") {
                let pts: Vec<Vec3> = pl.get(1).map(|v| v.refs()).unwrap_or_default().into_iter().filter_map(|r| f.point(r)).collect();
                for i in 0..pts.len() {
                    let (p, q) = (pts[i], pts[(i + 1) % pts.len()]);
                    let ei = self.poly_edge(p, q);
                    edges.push((ei, true));
                }
            } else {
                return self.skip(le.kind().to_string(), id);
            }
            if !edges.is_empty() {
                bounds.push(Bound { edges, outer, orientation });
            }
        }
        let natural = if bounds.is_empty() {
            match self.natural_loop(&surface) {
                Some(l) => Some(l),
                None => return self.skip(format!("unbounded {}", surface.kind()), id),
            }
        } else {
            None
        };
        self.m.faces.push(Face { step_id: id, surface_id: sid, surface, same_sense, bounds, shell, natural });
        Ok(())
    }

    /// Straight edge of a POLY_LOOP (vertices shared by position).
    fn poly_edge(&mut self, p: Vec3, q: Vec3) -> usize {
        let mut key = |p: Vec3| -> u32 {
            if let Some(i) = self.m.vertices.iter().position(|v| (*v - p).norm() <= 1e-12 * self.m.diag) {
                return i as u32;
            }
            self.push_vertex(p)
        };
        let (a, b) = (key(p), key(q));
        if let Some(i) = self.m.edges.iter().position(|e| e.v == [a, b] || e.v == [b, a]) {
            return i;
        }
        let curve = Curve::Line { p, d: q - p };
        let c = &curve;
        let ts = sample_curve(&|t| c.eval(t), 0.0, 1.0, 1, &self.m.tol);
        let last = ts.len() - 1;
        let samples = ts.iter().enumerate().map(|(k, &t)| (t, if k == 0 { a } else if k == last { b } else { self.push_vertex(c.eval(t)) })).collect();
        self.m.edges.push(Edge { step_id: 0, curve, v: [a, b], samples, pcurves: vec![] });
        self.m.edges.len() - 1
    }

    /// Boundary of the whole (bounded) parameter domain: the four sides,
    /// with periodic sides sharing vertices and collapsed sides (poles)
    /// reduced to one vertex.
    fn natural_loop(&mut self, s: &Surface) -> Option<ParamLoop> {
        let [(u0, u1), (v0, v1)] = s.natural_domain()?;
        let per = s.periods();
        let corners = [[u0, v0], [u1, v0], [u1, v1], [u0, v1]];
        let side = |k: usize| (corners[k], corners[(k + 1) % 4]);
        let tol = self.m.tol;
        let diag = self.m.diag;
        let sample_side = |a: [f64; 2], b: [f64; 2]| -> (Vec<f64>, bool) {
            let f = |t: f64| s.eval(a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t);
            let ts = sample_curve(&f, 0.0, 1.0, 4, &tol);
            let p0 = f(0.0);
            let degenerate = ts.iter().all(|&t| (f(t) - p0).norm() <= 1e-9 * diag) && (0..=8).all(|i| (f(i as f64 / 8.0) - p0).norm() <= 1e-9 * diag);
            (ts, degenerate)
        };
        // union-find over the four corners
        let mut parent = [0usize, 1, 2, 3];
        fn find(p: &mut [usize; 4], i: usize) -> usize {
            let mut i = i;
            while p[i] != i {
                i = p[i];
            }
            i
        }
        let join = |p: &mut [usize; 4], a: usize, b: usize| {
            let (ra, rb) = (find(p, a), find(p, b));
            if ra != rb {
                p[ra.max(rb)] = ra.min(rb);
            }
        };
        let sides: Vec<(Vec<f64>, bool)> = (0..4).map(|k| sample_side(side(k).0, side(k).1)).collect();
        if per[0].is_some() {
            join(&mut parent, 0, 1);
            join(&mut parent, 3, 2);
        }
        if per[1].is_some() {
            join(&mut parent, 0, 3);
            join(&mut parent, 1, 2);
        }
        for k in 0..4 {
            if sides[k].1 {
                join(&mut parent, k, (k + 1) % 4);
            }
        }
        let mut cid = [u32::MAX; 4];
        for k in 0..4 {
            let r = find(&mut parent, k);
            if cid[r] == u32::MAX {
                let c = corners[r];
                cid[r] = self.push_vertex(s.eval(c[0], c[1]));
            }
            cid[k] = cid[r];
        }
        // interior samples; the side opposite a periodic one reuses its vertices
        let mut interior: Vec<Vec<u32>> = vec![vec![]; 4];
        let mut params: Vec<Vec<f64>> = vec![vec![]; 4];
        for k in 0..4 {
            let partner = match k {
                2 if per[1].is_some() => Some(0),
                3 if per[0].is_some() => Some(1),
                _ => None,
            };
            if let Some(pk) = partner {
                // same curve, traversed the other way
                params[k] = params[pk].iter().rev().map(|t| 1.0 - t).collect();
                interior[k] = interior[pk].iter().rev().copied().collect();
                continue;
            }
            let (ts, degenerate) = &sides[k];
            let (a, b) = side(k);
            let inner: Vec<f64> = ts[1..ts.len() - 1].to_vec();
            interior[k] = inner
                .iter()
                .map(|&t| if *degenerate { cid[k] } else { self.push_vertex(s.eval(a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)) })
                .collect();
            params[k] = inner;
        }
        let mut lp = ParamLoop::default();
        for k in 0..4 {
            let (a, b) = side(k);
            lp.uv.push(a);
            lp.ids.push(cid[k]);
            for (t, id) in params[k].iter().zip(&interior[k]) {
                lp.uv.push([a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t]);
                lp.ids.push(*id);
            }
        }
        Some(lp)
    }
}

/// Bounding-box diagonal of a STEP file's geometry (points, and the extent
/// of circles and elementary surfaces).
fn file_diag(f: &StepFile) -> f64 {
    let mut lo = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = -lo;
    let mut add = |p: Vec3, r: f64| {
        lo = Vec3::new(lo.x.min(p.x - r), lo.y.min(p.y - r), lo.z.min(p.z - r));
        hi = Vec3::new(hi.x.max(p.x + r), hi.y.max(p.y + r), hi.z.max(p.z + r));
    };
    let mut in_2d = std::collections::HashSet::new();
    for id in f.ids_of("DEFINITIONAL_REPRESENTATION") {
        // points of parameter-space curves are not geometry
        let mut stack: Vec<u64> = f.part(id, "DEFINITIONAL_REPRESENTATION").and_then(|a| a.get(1)).map(|v| v.refs()).unwrap_or_default();
        let mut guard = 0;
        while let Some(r) = stack.pop() {
            guard += 1;
            if guard > 10_000 || !in_2d.insert(r) {
                continue;
            }
            if let Some(e) = f.get(r) {
                for (_, args) in &e.parts {
                    for a in args {
                        if let Some(x) = a.rf() {
                            stack.push(x);
                        }
                        stack.extend(a.refs());
                    }
                }
            }
        }
    }
    for id in f.ids_of("CARTESIAN_POINT") {
        if in_2d.contains(&id) {
            continue;
        }
        let dims = f.part(id, "CARTESIAN_POINT").and_then(|a| a.get(1)).and_then(|v| v.list()).map(|l| l.len()).unwrap_or(0);
        if dims == 3 {
            if let Some(p) = f.point(id) {
                add(p, 0.0);
            }
        }
    }
    for (kind, k) in [("CIRCLE", 2), ("ELLIPSE", 2), ("CYLINDRICAL_SURFACE", 2), ("CONICAL_SURFACE", 2), ("SPHERICAL_SURFACE", 2), ("TOROIDAL_SURFACE", 2)] {
        for id in f.ids_of(kind) {
            let a = f.part(id, kind).unwrap();
            let r = a.get(k).and_then(|v| v.num()).unwrap_or(0.0) * f.length_scale + if kind == "TOROIDAL_SURFACE" { a.get(3).and_then(|v| v.num()).unwrap_or(0.0) * f.length_scale } else { 0.0 };
            if let Some(fr) = a.get(1).and_then(|v| v.rf()).and_then(|r| f.axis2(r)) {
                add(fr.o, r.abs());
            }
        }
    }
    let d = (hi - lo).norm();
    if d.is_finite() && d > 0.0 {
        d
    } else {
        1.0
    }
}

impl Model {
    pub fn from_step(text: &str, opts: &StepOptions) -> Result<Model, String> {
        Model::from_file(&StepFile::parse(text)?, opts)
    }

    pub fn from_file(f: &StepFile, opts: &StepOptions) -> Result<Model, String> {
        let diag = file_diag(f);
        let tol = Tol {
            chord: if opts.chord_tol > 0.0 { opts.chord_tol } else { opts.rel_chord_tol.max(1e-9) * diag },
            angle: opts.angle_tol_deg.clamp(1.0, 90.0).to_radians(),
            max_edge: if opts.max_edge > 0.0 { opts.max_edge } else { opts.rel_max_edge.max(1e-6) * diag },
        };
        let m = Model { vertices: vec![], edges: vec![], faces: vec![], notes: vec![], unsupported: BTreeMap::new(), length_scale: f.length_scale, diag, tol, n_shells: 0, bare_surfaces: false };
        let mut b = Builder { f, m, vmap: HashMap::new(), emap: HashMap::new(), smap: HashMap::new() };
        // shells and their faces; faces outside any shell form one more body
        let mut in_shell = std::collections::HashSet::new();
        let mut shells: Vec<Vec<u64>> = Vec::new();
        for kind in ["CLOSED_SHELL", "OPEN_SHELL", "CONNECTED_FACE_SET"] {
            for id in f.ids_of(kind) {
                let faces = f.part(id, kind).and_then(|a| a.get(1)).map(|v| v.refs()).unwrap_or_default();
                in_shell.extend(faces.iter().copied());
                shells.push(faces);
            }
        }
        let loose: Vec<u64> = f.ids_of("ADVANCED_FACE").into_iter().chain(f.ids_of("FACE_SURFACE")).filter(|i| !in_shell.contains(i)).collect();
        if !loose.is_empty() {
            shells.push(loose);
        }
        for (k, faces) in shells.iter().enumerate() {
            for &fid in faces {
                b.face(fid, k)?;
            }
        }
        b.m.n_shells = shells.len();
        if shells.is_empty() {
            // no topology: bare surfaces (each surface is a face on its whole domain)
            b.m.bare_surfaces = true;
            let mut ids: Vec<u64> = f.ids_of("GEOMETRIC_SET").into_iter().chain(f.ids_of("GEOMETRIC_CURVE_SET")).flat_map(|g| f.part(g, f.kind(g)).and_then(|a| a.get(1)).map(|v| v.refs()).unwrap_or_default()).collect();
            if ids.is_empty() {
                ids = f.ids_of("B_SPLINE_SURFACE_WITH_KNOTS");
            }
            let first_bare_vertex = b.m.vertices.len();
            for (k, sid) in ids.into_iter().enumerate() {
                let kind = f.kind(sid).to_string();
                if kind.contains("CURVE") || kind == "CARTESIAN_POINT" {
                    continue;
                }
                match b.surface(sid) {
                    Ok(s) => match b.natural_loop(&s) {
                        Some(l) => b.m.faces.push(Face { step_id: sid, surface_id: sid, surface: s, same_sense: true, bounds: vec![], shell: k, natural: Some(l) }),
                        None => b.skip(format!("unbounded {}", s.kind()), sid)?,
                    },
                    Err(kd) => b.skip(kd, sid)?,
                }
            }
            b.m.n_shells = b.m.faces.len();
            b.m.weld_natural(first_bare_vertex);
        }
        let mut m = b.m;
        if !m.unsupported.is_empty() {
            let list: Vec<String> = m.unsupported.iter().map(|(k, n)| format!("{n} x {k}")).collect();
            let msg = format!("unsupported STEP entities: {}", list.join(", "));
            if opts.strict {
                return Err(msg);
            }
            m.notes.insert(0, msg);
        }
        if m.faces.is_empty() {
            return Err("no supported faces or surfaces in the STEP file".into());
        }
        Ok(m)
    }

    /// Merge the boundary vertices (from `first` on) of bare surfaces that
    /// coincide, so patches that share a boundary share its vertices.
    fn weld_natural(&mut self, first: usize) {
        let tol = 1e-7 * self.diag;
        let key = |p: Vec3| ((p.x / tol).floor() as i64, (p.y / tol).floor() as i64, (p.z / tol).floor() as i64);
        let mut grid: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
        let mut map: Vec<u32> = (0..self.vertices.len() as u32).collect();
        for i in first..self.vertices.len() {
            let p = self.vertices[i];
            let k = key(p);
            let mut found = None;
            'search: for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(c) = grid.get(&(k.0 + dx, k.1 + dy, k.2 + dz)) {
                            for &j in c {
                                if (self.vertices[j as usize] - p).norm() <= tol {
                                    found = Some(j);
                                    break 'search;
                                }
                            }
                        }
                    }
                }
            }
            match found {
                Some(j) => map[i] = j,
                None => grid.entry(k).or_default().push(i as u32),
            }
        }
        for f in &mut self.faces {
            if let Some(l) = &mut f.natural {
                for id in &mut l.ids {
                    *id = map[*id as usize];
                }
            }
        }
    }

    /// Kinds of the supported surfaces, with counts.
    pub fn surface_kinds(&self) -> BTreeMap<&'static str, usize> {
        let mut m = BTreeMap::new();
        for f in &self.faces {
            *m.entry(f.surface.kind()).or_default() += 1;
        }
        m
    }

    /// Kinds of the edge curves, with counts.
    pub fn curve_kinds(&self) -> BTreeMap<&'static str, usize> {
        let mut m = BTreeMap::new();
        for e in &self.edges {
            *m.entry(e.curve.kind()).or_default() += 1;
        }
        m
    }

    /// Position of a domain vertex id.
    pub fn vertex(&self, d: &Domain, id: u32) -> Vec3 {
        let n = self.vertices.len();
        if (id as usize) < n {
            self.vertices[id as usize]
        } else {
            d.extra[id as usize - n]
        }
    }

    /// The trimmed parameter-space domain of face `fi`.
    pub fn face_domain(&self, fi: usize) -> Result<Domain, String> {
        let face = &self.faces[fi];
        let mut dom = Domain::default();
        if let Some(l) = &face.natural {
            dom.loops.push(l.clone());
            return Ok(dom);
        }
        let s = &face.surface;
        let per = s.periods();
        let mut uses: HashMap<usize, usize> = HashMap::new();
        for b in &face.bounds {
            for &(e, _) in &b.edges {
                *uses.entry(e).or_default() += 1;
            }
        }
        let face_sign = if face.same_sense { 1.0 } else { -1.0 };
        let mut loops: Vec<(Vec<LoopPt>, [i64; 2], f64)> = Vec::new();
        for b in &face.bounds {
            let pts = self.raw_loop(face, b, &uses, &mut dom.used_pcurves);
            if pts.is_empty() {
                continue;
            }
            let outer = b.outer || face.bounds.len() == 1;
            let expected = if outer || face.bounds.iter().any(|x| x.outer) { face_sign * if b.orientation { 1.0 } else { -1.0 } * if outer { 1.0 } else { -1.0 } } else { 0.0 };
            let (fin, w) = unwrap_loop(s, &pts, b, &uses, expected);
            loops.push((fin, w, expected));
        }
        // winding loops need closing: two opposite ones get an artificial seam, a single one closes through a pole
        let winding: Vec<usize> = (0..loops.len()).filter(|&i| loops[i].1 != [0, 0]).collect();
        let mut plain: Vec<ParamLoop> = loops.iter().enumerate().filter(|(i, _)| !winding.contains(i)).map(|(_, l)| to_param(&l.0)).collect();
        match winding.len() {
            0 => {}
            2 => {
                let (a, b) = (&loops[winding[0]], &loops[winding[1]]);
                let d = if a.1[0] != 0 { 0 } else { 1 };
                if a.1[1 - d] != 0 || b.1[1 - d] != 0 || a.1[d] != -b.1[d] || a.1[d].abs() != 1 {
                    return Err(format!("face #{}: loops wind around the surface in an unsupported way", face.step_id));
                }
                let l = self.join_seam(s, &mut dom, &a.0, &b.0, d, a.1[d], per[d].unwrap());
                plain.push(l);
                dom.artificial_seams += 1;
            }
            1 => {
                let a = &loops[winding[0]];
                let d = if a.1[0] != 0 { 0 } else { 1 };
                if a.1[1 - d] != 0 || a.1[d].abs() != 1 {
                    return Err(format!("face #{}: a loop winds around the surface in an unsupported way", face.step_id));
                }
                let l = self.close_pole(s, &mut dom, &a.0, d, a.1[d], a.2, per[d].unwrap())?;
                plain.push(l);
                dom.pole_closures += 1;
            }
            _ => return Err(format!("face #{}: {} loops wind around the surface (unsupported)", face.step_id, winding.len())),
        }
        dom.loops = plain;
        Ok(dom)
    }

    fn raw_loop(&self, face: &Face, b: &Bound, uses: &HashMap<usize, usize>, used_pc: &mut usize) -> Vec<LoopPt> {
        let s = &face.surface;
        let mut pts = Vec::new();
        let mut prev: Option<[f64; 2]> = None;
        for (k, &(ei, fwd)) in b.edges.iter().enumerate() {
            let e = &self.edges[ei];
            let seam = uses[&ei] > 1;
            let pc = if seam { None } else { e.pcurves.iter().find(|(sid, _)| *sid == face.surface_id).map(|(_, c)| c) };
            let n = e.samples.len();
            for j in 0..n.saturating_sub(1) {
                let (t, id) = if fwd { e.samples[j] } else { e.samples[n - 1 - j] };
                let p = self.vertices[id as usize];
                let mut uv = None;
                if let Some(c) = pc {
                    let q = c.eval(t);
                    let cand = s.newton(p, [q.x, q.y]);
                    if (s.eval(cand[0], cand[1]) - p).norm() <= 1e-6 * self.diag {
                        uv = Some(cand);
                        *used_pc += 1;
                    }
                }
                let uv = uv.unwrap_or_else(|| s.invert(p, prev));
                pts.push(LoopPt { uv, id, use_k: k, sing: s.singular(uv[0], uv[1]) });
                prev = Some(uv);
            }
        }
        pts
    }

    /// Join two loops winding in direction `d` (`wa` = winding of `a`) by a
    /// seam, giving one simple polygon.
    fn join_seam(&self, s: &Surface, dom: &mut Domain, a: &[LoopPt], b: &[LoopPt], d: usize, wa: i64, per: f64) -> ParamLoop {
        let wb = -wa;
        let shift = |p: [f64; 2], k: f64| {
            let mut q = p;
            q[d] += k * per;
            q
        };
        // start of a: its first non-singular point; partner in b: closest in the seam direction
        let ia = a.iter().position(|p| !p.sing.iter().any(|x| *x)).unwrap_or(0);
        let pa = a[ia].uv;
        let jb = (0..b.len()).filter(|&j| !b[j].sing.iter().any(|x| *x)).min_by(|&i, &j| wrap(b[i].uv[d] - pa[d], per).abs().total_cmp(&wrap(b[j].uv[d] - pa[d], per).abs())).unwrap_or(0);
        let mut out = ParamLoop::default();
        let n = a.len();
        for k in 0..=n {
            let i = (ia + k) % n;
            let p = if ia + k >= n { shift(a[i].uv, wa as f64) } else { a[i].uv };
            out.uv.push(p);
            out.ids.push(a[i].id);
        }
        let a_end = *out.uv.last().unwrap();
        let m = b.len();
        let mut bp: Vec<([f64; 2], u32)> = (0..=m)
            .map(|k| {
                let j = (jb + k) % m;
                (if jb + k >= m { shift(b[j].uv, wb as f64) } else { b[j].uv }, b[j].id)
            })
            .collect();
        let off = ((a_end[d] - bp[0].0[d]) / per).round();
        for q in &mut bp {
            q.0 = shift(q.0, off);
        }
        // seam from the end of a to the start of b, sampled once and used twice
        let (s0, s1) = (a_end, bp[0].0);
        let seam = self.param_line(s, dom, s0, s1);
        for (uv, id) in &seam {
            out.uv.push(*uv);
            out.ids.push(*id);
        }
        for (uv, id) in &bp {
            out.uv.push(*uv);
            out.ids.push(*id);
        }
        for (uv, id) in seam.iter().rev() {
            out.uv.push(shift(*uv, -(wa as f64)));
            out.ids.push(*id);
        }
        out
    }

    /// Close a loop winding once in direction `d` through the surface pole
    /// (in the other direction) on the face's side.
    fn close_pole(&self, s: &Surface, dom: &mut Domain, a: &[LoopPt], d: usize, wa: i64, expected: f64, per: f64) -> Result<ParamLoop, String> {
        let o = 1 - d;
        let (vmin, vmax) = a.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| (lo.min(p.uv[o]), hi.max(p.uv[o])));
        let poles = s.poles(o);
        // the face lies left of the loop when it is counter-clockwise: +v of a +u loop, -u of a +v loop
        let side = wa as f64 * expected * if d == 0 { 1.0 } else { -1.0 };
        let eps = 1e-9 * (1.0 + vmax.abs() + vmin.abs());
        let up = poles.iter().copied().filter(|&v| v >= vmax - eps).min_by(|x, y| x.total_cmp(y));
        let down = poles.iter().copied().filter(|&v| v <= vmin + eps).max_by(|x, y| x.total_cmp(y));
        let vp = if side > 0.0 {
            up
        } else if side < 0.0 {
            down
        } else {
            up.or(down)
        }
        .ok_or_else(|| format!("{}: a loop winds around the surface and there is no pole to close it", s.kind()))?;
        let ia = a.iter().position(|p| !p.sing.iter().any(|x| *x)).unwrap_or(0);
        let n = a.len();
        let mut out = ParamLoop::default();
        for k in 0..=n {
            let i = (ia + k) % n;
            let mut p = a[i].uv;
            if ia + k >= n {
                p[d] += wa as f64 * per;
            }
            out.uv.push(p);
            out.ids.push(a[i].id);
        }
        let back = |mut q: [f64; 2]| {
            q[d] -= wa as f64 * per;
            q
        };
        let a_end = *out.uv.last().unwrap();
        let mut top = a_end;
        top[o] = vp;
        let line = self.param_line(s, dom, a_end, top);
        let pole = self.vertices.len() as u32 + dom.extra.len() as u32;
        dom.extra.push(s.eval(top[0], top[1]));
        for (uv, id) in &line {
            out.uv.push(*uv);
            out.ids.push(*id);
        }
        out.uv.push(top);
        out.ids.push(pole);
        out.uv.push(back(top));
        out.ids.push(pole);
        for (uv, id) in line.iter().rev() {
            out.uv.push(back(*uv));
            out.ids.push(*id);
        }
        Ok(out)
    }

    /// Interior samples of the straight parameter line p0 -> p1 (new
    /// face-local vertices).
    fn param_line(&self, s: &Surface, dom: &mut Domain, p0: [f64; 2], p1: [f64; 2]) -> Vec<([f64; 2], u32)> {
        let at = |t: f64| [p0[0] + (p1[0] - p0[0]) * t, p0[1] + (p1[1] - p0[1]) * t];
        let ts = sample_curve(&|t| {
            let q = at(t);
            s.eval(q[0], q[1])
        }, 0.0, 1.0, 1, &self.tol);
        let mut out = Vec::new();
        for &t in &ts[1..ts.len() - 1] {
            let q = at(t);
            let id = self.vertices.len() as u32 + dom.extra.len() as u32;
            dom.extra.push(s.eval(q[0], q[1]));
            out.push((q, id));
        }
        out
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct LoopPt {
    pub uv: [f64; 2],
    pub id: u32,
    /// index of the oriented edge (in its bound) this point comes from
    pub use_k: usize,
    pub sing: [bool; 2],
}

fn to_param(pts: &[LoopPt]) -> ParamLoop {
    ParamLoop { uv: pts.iter().map(|p| p.uv).collect(), ids: pts.iter().map(|p| p.id).collect() }
}

/// Unwrap a loop's periodic parameters (see the module notes).  Returns the
/// final points (singular points duplicated) and the loop's winding numbers.
fn unwrap_loop(s: &Surface, pts: &[LoopPt], b: &Bound, uses: &HashMap<usize, usize>, expected: f64) -> (Vec<LoopPt>, [i64; 2]) {
    let per = s.periods();
    // seams: edges used twice in this loop -> (first use, second use, direction)
    let mut seams: Vec<(usize, usize, usize)> = Vec::new();
    for (k, &(e, _)) in b.edges.iter().enumerate() {
        if uses[&e] < 2 {
            continue;
        }
        if let Some(k2) = (k + 1..b.edges.len()).find(|&j| b.edges[j].0 == e) {
            // seam direction: the parameter that is constant (mod period) along the edge
            let vals: Vec<[f64; 2]> = pts.iter().filter(|p| p.use_k == k && !p.sing.iter().any(|x| *x)).map(|p| p.uv).collect();
            let constant = |d: usize| per[d].map(|pp| vals.windows(2).all(|w| wrap(w[1][d] - w[0][d], pp).abs() < 1e-7 * pp)).unwrap_or(false);
            let dir = if constant(0) {
                0
            } else if constant(1) {
                1
            } else {
                continue;
            };
            seams.push((k, k2, dir));
        }
    }
    seams.truncate(4);
    let combos = 1usize << seams.len();
    let mut best: Option<((usize, u8), Vec<LoopPt>, [i64; 2])> = None;
    for c in 0..combos {
        let signs: Vec<f64> = (0..seams.len()).map(|j| if c >> j & 1 == 1 { -1.0 } else { 1.0 }).collect();
        let vals = unwrap_with(pts, &per, &seams, &signs);
        let fin = duplicate_singular(pts, &vals);
        let (jumps, w) = jumps_and_winding(&fin, &per);
        let area = to_param(&fin).signed_area();
        let mismatch = if expected != 0.0 && area * expected < 0.0 { 1u8 } else { 0 };
        // a winding loop counts its closure as a jump only if it has seams (seamed loops must close)
        let bad = if seams.is_empty() { jumps.saturating_sub((w[0] != 0) as usize + (w[1] != 0) as usize) } else { jumps };
        let score = (bad, mismatch);
        if best.as_ref().map(|b| score < b.0).unwrap_or(true) {
            best = Some((score, fin, w));
        }
        if seams.is_empty() {
            break;
        }
    }
    let (_, fin, w) = best.unwrap();
    (fin, w)
}

fn unwrap_with(pts: &[LoopPt], per: &[Option<f64>; 2], seams: &[(usize, usize, usize)], signs: &[f64]) -> Vec<[f64; 2]> {
    let mut out: Vec<[f64; 2]> = pts.iter().map(|p| p.uv).collect();
    for d in 0..2 {
        let Some(pp) = per[d] else { continue };
        let mut last: Option<f64> = None;
        let mut first_vals: Vec<HashMap<u32, f64>> = vec![HashMap::new(); seams.len()];
        for i in 0..pts.len() {
            let p = &pts[i];
            if p.sing[d] {
                if let Some(l) = last {
                    out[i][d] = l;
                }
                continue;
            }
            let raw = p.uv[d];
            let cont = match last {
                Some(l) => raw + pp * ((l - raw) / pp).round(),
                None => raw,
            };
            let mut val = cont;
            for (j, &(k1, k2, dir)) in seams.iter().enumerate() {
                if dir != d {
                    continue;
                }
                if p.use_k == k2 {
                    if let Some(f) = first_vals[j].get(&p.id) {
                        val = f + signs[j] * pp;
                    }
                }
                if p.use_k == k1 {
                    first_vals[j].insert(p.id, cont);
                }
            }
            out[i][d] = val;
            last = Some(val);
        }
    }
    out
}

/// Give singular points the free parameter of their neighbours, duplicating
/// them where the previous and next neighbours disagree.
fn duplicate_singular(pts: &[LoopPt], vals: &[[f64; 2]]) -> Vec<LoopPt> {
    let n = pts.len();
    let regular = |i: usize, d: usize| !pts[i].sing[d];
    let mut out = Vec::with_capacity(n + 4);
    for i in 0..n {
        let p = pts[i];
        if !p.sing[0] && !p.sing[1] {
            out.push(LoopPt { uv: vals[i], ..p });
            continue;
        }
        let mut a = vals[i];
        let mut b = vals[i];
        for d in 0..2 {
            if !p.sing[d] {
                continue;
            }
            if let Some(j) = (1..n).map(|k| (i + n - k) % n).find(|&j| regular(j, d)) {
                a[d] = vals[j][d];
            }
            if let Some(j) = (1..n).map(|k| (i + k) % n).find(|&j| regular(j, d)) {
                b[d] = vals[j][d];
            }
        }
        out.push(LoopPt { uv: a, ..p });
        if (b[0] - a[0]).abs() > 1e-12 || (b[1] - a[1]).abs() > 1e-12 {
            out.push(LoopPt { uv: b, ..p });
        }
    }
    out
}

/// Jumps of more than half a period between consecutive regular points
/// (including the closing pair), and the winding numbers.
fn jumps_and_winding(fin: &[LoopPt], per: &[Option<f64>; 2]) -> (usize, [i64; 2]) {
    let n = fin.len();
    let mut jumps = 0;
    let mut w = [0i64; 2];
    for d in 0..2 {
        let Some(pp) = per[d] else { continue };
        let mut acc = 0.0;
        for i in 0..n {
            let (a, b) = (&fin[i], &fin[(i + 1) % n]);
            if a.sing[d] || b.sing[d] {
                continue;
            }
            let delta = b.uv[d] - a.uv[d];
            let wr = wrap(delta, pp);
            if (wr - delta).abs() > 0.5 * pp {
                jumps += 1;
            }
            acc += wr - delta;
        }
        w[d] = (acc / pp).round() as i64;
    }
    (jumps, w)
}
