//! Panel meshes.
//!
//! * Thin meshes put one layer of quads on the camber surface of every
//!   wing (vortex-lattice methods); fuselages are ignored.
//! * Thick meshes wrap the upper and lower surfaces, close the wing tips
//!   (and open roots) with patches, and include fuselages (panel method).
//!
//! Panels of a strip are ordered from the leading edge to the trailing edge
//! for thin meshes, and from the lower trailing edge around the leading edge
//! to the upper trailing edge for thick meshes.  Node order is
//! [LA, TA, TB, LB] (leading/trailing, side A/B); the normal is
//! (TB-LA) x (LB-TA) and points outwards (or up for thin surfaces).
//!
//! Every body in `Plane::fuses` is meshed (its panels carry `body_id`).  A
//! two-sided wing root is trimmed against the body it intersects (bounding-box
//! prefilter, then the section intersection), which gives that wing-body pair
//! its carry-through segment.  Bodies that interpenetrate each other are an
//! error (BVH over the body triangles).

use crate::foildb::{FoilDb, FoilShape};
use crate::geom::{Fuse, Plane, WingKind};
use crate::PlaneError;
use foiley_core::interp::distribution;
use foiley_core::{Spacing, Vec3, DEG};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::Arc;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PanelKind {
    Mid,
    Bot,
    Top,
    Tip,
    Fuse,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Panel {
    /// [LA, TA, TB, LB]
    pub nodes: [Vec3; 4],
    pub normal: Vec3,
    pub area: f64,
    /// collocation point (centroid)
    pub centre: Vec3,
    /// unit tangent vectors of the local frame (l x m = normal)
    pub l: Vec3,
    pub m: Vec3,
    /// characteristic size (max diagonal)
    pub size: f64,
    /// VLM control point (3/4 chord)
    pub ctrl: Vec3,
    /// VLM bound vortex end points (1/4 chord), side A and B
    pub va: Vec3,
    pub vb: Vec3,
    pub kind: PanelKind,
    /// wing index, or usize::MAX for fuselage panels
    pub wing: usize,
    /// body (index into `Plane::fuses`) of a fuselage panel; None for wing panels
    #[serde(default)]
    pub body_id: Option<usize>,
    pub strip: Option<usize>,
    pub trailing: bool,
    pub leading: bool,
    pub flap: Option<usize>,
}

impl Panel {
    fn new(nodes: [Vec3; 4], kind: PanelKind, wing: usize, strip: Option<usize>) -> Panel {
        let [la, ta, tb, lb] = nodes;
        let nn = (tb - la).cross(lb - ta);
        let area = 0.5 * nn.norm();
        let normal = nn.normalized();
        let centre = (la + ta + tb + lb) * 0.25;
        // m: from centre towards the midpoint of side B ("spanwise"), l = m x n
        let mut m = ((lb + tb) * 0.5 - centre).normalized();
        m = (m - normal * normal.dot(m)).normalized();
        let l = m.cross(normal);
        let size = (tb - la).norm().max((lb - ta).norm());
        let va = la + (ta - la) * 0.25;
        let vb = lb + (tb - lb) * 0.25;
        let ctrl = ((la + (ta - la) * 0.75) + (lb + (tb - lb) * 0.75)) * 0.5;
        Panel { nodes, normal, area, centre, l, m, size, ctrl, va, vb, kind, wing, body_id: None, strip, trailing: false, leading: false, flap: None }
    }

    /// A triangular panel (stored with the last node repeated).
    pub fn tri(a: Vec3, b: Vec3, c: Vec3, kind: PanelKind, wing: usize, strip: Option<usize>) -> Panel {
        let nn = (b - a).cross(c - a);
        let area = 0.5 * nn.norm();
        let normal = nn.normalized();
        let centre = (a + b + c) / 3.0;
        let mut m = (b - a).normalized();
        m = (m - normal * normal.dot(m)).normalized();
        let l = m.cross(normal);
        let size = (b - a).norm().max((c - b).norm()).max((a - c).norm());
        Panel { nodes: [a, b, c, c], normal, area, centre, l, m, size, ctrl: centre, va: a, vb: b, kind, wing, body_id: None, strip, trailing: false, leading: false, flap: None }
    }

    pub fn is_tri(&self) -> bool {
        (self.nodes[3] - self.nodes[2]).norm2() == 0.0
    }
}

/// A spanwise strip of a wing (one chordwise row of panels).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Strip {
    pub wing: usize,
    pub left: bool,
    /// surface index (panel between sections `surf` and `surf+1`)
    pub surf: usize,
    /// spanwise fraction of the strip centre within the surface
    pub tau: f64,
    pub foil_a: String,
    pub foil_b: String,
    pub chord: f64,
    pub area: f64,
    pub width: f64,
    /// quarter-chord point
    pub c4: Vec3,
    pub le: Vec3,
    /// trailing-edge end points (sides A and B)
    pub te_a: Vec3,
    pub te_b: Vec3,
    /// mid-surface normal
    pub normal: Vec3,
    /// signed developed span position of the strip centre (negative = left half)
    pub span_pos: f64,
    /// geometric twist, degrees
    pub twist: f64,
    /// panels of the strip (see module docs for the order)
    pub panels: Vec<usize>,
    /// thick meshes: panels touching the trailing edge on the lower and upper side
    #[serde(default)]
    pub te_low: usize,
    #[serde(default)]
    pub te_up: usize,
    /// second triangle of a TE quad (uniform doublets use the mean; linear
    /// doublets take the trace from the triangle that owns the TE segment)
    #[serde(default)]
    pub te_low2: Option<usize>,
    #[serde(default)]
    pub te_up2: Option<usize>,
}

impl Strip {
    /// For thin strips the trailing panel; for thick the lower trailing panel.
    pub fn te_panel_low(&self) -> usize {
        if self.panels.is_empty() {
            0
        } else {
            self.panels[0]
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WingMesh {
    pub name: String,
    pub kind: WingKind,
    pub strips: std::ops::Range<usize>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FlapInfo {
    pub wing: usize,
    /// index of the flap within its wing (left tip to right tip)
    pub index: usize,
    pub name: String,
    pub hinge_a: Vec3,
    pub hinge_b: Vec3,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Mesh {
    pub panels: Vec<Panel>,
    pub strips: Vec<Strip>,
    pub wings: Vec<WingMesh>,
    pub flaps: Vec<FlapInfo>,
    pub thick: bool,
    /// neighbour panels used for surface-gradient reconstruction
    pub neighbours: Vec<Vec<usize>>,
    /// number of wing panels (fuse panels follow)
    pub n_wing_panels: usize,
    /// wing roots trimmed on a fuselage: (left root strip, right root strip);
    /// their circulation is carried across the fuselage in the Trefftz plane
    #[serde(default)]
    pub carry: Vec<(usize, usize)>,
}

/// Options for mesh generation.
#[derive(Clone, Debug, Default)]
pub struct MeshOptions {
    pub thick: bool,
    pub include_fuse: bool,
    /// flap deflections, degrees (+ = trailing edge down), keyed by (wing, flap index)
    pub flaps: HashMap<(usize, usize), f64>,
    /// extra tilt per wing (control polars), degrees
    pub tilts: HashMap<usize, f64>,
    /// use a flat plate for foils that are not loaded instead of failing
    pub allow_missing_foils: bool,
    /// split every panel into triangles (triangle panel method)
    pub triangulate: bool,
    /// display only: add fuselage panels to thin meshes too
    pub fuse_in_thin: bool,
    /// thin meshes: deflect flaps by turning the flap panels' normals about
    /// the hinge instead of moving their nodes (the small-disturbance model
    /// of AVL).  The vortex lattices use it: with the nodes moved, the ring
    /// sides of a deflected flap and of the undeflected wing beside it,
    /// which carry the chordwise running sum of circulation, separate at the
    /// flap's side edge into a close vortex pair whose strength grows with
    /// the lift, so VLM2's flap increments fell with angle of attack.
    pub flap_normals: bool,
}

struct SurfGrid {
    /// chord fractions (LE = 0 .. TE = 1)
    xi: Vec<f64>,
    /// spanwise fractions 0..1 (side A -> side B)
    ts: Vec<f64>,
    /// [k][l] node grids
    mid: Vec<Vec<Vec3>>,
    low: Vec<Vec<Vec3>>,
    up: Vec<Vec<Vec3>>,
    flap_from: Option<usize>,
}

fn rotate_about(p: Vec3, a: Vec3, axis: Vec3, ang: f64) -> Vec3 {
    a + (p - a).rot_axis(axis, ang)
}

/// Build the mesh of a plane.
pub fn build(plane: &Plane, db: &FoilDb, opt: &MeshOptions) -> Result<Mesh, PlaneError> {
    plane.validate().map_err(PlaneError::Geometry)?;
    let mut mesh = Mesh { panels: vec![], strips: vec![], wings: vec![], flaps: vec![], thick: opt.thick, neighbours: vec![], n_wing_panels: 0, carry: vec![] };
    let shape = |name: &str| -> Result<Arc<FoilShape>, PlaneError> {
        match db.shape(name) {
            Some(s) => Ok(s),
            None if opt.allow_missing_foils => Ok(Arc::new(FoilShape::flat_plate())),
            None => Err(PlaneError::MissingFoil(name.to_string())),
        }
    };
    // bodies that trim wing roots (thick analysis meshes) and bodies that are meshed
    let trim_bodies: Vec<(&Fuse, (Vec3, Vec3))> = if opt.thick && opt.include_fuse { plane.fuses.iter().map(|f| (f, body_bbox(f))).collect() } else { vec![] };
    let mesh_bodies: &[Fuse] = if opt.include_fuse && (opt.thick || opt.fuse_in_thin) { &plane.fuses } else { &[] };

    for (iw, wing0) in plane.wings.iter().enumerate() {
        let mut wing = wing0.clone();
        if let Some(t) = opt.tilts.get(&iw) {
            wing.ry += t;
        }
        let first_strip = mesh.strips.len();
        let nsurf = wing.sections.len() - 1;
        let sides: Vec<bool> = if wing.two_sided { vec![true, false] } else { vec![false] };
        let mut flap_counter = 0usize;
        let mut root_trimmed = false;
        for &left in &sides {
            let frames = wing.frames(left);
            // surfaces from tip to root for the left side, root to tip for the right
            let order: Vec<usize> = if left { (0..nsurf).rev().collect() } else { (0..nsurf).collect() };
            for &is in &order {
                let s0 = &wing.sections[is];
                let s1 = &wing.sections[is + 1];
                if (s1.y - s0.y).abs() < 1e-4 {
                    continue;
                }
                let (fa_name, fb_name) = (wing.foil(is, left).to_string(), wing.foil(is + 1, left).to_string());
                let (sha, shb) = (shape(&fa_name)?, shape(&fb_name)?);
                let (fi, fo) = (&frames[is], &frames[is + 1]);
                // chordwise distribution, split at the hinge if both foils have one
                let hinge = match (sha.hinge, shb.hinge) {
                    (Some(a), Some(b)) => Some((0.5 * (a.0 + b.0), a.1, b.1)),
                    _ => None,
                };
                let nx = s0.nx.max(1);
                let (xi, flap_from) = match hinge {
                    Some((xh, _, _)) if nx >= 2 && xh > 0.02 && xh < 0.98 => {
                        let nf = (((1.0 - xh) * nx as f64).round() as usize).clamp(1, nx - 1);
                        let nl = nx - nf;
                        let mut v: Vec<f64> = distribution(s0.xdist, nl).into_iter().map(|t| t * xh).collect();
                        v.pop();
                        v.extend(distribution(s0.xdist, nf).into_iter().map(|t| xh + t * (1.0 - xh)));
                        (v, Some(nl))
                    }
                    _ => (distribution(s0.xdist, nx), None),
                };
                let point = |t: f64, x: f64, which: u8| -> Vec3 {
                    let ea = match which {
                        0 => sha.camber_at(x),
                        1 => sha.lower_at(x),
                        _ => sha.upper_at(x),
                    };
                    let eb = match which {
                        0 => shb.camber_at(x),
                        1 => shb.lower_at(x),
                        _ => shb.upper_at(x),
                    };
                    fi.point(x, ea).lerp(fo.point(x, eb), t)
                };
                // wing root trimmed on the surface of the body it intersects, node by node
                let mut troot = vec![0.0; xi.len()];
                let root_body = if wing.two_sided && is == 0 {
                    // bounding-box prefilter, then the section test on the root nodes
                    trim_bodies
                        .iter()
                        .find(|(f, bb)| {
                            xi.iter().any(|&x| {
                                let p = point(0.0, x, 0);
                                in_box(p, *bb) && inside_body(f, p)
                            })
                        })
                        .map(|(f, _)| *f)
                } else {
                    None
                };
                if let Some(f) = root_body {
                    let inside = |p: Vec3| inside_body(f, p);
                    for (l, &x) in xi.iter().enumerate() {
                        if !inside(point(0.0, x, 0)) {
                            continue;
                        }
                        let (mut lo, mut hi) = (0.0, 0.9);
                        if inside(point(hi, x, 0)) {
                            troot[l] = hi;
                            continue;
                        }
                        for _ in 0..40 {
                            let m = 0.5 * (lo + hi);
                            if inside(point(m, x, 0)) {
                                lo = m;
                            } else {
                                hi = m;
                            }
                        }
                        troot[l] = hi;
                    }
                }
                let t_in = troot.iter().sum::<f64>() / troot.len().max(1) as f64;
                if t_in > 1e-9 {
                    root_trimmed = true;
                }
                let rs: Vec<f64> = distribution(s0.ydist, s0.ny.max(1));
                // orientation: side A = left in plane axes
                let rs_or: Vec<f64> = if left { rs.iter().rev().cloned().collect() } else { rs.clone() };
                let ts: Vec<f64> = rs_or.iter().map(|r| t_in + r * (1.0 - t_in)).collect();
                let grid_of = |which: u8| -> Vec<Vec<Vec3>> { rs_or.iter().map(|&r| xi.iter().enumerate().map(|(l, &x)| point(troot[l] + r * (1.0 - troot[l]), x, which)).collect()).collect() };
                let mut g = SurfGrid { mid: grid_of(0), low: grid_of(1), up: grid_of(2), xi: xi.clone(), ts: ts.clone(), flap_from };
                // close the trailing edge
                let nl = g.xi.len() - 1;
                for k in 0..g.ts.len() {
                    let te = (g.low[k][nl] + g.up[k][nl]) * 0.5;
                    g.low[k][nl] = te;
                    g.up[k][nl] = te;
                }
                // flap
                let mut flap_id = None;
                let mut normal_turn: Option<(Vec3, f64)> = None;
                if let (Some((xh, ya, yb)), Some(_)) = (hinge, flap_from) {
                    let ha = fi.point(xh, sha.lower_at(xh) + ya * (sha.upper_at(xh) - sha.lower_at(xh)));
                    let hb = fo.point(xh, shb.lower_at(xh) + yb * (shb.upper_at(xh) - shb.lower_at(xh)));
                    // hinge points ordered A (left) -> B (right)
                    let (h_a, h_b) = if left { (hb, ha) } else { (ha, hb) };
                    let id = mesh.flaps.len();
                    let idx = flap_counter;
                    flap_counter += 1;
                    let side = if left { "left" } else { "right" };
                    mesh.flaps.push(FlapInfo { wing: iw, index: idx, name: format!("{} flap {} ({side})", wing.name, idx + 1), hinge_a: h_a, hinge_b: h_b });
                    flap_id = Some(id);
                    if let Some(d) = opt.flaps.get(&(iw, idx)) {
                        if d.abs() > 1e-9 && opt.flap_normals && !opt.thick {
                            // vortex-lattice model: the sheet stays where it is and
                            // the flap panels' normals turn (see MeshOptions::flap_normals)
                            normal_turn = Some(((h_b - h_a).normalized(), d * DEG));
                        } else if d.abs() > 1e-9 {
                            let axis = (h_b - h_a).normalized();
                            let ang = d * DEG;
                            let ff = g.flap_from.unwrap();
                            for k in 0..g.ts.len() {
                                let hk = h_a.lerp(h_b, if left { 1.0 - g.ts[k] } else { g.ts[k] });
                                let _ = hk;
                                for l in ff + 1..g.xi.len() {
                                    g.mid[k][l] = rotate_about(g.mid[k][l], h_a, axis, ang);
                                    g.low[k][l] = rotate_about(g.low[k][l], h_a, axis, ang);
                                    g.up[k][l] = rotate_about(g.up[k][l], h_a, axis, ang);
                                }
                            }
                        }
                    }
                }
                let flap_from = g.flap_from;
                let is_flap = move |l: usize| flap_from.map(|f| l >= f).unwrap_or(false);
                // strips
                let ny = g.ts.len() - 1;
                for k in 0..ny {
                    let (t0, t1) = (g.ts[k], g.ts[k + 1]);
                    let tm = 0.5 * (t0 + t1);
                    let lam = |a: Vec3, b: Vec3| a.lerp(b, 0.5);
                    let le = lam(g.mid[k][0], g.mid[k + 1][0]);
                    let te_a = g.mid[k][nl];
                    let te_b = g.mid[k + 1][nl];
                    let te = lam(te_a, te_b);
                    let chord = (te - le).norm();
                    let width = {
                        let d = g.mid[k + 1][0] - g.mid[k][0];
                        // width perpendicular to the chord
                        let cx = (te - le).normalized();
                        (d - cx * d.dot(cx)).norm()
                    };
                    let nrm = (g.mid[k + 1][nl] - g.mid[k][0]).cross(g.mid[k + 1][0] - g.mid[k][nl]).normalized();
                    let span_dev = s0.y + tm * (s1.y - s0.y);
                    let twist = s0.twist + tm * (s1.twist - s0.twist) + wing.ry;
                    let strip_id = mesh.strips.len();
                    let mut st = Strip {
                        wing: iw,
                        left,
                        surf: is,
                        tau: tm,
                        foil_a: fa_name.clone(),
                        foil_b: fb_name.clone(),
                        chord,
                        area: chord * width,
                        width,
                        c4: le.lerp(te, 0.25),
                        le,
                        te_a,
                        te_b,
                        normal: nrm,
                        span_pos: if left { -span_dev } else { span_dev },
                        twist,
                        panels: vec![],
                        te_low: 0,
                        te_up: 0,
                        te_low2: None,
                        te_up2: None,
                    };
                    let push = |mesh: &mut Mesh, nodes: [Vec3; 4], kind: PanelKind, lead: bool, trail: bool, flap: bool| -> usize {
                        let mut p = Panel::new(nodes, kind, iw, Some(strip_id));
                        p.leading = lead;
                        p.trailing = trail;
                        if flap {
                            p.flap = flap_id;
                            if let Some((axis, ang)) = normal_turn {
                                p.normal = p.normal.rot_axis(axis, ang);
                                p.l = p.m.cross(p.normal);
                            }
                        }
                        mesh.panels.push(p);
                        mesh.panels.len() - 1
                    };
                    if !opt.thick {
                        for l in 0..nl {
                            let nodes = [g.mid[k][l], g.mid[k][l + 1], g.mid[k + 1][l + 1], g.mid[k + 1][l]];
                            let id = push(&mut mesh, nodes, PanelKind::Mid, l == 0, l == nl - 1, is_flap(l));
                            st.panels.push(id);
                        }
                    } else {
                        // lower surface TE -> LE (A/B swapped so the normal points down/out)
                        for l in (0..nl).rev() {
                            let nodes = [g.low[k + 1][l], g.low[k + 1][l + 1], g.low[k][l + 1], g.low[k][l]];
                            let id = push(&mut mesh, nodes, PanelKind::Bot, l == 0, l == nl - 1, is_flap(l));
                            st.panels.push(id);
                        }
                        for l in 0..nl {
                            let nodes = [g.up[k][l], g.up[k][l + 1], g.up[k + 1][l + 1], g.up[k + 1][l]];
                            let id = push(&mut mesh, nodes, PanelKind::Top, l == 0, l == nl - 1, is_flap(l));
                            st.panels.push(id);
                        }
                    }
                    st.te_low = st.panels.first().copied().unwrap_or(0);
                    st.te_up = st.panels.last().copied().unwrap_or(0);
                    mesh.strips.push(st);
                }
                // closing patches for thick meshes: tips and open roots
                if opt.thick {
                    let is_tip = is == nsurf - 1;
                    let is_root = is == 0 && (!wing.two_sided || t_in > 0.0 || s0.y > 1e-6);
                    let root_open_to_fuse = t_in > 1e-9;
                    let mut caps: Vec<usize> = Vec::new();
                    // grid column index of the outer / inner side
                    let (outer, inner) = if left { (0, g.ts.len() - 1) } else { (g.ts.len() - 1, 0) };
                    if is_tip {
                        caps.push(outer);
                    }
                    if is_root && !root_open_to_fuse && (wing.close_inner || !wing.two_sided || s0.y > 1e-6) {
                        caps.push(inner);
                    }
                    for c in caps {
                        let ref_n = mesh.strips.last().map(|s| s.normal).unwrap_or(Vec3::Z);
                        let out_dir = {
                            let other = if c == 0 { 1 } else { c - 1 };
                            (g.mid[c][0] - g.mid[other][0]).normalized()
                        };
                        let _ = ref_n;
                        for l in 0..nl {
                            let mut nodes = [g.low[c][l], g.low[c][l + 1], g.up[c][l + 1], g.up[c][l]];
                            let mut p = Panel::new(nodes, PanelKind::Tip, iw, None);
                            if p.area < 1e-14 {
                                continue;
                            }
                            if p.normal.dot(out_dir) < 0.0 {
                                nodes = [g.up[c][l], g.up[c][l + 1], g.low[c][l + 1], g.low[c][l]];
                                p = Panel::new(nodes, PanelKind::Tip, iw, None);
                            }
                            p.flap = if is_flap(l) { flap_id } else { None };
                            mesh.panels.push(p);
                        }
                    }
                }
            }
        }
        if root_trimmed && wing.two_sided {
            // left side strips run tip -> root, right side root -> tip
            let ls = (first_strip..mesh.strips.len()).filter(|&i| mesh.strips[i].left && mesh.strips[i].surf == 0).last();
            let rs = (first_strip..mesh.strips.len()).find(|&i| !mesh.strips[i].left && mesh.strips[i].surf == 0);
            if let (Some(l), Some(r)) = (ls, rs) {
                mesh.carry.push((l, r));
            }
        }
        mesh.wings.push(WingMesh { name: wing.name.clone(), kind: wing.kind, strips: first_strip..mesh.strips.len() });
    }
    mesh.n_wing_panels = mesh.panels.len();

    // fuselages: every body, panels tagged with its index
    let mut body_ranges: Vec<std::ops::Range<usize>> = Vec::with_capacity(mesh_bodies.len());
    for (ib, f) in mesh_bodies.iter().enumerate() {
        let body_start = mesh.panels.len();
        f.validate().map_err(PlaneError::Geometry)?;
        match (f.grid(f.nx, f.nh, false), f.grid(f.nx, f.nh, true)) {
            (Some(gr), Some(gl)) => {
                for g in [&gr, &gl] {
                    for i in 0..g.len() - 1 {
                        // section centre for the outward check
                        let c: Vec3 = (gr[i].iter().chain(gl[i].iter()).fold(Vec3::ZERO, |a, p| a + *p) + gr[i + 1].iter().chain(gl[i + 1].iter()).fold(Vec3::ZERO, |a, p| a + *p)) / (2.0 * (gr[i].len() + gl[i].len()) as f64);
                        for k in 0..g[i].len() - 1 {
                            let (p00, p01, p10, p11) = (g[i][k], g[i][k + 1], g[i + 1][k], g[i + 1][k + 1]);
                            let mut p = Panel::new([p00, p10, p11, p01], PanelKind::Fuse, usize::MAX, None);
                            if p.area < 1e-14 {
                                // collapsed nose/tail: keep as a triangle if possible
                                let t = if (p10 - p11).norm() < 1e-12 { Panel::tri(p00, p10, p01, PanelKind::Fuse, usize::MAX, None) } else { Panel::tri(p00, p11, p01, PanelKind::Fuse, usize::MAX, None) };
                                if t.area < 1e-14 {
                                    continue;
                                }
                                p = t;
                            }
                            if p.normal.dot(p.centre - c) < 0.0 {
                                let [a, b, cc, d] = p.nodes;
                                p = if p.is_tri() { Panel::tri(a, cc, b, PanelKind::Fuse, usize::MAX, None) } else { Panel::new([d, cc, b, a], PanelKind::Fuse, usize::MAX, None) };
                            }
                            mesh.panels.push(p);
                        }
                    }
                }
            }
            _ => {
                for t in f.world_triangles(0, 0) {
                    let p = Panel::tri(t[0], t[1], t[2], PanelKind::Fuse, usize::MAX, None);
                    if p.area > 1e-14 {
                        mesh.panels.push(p);
                    }
                }
            }
        }
        for p in &mut mesh.panels[body_start..] {
            p.body_id = Some(ib);
        }
        body_ranges.push(body_start..mesh.panels.len());
    }
    if opt.thick && body_ranges.len() >= 2 {
        check_bodies(&mesh.panels, &body_ranges, mesh_bodies)?;
    }

    if opt.triangulate {
        triangulate(&mut mesh);
    }

    // neighbours: panels sharing an edge (not across the trailing edge)
    let n = mesh.panels.len();
    let mut nb: Vec<Vec<usize>> = vec![Vec::new(); n];
    {
        use std::collections::HashMap;
        let size = mesh.panels.iter().map(|p| p.size).fold(0.0, f64::max).max(1e-9);
        let tol = 1e-7 * size;
        let key = |p: Vec3| ((p.x / tol).round() as i64, (p.y / tol).round() as i64, (p.z / tol).round() as i64);
        let mut edges: HashMap<((i64, i64, i64), (i64, i64, i64)), Vec<usize>> = HashMap::new();
        for (i, p) in mesh.panels.iter().enumerate() {
            for k in 0..4 {
                let (a, b) = (key(p.nodes[k]), key(p.nodes[(k + 1) % 4]));
                if a == b {
                    continue;
                }
                let e = if a < b { (a, b) } else { (b, a) };
                edges.entry(e).or_default().push(i);
            }
        }
        for list in edges.values() {
            for &i in list {
                for &j in list {
                    if i == j || nb[i].contains(&j) {
                        continue;
                    }
                    let (pi, pj) = (&mesh.panels[i], &mesh.panels[j]);
                    let across_te = pi.trailing && pj.trailing && pi.strip == pj.strip && pi.kind != pj.kind && matches!(pi.kind, PanelKind::Top | PanelKind::Bot) && matches!(pj.kind, PanelKind::Top | PanelKind::Bot);
                    if !across_te {
                        nb[i].push(j);
                    }
                }
            }
        }
    }
    // sorted, so that results do not depend on the HashMap's iteration order
    for list in &mut nb {
        list.sort_unstable();
    }
    mesh.neighbours = nb;
    if mesh.panels.iter().any(|p| !p.area.is_finite() || !p.centre.is_finite()) {
        return Err(PlaneError::Geometry("mesh contains degenerate panels".into()));
    }
    Ok(mesh)
}

/// Split every quadrilateral into two triangles (LA,TA,TB) + (LA,TB,LB);
/// the first one keeps the trailing-edge segment.
fn triangulate(mesh: &mut Mesh) {
    let old = std::mem::take(&mut mesh.panels);
    let mut map: Vec<(usize, Option<usize>)> = Vec::with_capacity(old.len());
    for p in old {
        if p.is_tri() {
            map.push((mesh.panels.len(), None));
            mesh.panels.push(p);
            continue;
        }
        // mirror-symmetric split (left/right and upper/lower): the diagonal
        // goes through the node with the smallest |y| (ties: smallest x, then z)
        let tol = 1e-9 * p.size.max(1e-12);
        let key = |v: Vec3| ((v.y.abs() / tol).round() as i64, (v.x / tol).round() as i64, (v.z / tol).round() as i64);
        let i0 = (0..4).min_by_key(|&i| key(p.nodes[i])).unwrap();
        let n = |k: usize| p.nodes[(i0 + k) % 4];
        let ta = Panel::tri(n(0), n(1), n(2), p.kind, p.wing, p.strip);
        let tb = Panel::tri(n(0), n(2), n(3), p.kind, p.wing, p.strip);
        // the first triangle keeps the trailing-edge segment (nodes 1-2 of the quad)
        let has_te = |k0: usize, k1: usize, k2: usize| [k0, k1, k2].iter().filter(|&&k| (i0 + k) % 4 == 1 || (i0 + k) % 4 == 2).count() == 2;
        let (mut t1, mut t2) = if has_te(0, 1, 2) { (ta, tb) } else { (tb, ta) };
        for t in [&mut t1, &mut t2] {
            t.flap = p.flap;
            t.leading = p.leading;
            t.body_id = p.body_id;
        }
        t1.trailing = p.trailing;
        let i1 = mesh.panels.len();
        let ok1 = t1.area > 1e-14;
        let ok2 = t2.area > 1e-14;
        if ok1 {
            mesh.panels.push(t1);
        }
        let i2 = mesh.panels.len();
        if ok2 {
            mesh.panels.push(t2);
        }
        map.push((if ok1 { i1 } else { i2 }, if ok1 && ok2 { Some(i2) } else { None }));
    }
    for s in &mut mesh.strips {
        let mut np = Vec::with_capacity(s.panels.len() * 2);
        for &q in &s.panels {
            np.push(map[q].0);
            if let Some(b) = map[q].1 {
                np.push(b);
            }
        }
        s.te_low2 = map[s.te_low].1;
        s.te_up2 = map[s.te_up].1;
        s.te_low = map[s.te_low].0;
        s.te_up = map[s.te_up].0;
        s.panels = np;
    }
    mesh.n_wing_panels = mesh.panels.iter().filter(|p| p.kind != PanelKind::Fuse).count();
}

/// Is `p` inside the cross section of body `f` at the plane x of `p`?
fn inside_body(f: &Fuse, p: Vec3) -> bool {
    match f.section_at(p.x) {
        Some((hw, hh, zc)) if hw > 0.0 && hh > 0.0 => (p.y / hw).powi(2) + ((p.z - zc) / hh).powi(2) < 1.0,
        _ => false,
    }
}

/// Conservative bounding box of a body (sampled surface, inflated by 10% of
/// its diagonal so that the spline / NURBS surface between samples is inside).
fn body_bbox(f: &Fuse) -> (Vec3, Vec3) {
    let (lo, hi) = pts_bbox(f.world_triangles(f.nx.max(2), f.nh.max(2)).iter().flatten());
    let pad = 0.1 * (hi - lo).norm() + 1e-9;
    (lo - Vec3::new(pad, pad, pad), hi + Vec3::new(pad, pad, pad))
}

fn pts_bbox<'a>(pts: impl Iterator<Item = &'a Vec3>) -> (Vec3, Vec3) {
    let mut lo = Vec3::new(f64::MAX, f64::MAX, f64::MAX);
    let mut hi = Vec3::new(f64::MIN, f64::MIN, f64::MIN);
    for p in pts {
        lo = Vec3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = Vec3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    }
    (lo, hi)
}

fn in_box(p: Vec3, (lo, hi): (Vec3, Vec3)) -> bool {
    p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y && p.z >= lo.z && p.z <= hi.z
}

fn boxes_overlap(a: (Vec3, Vec3), b: (Vec3, Vec3)) -> bool {
    a.0.x <= b.1.x && b.0.x <= a.1.x && a.0.y <= b.1.y && b.0.y <= a.1.y && a.0.z <= b.1.z && b.0.z <= a.1.z
}

/// Triangles of a panel range (quads split along LA-TB).
fn panel_tris(panels: &[Panel]) -> Vec<[Vec3; 3]> {
    let mut v = Vec::with_capacity(2 * panels.len());
    for p in panels {
        let [a, b, c, d] = p.nodes;
        v.push([a, b, c]);
        if !p.is_tri() {
            v.push([a, c, d]);
        }
    }
    v
}

/// Axis-aligned bounding-volume hierarchy over triangles.
struct Bvh {
    tris: Vec<[Vec3; 3]>,
    /// (box, first, count): count 0 = inner node with children `first` and
    /// `first + 1`, otherwise a leaf holding `tris[first..first + count]`
    nodes: Vec<((Vec3, Vec3), usize, usize)>,
}

impl Bvh {
    fn new(mut tris: Vec<[Vec3; 3]>) -> Bvh {
        let mut nodes = Vec::with_capacity(tris.len() / 2 + 1);
        nodes.push(((Vec3::ZERO, Vec3::ZERO), 0, 0));
        let n = tris.len();
        Self::build(&mut tris, &mut nodes, 0, 0, n);
        Bvh { tris, nodes }
    }

    fn build(tris: &mut [[Vec3; 3]], nodes: &mut Vec<((Vec3, Vec3), usize, usize)>, node: usize, start: usize, end: usize) {
        let bb = pts_bbox(tris[start..end].iter().flatten());
        if end - start <= 8 {
            nodes[node] = (bb, start, end - start);
            return;
        }
        let ext = bb.1 - bb.0;
        let axis = if ext.x >= ext.y && ext.x >= ext.z {
            0
        } else if ext.y >= ext.z {
            1
        } else {
            2
        };
        let key = |t: &[Vec3; 3]| {
            let c = t[0] + t[1] + t[2];
            [c.x, c.y, c.z][axis]
        };
        let mid = (start + end) / 2;
        tris[start..end].select_nth_unstable_by(mid - start, |a, b| key(a).total_cmp(&key(b)));
        let left = nodes.len();
        nodes.push(((Vec3::ZERO, Vec3::ZERO), 0, 0));
        nodes.push(((Vec3::ZERO, Vec3::ZERO), 0, 0));
        nodes[node] = (bb, left, 0);
        Self::build(tris, nodes, left, start, mid);
        Self::build(tris, nodes, left + 1, mid, end);
    }

    /// Does any triangle of the tree intersect `t` (box `tb`)?
    fn intersects(&self, t: &[Vec3; 3], tb: (Vec3, Vec3), stack: &mut Vec<usize>) -> bool {
        stack.clear();
        stack.push(0);
        while let Some(i) = stack.pop() {
            let (bb, first, count) = self.nodes[i];
            if !boxes_overlap(bb, tb) {
                continue;
            }
            if count == 0 {
                stack.push(first);
                stack.push(first + 1);
            } else if self.tris[first..first + count].iter().any(|u| tris_intersect(t, u)) {
                return true;
            }
        }
        false
    }
}

/// Segment p-q against triangle t (Moller-Trumbore); parallel segments do not hit.
fn segment_hits(p: Vec3, q: Vec3, t: &[Vec3; 3]) -> bool {
    let d = q - p;
    let (e1, e2) = (t[1] - t[0], t[2] - t[0]);
    let h = d.cross(e2);
    let det = e1.dot(h);
    if det.abs() <= 1e-14 * e1.norm() * e2.norm() * d.norm() {
        return false;
    }
    let s = p - t[0];
    let u = s.dot(h) / det;
    if !(0.0..=1.0).contains(&u) {
        return false;
    }
    let qv = s.cross(e1);
    let v = d.dot(qv) / det;
    if v < 0.0 || u + v > 1.0 {
        return false;
    }
    let tt = e2.dot(qv) / det;
    (0.0..=1.0).contains(&tt)
}

/// Two triangles intersect when an edge of one crosses the other (coplanar contact is ignored).
fn tris_intersect(a: &[Vec3; 3], b: &[Vec3; 3]) -> bool {
    (0..3).any(|k| segment_hits(a[k], a[(k + 1) % 3], b) || segment_hits(b[k], b[(k + 1) % 3], a))
}

/// Ray-parity test of a point against a (nearly) closed triangulated body.
fn point_in_tris(p: Vec3, tris: &[[Vec3; 3]], reach: f64) -> bool {
    // an oblique direction that does not run along mesh lines or out of an open nose / tail
    let dir = Vec3::new(0.2113, 0.9318, 0.2951).normalized();
    let q = p + dir * reach;
    tris.iter().filter(|t| segment_hits(p, q, t)).count() % 2 == 1
}

/// Body-body proximity and intersection check.  Bodies whose boxes overlap
/// are tested triangle against triangle through a BVH, and for containment
/// by ray parity.  Interpenetrating bodies are an error naming both.
fn check_bodies(panels: &[Panel], ranges: &[std::ops::Range<usize>], bodies: &[Fuse]) -> Result<(), PlaneError> {
    let tris: Vec<Vec<[Vec3; 3]>> = ranges.iter().map(|r| panel_tris(&panels[r.clone()])).collect();
    let boxes: Vec<(Vec3, Vec3)> = tris.iter().map(|t| pts_bbox(t.iter().flatten())).collect();
    let label = |i: usize| format!("'{}' (body {})", bodies[i].name, i + 1);
    let mut stack = Vec::new();
    for j in 1..tris.len() {
        let mut tree: Option<Bvh> = None;
        for i in 0..j {
            if tris[i].is_empty() || tris[j].is_empty() || !boxes_overlap(boxes[i], boxes[j]) {
                continue;
            }
            let tree = tree.get_or_insert_with(|| Bvh::new(tris[j].clone()));
            for t in &tris[i] {
                let tb = pts_bbox(t.iter());
                if boxes_overlap(tb, boxes[j]) && tree.intersects(t, tb, &mut stack) {
                    let c = (t[0] + t[1] + t[2]) / 3.0;
                    return Err(PlaneError::Geometry(format!(
                        "bodies {} and {} interpenetrate: their surfaces intersect near ({:.4}, {:.4}, {:.4}) m; move or resize one of them",
                        label(i),
                        label(j),
                        c.x,
                        c.y,
                        c.z
                    )));
                }
            }
            let reach = 2.0 * ((boxes[i].1 - boxes[i].0).norm() + (boxes[j].1 - boxes[j].0).norm());
            for (inner, outer) in [(i, j), (j, i)] {
                let (bi, bo) = (boxes[inner], boxes[outer]);
                let nested = in_box(bi.0, bo) && in_box(bi.1, bo);
                if nested && point_in_tris(tris[inner][0][0], &tris[outer], reach) {
                    return Err(PlaneError::Geometry(format!("body {} lies inside body {}: bodies must not overlap", label(inner), label(outer))));
                }
            }
        }
    }
    Ok(())
}

impl Mesh {
    pub fn total_wing_area(&self) -> f64 {
        self.strips.iter().map(|s| s.area).sum()
    }
    /// Strip index range of a wing.
    pub fn wing_strips(&self, iw: usize) -> std::ops::Range<usize> {
        self.wings.get(iw).map(|w| w.strips.clone()).unwrap_or(0..0)
    }
    /// Spacing helper re-export for UIs.
    pub fn spacing_labels() -> Vec<(&'static str, Spacing)> {
        Spacing::ALL.iter().map(|s| (s.label(), *s)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Wing, WingSection};
    use foiley_foil::Foil;

    fn db() -> FoilDb {
        let mut d = FoilDb::new();
        d.add_foil(&Foil::naca("0012").unwrap());
        let mut f = Foil::naca("2412").unwrap();
        f.name = "flapped".into();
        f.hinge = Some((0.75, 0.5));
        d.add_foil(&f);
        d
    }

    #[test]
    fn thin_mesh_rect() {
        let w = Wing::new("w", WingKind::Main, vec![WingSection::new(0.0, 0.2, 0.0, "NACA 0012"), WingSection::new(1.0, 0.2, 0.0, "NACA 0012")]);
        let p = Plane { name: "p".into(), description: String::new(), wings: vec![w], fuses: vec![], point_masses: vec![] };
        let m = build(&p, &db(), &MeshOptions::default()).unwrap();
        assert_eq!(m.strips.len(), 22);
        assert_eq!(m.panels.len(), 22 * 13);
        let a: f64 = m.panels.iter().map(|p| p.area).sum();
        assert!((a - 0.4).abs() < 1e-3, "{a}");
        assert!(m.panels.iter().all(|p| p.normal.z > 0.99));
        // strips ordered left to right
        assert!(m.strips[0].span_pos < 0.0 && m.strips[21].span_pos > 0.0);
        assert!((m.total_wing_area() - 0.4).abs() < 1e-3);
    }

    #[test]
    fn thick_mesh_closed() {
        let w = Wing::new("w", WingKind::Main, vec![WingSection::new(0.0, 0.2, 0.0, "NACA 0012"), WingSection::new(1.0, 0.2, 0.0, "NACA 0012")]);
        let p = Plane { name: "p".into(), description: String::new(), wings: vec![w], fuses: vec![], point_masses: vec![] };
        let m = build(&p, &db(), &MeshOptions { thick: true, ..Default::default() }).unwrap();
        // closed surface: sum of n*A ~ 0
        let s = m.panels.iter().fold(Vec3::ZERO, |a, p| a + p.normal * p.area);
        assert!(s.norm() < 1e-6, "{s:?}");
        // volume via divergence theorem ~ area(0012) * c^2 * span
        let v: f64 = m.panels.iter().map(|p| p.centre.dot(p.normal) * p.area).sum::<f64>() / 3.0;
        assert!((v - 0.0822 * 0.04 * 2.0).abs() / v < 0.03, "{v}");
    }

    #[test]
    fn flaps_deflect() {
        let mut w = Wing::new("w", WingKind::Main, vec![WingSection::new(0.0, 0.2, 0.0, "flapped"), WingSection::new(1.0, 0.2, 0.0, "flapped")]);
        w.sections[0].ny = 4;
        let p = Plane { name: "p".into(), description: String::new(), wings: vec![w], fuses: vec![], point_masses: vec![] };
        let mut o = MeshOptions::default();
        let m0 = build(&p, &db(), &o).unwrap();
        assert_eq!(m0.flaps.len(), 2);
        o.flaps.insert((0, 1), 10.0);
        let m1 = build(&p, &db(), &o).unwrap();
        let s = m1.strips.len() - 1;
        assert!(m1.strips[s].te_b.z < m0.strips[s].te_b.z - 0.005, "right flap TE should go down");
        assert!((m1.strips[0].te_a.z - m0.strips[0].te_a.z).abs() < 1e-12, "left flap untouched");
    }

    #[test]
    fn fuse_mesh() {
        let mut p = Plane::default_plane("NACA 0012", "NACA 0012");
        p.fuses.push(crate::geom::Fuse::default_pod());
        let m = build(&p, &db(), &MeshOptions { thick: true, include_fuse: true, ..Default::default() }).unwrap();
        assert!(m.panels.len() > m.n_wing_panels);
        let fuse_out = m.panels[m.n_wing_panels..].iter().fold(Vec3::ZERO, |a, p| a + p.normal * p.area);
        // open nose/tail only: nearly closed
        assert!(fuse_out.norm() < 1e-3, "{fuse_out:?}");
    }
}
