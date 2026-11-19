//! Cached triangulated mesh geometry, perimeter edges, bounding sphere, and
//! exact 3D ray-vs-triangle picking shared by both GPU and CPU render paths.

use crate::view3d::camera::BoundingSphere;
use crate::widgets::colormap;
use eframe::egui::Color32;
use foiley_core::Vec3;
use foiley_plane::mesh::{Mesh, Panel, PanelKind};

/// Which wing or fuselage body is currently highlighted in Plane Design.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum SelectionHighlight {
    #[default]
    None,
    Wing(usize),
    Body(usize),
}

impl SelectionHighlight {
    pub fn matches_panel(self, p: &Panel) -> bool {
        match self {
            SelectionHighlight::None => false,
            SelectionHighlight::Wing(w) => p.body_id.is_none() && p.wing == w,
            SelectionHighlight::Body(b) => p.body_id == Some(b),
        }
    }
}

/// Per-panel scalar field (Cp or dCp loading) with validated finite display range.
pub struct ColorField<'a> {
    pub values: &'a [f32],
    pub min: f32,
    pub max: f32,
    pub label: String,
    pub revision: FieldRevision,
}

/// Identity of one saved result within a particular revision of the project.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FieldRevision {
    pub project: u64,
    pub revision: u64,
    pub polar: usize,
    pub point: usize,
}

/// Single render/pick triangle retaining its source panel index.
#[derive(Clone, Copy, Debug)]
pub struct MeshTri {
    pub a: Vec3,
    pub b: Vec3,
    pub c: Vec3,
    pub normal: [f32; 3],
    pub panel_idx: usize,
}

/// Chunked AABB for fast ray-picking rejection on large panel/CAD meshes.
#[derive(Clone, Copy, Debug)]
struct TriChunk {
    start: usize,
    end: usize,
    min: Vec3,
    max: Vec3,
}

/// Geometry derived once per `mesh_key`, reused across frames for rendering and picking.
pub struct MeshCache {
    pub key: u64,
    pub panel_count: usize,
    pub bounds: BoundingSphere,
    pub aabb_min: Vec3,
    pub aabb_max: Vec3,
    pub triangles: Vec<MeshTri>,
    /// Pairs of world-space endpoints for `GL_LINES` perimeter edges.
    pub edge_positions: Vec<[f32; 3]>,
    chunks: Vec<TriChunk>,
}

/// Result of a 3D ray pick against the mesh.
#[derive(Clone, Copy, Debug)]
pub struct PickHit {
    pub panel_idx: usize,
    #[allow(dead_code)]
    pub t: f64,
    pub point: Vec3,
}

impl MeshCache {
    pub fn build(key: u64, mesh: &Mesh) -> Self {
        let mut aabb_min = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut aabb_max = Vec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        let mut triangles = Vec::with_capacity(mesh.panels.len() * 2);
        let mut edge_positions = Vec::with_capacity(mesh.panels.len() * 8);

        for (panel_idx, p) in mesh.panels.iter().enumerate() {
            let is_tri = p.is_tri();
            let vcount = if is_tri { 3 } else { 4 };
            for &n in &p.nodes[..vcount] {
                if n.x.is_finite() && n.y.is_finite() && n.z.is_finite() {
                    aabb_min.x = aabb_min.x.min(n.x);
                    aabb_min.y = aabb_min.y.min(n.y);
                    aabb_min.z = aabb_min.z.min(n.z);
                    aabb_max.x = aabb_max.x.max(n.x);
                    aabb_max.y = aabb_max.y.max(n.y);
                    aabb_max.z = aabb_max.z.max(n.z);
                }
            }

            let normal = [p.normal.x as f32, p.normal.y as f32, p.normal.z as f32];
            let [n0, n1, n2, n3] = p.nodes;
            triangles.push(MeshTri {
                a: n0,
                b: n1,
                c: n2,
                normal,
                panel_idx,
            });
            if !is_tri {
                triangles.push(MeshTri {
                    a: n0,
                    b: n2,
                    c: n3,
                    normal,
                    panel_idx,
                });
            }

            let verts = &p.nodes[..vcount];
            for e in 0..vcount {
                let a = verts[e];
                let b = verts[(e + 1) % vcount];
                edge_positions.push([a.x as f32, a.y as f32, a.z as f32]);
                edge_positions.push([b.x as f32, b.y as f32, b.z as f32]);
            }
        }

        let bounds = if aabb_min.x.is_finite() && aabb_max.x.is_finite() {
            let center = (aabb_min + aabb_max) * 0.5;
            let mut r2 = 0.0_f64;
            for p in &mesh.panels {
                let vcount = if p.is_tri() { 3 } else { 4 };
                for &n in &p.nodes[..vcount] {
                    let d = n - center;
                    let d2 = d.norm2();
                    if d2.is_finite() && d2 > r2 {
                        r2 = d2;
                    }
                }
            }
            BoundingSphere {
                center,
                radius: r2.sqrt().max(1e-3),
            }
        } else {
            aabb_min = Vec3::ZERO;
            aabb_max = Vec3::ZERO;
            BoundingSphere::default()
        };

        // Build 64-triangle chunks for fast ray rejection.
        const CHUNK: usize = 64;
        let mut chunks = Vec::with_capacity(triangles.len().div_ceil(CHUNK));
        for start in (0..triangles.len()).step_by(CHUNK) {
            let end = (start + CHUNK).min(triangles.len());
            let mut cmin = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
            let mut cmax = Vec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
            for tri in &triangles[start..end] {
                for v in [tri.a, tri.b, tri.c] {
                    cmin.x = cmin.x.min(v.x);
                    cmin.y = cmin.y.min(v.y);
                    cmin.z = cmin.z.min(v.z);
                    cmax.x = cmax.x.max(v.x);
                    cmax.y = cmax.y.max(v.y);
                    cmax.z = cmax.z.max(v.z);
                }
            }
            chunks.push(TriChunk {
                start,
                end,
                min: cmin,
                max: cmax,
            });
        }

        Self {
            key,
            panel_count: mesh.panels.len(),
            bounds,
            aabb_min,
            aabb_max,
            triangles,
            edge_positions,
            chunks,
        }
    }

    /// Nearest-hit ray-vs-triangle pick clipped to `[t_min, t_max]` (matching the
    /// camera's near/far planes).
    pub fn pick_nearest(
        &self,
        ray_origin: Vec3,
        ray_dir: Vec3,
        t_min: f64,
        t_max: f64,
    ) -> Option<PickHit> {
        if self.triangles.is_empty() || t_max <= t_min {
            return None;
        }
        let inv_dir = Vec3::new(
            inv_safe(ray_dir.x),
            inv_safe(ray_dir.y),
            inv_safe(ray_dir.z),
        );
        if !ray_hits_aabb(
            ray_origin,
            inv_dir,
            self.aabb_min,
            self.aabb_max,
            t_min,
            t_max,
        ) {
            return None;
        }

        let mut best_t = t_max;
        let mut best_panel: Option<usize> = None;

        for ch in &self.chunks {
            if !ray_hits_aabb(ray_origin, inv_dir, ch.min, ch.max, t_min, best_t) {
                continue;
            }
            for tri in &self.triangles[ch.start..ch.end] {
                if let Some(t) =
                    ray_triangle(ray_origin, ray_dir, tri.a, tri.b, tri.c, t_min, best_t)
                {
                    if t < best_t {
                        best_t = t;
                        best_panel = Some(tri.panel_idx);
                    }
                }
            }
        }

        best_panel.map(|panel_idx| PickHit {
            panel_idx,
            t: best_t,
            point: ray_origin + ray_dir * best_t,
        })
    }
}

/// Base RGB color (0..255) for a panel before lighting.
pub fn panel_base_color(
    panel: &Panel,
    panel_idx: usize,
    field: Option<&ColorField<'_>>,
    highlight: SelectionHighlight,
) -> Color32 {
    if let Some(f) = field {
        if let Some(&val) = f.values.get(panel_idx) {
            if val.is_finite() {
                let denom = (f.max - f.min).abs();
                let t = if denom > 1e-9 {
                    (val - f.min) / (f.max - f.min)
                } else {
                    0.5
                };
                return colormap(t);
            } else {
                return Color32::from_rgb(115, 118, 128);
            }
        }
    }

    let base = match panel.kind {
        PanelKind::Fuse => Color32::from_rgb(172, 178, 194),
        PanelKind::Tip => Color32::from_rgb(142, 158, 194),
        _ if panel.flap.is_some() => Color32::from_rgb(218, 162, 88),
        _ => Color32::from_rgb(112, 158, 224),
    };

    if highlight.matches_panel(panel) {
        // Subtle crisp lift toward bright sky-cyan so the active wing/body stands out.
        let blend = |a: u8, b: u8| {
            ((a as f32 * 0.72 + b as f32 * 0.28)
                .round()
                .clamp(0.0, 255.0)) as u8
        };
        Color32::from_rgb(
            blend(base.r(), 135),
            blend(base.g(), 230),
            blend(base.b(), 255),
        )
    } else {
        base
    }
}

/// Perimeter edges for a single panel (used for crisp hover outline without
/// occluding Cp colors).
pub fn panel_perimeter_edges(panel: &Panel) -> Vec<(Vec3, Vec3)> {
    let n = if panel.is_tri() { 3 } else { 4 };
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        out.push((panel.nodes[i], panel.nodes[(i + 1) % n]));
    }
    out
}

fn inv_safe(v: f64) -> f64 {
    if v.abs() < 1e-18 {
        if v.is_sign_negative() {
            -1e18
        } else {
            1e18
        }
    } else {
        1.0 / v
    }
}

fn ray_hits_aabb(o: Vec3, inv_d: Vec3, min: Vec3, max: Vec3, t_min: f64, t_max: f64) -> bool {
    let tx1 = (min.x - o.x) * inv_d.x;
    let tx2 = (max.x - o.x) * inv_d.x;
    let mut tmin = t_min.max(tx1.min(tx2));
    let mut tmax = t_max.min(tx1.max(tx2));

    let ty1 = (min.y - o.y) * inv_d.y;
    let ty2 = (max.y - o.y) * inv_d.y;
    tmin = tmin.max(ty1.min(ty2));
    tmax = tmax.min(ty1.max(ty2));

    let tz1 = (min.z - o.z) * inv_d.z;
    let tz2 = (max.z - o.z) * inv_d.z;
    tmin = tmin.max(tz1.min(tz2));
    tmax = tmax.min(tz1.max(tz2));

    tmax >= tmin
}

/// Two-sided Möller–Trumbore ray-triangle intersection in `[t_min, t_max]`.
fn ray_triangle(
    orig: Vec3,
    dir: Vec3,
    v0: Vec3,
    v1: Vec3,
    v2: Vec3,
    t_min: f64,
    t_max: f64,
) -> Option<f64> {
    let e1 = v1 - v0;
    let e2 = v2 - v0;
    let h = dir.cross(e2);
    let a = e1.dot(h);
    if a.abs() < 1e-12 {
        return None;
    }
    let f = 1.0 / a;
    let s = orig - v0;
    let u = f * s.dot(h);
    if !(-1e-7..=1.0 + 1e-7).contains(&u) {
        return None;
    }
    let q = s.cross(e1);
    let v = f * dir.dot(q);
    if v < -1e-7 || u + v > 1.0 + 1e-7 {
        return None;
    }
    let t = f * e2.dot(q);
    if t >= t_min && t <= t_max {
        Some(t)
    } else {
        None
    }
}
