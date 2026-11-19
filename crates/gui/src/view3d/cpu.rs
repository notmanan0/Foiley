//! CPU `egui::Painter` fallback renderer for headless tests, legacy OpenGL
//! contexts without GLSL 140 / ES 300 or depth buffers, or forced CPU mode.
//!
//! Uses the exact same `Camera`, projection, lighting model, and `MeshCache`
//! as the hardware `glow` renderer, with eye-space near-plane clipping in
//! perspective mode and strict back-to-front chunk ordering.

use crate::view3d::camera::Camera;
use crate::view3d::pick::{panel_base_color, ColorField, SelectionHighlight};
use eframe::egui::{self, Color32, Painter, Pos2, Rect, Shape, Stroke};
use foiley_core::Vec3;
use foiley_plane::mesh::Mesh;

/// Draw reference grid, filled mesh panels, edges, 3D overlay lines, and hover outline
/// purely via `egui::Painter`.
#[allow(clippy::too_many_arguments)]
pub fn paint_cpu(
    painter: &Painter,
    rect: Rect,
    cam: &Camera,
    perspective: bool,
    near: f64,
    mesh: &Mesh,
    field: Option<&ColorField<'_>>,
    highlight: SelectionHighlight,
    show_edges: bool,
    edge_color: Color32,
    grid_lines: &[(Vec3, Vec3, Color32)],
    overlay_lines: &[(Vec3, Vec3, Color32)],
    hovered_panel: Option<usize>,
) {
    // 1) Reference grid (behind geometry)
    for &(a, b, col) in grid_lines {
        if let Some((sa, sb)) = project_segment(cam, a, b, rect, near, perspective) {
            painter.line_segment([sa, sb], Stroke::new(1.0_f32, col));
        }
    }

    let basis = cam.basis();
    let view_dir = -basis.forward;
    let key_dir = Vec3::new(-0.35, 0.45, 0.82).normalized();
    let fill_dir = Vec3::new(0.55, -0.50, 0.30).normalized();

    // Sort panels strictly back-to-front by eye-space centroid depth.
    let mut order: Vec<(f64, usize)> = mesh
        .panels
        .iter()
        .enumerate()
        .map(|(i, p)| (-cam.world_to_eye(p.centre).z, i))
        .collect();
    order.sort_by(|a, b| b.0.total_cmp(&a.0));

    let line_width = 1.0 / painter.ctx().pixels_per_point();
    let edge_stroke = Stroke::new(line_width, edge_color);

    let mut shapes: Vec<Shape> = Vec::with_capacity(order.len() / 2 + 8);
    let mut emesh = egui::Mesh::default();

    for &(_, idx) in &order {
        let p = &mesh.panels[idx];
        let vcount = if p.is_tri() { 3 } else { 4 };
        let eye_verts: Vec<Vec3> = p.nodes[..vcount]
            .iter()
            .map(|&n| cam.world_to_eye(n))
            .collect();
        let clipped = perspective.then(|| clip_poly_near(&eye_verts, near));
        let polygon = clipped.as_deref().unwrap_or(&eye_verts);
        if polygon.len() < 3 {
            continue;
        }
        let pts: Vec<Pos2> = polygon
            .iter()
            .map(|&e| project_eye(cam, e, rect, perspective))
            .collect();
        let area = pts[1..]
            .windows(2)
            .map(|edge| {
                let a = edge[0] - pts[0];
                let b = edge[1] - pts[0];
                a.x * b.y - a.y * b.x
            })
            .sum::<f32>();
        if area.abs() < 1e-4 || !area.is_finite() {
            continue;
        }

        let base = panel_base_color(p, idx, field, highlight);
        let col = shade_color(base, p.normal, view_dir, key_dir, fill_dir);

        let hovered = hovered_panel == Some(idx);
        if show_edges || hovered {
            // Draw filled face then its own perimeter edge in painter order so back
            // edges stay hidden behind nearer panels.
            let mut face_mesh = egui::Mesh::default();
            for &q in &pts {
                face_mesh.colored_vertex(q, col);
            }
            for k in 1..(pts.len() as u32 - 1) {
                face_mesh.add_triangle(0, k, k + 1);
            }
            shapes.push(Shape::mesh(face_mesh));
            let stroke = if hovered {
                Stroke::new(line_width, Color32::from_rgb(255, 240, 120))
            } else {
                edge_stroke
            };
            // Clip the original perimeter. The new edge made by clipping the
            // filled polygon is not a mesh edge.
            for i in 0..vcount {
                let a = eye_verts[i];
                let b = eye_verts[(i + 1) % vcount];
                if let Some(points) = project_eye_segment(cam, a, b, rect, near, perspective) {
                    if points[0].distance_sq(points[1]) > 0.01 * line_width * line_width {
                        shapes.push(Shape::line_segment(points, stroke));
                    }
                }
            }
        } else {
            if emesh.vertices.len() + pts.len() > 60_000 {
                shapes.push(Shape::mesh(std::mem::take(&mut emesh)));
            }
            let i0 = emesh.vertices.len() as u32;
            for &q in &pts {
                emesh.colored_vertex(q, col);
            }
            for k in 1..(pts.len() as u32 - 1) {
                emesh.add_triangle(i0, i0 + k, i0 + k + 1);
            }
        }
    }

    if !emesh.vertices.is_empty() {
        shapes.push(Shape::mesh(emesh));
    }
    painter.extend(shapes);

    // 2) 3D overlay lines (wake, lift vectors)
    for &(a, b, col) in overlay_lines {
        if let Some((sa, sb)) = project_segment(cam, a, b, rect, near, perspective) {
            painter.line_segment([sa, sb], Stroke::new(line_width, col));
        }
    }
}

fn shade_color(
    base: Color32,
    normal: Vec3,
    view_dir: Vec3,
    key_dir: Vec3,
    fill_dir: Vec3,
) -> Color32 {
    let mut n = normal.normalized();
    if n.dot(view_dir) < 0.0 {
        n = -n;
    }
    let diff_key = n.dot(key_dir).max(0.0);
    let diff_fill = n.dot(fill_dir).max(0.0);
    let hemi = 0.5 + 0.5 * n.z;
    let ambient = 0.32 + 0.20 * hemi;
    let shade = (ambient + 0.44 * diff_key + 0.14 * diff_fill) as f32;
    let half_vec = (key_dir + view_dir).normalized();
    let spec = (n.dot(half_vec).max(0.0).powi(28) * 0.08 * 255.0) as f32;
    let ch = |c: u8| ((c as f32 * shade + spec).round().clamp(0.0, 255.0)) as u8;
    Color32::from_rgb(ch(base.r()), ch(base.g()), ch(base.b()))
}

fn project_eye(cam: &Camera, e: Vec3, rect: Rect, perspective: bool) -> Pos2 {
    let w = (rect.width() as f64).max(1.0);
    let h = (rect.height() as f64).max(1.0);
    let aspect = w / h;
    let (ndc_x, ndc_y) = if perspective {
        let depth = (-e.z).max(1e-6);
        let f = 1.0 / (0.5 * cam.fov_y_rad).tan();
        ((e.x * f / aspect) / depth, (e.y * f) / depth)
    } else {
        let half_h = cam.target_half_height().max(1e-6);
        let half_w = half_h * aspect;
        (e.x / half_w, e.y / half_h)
    };
    Pos2::new(
        (rect.left() as f64 + (ndc_x * 0.5 + 0.5) * w) as f32,
        (rect.top() as f64 + (0.5 - ndc_y * 0.5) * h) as f32,
    )
}

/// Sutherland-Hodgman polygon clip against the eye-space near plane `-e.z >= near`.
pub fn clip_poly_near(verts: &[Vec3], near: f64) -> Vec<Vec3> {
    if verts.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(verts.len() + 2);
    let mut prev = *verts.last().unwrap();
    let mut prev_d = -prev.z - near;
    for &cur in verts {
        let cur_d = -cur.z - near;
        let prev_in = prev_d >= 0.0;
        let cur_in = cur_d >= 0.0;
        if prev_in != cur_in {
            let denom = prev_d - cur_d;
            if denom.abs() > 1e-12 {
                let t = (prev_d / denom).clamp(0.0, 1.0);
                out.push(prev + (cur - prev) * t);
            }
        }
        if cur_in {
            out.push(cur);
        }
        prev = cur;
        prev_d = cur_d;
    }
    out
}

/// Project a 3D line segment with exact near-plane clipping in perspective mode.
pub fn project_segment(
    cam: &Camera,
    a: Vec3,
    b: Vec3,
    rect: Rect,
    near: f64,
    perspective: bool,
) -> Option<(Pos2, Pos2)> {
    project_eye_segment(
        cam,
        cam.world_to_eye(a),
        cam.world_to_eye(b),
        rect,
        near,
        perspective,
    )
    .map(|[a, b]| (a, b))
}

fn project_eye_segment(
    cam: &Camera,
    mut ea: Vec3,
    mut eb: Vec3,
    rect: Rect,
    near: f64,
    perspective: bool,
) -> Option<[Pos2; 2]> {
    if perspective {
        let da = -ea.z - near;
        let db = -eb.z - near;
        if da < 0.0 && db < 0.0 {
            return None;
        }
        if da < 0.0 {
            let t = da / (da - db);
            ea = ea + (eb - ea) * t;
        } else if db < 0.0 {
            let t = da / (da - db);
            eb = ea + (eb - ea) * t;
        }
    }
    Some([
        project_eye(cam, ea, rect, perspective),
        project_eye(cam, eb, rect, perspective),
    ])
}
