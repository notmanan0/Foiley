//! Lightweight 3D viewer: orthographic projection, painter's-algorithm
//! filled panels with optional colour mapping, drag to rotate, right-drag
//! to pan, wheel to zoom, double-click to reset.

use crate::widgets::colormap;
use eframe::egui::{self, Color32, Pos2, Rect, Sense, Shape, Stroke, Ui, Vec2};
use foiley_core::Vec3;
use foiley_plane::mesh::{Mesh, PanelKind};

#[derive(Clone)]
pub struct Camera {
    pub yaw: f32,
    pub pitch: f32,
    pub zoom: f32,
    pub pan: Vec2,
}

impl Default for Camera {
    fn default() -> Self {
        Camera { yaw: -2.3, pitch: 0.55, zoom: 1.0, pan: Vec2::ZERO }
    }
}

pub struct Overlay {
    pub lines: Vec<(Vec3, Vec3, Color32)>,
    pub points: Vec<(Vec3, Color32, String)>,
}

/// Per-panel scalar to colour, with its range and a legend label.
pub struct ColorField<'a> {
    pub values: &'a [f32],
    pub min: f32,
    pub max: f32,
    pub label: String,
}

impl Camera {
    fn rotate(&self, p: Vec3) -> (f32, f32, f32) {
        // world: x aft, y right, z up -> view: u right, v up, d depth
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();
        let x = p.x as f32;
        let y = p.y as f32;
        let z = p.z as f32;
        let x1 = cy * x - sy * y;
        let y1 = sy * x + cy * y;
        // pitch about the screen-horizontal axis
        let u = y1;
        let v = cp * z - sp * x1;
        let d = sp * z + cp * x1;
        (u, v, d)
    }
}

/// Draw a mesh; returns the hovered panel index.
pub fn show(ui: &mut Ui, cam: &mut Camera, mesh: Option<&Mesh>, field: Option<&ColorField>, overlay: &Overlay, show_edges: bool) -> Option<usize> {
    let size = ui.available_size().max(Vec2::new(200.0, 200.0));
    let (resp, painter) = ui.allocate_painter(size, Sense::click_and_drag());
    let rect = resp.rect;
    painter.rect_filled(rect, 4.0, ui.visuals().extreme_bg_color);
    if resp.dragged_by(egui::PointerButton::Primary) {
        let d = resp.drag_delta();
        cam.yaw += d.x * 0.008;
        cam.pitch = (cam.pitch + d.y * 0.008).clamp(-1.55, 1.55);
    }
    if resp.dragged_by(egui::PointerButton::Secondary) || resp.dragged_by(egui::PointerButton::Middle) {
        cam.pan += resp.drag_delta();
    }
    if resp.hovered() {
        let s = ui.input(|i| i.smooth_scroll_delta.y);
        if s != 0.0 {
            cam.zoom = (cam.zoom * (1.0 + s * 0.0015)).clamp(0.05, 50.0);
        }
    }
    if resp.double_clicked() {
        *cam = Camera::default();
    }
    let Some(mesh) = mesh else {
        painter.text(rect.center(), egui::Align2::CENTER_CENTER, "building mesh...", egui::FontId::proportional(14.0), ui.visuals().weak_text_color());
        return None;
    };
    // bounding box in view space for auto-fit
    let mut lo = (f32::MAX, f32::MAX);
    let mut hi = (f32::MIN, f32::MIN);
    for p in &mesh.panels {
        for n in &p.nodes {
            let (u, v, _) = cam.rotate(*n);
            lo = (lo.0.min(u), lo.1.min(v));
            hi = (hi.0.max(u), hi.1.max(v));
        }
    }
    let ext = ((hi.0 - lo.0).max(hi.1 - lo.1)).max(1e-6);
    let scale = 0.85 * rect.width().min(rect.height()) / ext * cam.zoom;
    let c = ((lo.0 + hi.0) * 0.5, (lo.1 + hi.1) * 0.5);
    let to_screen = |p: Vec3| -> (Pos2, f32) {
        let (u, v, d) = cam.rotate(p);
        (Pos2::new(rect.center().x + (u - c.0) * scale + cam.pan.x, rect.center().y - (v - c.1) * scale + cam.pan.y), d)
    };
    // light direction in view space
    let light = (0.3f32, 0.5f32, 0.8f32);
    let ln = (light.0 * light.0 + light.1 * light.1 + light.2 * light.2).sqrt();
    let mut order: Vec<(f32, usize)> = mesh.panels.iter().enumerate().map(|(i, p)| (to_screen(p.centre).1, i)).collect();
    order.sort_by(|a, b| b.0.total_cmp(&a.0));
    let mut shapes = Vec::with_capacity(order.len() * 2);
    let mut emesh = egui::Mesh::default();
    let mouse = resp.hover_pos();
    let mut hovered = None;
    let edge = Stroke::new(0.5, Color32::from_black_alpha(90));
    for (_, i) in &order {
        let p = &mesh.panels[*i];
        let pts: Vec<Pos2> = p.nodes.iter().map(|n| to_screen(*n).0).collect();
        let base = match (field, p.kind) {
            (Some(f), _) if *i < f.values.len() => {
                let t = (f.values[*i] - f.min) / (f.max - f.min).max(1e-9);
                colormap(t)
            }
            (_, PanelKind::Fuse) => Color32::from_rgb(170, 170, 185),
            (_, PanelKind::Tip) => Color32::from_rgb(150, 160, 190),
            _ if p.flap.is_some() => Color32::from_rgb(210, 160, 90),
            _ => Color32::from_rgb(120, 160, 220),
        };
        // simple Lambert shading using the view-space normal
        let (nu, nv, nd) = {
            let (a, b, cc) = cam.rotate(p.normal);
            let (o0, o1, o2) = cam.rotate(Vec3::ZERO);
            (a - o0, b - o1, cc - o2)
        };
        let lam = ((nu * light.0 + nv * light.1 - nd * light.2) / ln).abs();
        let k = 0.55 + 0.45 * lam;
        let col = Color32::from_rgb((base.r() as f32 * k) as u8, (base.g() as f32 * k) as u8, (base.b() as f32 * k) as u8);
        let i0 = emesh.vertices.len() as u32;
        for q in &pts {
            emesh.colored_vertex(*q, col);
        }
        emesh.add_triangle(i0, i0 + 1, i0 + 2);
        emesh.add_triangle(i0, i0 + 2, i0 + 3);
        if show_edges && mesh.panels.len() < 20000 {
            shapes.push(Shape::closed_line(pts.clone(), edge));
        }
        if let Some(m) = mouse {
            if point_in_quad(m, &pts) {
                hovered = Some(*i);
            }
        }
        if emesh.vertices.len() > 60000 {
            shapes.push(Shape::mesh(std::mem::take(&mut emesh)));
        }
    }
    // keep painter order: meshes first in chunks, then edges on top
    let mut all = Vec::with_capacity(shapes.len() + 1);
    all.push(Shape::mesh(emesh));
    all.extend(shapes);
    painter.extend(all);
    for (a, b, col) in &overlay.lines {
        painter.line_segment([to_screen(*a).0, to_screen(*b).0], Stroke::new(1.5, *col));
    }
    for (p, col, label) in &overlay.points {
        let s = to_screen(*p).0;
        painter.circle_filled(s, 4.0, *col);
        if !label.is_empty() {
            painter.text(s + Vec2::new(6.0, -6.0), egui::Align2::LEFT_BOTTOM, label, egui::FontId::proportional(12.0), *col);
        }
    }
    // axes triad
    let o = Pos2::new(rect.left() + 40.0, rect.bottom() - 40.0);
    for (d, name, col) in [(Vec3::X, "x", Color32::RED), (Vec3::Y, "y", Color32::GREEN), (Vec3::Z, "z", Color32::from_rgb(80, 120, 255))] {
        let (u, v, _) = cam.rotate(d);
        let e = o + Vec2::new(u, -v) * 25.0;
        painter.line_segment([o, e], Stroke::new(2.0, col));
        painter.text(e, egui::Align2::CENTER_CENTER, name, egui::FontId::proportional(12.0), col);
    }
    if let Some(f) = field {
        // legend bar
        let r = Rect::from_min_size(Pos2::new(rect.right() - 30.0, rect.top() + 30.0), Vec2::new(14.0, 160.0));
        let n = 40;
        for k in 0..n {
            let t0 = k as f32 / n as f32;
            let y0 = r.bottom() - t0 * r.height();
            let y1 = r.bottom() - (k + 1) as f32 / n as f32 * r.height();
            painter.rect_filled(Rect::from_min_max(Pos2::new(r.left(), y1), Pos2::new(r.right(), y0)), 0.0, colormap(t0));
        }
        let tc = ui.visuals().text_color();
        painter.text(Pos2::new(r.left() - 4.0, r.top()), egui::Align2::RIGHT_CENTER, format!("{:.3}", f.max), egui::FontId::proportional(11.0), tc);
        painter.text(Pos2::new(r.left() - 4.0, r.bottom()), egui::Align2::RIGHT_CENTER, format!("{:.3}", f.min), egui::FontId::proportional(11.0), tc);
        painter.text(Pos2::new(r.right(), r.top() - 8.0), egui::Align2::RIGHT_BOTTOM, &f.label, egui::FontId::proportional(12.0), tc);
    }
    painter.text(rect.left_bottom() + Vec2::new(80.0, -8.0), egui::Align2::LEFT_BOTTOM, format!("{} panels - drag: rotate, right-drag: pan, wheel: zoom, double-click: reset", mesh.panels.len()), egui::FontId::proportional(11.0), ui.visuals().weak_text_color());
    hovered
}

fn point_in_quad(p: Pos2, q: &[Pos2]) -> bool {
    let mut inside = false;
    let n = q.len();
    let mut j = n - 1;
    for i in 0..n {
        if ((q[i].y > p.y) != (q[j].y > p.y)) && (p.x < (q[j].x - q[i].x) * (p.y - q[i].y) / (q[j].y - q[i].y + 1e-12) + q[i].x) {
            inside = !inside;
        }
        j = i;
    }
    inside
}
