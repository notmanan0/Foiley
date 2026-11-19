//! Interactive 3D aircraft viewport with world-space orbit camera,
//! hardware Z-buffered `eframe::glow` rendering (with automatic CPU fallback),
//! exact 3D ray picking, clickable orientation gizmo, and aero view presets.

pub mod camera;
pub mod cpu;
pub mod gl;
pub mod pick;

pub use camera::{BoundingSphere, Camera, NavPreset, OrbitMode, ViewPrefs, ViewPreset};
pub use gl::{GlFrame, GlRenderer, MeshBufferKey};
pub use pick::{ColorField, FieldRevision, MeshCache, PickHit, SelectionHighlight};

use crate::widgets::colormap;
use eframe::egui::{
    self, Color32, CornerRadius, FontId, Id, Key, Painter, PointerButton, Pos2, Rect, RichText,
    Sense, Stroke, Ui, UiBuilder, Vec2,
};
use eframe::egui_glow::CallbackFn;
use foiley_core::Vec3;
use foiley_plane::mesh::Mesh;
use std::sync::{Arc, Mutex};

/// 3D overlay primitives drawn with the mesh (lines in 3D, labeled markers on top).
#[derive(Clone, Debug, Default)]
pub struct Overlay {
    pub lines: Vec<(Vec3, Vec3, Color32)>,
    pub points: Vec<(Vec3, Color32, String)>,
}

/// Inputs for one 3D viewport frame.
pub struct ViewportInput<'a> {
    /// Index of the active plane, independent of its editable name.
    pub plane_id: usize,
    /// Active mesh and its content key.
    pub mesh: Option<(u64, &'a Mesh)>,
    /// Optional validated per-panel scalar field (Cp / dCp).
    pub field: Option<&'a ColorField<'a>>,
    /// 3D lines (wake, lift vectors) and labeled 3D points (CoG, masses, NP).
    pub overlay: &'a Overlay,
    /// Show panel perimeter edges.
    pub show_edges: bool,
    /// Active wing or fuselage body highlight in Plane Design.
    pub highlight: SelectionHighlight,
    /// Optional zero-lag hover readout formatter for the currently hovered panel index.
    pub hover_label: Option<&'a dyn Fn(usize) -> String>,
    /// Optional non-blocking notice banner (e.g. when analysis polar doesn't match mesh).
    pub notice: Option<&'a str>,
}

/// Pointer picking events reported back to the caller.
#[derive(Clone, Copy, Debug, Default)]
pub struct ViewportEvent {
    #[allow(dead_code)]
    pub hovered_panel: Option<usize>,
    pub clicked_panel: Option<usize>,
}

#[cfg(feature = "snap")]
#[derive(Clone, serde::Serialize)]
pub(crate) struct ViewportSnapshot {
    pub rect: [f32; 4],
    pub target: [f64; 3],
    pub distance: f64,
    pub yaw: f64,
    pub pitch: f64,
    pub orbit_mode: &'static str,
    pub perspective: bool,
    pub grid: bool,
    pub animating: bool,
    pub hovered_panel: Option<usize>,
    pub clicked_panel: Option<usize>,
    pub renderer: &'static str,
    pub mesh_key: Option<u64>,
    pub field: Option<String>,
    pub notice: Option<String>,
    pub overlay_labels: Vec<String>,
}

#[cfg(feature = "snap")]
pub(crate) fn snapshot_id(salt: &str) -> Id {
    Id::new(("foiley_viewport_snapshot", salt))
}

#[cfg(feature = "snap")]
#[allow(clippy::too_many_arguments)]
fn record_snapshot(
    ui: &Ui,
    salt: &str,
    rect: Rect,
    cam: &Camera,
    prefs: &ViewPrefs,
    input: &ViewportInput<'_>,
    event: ViewportEvent,
    renderer: &'static str,
) {
    ui.ctx().data_mut(|data| {
        if data.get_temp::<bool>(Id::new("foiley_snap")) == Some(true) {
            data.insert_temp(
                snapshot_id(salt),
                ViewportSnapshot {
                    rect: [rect.min.x, rect.min.y, rect.max.x, rect.max.y],
                    target: [cam.target.x, cam.target.y, cam.target.z],
                    distance: cam.distance,
                    yaw: cam.yaw,
                    pitch: cam.pitch,
                    orbit_mode: prefs.orbit_mode.as_str(),
                    perspective: prefs.perspective,
                    grid: prefs.show_grid,
                    animating: cam.is_animating(),
                    hovered_panel: event.hovered_panel,
                    clicked_panel: event.clicked_panel,
                    renderer,
                    mesh_key: input.mesh.map(|(key, _)| key),
                    field: input.field.map(|field| field.label.clone()),
                    notice: input.notice.map(str::to_owned),
                    overlay_labels: input
                        .overlay
                        .points
                        .iter()
                        .map(|(_, _, label)| label.clone())
                        .collect(),
                },
            );
        }
    });
}

/// Render the interactive 3D viewport and return hover/click picking results.
#[allow(clippy::too_many_arguments)]
pub fn show(
    ui: &mut Ui,
    id_salt: &str,
    cam: &mut Camera,
    prefs: &mut ViewPrefs,
    cache: &mut Option<MeshCache>,
    gl_renderer: &Arc<Mutex<GlRenderer>>,
    input: &ViewportInput<'_>,
) -> ViewportEvent {
    let available = ui.available_rect_before_wrap();
    let visible = available
        .intersect(ui.clip_rect())
        .intersect(ui.ctx().content_rect());
    let size = (visible.max - available.min).max(Vec2::splat(1.0));
    let (rect, _) = ui.allocate_exact_size(size, Sense::hover());
    let dark = ui.visuals().dark_mode;
    let now_sec = ui.input(|i| i.time);

    // Step active smooth camera animation.
    if cam.step_anim(now_sec) {
        ui.ctx().request_repaint();
    }

    // Rebuild cached triangles/edges/bounding sphere only when mesh_key changes.
    if let Some((key, mesh)) = input.mesh {
        let needs_rebuild = cache.as_ref().map(|c| c.key != key).unwrap_or(true);
        if needs_rebuild {
            *cache = Some(MeshCache::build(key, mesh));
        }
    } else {
        *cache = None;
    }

    // Track fit intent independently of mesh edit revisions: auto-fit on initial
    // valid mesh or actual plane change, never on minor geometry tweaks.
    let bounds = cache.as_ref().map(|c| c.bounds);
    if let Some(sphere) = bounds {
        let plane_changed = cam.fitted_plane != Some(input.plane_id);
        if cam.needs_fit || plane_changed {
            cam.fit_sphere(sphere, rect.size(), prefs.perspective, None);
            cam.fitted_plane = Some(input.plane_id);
        }
    }

    // Paint subtle slate gradient backdrop.
    let painter = ui.painter_at(rect);
    draw_backdrop(&painter, rect, dark);

    let Some((mesh_key, mesh)) = input.mesh else {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "building mesh...",
            FontId::proportional(14.0),
            ui.visuals().weak_text_color(),
        );
        #[cfg(feature = "snap")]
        record_snapshot(
            ui,
            id_salt,
            rect,
            cam,
            prefs,
            input,
            ViewportEvent::default(),
            "cpu",
        );
        return ViewportEvent::default();
    };
    let Some(mesh_cache) = cache.as_ref() else {
        return ViewportEvent::default();
    };

    // Canvas interaction area (registered before floating chrome so HUD/gizmo win
    // pointer priority when overlapping).
    let canvas_id = ui.make_persistent_id(Id::new("view3d_canvas").with(id_salt));
    let canvas_resp = ui.interact(rect, canvas_id, Sense::click_and_drag());

    // Use exact HUD rect recorded from the previous frame (with sensible fallback on frame 1).
    let hud_rect_id = ui.make_persistent_id(Id::new("view3d_hud_rect").with(id_salt));
    let hud_rect = ui
        .ctx()
        .memory(|m| m.data.get_temp::<Rect>(hud_rect_id))
        .unwrap_or_else(|| {
            Rect::from_min_max(
                Pos2::new(rect.right() - 360.0, rect.top() + 4.0),
                Pos2::new(rect.right() - 4.0, rect.top() + 38.0),
            )
        });
    let gizmo_center = Pos2::new(rect.left() + 50.0, rect.bottom() - 62.0);
    let gizmo_radius = 30.0_f32;
    let gizmo_rect = Rect::from_center_size(gizmo_center, Vec2::splat(gizmo_radius * 2.0 + 22.0));

    let hover_pos = canvas_resp.hover_pos();
    let over_chrome = !canvas_resp.hovered()
        || hover_pos.is_some_and(|p| hud_rect.contains(p) || gizmo_rect.contains(p));

    // Compute dynamic clip planes for both rendering and picking.
    let grid_spec = if prefs.show_grid {
        Some(compute_grid_spec(mesh_cache.bounds))
    } else {
        None
    };
    let extra_extent = compute_extra_extent(mesh_cache.bounds, grid_spec, &input.overlay.lines);
    let (near, far) = cam.clip_planes(Some(mesh_cache.bounds), extra_extent, prefs.perspective);

    // Exact 3D ray picking when pointer is over the canvas and not dragging or over HUD/gizmo.
    let is_dragging = canvas_resp.dragged();
    let pick_hit: Option<PickHit> = if !is_dragging && !over_chrome {
        hover_pos.and_then(|pos| pick_at(mesh_cache, cam, rect, prefs.perspective, pos, near, far))
    } else {
        None
    };

    // Handle pointer drags (orbit / pan / Ctrl+MMB zoom).
    let modifiers = ui.input(|i| i.modifiers);
    if canvas_resp.dragged() {
        let delta = canvas_resp.drag_delta();
        let primary = canvas_resp.dragged_by(PointerButton::Primary);
        let middle = canvas_resp.dragged_by(PointerButton::Middle);
        let secondary = canvas_resp.dragged_by(PointerButton::Secondary);

        if middle && (modifiers.ctrl || modifiers.command) {
            if delta.y != 0.0 {
                let factor = (f64::from(delta.y) * 0.008).exp();
                cam.zoom_anchored(
                    factor,
                    rect.center(),
                    rect,
                    prefs.perspective,
                    cam.target,
                    mesh_cache.bounds.radius,
                );
            }
        } else if secondary
            || (middle && modifiers.shift)
            || (prefs.nav_preset == NavPreset::Standard && primary && modifiers.shift)
        {
            cam.pan(delta, rect.size());
        } else if middle
            || (prefs.nav_preset == NavPreset::Standard
                && primary
                && !modifiers.shift
                && !modifiers.ctrl
                && !modifiers.command)
        {
            cam.orbit(delta, prefs.orbit_mode);
        }
    }

    // Handle wheel & pinch zoom (deduplicating Ctrl+scroll so it never double-fires).
    if !over_chrome && hover_pos.is_some() && canvas_resp.hovered() {
        let (zoom_delta, scroll_y) = ui.input(|i| (i.zoom_delta(), i.smooth_scroll_delta.y));
        let zoom_factor = if (zoom_delta - 1.0).abs() > 1e-4 {
            Some(1.0 / zoom_delta as f64)
        } else if scroll_y.abs() > 1e-3 && !modifiers.ctrl && !modifiers.command {
            Some((-scroll_y as f64 * 0.0022).exp())
        } else {
            None
        };
        if let Some(factor) = zoom_factor {
            let cursor = hover_pos.unwrap_or_else(|| rect.center());
            let (ro, rd) = cam.screen_ray(cursor, rect, prefs.perspective);
            let basis = cam.basis();
            let anchor = pick_hit
                .map(|h| h.point)
                .or_else(|| camera::intersect_plane(ro, rd, cam.target, basis.forward))
                .unwrap_or(cam.target);
            cam.zoom_anchored(
                factor,
                cursor,
                rect,
                prefs.perspective,
                anchor,
                mesh_cache.bounds.radius,
            );
        }
    }

    // Double-click MMB (or double-click LMB on empty space) fits the aircraft.
    if !over_chrome
        && (canvas_resp.double_clicked_by(PointerButton::Middle)
            || (canvas_resp.double_clicked_by(PointerButton::Primary) && pick_hit.is_none()))
    {
        cam.fit_sphere(
            mesh_cache.bounds,
            rect.size(),
            prefs.perspective,
            Some(now_sec),
        );
        ui.ctx().request_repaint();
    }

    // Keyboard shortcuts when viewport is hovered or focused and no text input owns keyboard.
    if canvas_resp.has_focus() || (canvas_resp.hovered() && !ui.ctx().wants_keyboard_input()) {
        let events = ui.input(|i| i.events.clone());
        for ev in events {
            if let egui::Event::Key {
                mut key,
                physical_key,
                pressed: true,
                repeat,
                modifiers: km,
            } = ev
            {
                let alt_only = km.alt;
                if alt_only {
                    continue;
                }
                let opp = km.shift || km.ctrl || km.command;
                if km.shift {
                    if let Some(digit @ (Key::Num1 | Key::Num3 | Key::Num7)) = physical_key {
                        key = digit;
                    }
                }
                match key {
                    Key::Num1 => {
                        let p = if opp {
                            ViewPreset::Rear
                        } else {
                            ViewPreset::Front
                        };
                        cam.set_preset(p, None, Some(now_sec));
                        ui.ctx().request_repaint();
                    }
                    Key::Num3 => {
                        let p = if opp {
                            ViewPreset::Left
                        } else {
                            ViewPreset::Right
                        };
                        cam.set_preset(p, None, Some(now_sec));
                        ui.ctx().request_repaint();
                    }
                    Key::Num7 => {
                        let p = if opp {
                            ViewPreset::Bottom
                        } else {
                            ViewPreset::Top
                        };
                        cam.set_preset(p, None, Some(now_sec));
                        ui.ctx().request_repaint();
                    }
                    Key::Num0 | Key::R if !opp && !repeat => {
                        cam.set_preset(
                            ViewPreset::Iso,
                            Some((mesh_cache.bounds, rect.size(), prefs.perspective)),
                            Some(now_sec),
                        );
                        ui.ctx().request_repaint();
                    }
                    Key::Num5 if !opp && !repeat => {
                        prefs.perspective = !prefs.perspective;
                    }
                    Key::F if !opp && !repeat => {
                        cam.fit_sphere(
                            mesh_cache.bounds,
                            rect.size(),
                            prefs.perspective,
                            Some(now_sec),
                        );
                        ui.ctx().request_repaint();
                    }
                    Key::G if !opp && !repeat => {
                        prefs.show_grid = !prefs.show_grid;
                    }
                    _ => {}
                }
            }
        }
    }

    // Recompute clip planes after any camera updates this frame.
    let (near, far) = cam.clip_planes(Some(mesh_cache.bounds), extra_extent, prefs.perspective);
    let aspect = (rect.width() as f64 / (rect.height() as f64).max(1.0)).max(1e-4);

    // Build grid line geometry if enabled.
    let grid_segments = grid_spec
        .map(|g| build_grid_lines(g, dark))
        .unwrap_or_default();

    // Render 3D scene via hardware glow Z-buffer (or automatic CPU fallback).
    let hovered_panel = pick_hit.map(|h| h.panel_idx);
    let edge_color32 = if dark {
        Color32::from_rgba_unmultiplied(12, 16, 24, 145)
    } else {
        Color32::from_rgba_unmultiplied(20, 26, 36, 130)
    };

    let use_gl = {
        let g = gl_renderer.lock().unwrap();
        g.can_use_gl()
    };

    if use_gl {
        let field_rev = input.field.map(|f| f.revision);
        let mesh_buf_key = MeshBufferKey {
            mesh_key,
            field_rev,
            selection: input.highlight,
        };

        let (need_tri_upload, need_edge_upload) = {
            let g = gl_renderer.lock().unwrap();
            (
                !g.has_mesh_buffer(mesh_buf_key),
                !g.has_edge_buffer(mesh_key),
            )
        };

        let tri_upload = if need_tri_upload {
            Some(build_tri_vertex_buffer(
                mesh,
                mesh_cache,
                input.field,
                input.highlight,
            ))
        } else {
            None
        };
        let edge_upload = if need_edge_upload {
            Some(flatten_vec3(&mesh_cache.edge_positions))
        } else {
            None
        };

        let basis = cam.basis();
        let view_dir = [
            (-basis.forward.x) as f32,
            (-basis.forward.y) as f32,
            (-basis.forward.z) as f32,
        ];
        let view_proj = cam.view_proj_matrix(aspect, near, far, prefs.perspective);

        let mut hover_lines = Vec::new();
        if let Some(hi) = hovered_panel {
            if let Some(p) = mesh.panels.get(hi) {
                let hc = rgba_f32(Color32::from_rgb(255, 240, 120));
                for (a, b) in pick::panel_perimeter_edges(p) {
                    push_colored_segment(&mut hover_lines, a, b, hc);
                }
            }
        }

        let gl_frame = GlFrame {
            view_proj,
            view_dir,
            mesh_buf_key,
            tri_upload,
            tri_vertex_count: (mesh_cache.triangles.len() * 3) as i32,
            edge_mesh_key: mesh_key,
            edge_upload,
            edge_vertex_count: mesh_cache.edge_positions.len() as i32,
            show_edges: input.show_edges,
            edge_color: rgba_f32(edge_color32),
            grid_lines: pack_colored_segments(&grid_segments),
            overlay_lines: pack_colored_segments(&input.overlay.lines),
            hover_lines,
        };

        let renderer_arc = Arc::clone(gl_renderer);
        let ctx = ui.ctx().clone();
        let cb = egui::PaintCallback {
            rect,
            callback: Arc::new(CallbackFn::new(move |info, painter| {
                let mut guard = renderer_arc.lock().unwrap();
                if !guard.paint(painter, &info, &gl_frame) {
                    ctx.request_repaint();
                }
            })),
        };
        painter.add(cb);
    } else {
        cpu::paint_cpu(
            &painter,
            rect,
            cam,
            prefs.perspective,
            near,
            mesh,
            input.field,
            input.highlight,
            input.show_edges,
            edge_color32,
            &grid_segments,
            &input.overlay.lines,
            hovered_panel,
        );
    }

    // 2D egui overlays on top of the 3D scene (CoG, ballast, NP pills).
    draw_overlay_points(
        &painter,
        rect,
        cam,
        near,
        prefs.perspective,
        &input.overlay.points,
    );

    // Optional Cp / dCp legend bar.
    if let Some(f) = input.field {
        draw_legend(&painter, rect, f, ui.visuals().text_color(), dark);
    }

    // Optional notice banner (top-left) when results don't match active mesh.
    if let Some(notice) = input.notice {
        draw_pill_badge(
            &painter,
            rect.left_top() + Vec2::new(10.0, 43.0),
            egui::Align2::LEFT_TOP,
            notice,
            Color32::from_rgb(245, 195, 95),
            Color32::from_rgba_unmultiplied(45, 32, 12, 210),
        );
    }

    // Optional hover readout pill (bottom-right above hint).
    if let (Some(hi), Some(fmt)) = (hovered_panel, input.hover_label) {
        let readout = fmt(hi);
        if !readout.is_empty() {
            draw_pill_badge(
                &painter,
                Pos2::new(rect.right() - 10.0, rect.bottom() - 28.0),
                egui::Align2::RIGHT_BOTTOM,
                &readout,
                Color32::from_rgb(235, 242, 252),
                Color32::from_rgba_unmultiplied(16, 22, 32, 220),
            );
        }
    }

    // Interactive 3D Orientation Gizmo (bottom-left).
    let gizmo_clicked =
        draw_orientation_gizmo(ui, &painter, gizmo_center, gizmo_radius, cam, id_salt);
    if let Some(preset) = gizmo_clicked {
        cam.set_preset(preset, None, Some(now_sec));
        ui.ctx().request_repaint();
    }

    // Floating compact HUD toolbar (top-right).
    let painted_hud = draw_hud_toolbar(ui, rect, cam, prefs, mesh_cache.bounds, now_sec, id_salt);
    ui.ctx()
        .memory_mut(|m| m.data.insert_temp(hud_rect_id, painted_hud.expand(2.0)));

    // Bottom status bar (panel count + grid spacing + active navigation hint).
    let mut status = format!("{} panels", mesh_cache.panel_count);
    if let Some(g) = grid_spec {
        status.push_str(&format!(" | grid {:.2} m", g.step));
    }
    if !use_gl {
        status.push_str(" | CPU");
    }
    if rect.width() >= 660.0 {
        status.push_str("  -  ");
        status.push_str(prefs.nav_preset.hint());
    }
    painter.text(
        rect.left_bottom() + Vec2::new(10.0, -7.0),
        egui::Align2::LEFT_BOTTOM,
        status,
        FontId::proportional(11.0),
        if dark {
            Color32::from_rgb(155, 164, 178)
        } else {
            Color32::from_rgb(75, 82, 94)
        },
    );

    // Click-to-select fires only on clean LMB click over canvas (not HUD/gizmo).
    let clicked_panel = if !over_chrome
        && gizmo_clicked.is_none()
        && canvas_resp.clicked_by(PointerButton::Primary)
        && !modifiers.shift
        && !modifiers.ctrl
        && !modifiers.command
    {
        hovered_panel
    } else {
        None
    };

    let event = ViewportEvent {
        hovered_panel,
        clicked_panel,
    };
    #[cfg(feature = "snap")]
    record_snapshot(
        ui,
        id_salt,
        rect,
        cam,
        prefs,
        input,
        event,
        if use_gl { "gl" } else { "cpu" },
    );
    event
}

fn pick_at(
    cache: &MeshCache,
    cam: &Camera,
    rect: Rect,
    perspective: bool,
    pos: Pos2,
    near: f64,
    far: f64,
) -> Option<PickHit> {
    let (origin, direction) = cam.screen_ray(pos, rect, perspective);
    let depth_per_unit = direction.dot(cam.basis().forward);
    cache.pick_nearest(
        origin,
        direction,
        near / depth_per_unit,
        far / depth_per_unit,
    )
}

fn draw_backdrop(painter: &Painter, rect: Rect, dark: bool) {
    let (top, bot) = if dark {
        (Color32::from_rgb(18, 22, 30), Color32::from_rgb(10, 13, 19))
    } else {
        (
            Color32::from_rgb(232, 237, 244),
            Color32::from_rgb(214, 221, 232),
        )
    };
    let mut m = egui::Mesh::default();
    m.colored_vertex(rect.left_top(), top);
    m.colored_vertex(rect.right_top(), top);
    m.colored_vertex(rect.right_bottom(), bot);
    m.colored_vertex(rect.left_bottom(), bot);
    m.add_triangle(0, 1, 2);
    m.add_triangle(0, 2, 3);
    painter.add(egui::Shape::mesh(m));
    painter.rect_stroke(
        rect,
        CornerRadius::same(4),
        Stroke::new(
            1.0,
            if dark {
                Color32::from_rgb(42, 49, 62)
            } else {
                Color32::from_rgb(180, 190, 204)
            },
        ),
        egui::StrokeKind::Inside,
    );
}

#[derive(Clone, Copy)]
struct GridSpec {
    step: f64,
    half_steps: i32,
    origin_x: f64,
    z: f64,
}

fn compute_grid_spec(bounds: BoundingSphere) -> GridSpec {
    let target_step = (bounds.radius * 0.35).max(0.05);
    let exp = target_step.log10().floor();
    let base = 10.0_f64.powf(exp);
    let norm = target_step / base;
    let step = if norm < 1.5 {
        base
    } else if norm < 3.5 {
        2.0 * base
    } else if norm < 7.5 {
        5.0 * base
    } else {
        10.0 * base
    };
    let half_steps = ((bounds.radius * 2.2 / step).ceil() as i32).clamp(4, 14);
    let origin_x = (bounds.center.x / step).round() * step;
    GridSpec {
        step,
        half_steps,
        origin_x,
        z: 0.0,
    }
}

fn compute_extra_extent(
    bounds: BoundingSphere,
    grid: Option<GridSpec>,
    overlay_lines: &[(Vec3, Vec3, Color32)],
) -> f64 {
    let mut extra = 0.0_f64;
    if let Some(g) = grid {
        let grid_r = g.step * g.half_steps as f64 * 1.45;
        extra = extra.max((grid_r - bounds.radius).max(0.0));
    }
    // Include overlay extents (e.g. wakes) in clip planes so they aren't cut off,
    // without polluting Fit's bounding sphere.
    for &(a, b, _) in overlay_lines {
        for p in [a, b] {
            let d = (p - bounds.center).norm();
            if d.is_finite() {
                extra = extra.max((d - bounds.radius).max(0.0));
            }
        }
    }
    extra
}

fn build_grid_lines(g: GridSpec, dark: bool) -> Vec<(Vec3, Vec3, Color32)> {
    let n = g.half_steps;
    let min_x = g.origin_x - n as f64 * g.step;
    let max_x = g.origin_x + n as f64 * g.step;
    let min_y = -(n as f64) * g.step;
    let max_y = (n as f64) * g.step;

    let minor = if dark {
        Color32::from_rgba_unmultiplied(120, 138, 165, 36)
    } else {
        Color32::from_rgba_unmultiplied(70, 85, 110, 40)
    };
    let major = if dark {
        Color32::from_rgba_unmultiplied(145, 165, 195, 68)
    } else {
        Color32::from_rgba_unmultiplied(60, 75, 100, 72)
    };
    let x_axis = Color32::from_rgba_unmultiplied(215, 85, 85, 110);
    let y_axis = Color32::from_rgba_unmultiplied(85, 195, 115, 110);

    let mut out = Vec::with_capacity((2 * n as usize + 1) * 2);
    for i in -n..=n {
        let x = g.origin_x + i as f64 * g.step;
        let is_zero = x.abs() < g.step * 1e-4;
        let col = if is_zero {
            y_axis
        } else if i % 5 == 0 {
            major
        } else {
            minor
        };
        out.push((Vec3::new(x, min_y, g.z), Vec3::new(x, max_y, g.z), col));

        let y = i as f64 * g.step;
        let col_y = if i == 0 {
            x_axis
        } else if i % 5 == 0 {
            major
        } else {
            minor
        };
        out.push((Vec3::new(min_x, y, g.z), Vec3::new(max_x, y, g.z), col_y));
    }
    out
}

fn draw_overlay_points(
    painter: &Painter,
    rect: Rect,
    cam: &Camera,
    near: f64,
    perspective: bool,
    points: &[(Vec3, Color32, String)],
) {
    for (p, col, label) in points {
        let Some((s, _)) = cam.project_point(*p, rect, near, perspective) else {
            continue;
        };
        if !rect.expand(40.0).contains(s) {
            continue;
        }
        painter.circle_filled(s, 5.2, Color32::from_rgba_unmultiplied(10, 14, 20, 210));
        painter.circle_filled(s, 3.6, *col);
        if !label.is_empty() {
            draw_pill_badge(
                painter,
                s + Vec2::new(8.0, -6.0),
                egui::Align2::LEFT_BOTTOM,
                label,
                *col,
                Color32::from_rgba_unmultiplied(14, 18, 26, 205),
            );
        }
    }
}

fn draw_pill_badge(
    painter: &Painter,
    pos: Pos2,
    align: egui::Align2,
    text: &str,
    text_col: Color32,
    bg: Color32,
) {
    let font = FontId::proportional(11.5);
    let pad = Vec2::new(6.0, 2.5);
    let width = (painter.clip_rect().width() - pad.x * 2.0 - 24.0).max(32.0);
    let galley = painter.layout(text.to_string(), font, text_col, width);
    let pill_size = galley.size() + pad * 2.0;
    let text_rect = align.anchor_size(pos, pill_size);
    painter.rect_filled(text_rect, CornerRadius::same(4), bg);
    painter.rect_stroke(
        text_rect,
        CornerRadius::same(4),
        Stroke::new(
            0.8,
            Color32::from_rgba_unmultiplied(text_col.r(), text_col.g(), text_col.b(), 110),
        ),
        egui::StrokeKind::Inside,
    );
    painter.galley(text_rect.min + pad, galley, text_col);
}

fn draw_legend(painter: &Painter, rect: Rect, f: &ColorField<'_>, text_col: Color32, dark: bool) {
    let r = Rect::from_min_size(
        Pos2::new(rect.right() - 28.0, rect.top() + 52.0),
        Vec2::new(14.0, 150.0),
    );
    let bg = r
        .expand2(Vec2::new(56.0, 20.0))
        .translate(Vec2::new(-24.0, 0.0));
    painter.rect_filled(
        bg,
        CornerRadius::same(5),
        if dark {
            Color32::from_rgba_unmultiplied(14, 18, 26, 185)
        } else {
            Color32::from_rgba_unmultiplied(244, 247, 252, 205)
        },
    );
    let n = 40;
    for k in 0..n {
        let t0 = k as f32 / n as f32;
        let y0 = r.bottom() - t0 * r.height();
        let y1 = r.bottom() - (k + 1) as f32 / n as f32 * r.height();
        painter.rect_filled(
            Rect::from_min_max(Pos2::new(r.left(), y1), Pos2::new(r.right(), y0)),
            0.0,
            colormap(t0),
        );
    }
    painter.rect_stroke(
        r,
        0.0,
        Stroke::new(0.8, Color32::from_white_alpha(60)),
        egui::StrokeKind::Inside,
    );
    painter.text(
        Pos2::new(r.left() - 5.0, r.top()),
        egui::Align2::RIGHT_CENTER,
        format!("{:.3}", f.max),
        FontId::proportional(11.0),
        text_col,
    );
    painter.text(
        Pos2::new(r.left() - 5.0, r.center().y),
        egui::Align2::RIGHT_CENTER,
        format!("{:.3}", 0.5 * (f.min + f.max)),
        FontId::proportional(10.5),
        Color32::from_rgba_unmultiplied(text_col.r(), text_col.g(), text_col.b(), 170),
    );
    painter.text(
        Pos2::new(r.left() - 5.0, r.bottom()),
        egui::Align2::RIGHT_CENTER,
        format!("{:.3}", f.min),
        FontId::proportional(11.0),
        text_col,
    );
    painter.text(
        Pos2::new(r.right(), r.top() - 6.0),
        egui::Align2::RIGHT_BOTTOM,
        &f.label,
        FontId::proportional(11.5),
        text_col,
    );
}

fn draw_orientation_gizmo(
    ui: &mut Ui,
    painter: &Painter,
    center: Pos2,
    radius: f32,
    cam: &Camera,
    id_salt: &str,
) -> Option<ViewPreset> {
    let gizmo_rect = Rect::from_center_size(center, Vec2::splat(radius * 2.0 + 22.0));
    let id = ui.make_persistent_id(Id::new("view3d_gizmo").with(id_salt));
    let resp = ui.interact(gizmo_rect, id, Sense::click());
    let hovered_disc = resp.hovered();

    // Subtle circular backing disc.
    painter.circle_filled(
        center,
        radius + 9.0,
        Color32::from_rgba_unmultiplied(20, 26, 36, if hovered_disc { 175 } else { 110 }),
    );
    if hovered_disc {
        painter.circle_stroke(
            center,
            radius + 9.0,
            Stroke::new(1.0, Color32::from_white_alpha(45)),
        );
    }

    let b = cam.basis();
    struct AxisNode {
        dir: Vec3,
        label: &'static str,
        pos_axis: bool,
        col: Color32,
        preset: ViewPreset,
        screen: Pos2,
        depth: f64,
    }

    let defs = [
        (
            Vec3::X,
            "+X",
            true,
            Color32::from_rgb(232, 86, 86),
            ViewPreset::Rear,
        ),
        (
            -Vec3::X,
            "-X",
            false,
            Color32::from_rgb(175, 78, 78),
            ViewPreset::Front,
        ),
        (
            Vec3::Y,
            "+Y",
            true,
            Color32::from_rgb(86, 206, 114),
            ViewPreset::Right,
        ),
        (
            -Vec3::Y,
            "-Y",
            false,
            Color32::from_rgb(76, 155, 96),
            ViewPreset::Left,
        ),
        (
            Vec3::Z,
            "+Z",
            true,
            Color32::from_rgb(92, 148, 245),
            ViewPreset::Top,
        ),
        (
            -Vec3::Z,
            "-Z",
            false,
            Color32::from_rgb(76, 112, 185),
            ViewPreset::Bottom,
        ),
    ];

    let mut nodes: Vec<AxisNode> = defs
        .into_iter()
        .map(|(dir, label, pos_axis, col, preset)| {
            let sx = dir.dot(b.right) as f32;
            let sy = dir.dot(b.up) as f32;
            let depth = -dir.dot(b.forward); // larger = closer to viewer
            AxisNode {
                dir,
                label,
                pos_axis,
                col,
                preset,
                screen: center + Vec2::new(sx, -sy) * radius,
                depth,
            }
        })
        .collect();

    // Sort back-to-front so front axes paint over rear axes.
    nodes.sort_by(|a, b| a.depth.total_cmp(&b.depth));

    // Find front-most hovered axis cap.
    let hovered_node: Option<usize> = resp.hover_pos().and_then(|hp| {
        nodes
            .iter()
            .enumerate()
            .rev()
            .find(|(_, n)| {
                let r = if n.pos_axis { 9.5 } else { 7.5 };
                hp.distance(n.screen) <= r + 2.5
            })
            .map(|(i, _)| i)
    });

    for (i, n) in nodes.iter().enumerate() {
        let is_hov = hovered_node == Some(i);
        let alpha = if n.depth >= -0.05 { 245 } else { 145 };
        let cap_col = Color32::from_rgba_unmultiplied(n.col.r(), n.col.g(), n.col.b(), alpha);
        if n.pos_axis {
            painter.line_segment([center, n.screen], Stroke::new(2.0, cap_col));
            let cap_r = if is_hov { 9.5 } else { 8.0 };
            painter.circle_filled(n.screen, cap_r, cap_col);
            if is_hov {
                painter.circle_stroke(n.screen, cap_r + 1.0, Stroke::new(1.5, Color32::WHITE));
            }
            painter.text(
                n.screen,
                egui::Align2::CENTER_CENTER,
                &n.label[1..],
                FontId::proportional(10.5),
                Color32::from_rgb(16, 20, 26),
            );
        } else {
            if hovered_disc || is_hov {
                painter.line_segment(
                    [center, n.screen],
                    Stroke::new(
                        1.0,
                        Color32::from_rgba_unmultiplied(n.col.r(), n.col.g(), n.col.b(), 80),
                    ),
                );
            }
            let cap_r = if is_hov { 8.0 } else { 5.5 };
            painter.circle_filled(
                n.screen,
                cap_r,
                Color32::from_rgba_unmultiplied(24, 30, 42, if is_hov { 235 } else { 170 }),
            );
            painter.circle_stroke(
                n.screen,
                cap_r,
                Stroke::new(1.4, if is_hov { Color32::WHITE } else { cap_col }),
            );
            if is_hov {
                painter.text(
                    n.screen,
                    egui::Align2::CENTER_CENTER,
                    n.label,
                    FontId::proportional(9.0),
                    Color32::WHITE,
                );
            }
        }
        let _ = n.dir;
    }

    if let Some(i) = hovered_node {
        let preset = nodes[i].preset;
        resp.clone()
            .on_hover_text(format!("Click: {}", preset.label()));
        if resp.clicked() {
            return Some(preset);
        }
    } else if hovered_disc {
        resp.on_hover_text("Click an axis to snap view (1: Front, 3: Side, 7: Top, 0: Iso)");
    }

    None
}

fn draw_hud_toolbar(
    ui: &mut Ui,
    rect: Rect,
    cam: &mut Camera,
    prefs: &mut ViewPrefs,
    bounds: BoundingSphere,
    now_sec: f64,
    id_salt: &str,
) -> Rect {
    let bar_rect = Rect::from_min_max(
        Pos2::new(rect.left() + 8.0, rect.top() + 8.0),
        Pos2::new(rect.right() - 8.0, rect.top() + 36.0),
    );
    let mut child = ui.new_child(
        UiBuilder::new()
            .id_salt(Id::new("view3d_hud").with(id_salt))
            .max_rect(bar_rect)
            .layout(egui::Layout::right_to_left(egui::Align::Center)),
    );
    child.spacing_mut().item_spacing.x = 4.0;

    let inner = egui::Frame::new()
        .fill(Color32::from_rgba_unmultiplied(16, 20, 28, 200))
        .stroke(Stroke::new(0.8, Color32::from_white_alpha(35)))
        .corner_radius(CornerRadius::same(5))
        .inner_margin(egui::Margin::symmetric(6, 3))
        .show(&mut child, |ui| {
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 4.0;
                // Projection toggle
                let proj_label = if prefs.perspective { "Persp" } else { "Ortho" };
                if ui
                    .small_button(proj_label)
                    .on_hover_text("Toggle Perspective / Orthographic (5)")
                    .clicked()
                {
                    prefs.perspective = !prefs.perspective;
                }
                ui.separator();
                if ui
                    .small_button("Iso")
                    .on_hover_text("Isometric view (0 / R)")
                    .clicked()
                {
                    cam.set_preset(ViewPreset::Iso, None, Some(now_sec));
                    ui.ctx().request_repaint();
                }
                if ui
                    .small_button("Top")
                    .on_hover_text("Top view (7, Shift+7 Bottom)")
                    .clicked()
                {
                    cam.set_preset(ViewPreset::Top, None, Some(now_sec));
                    ui.ctx().request_repaint();
                }
                if ui
                    .small_button("Side")
                    .on_hover_text("Starboard side view (3, Shift+3 Port)")
                    .clicked()
                {
                    cam.set_preset(ViewPreset::Right, None, Some(now_sec));
                    ui.ctx().request_repaint();
                }
                if ui
                    .small_button("Front")
                    .on_hover_text("Nose-on front view (1, Shift+1 Rear)")
                    .clicked()
                {
                    cam.set_preset(ViewPreset::Front, None, Some(now_sec));
                    ui.ctx().request_repaint();
                }
                ui.separator();
                if ui
                    .small_button("Fit")
                    .on_hover_text("Frame aircraft in view (F or double-click MMB)")
                    .clicked()
                {
                    cam.fit_sphere(bounds, rect.size(), prefs.perspective, Some(now_sec));
                    ui.ctx().request_repaint();
                }
                let grid_btn = RichText::new("Grid").color(if prefs.show_grid {
                    Color32::from_rgb(130, 210, 255)
                } else {
                    ui.visuals().weak_text_color()
                });
                if ui
                    .small_button(grid_btn)
                    .on_hover_text("Toggle reference grid (G)")
                    .clicked()
                {
                    prefs.show_grid = !prefs.show_grid;
                }
                ui.menu_button("Nav", |ui| {
                    ui.label(RichText::new("Mouse controls").strong().small());
                    ui.label(prefs.nav_preset.hint());
                    for preset in [NavPreset::Standard, NavPreset::CadBlender] {
                        if ui
                            .selectable_value(&mut prefs.nav_preset, preset, preset.label())
                            .clicked()
                        {
                            ui.close();
                        }
                    }
                    ui.separator();
                    ui.label(RichText::new("Orbit mode").strong().small());
                    for mode in [OrbitMode::Turntable, OrbitMode::Trackball] {
                        if ui
                            .selectable_value(&mut prefs.orbit_mode, mode, mode.label())
                            .clicked()
                        {
                            cam.sync_mode(mode);
                            ui.close();
                        }
                    }
                    ui.separator();
                    if ui.button("Reset view (R)").clicked() {
                        cam.set_preset(
                            ViewPreset::Iso,
                            Some((bounds, rect.size(), prefs.perspective)),
                            Some(now_sec),
                        );
                        ui.ctx().request_repaint();
                        ui.close();
                    }
                });
            });
        });
    inner.response.rect
}

fn build_tri_vertex_buffer(
    mesh: &Mesh,
    cache: &MeshCache,
    field: Option<&ColorField<'_>>,
    highlight: SelectionHighlight,
) -> Vec<f32> {
    let mut out = Vec::with_capacity(cache.triangles.len() * 3 * 9);
    for tri in &cache.triangles {
        let col =
            pick::panel_base_color(&mesh.panels[tri.panel_idx], tri.panel_idx, field, highlight);
        let c = [
            col.r() as f32 / 255.0,
            col.g() as f32 / 255.0,
            col.b() as f32 / 255.0,
        ];
        for v in [tri.a, tri.b, tri.c] {
            out.extend_from_slice(&[
                v.x as f32,
                v.y as f32,
                v.z as f32,
                tri.normal[0],
                tri.normal[1],
                tri.normal[2],
                c[0],
                c[1],
                c[2],
            ]);
        }
    }
    out
}

fn flatten_vec3(pts: &[[f32; 3]]) -> Vec<f32> {
    let mut out = Vec::with_capacity(pts.len() * 3);
    for p in pts {
        out.extend_from_slice(p);
    }
    out
}

fn rgba_f32(c: Color32) -> [f32; 4] {
    [
        c.r() as f32 / 255.0,
        c.g() as f32 / 255.0,
        c.b() as f32 / 255.0,
        c.a() as f32 / 255.0,
    ]
}

fn push_colored_segment(out: &mut Vec<f32>, a: Vec3, b: Vec3, col: [f32; 4]) {
    out.extend_from_slice(&[
        a.x as f32, a.y as f32, a.z as f32, col[0], col[1], col[2], col[3],
    ]);
    out.extend_from_slice(&[
        b.x as f32, b.y as f32, b.z as f32, col[0], col[1], col[2], col[3],
    ]);
}

fn pack_colored_segments(lines: &[(Vec3, Vec3, Color32)]) -> Vec<f32> {
    let mut out = Vec::with_capacity(lines.len() * 14);
    for &(a, b, col) in lines {
        push_colored_segment(&mut out, a, b, rgba_f32(col));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use foiley_plane::mesh::{Mesh, Panel, PanelKind};

    fn make_quad(z: f64, wing: usize) -> Panel {
        let nodes = [
            Vec3::new(-0.5, -0.5, z),
            Vec3::new(0.5, -0.5, z),
            Vec3::new(0.5, 0.5, z),
            Vec3::new(-0.5, 0.5, z),
        ];
        let centre = (nodes[0] + nodes[1] + nodes[2] + nodes[3]) * 0.25;
        Panel {
            nodes,
            normal: Vec3::Z,
            area: 1.0,
            centre,
            l: Vec3::X,
            m: Vec3::Y,
            size: 2.0_f64.sqrt(),
            ctrl: centre,
            va: nodes[0],
            vb: nodes[3],
            kind: PanelKind::Mid,
            wing,
            body_id: None,
            strip: None,
            trailing: false,
            leading: false,
            flap: None,
        }
    }

    #[test]
    fn view_presets_match_foiley_axes_without_pole_singularity() {
        let mut cam = Camera::default();
        let all = [
            (ViewPreset::Front, Vec3::X, Vec3::Z),
            (ViewPreset::Rear, -Vec3::X, Vec3::Z),
            (ViewPreset::Right, -Vec3::Y, Vec3::Z),
            (ViewPreset::Left, Vec3::Y, Vec3::Z),
            (ViewPreset::Top, -Vec3::Z, -Vec3::X),
            (ViewPreset::Bottom, Vec3::Z, -Vec3::X),
        ];
        for (preset, exp_fwd, exp_up) in all {
            cam.set_preset(preset, None, None);
            let b = cam.basis();
            assert!(
                (b.forward - exp_fwd).norm() < 1e-9,
                "{preset:?} forward mismatch: {:?}",
                b.forward
            );
            assert!(
                (b.up - exp_up).norm() < 1e-9,
                "{preset:?} up mismatch: {:?}",
                b.up
            );
            assert!(b.right.dot(b.up).abs() < 1e-9, "{preset:?} non-orthogonal");
            assert!(
                (b.right.cross(b.up) + b.forward).norm() < 1e-9,
                "{preset:?} not right-handed"
            );
        }
        // Top view should place starboard (+Y) on the right; Bottom places starboard on the left.
        cam.set_preset(ViewPreset::Top, None, None);
        assert!((cam.basis().right - Vec3::Y).norm() < 1e-9);
        cam.set_preset(ViewPreset::Bottom, None, None);
        assert!((cam.basis().right + Vec3::Y).norm() < 1e-9);
    }

    #[test]
    fn cursor_anchored_zoom_keeps_anchor_pinned_in_ortho_and_perspective() {
        let rect = Rect::from_min_size(Pos2::new(40.0, 20.0), Vec2::new(800.0, 600.0));
        let cursor = Pos2::new(620.0, 190.0);

        for perspective in [true, false] {
            let mut cam = Camera::default();
            cam.target = Vec3::new(0.25, -0.1, 0.05);
            cam.distance = 6.0;
            let (ro, rd) = cam.screen_ray(cursor, rect, perspective);
            // Pick a 3D anchor off the target plane along the cursor ray.
            let anchor = ro + rd * 5.2;
            let (near0, _) = cam.clip_planes(None, 0.0, perspective);
            let (p_before, _) = cam.project_point(anchor, rect, near0, perspective).unwrap();
            assert!(
                (p_before - cursor).length() < 1e-2,
                "initial projection mismatch: {p_before:?} vs {cursor:?}"
            );

            for factor in [0.65, 0.72, 1.45, 1.3] {
                cam.zoom_anchored(factor, cursor, rect, perspective, anchor, 1.5);
                let (near1, _) = cam.clip_planes(None, 0.0, perspective);
                let (p_after, _) = cam.project_point(anchor, rect, near1, perspective).unwrap();
                assert!(
                    (p_after - cursor).length() < 0.05,
                    "cursor anchor drifted (persp={perspective}, factor={factor}): {p_after:?} vs {cursor:?}"
                );
            }
        }
    }

    #[test]
    fn ortho_and_perspective_preserve_target_plane_scale_and_match_gl_matrix() {
        let rect = Rect::from_min_size(Pos2::new(10.0, 15.0), Vec2::new(960.0, 540.0));
        let aspect = rect.width() as f64 / rect.height() as f64;
        let cam = Camera::default();
        let b = cam.basis();
        let p_target_off = cam.target + b.right * 0.65 + b.up * 0.35;

        let (near_p, far_p) = cam.clip_planes(None, 0.0, true);
        let (near_o, far_o) = cam.clip_planes(None, 0.0, false);
        let (sp, _) = cam.project_point(p_target_off, rect, near_p, true).unwrap();
        let (so, _) = cam
            .project_point(p_target_off, rect, near_o, false)
            .unwrap();
        assert!(
            (sp - so).length() < 1e-3,
            "target-plane scale mismatch on toggle: {sp:?} vs {so:?}"
        );

        // Verify GL column-major view_proj_matrix agrees with CPU project_point.
        for (persp, near, far) in [(true, near_p, far_p), (false, near_o, far_o)] {
            let mvp = cam.view_proj_matrix(aspect, near, far, persp);
            let test_pt = Vec3::new(0.4, -0.7, 0.25);
            let v = [test_pt.x as f32, test_pt.y as f32, test_pt.z as f32, 1.0];
            let clip = [
                mvp[0] * v[0] + mvp[4] * v[1] + mvp[8] * v[2] + mvp[12] * v[3],
                mvp[1] * v[0] + mvp[5] * v[1] + mvp[9] * v[2] + mvp[13] * v[3],
                mvp[2] * v[0] + mvp[6] * v[1] + mvp[10] * v[2] + mvp[14] * v[3],
                mvp[3] * v[0] + mvp[7] * v[1] + mvp[11] * v[2] + mvp[15] * v[3],
            ];
            let ndc_x = clip[0] / clip[3];
            let ndc_y = clip[1] / clip[3];
            let gl_screen = Pos2::new(
                rect.left() + (ndc_x * 0.5 + 0.5) * rect.width(),
                rect.top() + (0.5 - ndc_y * 0.5) * rect.height(),
            );
            let (cpu_screen, _) = cam.project_point(test_pt, rect, near, persp).unwrap();
            assert!(
                (gl_screen - cpu_screen).length() < 1e-2,
                "GL vs CPU screen mismatch (persp={persp}): {gl_screen:?} vs {cpu_screen:?}"
            );
        }
    }

    #[test]
    fn turntable_and_trackball_orbit_same_direction_and_sync_cleanly() {
        let mut tt = Camera::default();
        let mut tb = Camera::default();
        tb.sync_mode(OrbitMode::Trackball);

        let eye0 = tt.eye();
        tt.orbit(Vec2::new(20.0, 15.0), OrbitMode::Turntable);
        tb.orbit(Vec2::new(20.0, 15.0), OrbitMode::Trackball);
        let d_tt = (tt.eye() - eye0).normalized();
        let d_tb = (tb.eye() - eye0).normalized();
        assert!(
            d_tt.dot(d_tb) > 0.98,
            "orbit direction mismatch: dot={}",
            d_tt.dot(d_tb)
        );

        tb.sync_mode(OrbitMode::Turntable);
        let b = tb.basis();
        assert!(b.right.x.is_finite() && b.up.z.is_finite() && b.forward.y.is_finite());
    }

    #[test]
    fn nearest_hit_picking_prefers_front_panel_and_handles_triangles() {
        let back_quad = make_quad(0.0, 0);
        let front_quad = make_quad(1.0, 1);
        let tri_panel = Panel::tri(
            Vec3::new(2.0, -0.5, 0.5),
            Vec3::new(3.0, -0.5, 0.5),
            Vec3::new(2.5, 0.5, 0.5),
            PanelKind::Fuse,
            usize::MAX,
            None,
        );
        let mesh = Mesh {
            panels: vec![back_quad, front_quad, tri_panel],
            strips: vec![],
            wings: vec![],
            flaps: vec![],
            thick: true,
            neighbours: vec![],
            n_wing_panels: 2,
            carry: vec![],
        };
        let cache = MeshCache::build(42, &mesh);

        // 2 quads (4 tris, 8 edges) + 1 tri (1 tri, 3 edges) = 5 tris, 11 edges (22 verts).
        assert_eq!(cache.triangles.len(), 5);
        assert_eq!(cache.edge_positions.len(), 22);

        // Ray straight down from z = 5 hits front_quad (panel 1) first, not back_quad (panel 0).
        let hit = cache
            .pick_nearest(Vec3::new(0.0, 0.0, 5.0), -Vec3::Z, 0.1, 100.0)
            .unwrap();
        assert_eq!(hit.panel_idx, 1);
        assert!((hit.t - 4.0).abs() < 1e-9);

        // Clipping near plane past z=1 (t_min = 4.2) skips panel 1 and picks panel 0.
        let hit_clipped = cache
            .pick_nearest(Vec3::new(0.0, 0.0, 5.0), -Vec3::Z, 4.2, 100.0)
            .unwrap();
        assert_eq!(hit_clipped.panel_idx, 0);

        // Ray over the triangular panel hits panel 2 cleanly.
        let hit_tri = cache
            .pick_nearest(Vec3::new(2.5, 0.0, 5.0), -Vec3::Z, 0.1, 100.0)
            .unwrap();
        assert_eq!(hit_tri.panel_idx, 2);
    }

    #[test]
    fn cpu_near_plane_clipping_keeps_partial_triangles_in_front() {
        // Two vertices in front (depth 2.0), one vertex behind eye (depth -1.0).
        let tri = [
            Vec3::new(-1.0, 0.0, -2.0),
            Vec3::new(1.0, 0.0, -2.0),
            Vec3::new(0.0, 1.0, 1.0),
        ];
        let clipped = cpu::clip_poly_near(&tri, 0.5);
        assert_eq!(clipped.len(), 4);
        for v in &clipped {
            assert!(
                -v.z >= 0.5 - 1e-9,
                "clipped vertex behind near plane: {v:?}"
            );
        }
    }

    #[test]
    fn fit_preserves_the_interrupted_preset_orientation() {
        let mut cam = Camera::default();
        let screen = Vec2::new(800.0, 600.0);
        cam.set_preset(ViewPreset::Top, None, Some(0.0));
        assert!(cam.step_anim(0.04));
        let before = cam.basis();
        cam.fit_sphere(BoundingSphere::default(), screen, true, Some(0.04));
        assert!(!cam.step_anim(1.0));
        let after = cam.basis();
        assert!((before.forward - after.forward).norm() < 1e-12);
        assert!((before.up - after.up).norm() < 1e-12);
    }

    #[test]
    fn off_axis_pick_respects_eye_space_near_clip() {
        let mut cam = Camera::default();
        cam.set_preset(ViewPreset::Top, None, None);
        cam.distance = 1.0;
        let rect = Rect::from_min_size(Pos2::ZERO, Vec2::new(800.0, 600.0));
        let pos = Pos2::new(700.0, 300.0);
        let mesh = Mesh {
            panels: vec![make_quad(0.02, 0), make_quad(-0.2, 1)],
            strips: vec![],
            wings: vec![],
            flaps: vec![],
            thick: false,
            neighbours: vec![],
            n_wing_panels: 2,
            carry: vec![],
        };
        let cache = MeshCache::build(0, &mesh);
        assert_eq!(
            pick_at(&cache, &cam, rect, true, pos, 1.0, 10.0)
                .unwrap()
                .panel_idx,
            1
        );
    }

    #[test]
    fn egui_viewport_mmb_orbit_shift_pan_cad_preset_and_hotkeys() {
        let ctx = egui::Context::default();
        let mesh = Mesh {
            panels: vec![make_quad(0.0, 0)],
            strips: vec![],
            wings: vec![],
            flaps: vec![],
            thick: false,
            neighbours: vec![],
            n_wing_panels: 1,
            carry: vec![],
        };
        let mut cam = Camera::default();
        let mut prefs = ViewPrefs::default();
        let mut cache: Option<MeshCache> = None;
        let gl = Arc::new(Mutex::new(GlRenderer::default()));
        let ov = Overlay::default();
        let screen = Rect::from_min_size(Pos2::ZERO, Vec2::new(800.0, 600.0));
        let canvas_id = std::cell::Cell::new(None);

        let run_frame = |raw: egui::RawInput,
                         cam: &mut Camera,
                         prefs: &mut ViewPrefs,
                         cache: &mut Option<MeshCache>|
         -> ViewportEvent {
            let mut ev_out = ViewportEvent::default();
            let _ = ctx.run(raw, |ctx| {
                egui::CentralPanel::default().show(ctx, |ui| {
                    canvas_id.set(Some(
                        ui.make_persistent_id(Id::new("view3d_canvas").with("test_vp")),
                    ));
                    let fmt = |idx: usize| format!("panel #{idx}");
                    let inp = ViewportInput {
                        plane_id: 0,
                        mesh: Some((1, &mesh)),
                        field: None,
                        overlay: &ov,
                        show_edges: true,
                        highlight: SelectionHighlight::Wing(0),
                        hover_label: Some(&fmt),
                        notice: None,
                    };
                    ev_out = show(ui, "test_vp", cam, prefs, cache, &gl, &inp);
                });
            });
            ev_out
        };

        // Frame 1: auto-fit on initial valid mesh.
        let mut raw = egui::RawInput {
            screen_rect: Some(screen),
            ..Default::default()
        };
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);
        assert_eq!(cam.fitted_plane, Some(0));
        assert!(!cam.needs_fit);

        // MMB drag orbits in both Standard and CadBlender presets.
        prefs.nav_preset = NavPreset::CadBlender;
        let yaw0 = cam.yaw;
        let pitch0 = cam.pitch;
        let center = Pos2::new(400.0, 300.0);
        raw.events = vec![
            egui::Event::PointerMoved(center),
            egui::Event::PointerButton {
                pos: center,
                button: PointerButton::Middle,
                pressed: true,
                modifiers: egui::Modifiers::NONE,
            },
        ];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);
        raw.events = vec![egui::Event::PointerMoved(center + Vec2::new(40.0, 20.0))];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);
        assert!((cam.yaw - yaw0).abs() > 0.1, "MMB drag should orbit yaw");
        assert!(
            (cam.pitch - pitch0).abs() > 0.05,
            "MMB drag should orbit pitch"
        );
        raw.events = vec![egui::Event::PointerButton {
            pos: center + Vec2::new(40.0, 20.0),
            button: PointerButton::Middle,
            pressed: false,
            modifiers: egui::Modifiers::NONE,
        }];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);

        // Ctrl+MMB zoom retains the pivot during off-center input.
        raw.modifiers = egui::Modifiers::CTRL;
        let offset = center + Vec2::new(130.0, 0.0);
        let zoom_target = cam.target;
        let zoom_distance = cam.distance;
        raw.events = vec![
            egui::Event::PointerMoved(offset),
            egui::Event::PointerButton {
                pos: offset,
                button: PointerButton::Middle,
                pressed: true,
                modifiers: raw.modifiers,
            },
        ];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);
        raw.events = vec![egui::Event::PointerMoved(offset + Vec2::new(30.0, 0.0))];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);
        assert_eq!(cam.distance, zoom_distance);
        assert_eq!(cam.target, zoom_target);
        let end = offset + Vec2::new(30.0, -35.0);
        raw.events = vec![egui::Event::PointerMoved(end)];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);
        assert!(cam.distance < zoom_distance);
        assert!((cam.target - zoom_target).norm() < 1e-12);
        raw.events = vec![egui::Event::PointerButton {
            pos: end,
            button: PointerButton::Middle,
            pressed: false,
            modifiers: raw.modifiers,
        }];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);

        // Shift + MMB drag pans target in world space without rotating.
        let yaw_before_pan = cam.yaw;
        let target_before_pan = cam.target;
        raw.modifiers = egui::Modifiers::SHIFT;
        raw.events = vec![
            egui::Event::PointerMoved(center),
            egui::Event::PointerButton {
                pos: center,
                button: PointerButton::Middle,
                pressed: true,
                modifiers: egui::Modifiers::SHIFT,
            },
        ];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);
        raw.events = vec![egui::Event::PointerMoved(center + Vec2::new(50.0, 0.0))];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);
        assert!(
            (cam.target - target_before_pan).norm() > 1e-3,
            "Shift+MMB drag should pan target"
        );
        assert!(
            (cam.yaw - yaw_before_pan).abs() < 1e-9,
            "Shift+MMB pan must not orbit"
        );
        raw.events = vec![egui::Event::PointerButton {
            pos: center + Vec2::new(50.0, 0.0),
            button: PointerButton::Middle,
            pressed: false,
            modifiers: egui::Modifiers::SHIFT,
        }];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);

        // In CadBlender mode, LMB drag must NOT orbit.
        raw.modifiers = egui::Modifiers::NONE;
        let yaw_cad = cam.yaw;
        raw.events = vec![
            egui::Event::PointerMoved(center),
            egui::Event::PointerButton {
                pos: center,
                button: PointerButton::Primary,
                pressed: true,
                modifiers: egui::Modifiers::NONE,
            },
        ];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);
        raw.events = vec![egui::Event::PointerMoved(center + Vec2::new(60.0, 30.0))];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);
        assert!(
            (cam.yaw - yaw_cad).abs() < 1e-9,
            "CadBlender mode must reserve LMB drag"
        );
        raw.events = vec![egui::Event::PointerButton {
            pos: center + Vec2::new(60.0, 30.0),
            button: PointerButton::Primary,
            pressed: false,
            modifiers: egui::Modifiers::NONE,
        }];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);

        // Key `5` toggles perspective and ignores repeat; `G` toggles grid.
        let p0 = prefs.perspective;
        let g0 = prefs.show_grid;
        raw.events = vec![
            egui::Event::PointerMoved(center),
            egui::Event::Key {
                key: Key::Num5,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: egui::Modifiers::NONE,
            },
            egui::Event::Key {
                key: Key::Num5,
                physical_key: None,
                pressed: true,
                repeat: true,
                modifiers: egui::Modifiers::NONE,
            },
            egui::Event::Key {
                key: Key::G,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: egui::Modifiers::NONE,
            },
        ];
        run_frame(raw.clone(), &mut cam, &mut prefs, &mut cache);
        assert_eq!(
            prefs.perspective, !p0,
            "Key 5 should toggle perspective once (ignoring repeat)"
        );
        assert_eq!(prefs.show_grid, !g0, "Key G should toggle reference grid");

        ctx.memory_mut(|memory| memory.request_focus(canvas_id.get().unwrap()));
        raw.events = vec![
            egui::Event::PointerGone,
            egui::Event::Key {
                key: Key::G,
                physical_key: None,
                pressed: false,
                repeat: false,
                modifiers: egui::Modifiers::NONE,
            },
            egui::Event::Key {
                key: Key::G,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: egui::Modifiers::NONE,
            },
        ];
        run_frame(raw, &mut cam, &mut prefs, &mut cache);
        assert_eq!(
            prefs.show_grid, g0,
            "The focused canvas should handle G without pointer hover"
        );
    }
}
