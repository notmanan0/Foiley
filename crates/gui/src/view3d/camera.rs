//! World-space 3D orbit camera, projection matrices, view presets, and
//! cursor-anchored zoom math.
//!
//! Foiley world coordinates (see `foiley_plane::geom`):
//! - `+X` aft (tail), `-X` forward (nose)
//! - `+Y` starboard (right wing), `-Y` port (left wing)
//! - `+Z` up

use eframe::egui::{Pos2, Rect, Vec2};
use foiley_core::Vec3;
use std::f64::consts::{FRAC_PI_2, PI};

/// Mouse / trackpad navigation scheme.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NavPreset {
    /// Trackpad & mouse friendly: LMB click selects, LMB/MMB drag orbits,
    /// Shift+LMB / Shift+MMB / RMB drag pans, wheel / Ctrl+MMB zooms.
    #[default]
    Standard,
    /// CAD / Blender style: MMB drag orbits, Shift+MMB / RMB drag pans,
    /// wheel / Ctrl+MMB zooms, LMB click selects (LMB drag reserved).
    CadBlender,
}

impl NavPreset {
    pub fn label(self) -> &'static str {
        match self {
            NavPreset::Standard => "Standard (LMB/MMB)",
            NavPreset::CadBlender => "CAD / Blender (MMB)",
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            NavPreset::Standard => "standard",
            NavPreset::CadBlender => "cad",
        }
    }

    pub fn from_str(s: &str) -> Self {
        match s {
            "cad" => NavPreset::CadBlender,
            _ => NavPreset::Standard,
        }
    }

    pub fn hint(self) -> &'static str {
        match self {
            NavPreset::Standard => "LMB/MMB drag: orbit | Shift+drag or RMB: pan | wheel / Ctrl+MMB: zoom | LMB click: select | F: fit",
            NavPreset::CadBlender => "MMB drag: orbit | Shift+MMB or RMB: pan | wheel / Ctrl+MMB: zoom | LMB click: select | F: fit",
        }
    }
}

/// Orbit constraint mode.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum OrbitMode {
    /// Z-up constrained turntable (keeps aircraft horizon level).
    #[default]
    Turntable,
    /// Unconstrained quaternion trackball.
    Trackball,
}

impl OrbitMode {
    pub fn label(self) -> &'static str {
        match self {
            OrbitMode::Turntable => "Turntable (Z-up)",
            OrbitMode::Trackball => "Trackball (Free)",
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            OrbitMode::Turntable => "turntable",
            OrbitMode::Trackball => "trackball",
        }
    }

    pub fn from_str(s: &str) -> Self {
        match s {
            "trackball" => OrbitMode::Trackball,
            _ => OrbitMode::Turntable,
        }
    }
}

/// Persisted viewport preferences (owned by `App`, survives project loads).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ViewPrefs {
    pub nav_preset: NavPreset,
    pub orbit_mode: OrbitMode,
    pub perspective: bool,
    pub show_grid: bool,
}

impl Default for ViewPrefs {
    fn default() -> Self {
        Self {
            nav_preset: NavPreset::Standard,
            orbit_mode: OrbitMode::Turntable,
            perspective: true,
            show_grid: true,
        }
    }
}

/// Standard aerodynamic camera viewpoints.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ViewPreset {
    Iso,
    /// Eye on `-X` looking aft (`+X`): nose-on view.
    Front,
    /// Eye on `+X` looking forward (`-X`): tail-on view.
    Rear,
    /// Eye on `+Y` looking port (`-Y`): starboard side view.
    Right,
    /// Eye on `-Y` looking starboard (`+Y`): port side view.
    Left,
    /// Eye on `+Z` looking down (`-Z`), nose (`-X`) pointing up on screen.
    Top,
    /// Eye on `-Z` looking up (`+Z`), nose (`-X`) pointing up on screen.
    Bottom,
}

impl ViewPreset {
    pub fn yaw_pitch(self) -> (f64, f64) {
        match self {
            ViewPreset::Iso => (3.0 * PI / 4.0, 0.52),
            ViewPreset::Front => (PI, 0.0),
            ViewPreset::Rear => (0.0, 0.0),
            ViewPreset::Right => (FRAC_PI_2, 0.0),
            ViewPreset::Left => (-FRAC_PI_2, 0.0),
            ViewPreset::Top => (0.0, FRAC_PI_2),
            ViewPreset::Bottom => (PI, -FRAC_PI_2),
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            ViewPreset::Iso => "Iso",
            ViewPreset::Front => "Front (-X)",
            ViewPreset::Rear => "Rear (+X)",
            ViewPreset::Right => "Right (+Y)",
            ViewPreset::Left => "Left (-Y)",
            ViewPreset::Top => "Top (+Z)",
            ViewPreset::Bottom => "Bottom (-Z)",
        }
    }
}

/// 3D bounding sphere of the active mesh (excludes wake lines).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundingSphere {
    pub center: Vec3,
    pub radius: f64,
}

impl Default for BoundingSphere {
    fn default() -> Self {
        Self {
            center: Vec3::ZERO,
            radius: 1.0,
        }
    }
}

/// Orthonormal right-handed camera basis in world space.
/// `right × up = -forward` (eye looks along `forward`).
#[derive(Clone, Copy, Debug)]
pub struct CameraBasis {
    pub right: Vec3,
    pub up: Vec3,
    pub forward: Vec3,
}

#[derive(Clone, Debug)]
struct CameraAnim {
    start_target: Vec3,
    end_target: Vec3,
    start_distance: f64,
    end_distance: f64,
    start_quat: [f64; 4],
    end_quat: [f64; 4],
    end_yaw: f64,
    end_pitch: f64,
    start_time: f64,
    duration: f64,
}

/// World-space 3D orbit camera.
#[derive(Clone, Debug)]
pub struct Camera {
    pub target: Vec3,
    pub distance: f64,
    pub yaw: f64,
    pub pitch: f64,
    /// Unit quaternion `[x, y, z, w]` mapping camera `(right, up, -forward)` into world space.
    pub quat: [f64; 4],
    pub fov_y_rad: f64,
    /// Plane identity that was last auto-fitted; independent of mesh edit revisions.
    pub fitted_plane: Option<usize>,
    /// Explicitly request a fit on the next valid mesh frame.
    pub needs_fit: bool,
    anim: Option<CameraAnim>,
}

impl Default for Camera {
    fn default() -> Self {
        let (yaw, pitch) = ViewPreset::Iso.yaw_pitch();
        let quat = quat_from_basis(turntable_basis(yaw, pitch));
        Self {
            target: Vec3::ZERO,
            distance: 4.0,
            yaw,
            pitch,
            quat,
            fov_y_rad: 35.0_f64.to_radians(),
            fitted_plane: None,
            needs_fit: true,
            anim: None,
        }
    }
}

impl Camera {
    /// Current orthonormal basis, including any rotation interpolated by an animation.
    pub fn basis(&self) -> CameraBasis {
        basis_from_quat(self.quat)
    }

    /// World-space camera eye position.
    pub fn eye(&self) -> Vec3 {
        let b = self.basis();
        self.target - b.forward * self.distance
    }

    /// Half-height of the view frustum at the `target` plane (shared between
    /// Perspective and Orthographic so toggling projections preserves scale).
    pub fn target_half_height(&self) -> f64 {
        self.distance * (0.5 * self.fov_y_rad).tan()
    }

    /// Sync internal quaternion/Euler state when switching orbit modes.
    pub fn sync_mode(&mut self, mode: OrbitMode) {
        self.anim = None;
        if mode == OrbitMode::Turntable {
            self.quat = quat_from_basis(turntable_basis(self.yaw, self.pitch));
        }
    }

    /// Cancel any active transition animation (e.g. on user pointer input).
    fn cancel_anim(&mut self) {
        self.anim = None;
    }

    /// Advance smooth camera transition. Returns `true` if still animating.
    pub fn step_anim(&mut self, now_sec: f64) -> bool {
        let Some(a) = self.anim.clone() else {
            return false;
        };
        let raw = ((now_sec - a.start_time) / a.duration.max(1e-4)).clamp(0.0, 1.0);
        let t = 1.0 - (1.0 - raw).powi(3); // cubic ease-out
        self.target = a.start_target + (a.end_target - a.start_target) * t;
        self.distance = a.start_distance + (a.end_distance - a.start_distance) * t;
        self.quat = quat_slerp(a.start_quat, a.end_quat, t);
        if raw >= 1.0 {
            self.target = a.end_target;
            self.distance = a.end_distance;
            self.yaw = a.end_yaw;
            self.pitch = a.end_pitch;
            self.quat = a.end_quat;
            self.anim = None;
            false
        } else {
            (self.yaw, self.pitch) = orbit_angles(self.basis());
            true
        }
    }

    /// Orbit by logical pointer delta `(dx, dy)`.
    pub fn orbit(&mut self, delta: Vec2, mode: OrbitMode) {
        self.cancel_anim();
        const SENS: f64 = 0.0075;
        let dx = delta.x as f64 * SENS;
        let dy = delta.y as f64 * SENS;
        match mode {
            OrbitMode::Turntable => {
                self.yaw = wrap_angle(self.yaw - dx);
                let max_pitch = FRAC_PI_2 - 0.004;
                self.pitch = (self.pitch + dy).clamp(-max_pitch, max_pitch);
                self.quat = quat_from_basis(turntable_basis(self.yaw, self.pitch));
            }
            OrbitMode::Trackball => {
                let b = basis_from_quat(self.quat);
                let q_yaw = quat_from_axis_angle(b.up, -dx);
                let q_pitch = quat_from_axis_angle(b.right, -dy);
                self.quat = quat_normalize(quat_mul(quat_mul(q_yaw, q_pitch), self.quat));
                (self.yaw, self.pitch) = orbit_angles(self.basis());
            }
        }
    }

    /// Pan target in world space by logical pointer delta `(dx, dy)`.
    pub fn pan(&mut self, delta: Vec2, viewport_size: Vec2) {
        self.cancel_anim();
        let h = (viewport_size.y as f64).max(1.0);
        let units_per_pt = 2.0 * self.target_half_height() / h;
        let b = self.basis();
        self.target = self.target - b.right * (delta.x as f64 * units_per_pt)
            + b.up * (delta.y as f64 * units_per_pt);
    }

    /// Fit camera distance & target to a 3D bounding sphere without altering orientation.
    pub fn fit_sphere(
        &mut self,
        sphere: BoundingSphere,
        viewport_size: Vec2,
        perspective: bool,
        now_sec: Option<f64>,
    ) {
        let radius = sphere.radius.max(1e-3);
        let dist = fit_distance(radius, viewport_size, self.fov_y_rad, perspective);
        if let Some(now) = now_sec {
            self.anim = Some(CameraAnim {
                start_target: self.target,
                end_target: sphere.center,
                start_distance: self.distance,
                end_distance: dist,
                start_quat: self.quat,
                end_quat: self.quat,
                end_yaw: self.yaw,
                end_pitch: self.pitch,
                start_time: now,
                duration: 0.18,
            });
        } else {
            self.anim = None;
            self.target = sphere.center;
            self.distance = dist;
        }
        self.needs_fit = false;
    }

    /// Smoothly snap orientation to a standard aero preset (optionally also fitting sphere).
    pub fn set_preset(
        &mut self,
        preset: ViewPreset,
        fit: Option<(BoundingSphere, Vec2, bool)>,
        now_sec: Option<f64>,
    ) {
        let (yaw, pitch) = preset.yaw_pitch();
        let end_quat = quat_from_basis(turntable_basis(yaw, pitch));
        let (end_target, end_distance) = match fit {
            Some((s, vp, persp)) => (
                s.center,
                fit_distance(s.radius.max(1e-3), vp, self.fov_y_rad, persp),
            ),
            None => (self.target, self.distance),
        };
        if let Some(now) = now_sec {
            self.anim = Some(CameraAnim {
                start_target: self.target,
                end_target,
                start_distance: self.distance,
                end_distance,
                start_quat: self.quat,
                end_quat,
                end_yaw: yaw,
                end_pitch: pitch,
                start_time: now,
                duration: 0.20,
            });
        } else {
            self.anim = None;
            self.target = end_target;
            self.distance = end_distance;
            self.yaw = yaw;
            self.pitch = pitch;
            self.quat = end_quat;
        }
        if fit.is_some() {
            self.needs_fit = false;
        }
    }

    /// Cursor-anchored exponential zoom. Keeps `anchor` pinned under `cursor`
    /// in both orthographic and perspective projections.
    pub fn zoom_anchored(
        &mut self,
        scale_factor: f64,
        cursor: Pos2,
        rect: Rect,
        perspective: bool,
        anchor: Vec3,
        bounds_radius: f64,
    ) {
        if !scale_factor.is_finite()
            || scale_factor <= 0.0
            || rect.width() < 1.0
            || rect.height() < 1.0
        {
            return;
        }
        self.cancel_anim();
        let r = bounds_radius.max(0.05);
        let min_dist = (r * 0.02).max(1e-3);
        let max_dist = (r * 80.0).max(50.0);

        let b = self.basis();
        let mut new_dist = (self.distance * scale_factor).clamp(min_dist, max_dist);
        if perspective {
            // Keep the anchor strictly in front of the near plane.
            // Anchor depth relative to target is `anchor_off = (anchor - target).dot(forward)`.
            // Depth from eye is `new_dist + anchor_off`.
            let anchor_off = (anchor - self.target).dot(b.forward);
            let min_eye_clear = (r * 0.015).max(5e-3);
            let min_for_anchor = min_eye_clear - anchor_off;
            if min_for_anchor > new_dist {
                new_dist = min_for_anchor.clamp(min_dist, max_dist);
            }
        }
        self.distance = new_dist;

        // Re-cast the cursor ray after changing distance and shift target so
        // the view-facing plane through `anchor` intersects `anchor` again.
        let (ray_o, ray_d) = self.screen_ray(cursor, rect, perspective);
        if let Some(hit_new) = intersect_plane(ray_o, ray_d, anchor, b.forward) {
            self.target += anchor - hit_new;
        }
    }

    /// World-space ray `(origin, unit_dir)` passing through logical screen point `pos`.
    pub fn screen_ray(&self, pos: Pos2, rect: Rect, perspective: bool) -> (Vec3, Vec3) {
        let w = (rect.width() as f64).max(1.0);
        let h = (rect.height() as f64).max(1.0);
        let ndc_x = ((pos.x as f64 - rect.left() as f64) / w) * 2.0 - 1.0;
        let ndc_y = 1.0 - ((pos.y as f64 - rect.top() as f64) / h) * 2.0;
        let aspect = w / h;
        let b = self.basis();
        let half_h = self.target_half_height();
        let half_w = half_h * aspect;
        let eye = self.target - b.forward * self.distance;
        if perspective {
            let point_on_target =
                self.target + b.right * (ndc_x * half_w) + b.up * (ndc_y * half_h);
            let dir = (point_on_target - eye).normalized();
            (eye, dir)
        } else {
            let origin = eye + b.right * (ndc_x * half_w) + b.up * (ndc_y * half_h);
            (origin, b.forward)
        }
    }

    /// Compute dynamic eye-space near and far clip distances from scene bounds
    /// after pan/orbit so 24-bit Z precision stays tight.
    pub fn clip_planes(
        &self,
        bounds: Option<BoundingSphere>,
        extra_extent: f64,
        perspective: bool,
    ) -> (f64, f64) {
        let b = self.basis();
        let eye = self.target - b.forward * self.distance;
        let (center_depth, r) = match bounds {
            Some(s) => {
                let d = (s.center - eye).dot(b.forward);
                (d, s.radius.max(0.05) + extra_extent.max(0.0))
            }
            None => (
                self.distance,
                (self.distance * 2.0).max(1.0) + extra_extent.max(0.0),
            ),
        };
        let pad = (r * 0.25).max(0.05);
        // Ensure target plane is always comfortably inside clip range too.
        let z_min = (center_depth - r - pad).min(self.distance * 0.5);
        let z_max = (center_depth + r + pad).max(self.distance * 1.5);
        if perspective {
            let min_near = (self.distance * 1e-3).clamp(1e-4, 0.25);
            let near = z_min.max(min_near);
            let far = z_max.max(near + (r * 4.0).max(2.0));
            (near, far)
        } else {
            let span = (r * 4.0).max(self.distance * 2.0).max(2.0);
            let near = z_min.min(self.distance - span);
            let far = z_max.max(self.distance + span);
            (near, far.max(near + 1.0))
        }
    }

    /// OpenGL column-major 4x4 view matrix (world -> eye, looking down `-Z`).
    pub fn view_matrix(&self) -> [f32; 16] {
        let b = self.basis();
        let eye = self.target - b.forward * self.distance;
        let rx = b.right.x;
        let ry = b.right.y;
        let rz = b.right.z;
        let ux = b.up.x;
        let uy = b.up.y;
        let uz = b.up.z;
        let zx = -b.forward.x;
        let zy = -b.forward.y;
        let zz = -b.forward.z;
        let tx = -b.right.dot(eye);
        let ty = -b.up.dot(eye);
        let tz = b.forward.dot(eye); // -(-forward).dot(eye)
        [
            rx as f32, ux as f32, zx as f32, 0.0, ry as f32, uy as f32, zy as f32, 0.0, rz as f32,
            uz as f32, zz as f32, 0.0, tx as f32, ty as f32, tz as f32, 1.0,
        ]
    }

    /// OpenGL column-major 4x4 projection matrix (eye -> clip, NDC z in `[-1, 1]`).
    pub fn proj_matrix(&self, aspect: f64, near: f64, far: f64, perspective: bool) -> [f32; 16] {
        let aspect = aspect.max(1e-4);
        let nf = (far - near).max(1e-6);
        if perspective {
            let f = 1.0 / (0.5 * self.fov_y_rad).tan();
            let sx = f / aspect;
            let sy = f;
            let sz = -(far + near) / nf;
            let tz = -(2.0 * far * near) / nf;
            [
                sx as f32, 0.0, 0.0, 0.0, 0.0, sy as f32, 0.0, 0.0, 0.0, 0.0, sz as f32, -1.0, 0.0,
                0.0, tz as f32, 0.0,
            ]
        } else {
            let half_h = self.target_half_height().max(1e-6);
            let half_w = half_h * aspect;
            let sx = 1.0 / half_w;
            let sy = 1.0 / half_h;
            let sz = -2.0 / nf;
            let tz = -(far + near) / nf;
            [
                sx as f32, 0.0, 0.0, 0.0, 0.0, sy as f32, 0.0, 0.0, 0.0, 0.0, sz as f32, 0.0, 0.0,
                0.0, tz as f32, 1.0,
            ]
        }
    }

    /// Combined column-major `proj * view` matrix.
    pub fn view_proj_matrix(
        &self,
        aspect: f64,
        near: f64,
        far: f64,
        perspective: bool,
    ) -> [f32; 16] {
        mat4_mul(
            &self.proj_matrix(aspect, near, far, perspective),
            &self.view_matrix(),
        )
    }

    /// Transform a world point into camera eye coordinates `(ex, ey, ez)` where
    /// `+ex` is right, `+ey` is up, and `-ez` is depth in front of the eye.
    pub fn world_to_eye(&self, p: Vec3) -> Vec3 {
        let b = self.basis();
        let eye = self.target - b.forward * self.distance;
        let d = p - eye;
        Vec3::new(d.dot(b.right), d.dot(b.up), -d.dot(b.forward))
    }

    /// Project a world point to logical screen coordinates `(Pos2, depth)` where
    /// `depth` is positive distance in front of the eye. Returns `None` if clipped
    /// behind the near plane in perspective mode.
    pub fn project_point(
        &self,
        p: Vec3,
        rect: Rect,
        near: f64,
        perspective: bool,
    ) -> Option<(Pos2, f64)> {
        let e = self.world_to_eye(p);
        let depth = -e.z;
        let w = (rect.width() as f64).max(1.0);
        let h = (rect.height() as f64).max(1.0);
        let aspect = w / h;
        let (ndc_x, ndc_y) = if perspective {
            if depth <= near * 0.999 {
                return None;
            }
            let f = 1.0 / (0.5 * self.fov_y_rad).tan();
            ((e.x * f / aspect) / depth, (e.y * f) / depth)
        } else {
            let half_h = self.target_half_height().max(1e-6);
            let half_w = half_h * aspect;
            (e.x / half_w, e.y / half_h)
        };
        let sx = rect.left() as f64 + (ndc_x * 0.5 + 0.5) * w;
        let sy = rect.top() as f64 + (0.5 - ndc_y * 0.5) * h;
        Some((Pos2::new(sx as f32, sy as f32), depth))
    }

    #[cfg(feature = "snap")]
    pub fn is_animating(&self) -> bool {
        self.anim.is_some()
    }
}

/// Fit distance so a sphere of `radius` fits cleanly inside `viewport_size`.
pub fn fit_distance(radius: f64, viewport_size: Vec2, fov_y_rad: f64, perspective: bool) -> f64 {
    let w = (viewport_size.x as f64).max(1.0);
    let h = (viewport_size.y as f64).max(1.0);
    let aspect = (w / h).max(1e-3);
    let half_fov_y = 0.5 * fov_y_rad;
    let margin = 1.18;
    if perspective {
        let half_fov_x = (half_fov_y.tan() * aspect).atan();
        let narrow = half_fov_y.min(half_fov_x).max(1e-3);
        (radius / narrow.sin()) * margin
    } else {
        // ortho_half_height = distance * tan(half_fov_y)
        let min_half = (radius * margin) / aspect.min(1.0);
        min_half / half_fov_y.tan().max(1e-3)
    }
}

/// Intersect ray `(ray_o, ray_d)` with plane `(plane_p, plane_n)`.
pub fn intersect_plane(ray_o: Vec3, ray_d: Vec3, plane_p: Vec3, plane_n: Vec3) -> Option<Vec3> {
    let denom = ray_d.dot(plane_n);
    if denom.abs() < 1e-12 {
        return None;
    }
    let t = (plane_p - ray_o).dot(plane_n) / denom;
    Some(ray_o + ray_d * t)
}

/// Exact singularity-free turntable basis for any `(yaw, pitch)` including `±PI/2`.
pub fn turntable_basis(yaw: f64, pitch: f64) -> CameraBasis {
    let (sy, cy) = yaw.sin_cos();
    let (sp, cp) = pitch.sin_cos();
    let eye_dir = Vec3::new(cp * cy, cp * sy, sp);
    let forward = -eye_dir;
    let right = Vec3::new(-sy, cy, 0.0);
    let up = Vec3::new(-sp * cy, -sp * sy, cp);
    CameraBasis { right, up, forward }
}

fn orbit_angles(basis: CameraBasis) -> (f64, f64) {
    let back = -basis.forward;
    let pitch = back.z.clamp(-1.0, 1.0).asin();
    let yaw = if back.x.hypot(back.y) > 1e-5 {
        back.y.atan2(back.x)
    } else {
        (-basis.right.x).atan2(basis.right.y)
    };
    (yaw, pitch)
}

fn wrap_angle(mut a: f64) -> f64 {
    while a > PI {
        a -= 2.0 * PI;
    }
    while a < -PI {
        a += 2.0 * PI;
    }
    a
}

/// Multiply two column-major 4x4 matrices `a * b`.
pub fn mat4_mul(a: &[f32; 16], b: &[f32; 16]) -> [f32; 16] {
    let mut out = [0.0f32; 16];
    for col in 0..4 {
        for row in 0..4 {
            let mut s = 0.0f32;
            for k in 0..4 {
                s += a[k * 4 + row] * b[col * 4 + k];
            }
            out[col * 4 + row] = s;
        }
    }
    out
}

pub fn quat_normalize(q: [f64; 4]) -> [f64; 4] {
    let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if n < 1e-15 {
        [0.0, 0.0, 0.0, 1.0]
    } else {
        [q[0] / n, q[1] / n, q[2] / n, q[3] / n]
    }
}

pub fn quat_mul(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    [
        a[3] * b[0] + a[0] * b[3] + a[1] * b[2] - a[2] * b[1],
        a[3] * b[1] - a[0] * b[2] + a[1] * b[3] + a[2] * b[0],
        a[3] * b[2] + a[0] * b[1] - a[1] * b[0] + a[2] * b[3],
        a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2],
    ]
}

pub fn quat_from_axis_angle(axis: Vec3, angle: f64) -> [f64; 4] {
    let u = axis.normalized();
    let (s, c) = (0.5 * angle).sin_cos();
    [u.x * s, u.y * s, u.z * s, c]
}

pub fn quat_slerp(a: [f64; 4], mut b: [f64; 4], t: f64) -> [f64; 4] {
    let mut dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
    if dot < 0.0 {
        dot = -dot;
        b = [-b[0], -b[1], -b[2], -b[3]];
    }
    if dot > 0.9995 {
        return quat_normalize([
            a[0] + t * (b[0] - a[0]),
            a[1] + t * (b[1] - a[1]),
            a[2] + t * (b[2] - a[2]),
            a[3] + t * (b[3] - a[3]),
        ]);
    }
    let theta = dot.clamp(-1.0, 1.0).acos();
    let sin_t = theta.sin();
    let w0 = ((1.0 - t) * theta).sin() / sin_t;
    let w1 = (t * theta).sin() / sin_t;
    quat_normalize([
        w0 * a[0] + w1 * b[0],
        w0 * a[1] + w1 * b[1],
        w0 * a[2] + w1 * b[2],
        w0 * a[3] + w1 * b[3],
    ])
}

/// Convert camera basis `(right, up, -forward)` (rotation matrix columns) to unit quaternion.
pub fn quat_from_basis(b: CameraBasis) -> [f64; 4] {
    let z = -b.forward;
    let m00 = b.right.x;
    let m01 = b.up.x;
    let m02 = z.x;
    let m10 = b.right.y;
    let m11 = b.up.y;
    let m12 = z.y;
    let m20 = b.right.z;
    let m21 = b.up.z;
    let m22 = z.z;

    let tr = m00 + m11 + m22;
    let q = if tr > 0.0 {
        let s = (tr + 1.0).sqrt() * 2.0;
        [(m21 - m12) / s, (m02 - m20) / s, (m10 - m01) / s, 0.25 * s]
    } else if m00 > m11 && m00 > m22 {
        let s = (1.0 + m00 - m11 - m22).sqrt() * 2.0;
        [0.25 * s, (m01 + m10) / s, (m02 + m20) / s, (m21 - m12) / s]
    } else if m11 > m22 {
        let s = (1.0 + m11 - m00 - m22).sqrt() * 2.0;
        [(m01 + m10) / s, 0.25 * s, (m12 + m21) / s, (m02 - m20) / s]
    } else {
        let s = (1.0 + m22 - m00 - m11).sqrt() * 2.0;
        [(m02 + m20) / s, (m12 + m21) / s, 0.25 * s, (m10 - m01) / s]
    };
    quat_normalize(q)
}

/// Convert unit quaternion back to camera basis `(right, up, forward)`.
pub fn basis_from_quat(q: [f64; 4]) -> CameraBasis {
    let [x, y, z, w] = quat_normalize(q);
    let xx = x * x;
    let yy = y * y;
    let zz = z * z;
    let xy = x * y;
    let xz = x * z;
    let yz = y * z;
    let wx = w * x;
    let wy = w * y;
    let wz = w * z;

    let right = Vec3::new(1.0 - 2.0 * (yy + zz), 2.0 * (xy + wz), 2.0 * (xz - wy)).normalized();
    let up = Vec3::new(2.0 * (xy - wz), 1.0 - 2.0 * (xx + zz), 2.0 * (yz + wx)).normalized();
    let back = Vec3::new(2.0 * (xz + wy), 2.0 * (yz - wx), 1.0 - 2.0 * (xx + yy)).normalized();
    CameraBasis {
        right,
        up,
        forward: -back,
    }
}
