//! Hardware Z-buffered OpenGL 3D renderer built on `eframe::glow` and
//! `eframe::egui_glow::CallbackFn`.
//!
//! Designed to work on any desktop OpenGL 3.1+ (GLSL 140) or OpenGL ES 3.0+
//! (GLSL 300 es) GPU, including integrated GPUs and software GL fallbacks.
//! The viewport owns its framebuffer and depth buffer. If setup fails the
//! viewport uses the CPU painter.

use super::pick::{FieldRevision, SelectionHighlight};
use eframe::egui::{PaintCallbackInfo, Rect};
use eframe::egui_glow::Painter;
use eframe::egui_glow::ShaderVersion;
use eframe::glow::{self, HasContext};

/// Key identifying the uploaded triangle mesh vertex buffer contents.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MeshBufferKey {
    pub mesh_key: u64,
    pub field_rev: Option<FieldRevision>,
    pub selection: SelectionHighlight,
}

/// Owned per-frame draw payload moved into the `egui_glow::CallbackFn`.
pub struct GlFrame {
    pub view_proj: [f32; 16],
    pub view_dir: [f32; 3],
    pub mesh_buf_key: MeshBufferKey,
    /// Uploaded only when `mesh_buf_key` changed (`pos[3], normal[3], color[3]` per vertex).
    pub tri_upload: Option<Vec<f32>>,
    pub tri_vertex_count: i32,
    pub edge_mesh_key: u64,
    /// Uploaded only when `edge_mesh_key` changed (`pos[3]` per vertex).
    pub edge_upload: Option<Vec<f32>>,
    pub edge_vertex_count: i32,
    pub show_edges: bool,
    pub edge_color: [f32; 4],
    /// Dynamic colored lines (`pos[3], rgba[4]` per vertex): grid, wake, lift vectors, hover outline.
    pub grid_lines: Vec<f32>,
    pub overlay_lines: Vec<f32>,
    pub hover_lines: Vec<f32>,
}

struct GlPrograms {
    mesh_program: glow::Program,
    mesh_u_mvp: Option<glow::UniformLocation>,
    mesh_u_view_dir: Option<glow::UniformLocation>,
    mesh_vao: glow::VertexArray,
    mesh_vbo: glow::Buffer,
    mesh_key: Option<MeshBufferKey>,
    mesh_vertex_count: i32,

    edge_program: glow::Program,
    edge_u_mvp: Option<glow::UniformLocation>,
    edge_u_color: Option<glow::UniformLocation>,
    edge_vao: glow::VertexArray,
    edge_vbo: glow::Buffer,
    edge_key: Option<u64>,
    edge_vertex_count: i32,

    line_program: glow::Program,
    line_u_mvp: Option<glow::UniformLocation>,
    line_vao: glow::VertexArray,
    line_vbo: glow::Buffer,
    present_program: glow::Program,
    present_texture: Option<glow::UniformLocation>,
}

struct GlTarget {
    fbo: glow::Framebuffer,
    color: glow::Texture,
    depth: glow::Renderbuffer,
    size: [i32; 2],
}

impl GlTarget {
    unsafe fn new(gl: &glow::Context, size: [i32; 2]) -> Result<Self, String> {
        drain_gl_errors(gl);
        let fbo = gl.create_framebuffer()?;
        let color = gl
            .create_texture()
            .inspect_err(|_| gl.delete_framebuffer(fbo))?;
        let depth = gl.create_renderbuffer().inspect_err(|_| {
            gl.delete_texture(color);
            gl.delete_framebuffer(fbo);
        })?;
        let target = Self {
            fbo,
            color,
            depth,
            size,
        };
        let [w, h] = size;
        gl.bind_framebuffer(glow::FRAMEBUFFER, Some(fbo));
        gl.bind_texture(glow::TEXTURE_2D, Some(color));
        gl.tex_parameter_i32(
            glow::TEXTURE_2D,
            glow::TEXTURE_MIN_FILTER,
            glow::NEAREST as i32,
        );
        gl.tex_parameter_i32(
            glow::TEXTURE_2D,
            glow::TEXTURE_MAG_FILTER,
            glow::NEAREST as i32,
        );
        gl.tex_image_2d(
            glow::TEXTURE_2D,
            0,
            glow::RGBA8 as i32,
            w,
            h,
            0,
            glow::RGBA,
            glow::UNSIGNED_BYTE,
            glow::PixelUnpackData::Slice(None),
        );
        gl.framebuffer_texture_2d(
            glow::FRAMEBUFFER,
            glow::COLOR_ATTACHMENT0,
            glow::TEXTURE_2D,
            Some(color),
            0,
        );
        gl.bind_renderbuffer(glow::RENDERBUFFER, Some(depth));
        gl.framebuffer_renderbuffer(
            glow::FRAMEBUFFER,
            glow::DEPTH_ATTACHMENT,
            glow::RENDERBUFFER,
            Some(depth),
        );

        let mut complete = false;
        for format in [glow::DEPTH_COMPONENT24, glow::DEPTH_COMPONENT16] {
            gl.renderbuffer_storage(glow::RENDERBUFFER, format, w, h);
            let error = gl.get_error();
            if error == glow::NO_ERROR
                && gl.check_framebuffer_status(glow::FRAMEBUFFER) == glow::FRAMEBUFFER_COMPLETE
            {
                complete = true;
                break;
            }
            drain_gl_errors(gl);
        }
        gl.bind_renderbuffer(glow::RENDERBUFFER, None);
        gl.bind_texture(glow::TEXTURE_2D, None);
        if complete {
            Ok(target)
        } else {
            target.destroy(gl);
            Err(format!(
                "could not create an OpenGL viewport framebuffer ({w} × {h})"
            ))
        }
    }

    unsafe fn destroy(self, gl: &glow::Context) {
        gl.delete_framebuffer(self.fbo);
        gl.delete_renderbuffer(self.depth);
        gl.delete_texture(self.color);
    }
}

/// App-owned OpenGL state for the 3D viewport. Destroyed explicitly in `App::on_exit`.
pub struct GlRenderer {
    programs: Option<GlPrograms>,
    target: Option<GlTarget>,
    failed: Option<String>,
    pub force_cpu_fallback: bool,
}

impl Default for GlRenderer {
    fn default() -> Self {
        let force_cpu = std::env::var("FOILEY_VIEW3D_CPU")
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false);
        Self {
            programs: None,
            target: None,
            failed: None,
            force_cpu_fallback: force_cpu,
        }
    }
}

impl GlRenderer {
    /// Whether hardware GL rendering should be attempted.
    pub fn can_use_gl(&self) -> bool {
        !self.force_cpu_fallback && self.failed.is_none()
    }

    /// Last recorded GL initialization failure, if any.
    pub fn failure_reason(&self) -> Option<&str> {
        self.failed.as_deref()
    }

    /// Check whether the current cached GPU mesh buffer already matches `key`.
    pub fn has_mesh_buffer(&self, key: MeshBufferKey) -> bool {
        self.programs.as_ref().and_then(|p| p.mesh_key) == Some(key)
    }

    /// Check whether the current cached GPU edge buffer already matches `mesh_key`.
    pub fn has_edge_buffer(&self, mesh_key: u64) -> bool {
        self.programs.as_ref().and_then(|p| p.edge_key) == Some(mesh_key)
    }

    /// Compile shaders before the first viewport frame.
    pub fn init(&mut self, gl: &glow::Context) {
        self.ensure_init(gl);
    }

    /// Lazily compile shaders & allocate VAOs/VBOs on the GL context thread.
    fn ensure_init(&mut self, gl: &glow::Context) -> bool {
        if self.force_cpu_fallback || self.failed.is_some() {
            return false;
        }
        if self.programs.is_some() {
            return true;
        }
        match unsafe { init_gl(gl) } {
            Ok(p) => {
                self.programs = Some(p);
                true
            }
            Err(e) => {
                self.failed = Some(e);
                false
            }
        }
    }

    /// Execute one 3D viewport draw inside `egui_glow::CallbackFn`.
    pub fn paint(&mut self, painter: &Painter, info: &PaintCallbackInfo, frame: &GlFrame) -> bool {
        let gl = painter.gl();
        if !self.ensure_init(gl) {
            return false;
        }
        let Some((viewport, visible)) = viewport_pixels(info) else {
            return true;
        };
        let [x, y, w, h] = visible;
        let destination = painter.intermediate_fbo();

        unsafe {
            if self
                .target
                .as_ref()
                .is_none_or(|target| target.size != [w, h])
            {
                if let Some(old) = self.target.take() {
                    old.destroy(gl);
                }
                match GlTarget::new(gl, [w, h]) {
                    Ok(target) => self.target = Some(target),
                    Err(error) => {
                        gl.bind_framebuffer(glow::FRAMEBUFFER, destination);
                        self.failed = Some(error);
                        return false;
                    }
                }
            }
            let target = self.target.as_ref().unwrap();
            let progs = self.programs.as_mut().unwrap();

            // Upload cached mesh/edge buffers only when invalidated.
            if let Some(ref data) = frame.tri_upload {
                gl.bind_buffer(glow::ARRAY_BUFFER, Some(progs.mesh_vbo));
                let bytes = pack_f32(data);
                gl.buffer_data_u8_slice(glow::ARRAY_BUFFER, &bytes, glow::STATIC_DRAW);
                progs.mesh_key = Some(frame.mesh_buf_key);
                progs.mesh_vertex_count = frame.tri_vertex_count;
            }
            if let Some(ref data) = frame.edge_upload {
                gl.bind_buffer(glow::ARRAY_BUFFER, Some(progs.edge_vbo));
                let bytes = pack_f32(data);
                gl.buffer_data_u8_slice(glow::ARRAY_BUFFER, &bytes, glow::STATIC_DRAW);
                progs.edge_key = Some(frame.edge_mesh_key);
                progs.edge_vertex_count = frame.edge_vertex_count;
            }

            // Render with the canvas projection even when only part of the canvas is visible.
            gl.bind_framebuffer(glow::FRAMEBUFFER, Some(target.fbo));
            gl.viewport(viewport[0] - x, viewport[1] - y, viewport[2], viewport[3]);
            gl.enable(glow::SCISSOR_TEST);
            gl.scissor(0, 0, w, h);

            gl.enable(glow::DEPTH_TEST);
            gl.depth_mask(true);
            gl.clear_depth(1.0);
            gl.clear_color(0.0, 0.0, 0.0, 0.0);
            gl.clear(glow::COLOR_BUFFER_BIT | glow::DEPTH_BUFFER_BIT);

            gl.disable(glow::CULL_FACE);
            gl.enable(glow::BLEND);
            gl.blend_func(glow::ONE, glow::ONE_MINUS_SRC_ALPHA);
            gl.line_width(1.0);

            // 1) Reference grid (depth tested, writes depth = false so it stays unobtrusive)
            if !frame.grid_lines.is_empty() {
                gl.depth_func(glow::LEQUAL);
                gl.depth_mask(false);
                draw_colored_lines(gl, progs, &frame.view_proj, &frame.grid_lines);
            }

            // 2) Filled panels with polygon offset pushing filled surfaces slightly away
            //    so perimeter edges drawn next never z-fight or bleed through from behind.
            if progs.mesh_vertex_count > 0 {
                gl.depth_func(glow::LESS);
                gl.depth_mask(true);
                gl.enable(glow::POLYGON_OFFSET_FILL);
                gl.polygon_offset(1.0, 1.0);

                gl.use_program(Some(progs.mesh_program));
                gl.uniform_matrix_4_f32_slice(progs.mesh_u_mvp.as_ref(), false, &frame.view_proj);
                gl.uniform_3_f32_slice(progs.mesh_u_view_dir.as_ref(), &frame.view_dir);
                gl.bind_vertex_array(Some(progs.mesh_vao));
                gl.draw_arrays(glow::TRIANGLES, 0, progs.mesh_vertex_count);

                gl.disable(glow::POLYGON_OFFSET_FILL);
            }

            // 3) Panel perimeter edges (depth test LEQUAL, depth write OFF)
            if frame.show_edges && progs.edge_vertex_count > 0 {
                gl.depth_func(glow::LEQUAL);
                gl.depth_mask(false);
                gl.use_program(Some(progs.edge_program));
                gl.uniform_matrix_4_f32_slice(progs.edge_u_mvp.as_ref(), false, &frame.view_proj);
                gl.uniform_4_f32_slice(progs.edge_u_color.as_ref(), &frame.edge_color);
                gl.bind_vertex_array(Some(progs.edge_vao));
                gl.draw_arrays(glow::LINES, 0, progs.edge_vertex_count);
            }

            // 4) Depth-tested 3D overlay lines (wake, lift vectors)
            if !frame.overlay_lines.is_empty() {
                gl.depth_func(glow::LEQUAL);
                gl.depth_mask(false);
                draw_colored_lines(gl, progs, &frame.view_proj, &frame.overlay_lines);
            }

            // 5) Hovered panel perimeter outline (drawn cleanly on top of the hovered panel)
            if !frame.hover_lines.is_empty() {
                gl.depth_func(glow::LEQUAL);
                gl.depth_mask(false);
                draw_colored_lines(gl, progs, &frame.view_proj, &frame.hover_lines);
            }

            // Composite the premultiplied scene over the egui background.
            gl.bind_framebuffer(glow::FRAMEBUFFER, destination);
            gl.viewport(x, y, w, h);
            gl.scissor(x, y, w, h);
            gl.disable(glow::DEPTH_TEST);
            gl.use_program(Some(progs.present_program));
            gl.active_texture(glow::TEXTURE0);
            gl.bind_texture(glow::TEXTURE_2D, Some(target.color));
            gl.uniform_1_i32(progs.present_texture.as_ref(), 0);
            gl.bind_vertex_array(Some(progs.mesh_vao));
            gl.draw_arrays(glow::TRIANGLES, 0, 3);

            // egui_glow restores the rest of its painting state after the callback.
            gl.depth_mask(true);
            gl.bind_texture(glow::TEXTURE_2D, None);
            gl.bind_vertex_array(None);
            gl.bind_buffer(glow::ARRAY_BUFFER, None);
            gl.use_program(None);
        }
        true
    }

    /// Release all allocated OpenGL objects. Called from `App::on_exit(Some(gl))`.
    pub fn destroy(&mut self, gl: &glow::Context) {
        if let Some(target) = self.target.take() {
            unsafe { target.destroy(gl) };
        }
        if let Some(p) = self.programs.take() {
            unsafe {
                gl.delete_program(p.mesh_program);
                gl.delete_vertex_array(p.mesh_vao);
                gl.delete_buffer(p.mesh_vbo);
                gl.delete_program(p.edge_program);
                gl.delete_vertex_array(p.edge_vao);
                gl.delete_buffer(p.edge_vbo);
                gl.delete_program(p.line_program);
                gl.delete_vertex_array(p.line_vao);
                gl.delete_buffer(p.line_vbo);
                gl.delete_program(p.present_program);
            }
        }
    }
}

unsafe fn draw_colored_lines(
    gl: &glow::Context,
    progs: &GlPrograms,
    view_proj: &[f32; 16],
    data: &[f32],
) {
    let vert_count = (data.len() / 7) as i32;
    if vert_count <= 0 {
        return;
    }
    gl.use_program(Some(progs.line_program));
    gl.uniform_matrix_4_f32_slice(progs.line_u_mvp.as_ref(), false, view_proj);
    gl.bind_vertex_array(Some(progs.line_vao));
    gl.bind_buffer(glow::ARRAY_BUFFER, Some(progs.line_vbo));
    let bytes = pack_f32(data);
    gl.buffer_data_u8_slice(glow::ARRAY_BUFFER, &bytes, glow::STREAM_DRAW);
    gl.draw_arrays(glow::LINES, 0, vert_count);
}

fn viewport_pixels(info: &PaintCallbackInfo) -> Option<([i32; 4], [i32; 4])> {
    let height = info.screen_size_px[1] as i32;
    let ppp = info.pixels_per_point;
    let pixels = |r: Rect| {
        let x = (r.min.x * ppp).round() as i32;
        let y = height - (r.max.y * ppp).round() as i32;
        let w = (r.max.x * ppp).round() as i32 - x;
        let h = height - (r.min.y * ppp).round() as i32 - y;
        [x, y, w, h]
    };
    let screen = Rect::from_min_size(
        eframe::egui::Pos2::ZERO,
        eframe::egui::vec2(info.screen_size_px[0] as f32 / ppp, height as f32 / ppp),
    );
    let viewport = pixels(info.viewport);
    let visible = pixels(info.viewport.intersect(info.clip_rect).intersect(screen));
    (visible[2] > 0 && visible[3] > 0).then_some((viewport, visible))
}

fn drain_gl_errors(gl: &glow::Context) {
    unsafe {
        for _ in 0..8 {
            if gl.get_error() == glow::NO_ERROR {
                break;
            }
        }
    }
}

unsafe fn init_gl(gl: &glow::Context) -> Result<GlPrograms, String> {
    let sv = ShaderVersion::get(gl);
    if !sv.is_new_shader_interface() {
        return Err(format!(
            "OpenGL shader version {sv:?} does not support GLSL 140 / ES 300"
        ));
    }
    let header = sv.version_declaration();
    let frag_prec = if sv.is_embedded() {
        "precision mediump float;\n"
    } else {
        ""
    };

    let present_vs = format!(
        "{header}
        out vec2 v_uv;
        void main() {{
            vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
            v_uv = p;
            gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
        }}"
    );
    let present_fs = format!(
        "{header}{frag_prec}
        in vec2 v_uv;
        uniform sampler2D u_scene;
        out vec4 o_color;
        void main() {{ o_color = texture(u_scene, v_uv); }}"
    );
    let present_program = link_program(gl, &present_vs, &present_fs, &[])?;
    let present_texture = gl.get_uniform_location(present_program, "u_scene");

    // Lit panel mesh shader (two-sided hemisphere + specular lighting)
    let mesh_vs = format!(
        "{header}
        in vec3 a_pos;
        in vec3 a_normal;
        in vec3 a_color;
        uniform mat4 u_mvp;
        out vec3 v_normal;
        out vec3 v_color;
        void main() {{
            v_normal = a_normal;
            v_color = a_color;
            gl_Position = u_mvp * vec4(a_pos, 1.0);
        }}"
    );
    let mesh_fs = format!(
        "{header}{frag_prec}
        in vec3 v_normal;
        in vec3 v_color;
        uniform vec3 u_view_dir;
        out vec4 o_color;
        void main() {{
            vec3 n = normalize(v_normal);
            vec3 v = normalize(u_view_dir);
            if (dot(n, v) < 0.0) {{
                n = -n;
            }}
            vec3 key_dir = normalize(vec3(-0.35, 0.45, 0.82));
            vec3 fill_dir = normalize(vec3(0.55, -0.50, 0.30));
            float diff_key = max(dot(n, key_dir), 0.0);
            float diff_fill = max(dot(n, fill_dir), 0.0);
            float hemi = 0.5 + 0.5 * n.z;
            float ambient = mix(0.32, 0.52, hemi);
            float shade = ambient + 0.44 * diff_key + 0.14 * diff_fill;
            vec3 half_vec = normalize(key_dir + v);
            float spec = pow(max(dot(n, half_vec), 0.0), 28.0) * 0.08;
            vec3 col = clamp(v_color * shade + vec3(spec), 0.0, 1.0);
            o_color = vec4(col, 1.0);
        }}"
    );
    let mesh_program = link_program(
        gl,
        &mesh_vs,
        &mesh_fs,
        &[(0, "a_pos"), (1, "a_normal"), (2, "a_color")],
    )?;
    let mesh_u_mvp = gl.get_uniform_location(mesh_program, "u_mvp");
    let mesh_u_view_dir = gl.get_uniform_location(mesh_program, "u_view_dir");
    let mesh_vao = gl
        .create_vertex_array()
        .map_err(|e| format!("create mesh vao: {e}"))?;
    let mesh_vbo = gl
        .create_buffer()
        .map_err(|e| format!("create mesh vbo: {e}"))?;
    gl.bind_vertex_array(Some(mesh_vao));
    gl.bind_buffer(glow::ARRAY_BUFFER, Some(mesh_vbo));
    let stride = 9 * 4;
    gl.enable_vertex_attrib_array(0);
    gl.vertex_attrib_pointer_f32(0, 3, glow::FLOAT, false, stride, 0);
    gl.enable_vertex_attrib_array(1);
    gl.vertex_attrib_pointer_f32(1, 3, glow::FLOAT, false, stride, 3 * 4);
    gl.enable_vertex_attrib_array(2);
    gl.vertex_attrib_pointer_f32(2, 3, glow::FLOAT, false, stride, 6 * 4);

    // 2) Uniform-color edge shader
    let edge_vs = format!(
        "{header}
        in vec3 a_pos;
        uniform mat4 u_mvp;
        void main() {{
            gl_Position = u_mvp * vec4(a_pos, 1.0);
        }}"
    );
    let edge_fs = format!(
        "{header}{frag_prec}
        uniform vec4 u_color;
        out vec4 o_color;
        void main() {{
            o_color = u_color;
        }}"
    );
    let edge_program = link_program(gl, &edge_vs, &edge_fs, &[(0, "a_pos")])?;
    let edge_u_mvp = gl.get_uniform_location(edge_program, "u_mvp");
    let edge_u_color = gl.get_uniform_location(edge_program, "u_color");
    let edge_vao = gl
        .create_vertex_array()
        .map_err(|e| format!("create edge vao: {e}"))?;
    let edge_vbo = gl
        .create_buffer()
        .map_err(|e| format!("create edge vbo: {e}"))?;
    gl.bind_vertex_array(Some(edge_vao));
    gl.bind_buffer(glow::ARRAY_BUFFER, Some(edge_vbo));
    gl.enable_vertex_attrib_array(0);
    gl.vertex_attrib_pointer_f32(0, 3, glow::FLOAT, false, 3 * 4, 0);

    // 3) Per-vertex RGBA line shader (grid, wake, lift vectors, hover outline)
    let line_vs = format!(
        "{header}
        in vec3 a_pos;
        in vec4 a_color;
        uniform mat4 u_mvp;
        out vec4 v_color;
        void main() {{
            v_color = a_color;
            gl_Position = u_mvp * vec4(a_pos, 1.0);
        }}"
    );
    let line_fs = format!(
        "{header}{frag_prec}
        in vec4 v_color;
        out vec4 o_color;
        void main() {{
            o_color = v_color;
        }}"
    );
    let line_program = link_program(gl, &line_vs, &line_fs, &[(0, "a_pos"), (1, "a_color")])?;
    let line_u_mvp = gl.get_uniform_location(line_program, "u_mvp");
    let line_vao = gl
        .create_vertex_array()
        .map_err(|e| format!("create line vao: {e}"))?;
    let line_vbo = gl
        .create_buffer()
        .map_err(|e| format!("create line vbo: {e}"))?;
    gl.bind_vertex_array(Some(line_vao));
    gl.bind_buffer(glow::ARRAY_BUFFER, Some(line_vbo));
    let lstride = 7 * 4;
    gl.enable_vertex_attrib_array(0);
    gl.vertex_attrib_pointer_f32(0, 3, glow::FLOAT, false, lstride, 0);
    gl.enable_vertex_attrib_array(1);
    gl.vertex_attrib_pointer_f32(1, 4, glow::FLOAT, false, lstride, 3 * 4);

    gl.bind_vertex_array(None);
    gl.bind_buffer(glow::ARRAY_BUFFER, None);
    drain_gl_errors(gl);

    Ok(GlPrograms {
        mesh_program,
        mesh_u_mvp,
        mesh_u_view_dir,
        mesh_vao,
        mesh_vbo,
        mesh_key: None,
        mesh_vertex_count: 0,
        edge_program,
        edge_u_mvp,
        edge_u_color,
        edge_vao,
        edge_vbo,
        edge_key: None,
        edge_vertex_count: 0,
        line_program,
        line_u_mvp,
        line_vao,
        line_vbo,
        present_program,
        present_texture,
    })
}

unsafe fn link_program(
    gl: &glow::Context,
    vs_src: &str,
    fs_src: &str,
    attribs: &[(u32, &str)],
) -> Result<glow::Program, String> {
    let vs = compile_shader(gl, glow::VERTEX_SHADER, vs_src)?;
    let fs = match compile_shader(gl, glow::FRAGMENT_SHADER, fs_src) {
        Ok(s) => s,
        Err(e) => {
            gl.delete_shader(vs);
            return Err(e);
        }
    };
    let prog = gl
        .create_program()
        .map_err(|e| format!("create program: {e}"))?;
    gl.attach_shader(prog, vs);
    gl.attach_shader(prog, fs);
    for &(loc, name) in attribs {
        gl.bind_attrib_location(prog, loc, name);
    }
    gl.link_program(prog);
    gl.detach_shader(prog, vs);
    gl.detach_shader(prog, fs);
    gl.delete_shader(vs);
    gl.delete_shader(fs);
    if !gl.get_program_link_status(prog) {
        let log = gl.get_program_info_log(prog);
        gl.delete_program(prog);
        return Err(format!("GL program link failed: {log}"));
    }
    Ok(prog)
}

unsafe fn compile_shader(gl: &glow::Context, kind: u32, src: &str) -> Result<glow::Shader, String> {
    let sh = gl
        .create_shader(kind)
        .map_err(|e| format!("create shader: {e}"))?;
    gl.shader_source(sh, src);
    gl.compile_shader(sh);
    if !gl.get_shader_compile_status(sh) {
        let log = gl.get_shader_info_log(sh);
        gl.delete_shader(sh);
        return Err(format!("GL shader compile failed: {log}"));
    }
    Ok(sh)
}

fn pack_f32(data: &[f32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(std::mem::size_of_val(data));
    for &v in data {
        out.extend_from_slice(&v.to_ne_bytes());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use eframe::egui::pos2;

    #[test]
    fn clipped_gl_canvas_keeps_its_projection() {
        let info = PaintCallbackInfo {
            viewport: Rect::from_min_max(pos2(500.0, 200.0), pos2(1020.0, 670.0)),
            clip_rect: Rect::from_min_max(pos2(520.0, 0.0), pos2(860.0, 625.0)),
            pixels_per_point: 1.0,
            screen_size_px: [900, 600],
        };
        assert_eq!(
            viewport_pixels(&info),
            Some(([500, -70, 520, 470], [520, 0, 340, 400]))
        );
    }
}
