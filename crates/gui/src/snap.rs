//! Run Studio and rasterize its egui output without a desktop window.

mod output;
mod script;

use crate::view3d::{self, NavPreset, ViewportSnapshot};
use crate::{App, Launch, Module};
use eframe::egui::{self, Event, Key, Modifiers, PointerButton, Pos2, Rect, Vec2};
use egui_kittest::{kittest::Queryable, wgpu::WgpuTestRenderer, Harness};
use egui_wgpu::{
    wgpu, RenderState, RendererOptions, WgpuConfiguration, WgpuSetup, WgpuSetupExisting,
};
use foiley_app::JobStatus;
use foiley_core::progress::LogLevel;
use image::RgbaImage;
use script::{Instruction, Position, Step};
use serde::Serialize;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

pub use script::Crop;

#[derive(Clone, Copy, Debug, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Preset {
    Views,
    Orbit,
    Nav,
    Analysis,
    Narrow,
}

pub struct Options {
    pub launch: Launch,
    pub size: [u32; 2],
    pub scale: f32,
    pub light: bool,
    pub nav: Option<String>,
    pub perspective: Option<bool>,
    pub grid: Option<bool>,
    pub out: Option<PathBuf>,
    pub script: Option<PathBuf>,
    pub steps: Vec<String>,
    pub preset: Option<Preset>,
    pub gif: bool,
    pub crop: Crop,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            launch: Launch::default(),
            size: [1400, 900],
            scale: 1.0,
            light: false,
            nav: None,
            perspective: None,
            grid: None,
            out: None,
            script: None,
            steps: Vec::new(),
            preset: None,
            gif: false,
            crop: Crop::Full,
        }
    }
}

#[derive(Debug)]
pub enum SnapError {
    Input(String),
    Settle(String),
    Adapter(String),
}

impl SnapError {
    pub fn code(&self) -> i32 {
        match self {
            Self::Input(_) => 1,
            Self::Settle(_) => 2,
            Self::Adapter(_) => 3,
        }
    }

    fn output(error: impl std::fmt::Display) -> Self {
        Self::Input(format!("could not write snap output: {error}"))
    }
}

impl std::fmt::Display for SnapError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Input(reason) | Self::Settle(reason) | Self::Adapter(reason) => {
                f.write_str(reason)
            }
        }
    }
}

#[derive(Serialize)]
struct Capture {
    frame: u64,
    file: String,
    crop: Crop,
    viewport: Option<ViewportSnapshot>,
    plane: Option<String>,
    selected_point: Option<usize>,
    control: Option<f64>,
}

#[derive(Serialize)]
struct StepRecord {
    frame: u64,
    step: String,
}

struct Trace {
    name: String,
    captures: Vec<usize>,
}

struct Studio {
    app: App,
    ctx: egui::Context,
}

struct Runner<'a> {
    harness: Harness<'a, Option<Studio>>,
    out: PathBuf,
    frame: u64,
    steps: Vec<StepRecord>,
    captures: Vec<Capture>,
    files: Vec<PathBuf>,
    trace: Option<Trace>,
    crop: Crop,
    scale: f32,
    gif: bool,
    snaps: usize,
}

const FRAME_TIME: f64 = 1.0 / 60.0;

fn module(name: &str) -> Result<Module, String> {
    match name {
        "foil" => Ok(Module::FoilDesign),
        "2d" => Ok(Module::FoilAnalysis),
        "plane" => Ok(Module::PlaneDesign),
        "3d" => Ok(Module::PlaneAnalysis),
        "stability" => Ok(Module::Stability),
        _ => Err(format!("unknown module {name:?}")),
    }
}

fn module_name(module: Module) -> &'static str {
    match module {
        Module::FoilDesign => "foil",
        Module::FoilAnalysis => "2d",
        Module::PlaneDesign => "plane",
        Module::PlaneAnalysis => "3d",
        Module::Stability => "stability",
    }
}

fn renderer() -> Result<(WgpuTestRenderer, String), SnapError> {
    let mut setup = egui_kittest::wgpu::default_wgpu_setup();
    let WgpuSetup::CreateNew(config) = &setup else {
        unreachable!("default kittest setup")
    };
    let software = std::env::var_os("FOILEY_SNAP_SOFTWARE").is_some_and(|value| value == "1");
    let mut descriptor = config.instance_descriptor.clone();
    if cfg!(windows) && software {
        descriptor.backends = wgpu::Backends::DX12;
    }
    let instance = wgpu::Instance::new(&descriptor);
    if cfg!(windows) || software {
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            force_fallback_adapter: true,
            ..Default::default()
        }));
        match adapter {
            Ok(adapter) => {
                let (device, queue) = pollster::block_on(
                    adapter.request_device(&(config.device_descriptor)(&adapter)),
                )
                .map_err(|error| {
                    SnapError::Adapter(format!(
                        "could not open the software render adapter: {error}"
                    ))
                })?;
                setup = WgpuSetup::Existing(WgpuSetupExisting {
                    instance: instance.clone(),
                    adapter,
                    device,
                    queue,
                });
            }
            Err(error) if software => {
                return Err(SnapError::Adapter(format!(
                    "no software wgpu adapter: {error}"
                )))
            }
            Err(_) => {}
        }
    }
    let configuration = WgpuConfiguration {
        wgpu_setup: setup,
        ..Default::default()
    };
    let state = pollster::block_on(RenderState::create(
        &configuration,
        &instance,
        None,
        RendererOptions::PREDICTABLE,
    ))
    .map_err(|error| SnapError::Adapter(format!("no usable wgpu render adapter: {error}")))?;
    let info = state.adapter.get_info();
    let adapter = format!("{} ({:?}, {:?})", info.name, info.backend, info.device_type);
    Ok((WgpuTestRenderer::from_render_state(state), adapter))
}

/// Return absolute image paths in creation order, followed by the manifest.
pub fn run(mut options: Options) -> Result<Vec<PathBuf>, SnapError> {
    let started = Instant::now();
    let instructions = if let Some(path) = &options.script {
        let contents = std::fs::read_to_string(path)
            .map_err(|error| SnapError::Input(format!("{}: {error}", path.display())))?;
        script::parse(&contents)
    } else {
        script::parse(&options.steps.join("\n"))
    }
    .map_err(SnapError::Input)?;
    if options.script.is_some() && !options.steps.is_empty() {
        return Err(SnapError::Input("use --script or --do".into()));
    }
    options.launch.module = match options.preset {
        Some(Preset::Analysis) => {
            options.launch.run = true;
            Some(Module::PlaneAnalysis)
        }
        Some(Preset::Narrow) => {
            options.size = [900, 600];
            Some(Module::PlaneDesign)
        }
        Some(_) => Some(Module::PlaneDesign),
        None => options.launch.module,
    };
    if !options.scale.is_finite()
        || options.scale <= 0.0
        || options.size.iter().any(|&size| size == 0 || size > 8192)
    {
        return Err(SnapError::Input(
            "use a positive scale and a size from 1 to 8192 pixels per dimension".into(),
        ));
    }
    let nav = match options.nav.as_deref() {
        Some("standard") | None => NavPreset::Standard,
        Some("cad") => NavPreset::CadBlender,
        Some(nav) => {
            return Err(SnapError::Input(format!(
                "unknown navigation preset {nav:?}"
            )))
        }
    };
    let (renderer, adapter) = renderer()?;
    eprintln!("foiley snap uses the CPU 3D viewport with wgpu ({adapter}).");
    let out = match options.out.clone() {
        Some(out) => {
            std::fs::create_dir_all(&out).map_err(SnapError::output)?;
            out
        }
        None => new_output_dir(Path::new("target").join("snap"))?,
    };
    let out = std::path::absolute(out).map_err(SnapError::output)?;
    let arguments = serde_json::json!({
        "files": options.launch.files,
        "module": options.launch.module.map(module_name),
        "run": options.launch.run,
        "threads": options.launch.threads,
        "size": options.size,
        "scale": options.scale,
        "theme": if options.light { "light" } else { "dark" },
        "out": out,
        "script": options.script,
        "do": options.steps,
        "preset": options.preset,
        "gif": options.gif || matches!(options.preset, Some(Preset::Orbit)),
        "crop": options.crop,
        "nav": if matches!(options.preset, Some(Preset::Orbit | Preset::Nav)) { "cad" } else { nav.as_str() },
        "perspective": options.perspective.unwrap_or(true),
        "grid": options.grid.unwrap_or(true),
    });
    let mut launch = Some(options.launch);
    let theme = if options.light {
        egui::Theme::Light
    } else {
        egui::Theme::Dark
    };
    let size = egui::vec2(options.size[0] as f32, options.size[1] as f32) / options.scale;
    let mut frame = eframe::Frame::_new_kittest();
    let mut harness = Harness::builder()
        .with_size(size)
        .with_pixels_per_point(options.scale)
        .with_theme(theme)
        .with_step_dt(FRAME_TIME as f32)
        .with_max_steps(1)
        .renderer(renderer)
        .build_state(
            |ctx, state: &mut Option<Studio>| {
                let studio = state.get_or_insert_with(|| {
                    let mut app = App::new(
                        &eframe::CreationContext::_new_kittest(ctx.clone()),
                        launch.take().unwrap(),
                    );
                    app.headless = true;
                    app.dark = !options.light;
                    app.view3d_gl.lock().unwrap().force_cpu_fallback = true;
                    app.view3d_prefs.nav_preset = nav;
                    app.view3d_prefs.perspective = options.perspective.unwrap_or(true);
                    app.view3d_prefs.show_grid = options.grid.unwrap_or(true);
                    let mut visuals = if options.light {
                        egui::Visuals::light()
                    } else {
                        egui::Visuals::dark()
                    };
                    visuals.text_cursor.blink = false;
                    ctx.set_visuals(visuals);
                    ctx.data_mut(|data| data.insert_temp(egui::Id::new("foiley_snap"), true));
                    Studio {
                        app,
                        ctx: ctx.clone(),
                    }
                });
                ctx.data_mut(|data| {
                    for salt in ["plane_design", "plane_analysis"] {
                        data.remove::<ViewportSnapshot>(view3d::snapshot_id(salt));
                    }
                });
                eframe::App::update(&mut studio.app, ctx, &mut frame);
            },
            None,
        );
    // The harness creates and lays out the widgets before accepting input.
    let layout_frames = harness
        .state()
        .as_ref()
        .unwrap()
        .ctx
        .input(|i| (i.time / FRAME_TIME).round() as u64);
    harness.input_mut().focused = true;
    let mut runner = Runner {
        harness,
        out,
        frame: layout_frames,
        steps: Vec::new(),
        captures: Vec::new(),
        files: Vec::new(),
        trace: None,
        crop: options.crop,
        scale: options.scale,
        gif: options.gif,
        snaps: 0,
    };
    let result = runner.execute(options.preset, instructions);
    let studio = runner.studio();
    let warnings: Vec<_> = studio
        .app
        .log
        .iter()
        .filter(|line| matches!(line.level, LogLevel::Warn | LogLevel::Error))
        .map(|line| &line.text)
        .collect();
    let jobs: Vec<_> = studio
        .app
        .jobs
        .jobs()
        .iter()
        .map(|job| {
            serde_json::json!({
                "id": job.id, "title": job.title, "kind": job.kind,
                "status": format!("{:?}", job.status()), "seconds": job.elapsed().as_secs_f64()
            })
        })
        .collect();
    let manifest = serde_json::json!({
        "arguments": arguments, "adapter": adapter, "viewport_renderer": "cpu",
        "steps": runner.steps, "captures": runner.captures, "warnings": warnings,
        "jobs": jobs, "seconds": started.elapsed().as_secs_f64(),
        "error": result.as_ref().err().map(ToString::to_string),
    });
    let manifest_path = runner.out.join("manifest.json");
    let body = serde_json::to_vec_pretty(&manifest).map_err(SnapError::output)?;
    std::fs::write(&manifest_path, body).map_err(SnapError::output)?;
    result?;
    runner.files.push(manifest_path);
    Ok(runner.files)
}

/// A new directory `<parent>/<UTC timestamp>`.  It is created here, so two
/// runs started in the same second never share one: the second one gets
/// `<timestamp>-2`, and so on.
fn new_output_dir(parent: PathBuf) -> Result<PathBuf, SnapError> {
    let time = foiley_core::provenance::now_rfc3339();
    let stamp = format!(
        "{}-{}",
        time[..10].replace('-', ""),
        time[11..19].replace(':', "")
    );
    std::fs::create_dir_all(&parent).map_err(SnapError::output)?;
    let mut n = 1;
    loop {
        let dir = if n == 1 {
            parent.join(&stamp)
        } else {
            parent.join(format!("{stamp}-{n}"))
        };
        match std::fs::create_dir(&dir) {
            Ok(()) => return Ok(dir),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => n += 1,
            Err(error) => return Err(SnapError::output(format!("{}: {error}", dir.display()))),
        }
    }
}

impl Runner<'_> {
    fn studio(&self) -> &Studio {
        self.harness.state().as_ref().unwrap()
    }

    fn app_mut(&mut self) -> &mut App {
        &mut self.harness.state_mut().as_mut().unwrap().app
    }

    fn viewport(&self) -> Option<ViewportSnapshot> {
        let studio = self.studio();
        let salt = match studio.app.module {
            Module::PlaneDesign => "plane_design",
            Module::PlaneAnalysis => "plane_analysis",
            _ => return None,
        };
        studio
            .ctx
            .data(|data| data.get_temp(view3d::snapshot_id(salt)))
    }

    fn rect(&self) -> Option<Rect> {
        self.viewport().map(|viewport| {
            let [left, top, right, bottom] = viewport.rect;
            Rect::from_min_max(Pos2::new(left, top), Pos2::new(right, bottom))
        })
    }

    fn resolve(&self, pos: Position) -> Result<Pos2, SnapError> {
        pos.resolve(self.rect(), self.scale)
            .map_err(SnapError::Input)
    }

    fn step_frame(
        &mut self,
        modifiers: Modifiers,
        events: impl IntoIterator<Item = Event>,
    ) -> Result<(), SnapError> {
        self.frame += 1;
        let input = self.harness.input_mut();
        input.time = Some(self.frame as f64 * FRAME_TIME);
        input.modifiers = modifiers;
        input.events.extend(events);
        self.harness.step();
        if let Some(trace) = &mut self.trace {
            let file = format!("{}-{:03}.png", trace.name, trace.captures.len());
            trace.captures.push(self.captures.len());
            self.capture(file, self.crop)?;
        }
        Ok(())
    }

    fn settle(&mut self, seconds: f32) -> Result<(), SnapError> {
        let deadline = Instant::now() + Duration::from_secs_f32(seconds);
        let mut idle_frames = 0;
        loop {
            self.step_frame(Modifiers::NONE, [])?;
            let app = &self.studio().app;
            for job in app.jobs.jobs() {
                if let JobStatus::Failed(error) = job.status() {
                    return Err(SnapError::Settle(format!("{} failed: {error}", job.title)));
                }
            }
            let idle = app.jobs.running() == 0 && app.pending.is_empty() && !app.io_busy;
            if idle {
                let repaint = self
                    .harness
                    .output()
                    .viewport_output
                    .values()
                    .any(|viewport| viewport.repaint_delay != Duration::MAX);
                if !repaint {
                    return Ok(());
                }
                idle_frames += 1;
                if idle_frames == 600 {
                    return Err(SnapError::Settle(
                        "Studio still requests a repaint after 600 idle frames".into(),
                    ));
                }
            } else {
                if Instant::now() >= deadline {
                    let jobs = app
                        .jobs
                        .jobs()
                        .iter()
                        .filter(|job| !job.status().is_finished())
                        .map(|job| job.title.as_str())
                        .collect::<Vec<_>>();
                    return Err(SnapError::Settle(format!(
                        "settle timed out after {seconds}s; pending jobs: {}",
                        jobs.join(", ")
                    )));
                }
                std::thread::sleep(Duration::from_millis(10));
            }
        }
    }

    fn execute(
        &mut self,
        preset: Option<Preset>,
        instructions: Vec<Instruction>,
    ) -> Result<(), SnapError> {
        if let Some(error) = self
            .studio()
            .app
            .log
            .iter()
            .find(|line| matches!(line.level, LogLevel::Error))
        {
            return Err(SnapError::Input(error.text.clone()));
        }
        if let Some(preset) = preset {
            self.preset(preset)?;
        } else if instructions.is_empty() {
            return self.run_steps("settle; snap studio");
        }
        for instruction in instructions {
            self.instruction(instruction)?;
        }
        if self.trace.is_some() {
            return Err(SnapError::Input("trace is missing end".into()));
        }
        Ok(())
    }

    fn run_steps(&mut self, text: &str) -> Result<(), SnapError> {
        for instruction in script::parse(text).map_err(SnapError::Input)? {
            self.instruction(instruction)?;
        }
        Ok(())
    }

    fn instruction(&mut self, instruction: Instruction) -> Result<(), SnapError> {
        self.steps.push(StepRecord {
            frame: self.frame,
            step: instruction.source,
        });
        match instruction.step {
            Step::Settle(seconds) => self.settle(seconds)?,
            Step::Wait(frames) => {
                for _ in 0..frames {
                    self.step_frame(Modifiers::NONE, [])?;
                }
            }
            Step::Module(module) => {
                self.app_mut().module = module;
                self.step_frame(Modifiers::NONE, [])?;
            }
            Step::Snap(name, crop) => {
                self.snaps += 1;
                self.capture(
                    format!("{:03}-{name}.png", self.snaps),
                    crop.unwrap_or(self.crop),
                )?;
            }
            Step::Key(physical, mods) => {
                let key = if physical == Key::Num1 && mods.shift {
                    Key::Exclamationmark
                } else {
                    physical
                };
                let event = |pressed| Event::Key {
                    key,
                    physical_key: Some(physical),
                    pressed,
                    repeat: false,
                    modifiers: mods,
                };
                self.step_frame(mods, [event(true)])?;
                self.step_frame(Modifiers::NONE, [event(false)])?;
            }
            Step::Text(text) => self.step_frame(Modifiers::NONE, [Event::Text(text)])?,
            Step::Click(pos, button) => {
                let pos = self.resolve(pos)?;
                self.step_frame(
                    Modifiers::NONE,
                    [
                        Event::PointerMoved(pos),
                        mouse(pos, button, true, Modifiers::NONE),
                    ],
                )?;
                self.step_frame(
                    Modifiers::NONE,
                    [mouse(pos, button, false, Modifiers::NONE)],
                )?;
            }
            Step::Label(label) => {
                let matches: Vec<_> = self.harness.query_all_by_label(&label).collect();
                match matches.as_slice() {
                    [node] => node.click_accesskit(),
                    _ => {
                        return Err(SnapError::Input(format!(
                            "click {label:?} matched {} widgets",
                            matches.len()
                        )))
                    }
                };
                self.step_frame(Modifiers::NONE, [])?;
            }
            Step::Hover(pos) => self.hover(self.resolve(pos)?)?,
            Step::Drag(mods, button, from, delta, frames) => {
                let from = self.resolve(from)?;
                self.step_frame(
                    mods,
                    [Event::PointerMoved(from), mouse(from, button, true, mods)],
                )?;
                let to = self.motion(from, delta / self.scale, frames, mods)?;
                self.step_frame(Modifiers::NONE, [mouse(to, button, false, mods)])?;
            }
            Step::Scroll(dy, pos) => {
                if let Some(pos) = pos {
                    self.hover(self.resolve(pos)?)?;
                }
                self.step_frame(
                    Modifiers::NONE,
                    [Event::MouseWheel {
                        unit: egui::MouseWheelUnit::Point,
                        delta: egui::vec2(0.0, dy / self.scale),
                        modifiers: Modifiers::NONE,
                    }],
                )?;
            }
            Step::Pinch(factor, pos) => {
                if let Some(pos) = pos {
                    self.hover(self.resolve(pos)?)?;
                }
                self.step_frame(Modifiers::NONE, [Event::Zoom(factor)])?;
            }
            Step::Trace(name) => {
                if self.trace.is_some()
                    || self
                        .files
                        .iter()
                        .any(|file| file.ends_with(format!("{name}-sheet.png")))
                {
                    return Err(SnapError::Input(format!(
                        "trace {name:?} is nested or already used"
                    )));
                }
                self.trace = Some(Trace {
                    name,
                    captures: Vec::new(),
                });
            }
            Step::End => self.end_trace(self.gif)?,
        }
        Ok(())
    }

    fn hover(&mut self, pos: Pos2) -> Result<(), SnapError> {
        self.step_frame(Modifiers::NONE, [Event::PointerMoved(pos)])
    }

    fn motion(
        &mut self,
        from: Pos2,
        delta: Vec2,
        count: u32,
        mods: Modifiers,
    ) -> Result<Pos2, SnapError> {
        for index in 1..=count {
            let pos = from + delta * (index as f32 / count as f32);
            self.step_frame(mods, [Event::PointerMoved(pos)])?;
        }
        Ok(from + delta)
    }

    fn capture(&mut self, file: String, crop: Crop) -> Result<(), SnapError> {
        let image = self.rasterize(crop)?;
        let path = self.out.join(&file);
        image.save(&path).map_err(SnapError::output)?;
        let app = &self.studio().app;
        let polar = app.pui.sel_polar.as_ref().and_then(|(plane, name)| {
            app.project
                .plane_polars
                .iter()
                .find(|p| &p.plane == plane && &p.name == name)
        });
        let point = polar.and_then(|polar| polar.points.get(app.pui.sel_point));
        self.captures.push(Capture {
            frame: self.frame,
            file,
            crop,
            viewport: self.viewport(),
            plane: app
                .project
                .planes
                .get(app.sel_plane)
                .map(|p| p.name.clone()),
            selected_point: point.map(|_| app.pui.sel_point),
            control: point.map(|p| p.ctrl),
        });
        self.files.push(path);
        Ok(())
    }

    fn rasterize(&mut self, crop: Crop) -> Result<RgbaImage, SnapError> {
        let image = self
            .harness
            .render()
            .map_err(|error| SnapError::Adapter(format!("wgpu render failed: {error}")))?;
        if crop == Crop::Full {
            return Ok(image);
        }
        let rect = self
            .rect()
            .ok_or_else(|| SnapError::Input("no 3D viewport to crop".into()))?;
        let x = (rect.left() * self.scale)
            .round()
            .clamp(0.0, image.width() as f32) as u32;
        let y = (rect.top() * self.scale)
            .round()
            .clamp(0.0, image.height() as f32) as u32;
        let right = (rect.right() * self.scale)
            .round()
            .clamp(x as f32, image.width() as f32) as u32;
        let bottom = (rect.bottom() * self.scale)
            .round()
            .clamp(y as f32, image.height() as f32) as u32;
        if x == right || y == bottom {
            return Err(SnapError::Input(
                "the 3D viewport is outside the image".into(),
            ));
        }
        Ok(image::imageops::crop_imm(&image, x, y, right - x, bottom - y).to_image())
    }

    fn end_trace(&mut self, gif: bool) -> Result<(), SnapError> {
        let trace = self
            .trace
            .take()
            .ok_or_else(|| SnapError::Input("end without trace".into()))?;
        let frames: Vec<_> = trace
            .captures
            .into_iter()
            .map(|index| &self.captures[index])
            .collect();
        let file = format!("{}-sheet.png", trace.name);
        output::sheet(
            &self.out,
            &frames,
            &file,
            frames.len().clamp(1, 4) as u32,
            16,
        )?;
        self.files.push(self.out.join(file));
        if gif {
            let file = format!("{}.gif", trace.name);
            output::animation(&self.out, &frames, &file)?;
            self.files.push(self.out.join(file));
        }
        Ok(())
    }

    fn preset(&mut self, preset: Preset) -> Result<(), SnapError> {
        match preset {
            Preset::Views => {
                self.app_mut().view3d_prefs.perspective = true;
                self.app_mut().pui.view_thick = true;
                self.run_steps("settle; hover vp:0.5,0.5")?;
                for (suffix, toggle) in [("persp", ""), ("ortho", "key 5;")] {
                    self.run_steps(toggle)?;
                    for (key, name) in [("0", "iso"), ("7", "top"), ("1", "front"), ("3", "side")] {
                        self.run_steps(&format!(
                            "key {key}; settle; key f; settle; snap {name}-{suffix} viewport"
                        ))?;
                    }
                }
                let frames: Vec<_> = self.captures.iter().collect();
                output::sheet(&self.out, &frames, "views-sheet.png", 4, 8)?;
                self.files.push(self.out.join("views-sheet.png"));
            }
            Preset::Orbit => {
                self.app_mut().view3d_prefs.nav_preset = NavPreset::CadBlender;
                self.run_steps("settle")?;
                let from = self.resolve(Position::Viewport(0.5, 0.5))?;
                let button = PointerButton::Middle;
                self.step_frame(
                    Modifiers::NONE,
                    [
                        Event::PointerMoved(from),
                        mouse(from, button, true, Modifiers::NONE),
                    ],
                )?;
                self.run_steps("trace orbit")?;
                self.steps.push(StepRecord {
                    frame: self.frame,
                    step: "orbit 360 degrees, 24 middle-drag frames".into(),
                });
                let to = self.motion(
                    from,
                    Vec2::new(std::f32::consts::TAU / 0.0075, 0.0),
                    24,
                    Modifiers::NONE,
                )?;
                self.steps.push(StepRecord {
                    frame: self.frame,
                    step: "end".into(),
                });
                self.end_trace(true)?;
                self.step_frame(Modifiers::NONE, [mouse(to, button, false, Modifiers::NONE)])?;
            }
            Preset::Nav => {
                self.app_mut().view3d_prefs.nav_preset = NavPreset::CadBlender;
                self.run_steps(concat!(
                    "settle; hover vp:0.5,0.5; trace nav; ",
                    "drag shift+mmb vp:0.5,0.5 90 30; scroll 170 vp:0.65,0.55; ",
                    "wait 12; drag ctrl+mmb vp:0.65,0.55 0 -45; key f; wait 14; end"
                ))?;
            }
            Preset::Analysis => {
                self.run_steps("settle; snap analysis viewport")?;
                if self.viewport().is_none_or(|view| view.field.is_none()) {
                    return Err(SnapError::Settle(
                        "the analysis viewport has no Cp overlay".into(),
                    ));
                }
                let point = self.midspan_on_screen()?;
                self.steps.push(StepRecord {
                    frame: self.frame,
                    step: "hover mid-span panel".into(),
                });
                self.hover(point)?;
                self.run_steps("snap analysis-hover viewport")?;
                if self
                    .viewport()
                    .is_none_or(|view| view.hovered_panel.is_none())
                {
                    return Err(SnapError::Settle("no panel picked at mid-span".into()));
                }
            }
            Preset::Narrow => {
                self.app_mut().view3d_prefs.perspective = false;
                self.app_mut().view3d_prefs.show_grid = false;
                self.run_steps("settle; hover vp:0.5,0.5; key 7; settle")?;
                let center = self.resolve(Position::Viewport(0.5, 0.5))?;
                let from = self.midspan_on_screen()?;
                let delta = center - from;
                self.run_steps(&format!(
                    "drag rmb vp:0.5,0.5 {} {} 10; hover vp:0.5,0.5",
                    delta.x * self.scale,
                    delta.y * self.scale
                ))?;
                let hit = self.viewport().and_then(|view| view.hovered_panel);
                let image = self.rasterize(Crop::Full)?;
                let x = (center.x * self.scale) as u32;
                let y = (center.y * self.scale) as u32;
                let pixel = image.get_pixel_checked(x, y).ok_or_else(|| {
                    SnapError::Settle("narrow viewport center is outside the image".into())
                })?;
                if hit.is_none() || pixel[2].saturating_sub(pixel[0]) < 18 {
                    return Err(SnapError::Settle(
                        "the CPU viewport did not pick a panel drawn at the cursor".into(),
                    ));
                }
                self.run_steps("snap narrow viewport")?;
            }
        }
        Ok(())
    }

    fn midspan_on_screen(&self) -> Result<Pos2, SnapError> {
        let app = &self.studio().app;
        let rect = self
            .rect()
            .ok_or_else(|| SnapError::Input("no 3D viewport on screen".into()))?;
        let (_, mesh) = app
            .pui
            .mesh
            .as_ref()
            .ok_or_else(|| SnapError::Settle("no viewport mesh".into()))?;
        let plane = app
            .project
            .planes
            .get(app.sel_plane)
            .ok_or_else(|| SnapError::Input("the project has no plane".into()))?;
        // the mesh on screen may still be another plane's
        let (main, wing) = main_wing(plane)?;
        let goal = wing.planform_span() / 4.0;
        let cam = &app.pui.cam;
        let forward = cam.basis().forward;
        mesh.panels
            .iter()
            .filter_map(|panel| {
                let strip = mesh.strips.get(panel.strip?)?;
                if strip.wing != main {
                    return None;
                }
                let (pos, depth) =
                    cam.project_point(panel.centre, rect, 1e-4, app.view3d_prefs.perspective)?;
                let score = (strip.span_pos - goal).abs()
                    + if forward.dot(panel.normal) > 0.0 {
                        goal
                    } else {
                        0.0
                    };
                rect.shrink(12.0)
                    .contains(pos)
                    .then_some((score, depth, pos))
            })
            .min_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.total_cmp(&b.1)))
            .map(|(_, _, pos)| pos)
            .ok_or_else(|| SnapError::Settle("no visible mid-span panel".into()))
    }
}

/// The main wing of `plane` (else its first wing) and its index, for the
/// presets that aim at its mid-span.
fn main_wing(plane: &foiley_plane::Plane) -> Result<(usize, &foiley_plane::geom::Wing), SnapError> {
    let wings = &plane.wings;
    let main = wings
        .iter()
        .position(|wing| wing.kind == foiley_plane::geom::WingKind::Main)
        .unwrap_or(0);
    wings
        .get(main)
        .map(|wing| (main, wing))
        .ok_or_else(|| SnapError::Input("this preset needs a plane with a wing".into()))
}

fn mouse(pos: Pos2, button: PointerButton, pressed: bool, modifiers: Modifiers) -> Event {
    Event::PointerButton {
        pos,
        button,
        pressed,
        modifiers,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_plane_without_wings_has_no_mid_span_to_aim_at() {
        let mut plane = foiley_plane::Plane::default_plane("NACA 2412", "NACA 0009");
        assert_eq!(main_wing(&plane).unwrap().0, 0);
        plane.wings.rotate_left(1);
        assert_eq!(main_wing(&plane).unwrap().0, 2, "the main wing is found wherever it is");
        plane.wings.clear();
        assert!(matches!(main_wing(&plane), Err(SnapError::Input(reason)) if reason.contains("wing")));
    }

    #[test]
    fn runs_started_in_the_same_second_get_their_own_directories() {
        let parent = std::env::temp_dir().join(format!("foiley-snap-out-{}", std::process::id()));
        let dirs: Vec<_> = (0..3).map(|_| new_output_dir(parent.clone()).unwrap()).collect();
        assert!(dirs.iter().all(|dir| dir.is_dir() && dir.parent() == Some(parent.as_path())));
        assert!(dirs[0] != dirs[1] && dirs[1] != dirs[2] && dirs[0] != dirs[2], "{dirs:?}");
        std::fs::remove_dir_all(parent).unwrap();
    }
}
