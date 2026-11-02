//! Foiley Studio - graphical interface.
//!
//! The user interface thread never runs a solver: every analysis (and
//! meshing, saving, loading) is a background job.  The bottom panel shows
//! every job with a live progress bar, stage, point counts, elapsed time,
//! ETA and a Cancel button, so the window stays responsive at all times.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod foil_ui;
mod plane_ui;
mod view3d;
mod widgets;

use eframe::egui::{self, Color32, RichText};
use foiley_app::{JobId, JobKind, JobManager, JobOutput, JobStatus, Project};
use foiley_core::progress::LogLevel;
use foiley_foil::analysis::PointDetail;
use foiley_plane::{Plane, PlanePolarDef};
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::mpsc;
use std::time::{Duration, Instant};

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Module {
    FoilDesign,
    FoilAnalysis,
    PlaneDesign,
    PlaneAnalysis,
    Stability,
}

impl Module {
    const ALL: [(Module, &'static str, &'static str); 5] = [
        (Module::FoilDesign, "Foil design", "Create and modify airfoils: NACA, flaps, thickness, camber, LE radius, TE gap, paneling, interpolation"),
        (Module::FoilAnalysis, "2D analysis", "XFoil viscous/inviscid analysis: polars (types 1-4, 6), single points with Cp and boundary layer, batch polar meshes"),
        (Module::PlaneDesign, "Plane design", "Wings, sections, fuselage, masses; geometry, reference dimensions and inertia"),
        (Module::PlaneAnalysis, "3D analysis", "LLT, VLM and 3D panel analyses (polar types 1-8), span distributions, Cp"),
        (Module::Stability, "Stability", "Trimmed stability: derivatives, modes, root locus, time responses"),
    ];
}

pub enum Pending {
    FoilPolar,
    FoilPoint,
    Batch,
    PlanePolar,
    /// compute the missing 2D polars, then run this 3D analysis
    AutoPolarsThen(Box<(Plane, PlanePolarDef)>),
    Mesh(u64),
    Inverse,
    Optimise,
}

enum IoMsg {
    Loaded(Result<(Project, PathBuf), String>),
    Saved(Result<PathBuf, String>),
}

pub struct LogLine {
    pub level: LogLevel,
    pub text: String,
    pub job: String,
}

pub struct App {
    pub project: Project,
    pub path: Option<PathBuf>,
    pub dirty: bool,
    pub jobs: JobManager,
    pub pending: Vec<(JobId, Pending)>,
    pub module: Module,
    pub sel_foil: Option<String>,
    pub sel_plane: usize,
    pub log: Vec<LogLine>,
    pub log_min: LogLevel,
    undo: Vec<Project>,
    redo: Vec<Project>,
    pub toast: Option<(String, Instant, bool)>,
    pub fui: foil_ui::FoilUi,
    pub pui: plane_ui::PlaneUi,
    pub dark: bool,
    io_tx: mpsc::Sender<IoMsg>,
    io_rx: mpsc::Receiver<IoMsg>,
    io_busy: bool,
    autosave_at: Instant,
    pub details: HashMap<(String, String), Vec<PointDetail>>,
    pub point_detail: Option<PointDetail>,
    bottom_tab: usize,
}

impl App {
    fn new(cc: &eframe::CreationContext) -> Self {
        let mut jobs = JobManager::new(0);
        let ctx = cc.egui_ctx.clone();
        jobs.set_notify(move || ctx.request_repaint());
        let dark = cc.storage.and_then(|s| s.get_string("dark")).map(|v| v == "1").unwrap_or(true);
        cc.egui_ctx.set_visuals(if dark { egui::Visuals::dark() } else { egui::Visuals::light() });
        cc.egui_ctx.style_mut(|s| {
            s.spacing.item_spacing = egui::vec2(8.0, 5.0);
            s.spacing.interact_size.y = 22.0;
        });
        let (io_tx, io_rx) = mpsc::channel();
        let project = Project::starter();
        let sel_foil = project.foils.first().map(|f| f.name.clone());
        // command line: foiley-studio [project.foiley | plane.xml | foil.dat ...] [--module foil|2d|plane|3d|stability]
        let args: Vec<String> = std::env::args().skip(1).collect();
        let mut start_module = Module::FoilDesign;
        let mut run_now = false;
        let mut files = Vec::new();
        let mut it = args.iter();
        while let Some(a) = it.next() {
            if a == "--run" {
                run_now = true;
            } else if a == "--module" {
                start_module = match it.next().map(|s| s.as_str()) {
                    Some("2d") => Module::FoilAnalysis,
                    Some("plane") => Module::PlaneDesign,
                    Some("3d") => Module::PlaneAnalysis,
                    Some("stability") => Module::Stability,
                    _ => Module::FoilDesign,
                };
            } else {
                files.push(PathBuf::from(a));
            }
        }
        let mut app = App {
            project,
            path: None,
            dirty: false,
            jobs,
            pending: vec![],
            module: Module::FoilDesign,
            sel_foil,
            sel_plane: 0,
            log: vec![],
            log_min: LogLevel::Info,
            undo: vec![],
            redo: vec![],
            toast: None,
            fui: Default::default(),
            pui: Default::default(),
            dark,
            io_tx,
            io_rx,
            io_busy: false,
            autosave_at: Instant::now(),
            details: HashMap::new(),
            point_detail: None,
            bottom_tab: 0,
        };
        app.module = start_module;
        let _ = &mut run_now;
        for f in files {
            let ext = f.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default();
            if ext == "foiley" || ext == "json" {
                // load synchronously at start-up so the requested module shows the data
                match Project::load(&f) {
                    Ok(p) => {
                        app.project = p;
                        app.path = Some(f.clone());
                        app.sel_foil = app.project.foils.first().map(|f| f.name.clone());
                        for q in &app.project.foil_polars {
                            app.fui.visible.insert((q.foil.clone(), q.name.clone()));
                        }
                        for q in &app.project.plane_polars {
                            app.pui.visible.insert((q.plane.clone(), q.name.clone()));
                            app.pui.sel_polar = Some((q.plane.clone(), q.name.clone()));
                        }
                    }
                    Err(e) => app.notify(e, true),
                }
            } else {
                match app.project.import(&f) {
                    Ok(m) => app.notify(format!("imported {m}"), false),
                    Err(e) => app.notify(e, true),
                }
            }
        }
        if run_now && !app.project.planes.is_empty() {
            // start the current 3D analysis definition of the first plane (fine panel mesh)
            let mut plane = app.project.planes[0].clone();
            for w in &mut plane.wings {
                for s in &mut w.sections {
                    s.nx = s.nx.max(20);
                    s.ny = s.ny.max(30);
                }
            }
            app.pui.def.method = foiley_plane::analysis::Method::Panel;
            app.pui.def.viscous = false;
            app.module = Module::PlaneAnalysis;
            app.run_3d_public(&plane);
        }
        app
    }

    /// Record the project state before a modification (undo).
    pub fn snapshot(&mut self) {
        self.undo.push(self.project.clone());
        if self.undo.len() > 40 {
            self.undo.remove(0);
        }
        self.redo.clear();
        self.dirty = true;
    }
    fn undo(&mut self) {
        if let Some(p) = self.undo.pop() {
            self.redo.push(std::mem::replace(&mut self.project, p));
            self.dirty = true;
            self.notify("undone", false);
        }
    }
    fn redo(&mut self) {
        if let Some(p) = self.redo.pop() {
            self.undo.push(std::mem::replace(&mut self.project, p));
            self.dirty = true;
        }
    }

    pub fn notify(&mut self, msg: impl Into<String>, error: bool) {
        let m: String = msg.into();
        self.log.push(LogLine { level: if error { LogLevel::Error } else { LogLevel::Info }, text: m.clone(), job: String::new() });
        self.toast = Some((m, Instant::now(), error));
    }

    pub fn submit(&mut self, title: String, kind: JobKind, then: Pending) {
        let j = self.jobs.submit(title, kind);
        self.pending.push((j.id, then));
        self.bottom_tab = 0;
    }

    fn open_dialog(&mut self) {
        if let Some(p) = rfd::FileDialog::new().add_filter("Foiley project", &["foiley", "json"]).pick_file() {
            self.load_async(p);
        }
    }

    fn load_async(&mut self, p: PathBuf) {
        let tx = self.io_tx.clone();
        self.io_busy = true;
        std::thread::spawn(move || {
            let r = Project::load(&p).map(|pr| (pr, p));
            let _ = tx.send(IoMsg::Loaded(r));
        });
    }

    fn save(&mut self, as_new: bool) {
        let path = if as_new || self.path.is_none() {
            match rfd::FileDialog::new().add_filter("Foiley project", &["foiley"]).set_file_name(format!("{}.foiley", self.project.name)).save_file() {
                Some(p) => p,
                None => return,
            }
        } else {
            self.path.clone().unwrap()
        };
        let proj = self.project.clone();
        let tx = self.io_tx.clone();
        self.io_busy = true;
        std::thread::spawn(move || {
            let r = proj.save(&path).map(|_| path);
            let _ = tx.send(IoMsg::Saved(r));
        });
    }

    fn import_dialog(&mut self) {
        if let Some(files) = rfd::FileDialog::new()
            .add_filter("All supported", &["dat", "txt", "xml", "pol", "foiley", "json", "xfl", "fl5", "plr"])
            .add_filter("xflr5 / flow5 projects and polars", &["xfl", "fl5", "plr"])
            .add_filter("Airfoil coordinates", &["dat", "txt"])
            .add_filter("xflr5 / flow5 plane", &["xml"])
            .add_filter("XFoil polar", &["txt", "pol"])
            .pick_files()
        {
            self.snapshot();
            for f in files {
                match self.project.import(&f) {
                    Ok(m) => {
                        self.notify(format!("imported {m}"), false);
                        if f.extension().map(|e| e == "xml").unwrap_or(false) {
                            self.sel_plane = self.project.planes.len() - 1;
                            self.module = Module::PlaneDesign;
                        } else if let Some(last) = self.project.foils.last() {
                            self.sel_foil = Some(last.name.clone());
                        }
                    }
                    Err(e) => self.notify(format!("{}: {e}", f.display()), true),
                }
            }
        }
    }

    /// Collect finished jobs and log lines.
    fn process_jobs(&mut self) {
        for j in self.jobs.jobs() {
            for l in j.progress.drain_log() {
                self.log.push(LogLine { level: l.level, text: l.text, job: j.title.clone() });
            }
        }
        if self.log.len() > 5000 {
            self.log.drain(0..1000);
        }
        let mut done = Vec::new();
        for (i, (id, _)) in self.pending.iter().enumerate() {
            if let Some(j) = self.jobs.get(*id) {
                if j.status().is_finished() {
                    done.push(i);
                }
            }
        }
        for i in done.into_iter().rev() {
            let (id, what) = self.pending.remove(i);
            let Some(job) = self.jobs.get(id) else { continue };
            let st = job.status();
            let out = job.take_output();
            match (st, out) {
                (JobStatus::Done, Some(out)) | (JobStatus::Cancelled, Some(out)) => self.apply_output(&job.title, what, out),
                (JobStatus::Failed(e), _) => self.notify(format!("{} failed: {e}", job.title), true),
                (JobStatus::Cancelled, None) => self.notify(format!("{} cancelled", job.title), false),
                _ => {}
            }
        }
    }

    fn apply_output(&mut self, title: &str, what: Pending, out: JobOutput) {
        match (what, out) {
            (Pending::FoilPolar, JobOutput::FoilPolar { polar, details }) => {
                let key = (polar.foil.clone(), polar.name.clone());
                let msg = format!("{}: {} points{}", polar.name, polar.points.len(), if polar.failed.is_empty() { String::new() } else { format!(", {} not converged", polar.failed.len()) });
                if !details.is_empty() {
                    self.details.insert(key.clone(), details);
                }
                self.snapshot();
                self.fui.visible.insert(key);
                self.project.add_foil_polar(polar);
                self.notify(msg, false);
            }
            (Pending::FoilPoint, JobOutput::FoilPoint(d)) => match d {
                Some(d) => {
                    self.point_detail = Some(d);
                    self.notify(format!("{title}: converged"), false);
                }
                None => self.notify(format!("{title}: did not converge - try a different alpha, Re, or Ncrit"), true),
            },
            (Pending::Batch, JobOutput::Batch(pols)) => {
                self.snapshot();
                let n = pols.len();
                for p in pols {
                    self.fui.visible.insert((p.foil.clone(), p.name.clone()));
                    self.project.add_foil_polar(p);
                }
                self.notify(format!("{title}: {n} polars added"), false);
            }
            (Pending::AutoPolarsThen(b), JobOutput::Batch(pols)) => {
                self.snapshot();
                for p in pols {
                    self.project.add_foil_polar(p);
                }
                let (plane, def) = *b;
                let db = self.project.foil_db();
                let t = format!("{} - {}", plane.name, def.name);
                self.submit(t, JobKind::PlanePolar { plane, db, def }, Pending::PlanePolar);
            }
            (Pending::PlanePolar, JobOutput::PlanePolar(pol)) => {
                self.snapshot();
                let msg = format!("{}: {} points ({} unknowns){}", pol.name, pol.points.len(), pol.unknowns, if pol.failed.is_empty() { String::new() } else { format!(", {} failed", pol.failed.len()) });
                let key = (pol.plane.clone(), pol.name.clone());
                self.pui.visible.insert(key.clone());
                self.pui.sel_polar = Some(key);
                self.pui.sel_point = pol.points.len().saturating_sub(1).min(pol.points.len() / 2);
                self.project.add_plane_polar(pol);
                self.notify(msg, false);
            }
            (Pending::Inverse, JobOutput::Inverse(foil, d)) => {
                match foil {
                    Some(f) => {
                        self.fui.preview = Some(f);
                        self.notify(format!("{title}: new foil ready (orange preview) - CL {:.3}", d.cl), false);
                    }
                    None => self.fui.qdist = Some(d),
                }
            }
            (Pending::Optimise, JobOutput::Optimise(r)) => {
                let msg = match &r.best {
                    Some(b) => format!("{title}: best CL {:.4} CD {:.5} after {} evaluations", b.cl, b.cd, r.evaluations),
                    None => format!("{title}: no feasible design"),
                };
                self.fui.opt_result = Some(r);
                self.notify(msg, false);
            }
            (Pending::Mesh(key), JobOutput::Mesh(m)) => {
                self.pui.mesh = Some((key, m));
                self.pui.mesh_job = None;
            }
            _ => {}
        }
    }

    fn explorer(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.heading("Project");
            if self.dirty {
                ui.label(RichText::new("modified").small().italics());
            }
        });
        ui.add(egui::TextEdit::singleline(&mut self.project.name).hint_text("project name"));
        ui.separator();
        egui::ScrollArea::vertical().id_salt("explorer").show(ui, |ui| {
            egui::CollapsingHeader::new(RichText::new(format!("Foils ({})", self.project.foils.len())).strong()).default_open(true).show(ui, |ui| {
                let names: Vec<String> = self.project.foils.iter().map(|f| f.name.clone()).collect();
                for n in names {
                    let np = self.project.foil_polars.iter().filter(|p| p.foil == n).count();
                    let sel = self.sel_foil.as_deref() == Some(n.as_str());
                    let label = if np > 0 { format!("{n}  [{np}]") } else { n.clone() };
                    let r = ui.selectable_label(sel, label).on_hover_text(format!("{np} polar(s)"));
                    if r.clicked() {
                        self.sel_foil = Some(n.clone());
                        if !matches!(self.module, Module::FoilDesign | Module::FoilAnalysis) {
                            self.module = Module::FoilDesign;
                        }
                    }
                    r.context_menu(|ui| {
                        if ui.button("Duplicate").clicked() {
                            if let Some(f) = self.project.foil(&n).cloned() {
                                self.snapshot();
                                let mut g = f;
                                g.name = self.project.unique_foil_name(&g.name);
                                self.sel_foil = Some(g.name.clone());
                                self.project.foils.push(g);
                            }
                            ui.close();
                        }
                        if ui.button("Delete (with its polars)").clicked() {
                            self.snapshot();
                            self.project.remove_foil(&n);
                            ui.close();
                        }
                    });
                }
            });
            egui::CollapsingHeader::new(RichText::new(format!("Planes ({})", self.project.planes.len())).strong()).default_open(true).show(ui, |ui| {
                let mut del = None;
                for i in 0..self.project.planes.len() {
                    let name = self.project.planes[i].name.clone();
                    let np = self.project.plane_polars.iter().filter(|p| p.plane == name).count();
                    let r = ui.selectable_label(self.sel_plane == i && matches!(self.module, Module::PlaneDesign | Module::PlaneAnalysis | Module::Stability), if np > 0 { format!("{name}  [{np}]") } else { name.clone() });
                    if r.clicked() {
                        self.sel_plane = i;
                        if !matches!(self.module, Module::PlaneDesign | Module::PlaneAnalysis | Module::Stability) {
                            self.module = Module::PlaneDesign;
                        }
                    }
                    r.context_menu(|ui| {
                        if ui.button("Duplicate").clicked() {
                            let mut p = self.project.planes[i].clone();
                            p.name = format!("{} copy", p.name);
                            self.snapshot();
                            self.project.planes.push(p);
                            ui.close();
                        }
                        if ui.button("Delete (with its polars)").clicked() {
                            del = Some(i);
                            ui.close();
                        }
                    });
                }
                if let Some(i) = del {
                    self.snapshot();
                    let name = self.project.planes[i].name.clone();
                    self.project.planes.remove(i);
                    self.project.plane_polars.retain(|p| p.plane != name);
                    self.sel_plane = self.sel_plane.min(self.project.planes.len().saturating_sub(1));
                }
                if ui.button("+ New plane").on_hover_text("Default glider: wing, elevator, fin").clicked() {
                    self.snapshot();
                    let wf = self.project.foils.first().map(|f| f.name.clone()).unwrap_or_else(|| "NACA 2412".into());
                    let tf = self.project.foils.get(1).map(|f| f.name.clone()).unwrap_or_else(|| wf.clone());
                    let mut p = Plane::default_plane(&wf, &tf);
                    p.name = format!("Plane {}", self.project.planes.len() + 1);
                    self.project.planes.push(p);
                    self.sel_plane = self.project.planes.len() - 1;
                    self.module = Module::PlaneDesign;
                }
            });
            ui.separator();
            ui.label(RichText::new(format!("{} 2D polars, {} 3D polars", self.project.foil_polars.len(), self.project.plane_polars.len())).small());
        });
    }

    fn jobs_panel(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.bottom_tab, 0, format!("Jobs ({} running)", self.jobs.running()));
            ui.selectable_value(&mut self.bottom_tab, 1, format!("Log ({})", self.log.len()));
            ui.separator();
            if self.bottom_tab == 0 {
                if ui.button("Cancel all").clicked() {
                    self.jobs.cancel_all();
                }
                if ui.button("Clear finished").clicked() {
                    self.jobs.clear_finished();
                }
                ui.label(RichText::new(format!("{} compute threads (one core kept free for the interface)", self.jobs.threads)).small());
            } else {
                ui.label("show:");
                ui.selectable_value(&mut self.log_min, LogLevel::Debug, "all");
                ui.selectable_value(&mut self.log_min, LogLevel::Info, "info");
                ui.selectable_value(&mut self.log_min, LogLevel::Warn, "warnings");
                if ui.button("Clear").clicked() {
                    self.log.clear();
                }
            }
        });
        ui.separator();
        if self.bottom_tab == 0 {
            egui::ScrollArea::vertical().id_salt("jobs").auto_shrink([false, false]).show(ui, |ui| {
                // finished background meshing is not interesting: hide it
                let jobs: Vec<_> = self.jobs.jobs().iter().rev().filter(|j| !(j.kind == "mesh" && j.status() == JobStatus::Done)).cloned().collect();
                if jobs.is_empty() {
                    widgets::hint(ui, "No jobs yet. Analyses run here in the background; the window stays responsive and every job can be cancelled.");
                }
                for j in jobs {
                    ui.horizontal(|ui| {
                        let st = j.status();
                        let (ok, bad) = j.progress.counts();
                        let (txt, col) = match &st {
                            JobStatus::Running => ("running", Color32::from_rgb(90, 160, 255)),
                            JobStatus::Done => ("done", Color32::from_rgb(80, 190, 90)),
                            JobStatus::Cancelled => ("cancelled", Color32::GRAY),
                            JobStatus::Failed(_) => ("failed", Color32::from_rgb(230, 80, 70)),
                        };
                        ui.label(RichText::new(txt).color(col).strong());
                        ui.label(RichText::new(&j.title).strong());
                        ui.label(RichText::new(j.kind).small());
                        let f = if st == JobStatus::Done { 1.0 } else { j.progress.fraction() as f32 };
                        let mut stage = j.progress.stage_text();
                        if let JobStatus::Failed(e) = &st {
                            stage = e.clone();
                        }
                        ui.add(egui::ProgressBar::new(f).desired_width(260.0).text(format!("{:.0}%", f * 100.0)));
                        ui.label(widgets::fmt_secs(j.elapsed()));
                        if st == JobStatus::Running {
                            if let Some(eta) = j.progress.eta() {
                                ui.label(format!("ETA {}", widgets::fmt_secs(eta)));
                            }
                            if ui.button("Cancel").clicked() {
                                j.cancel();
                            }
                        }
                        if ok + bad > 0 {
                            ui.label(format!("{ok} ok"));
                            if bad > 0 {
                                ui.label(RichText::new(format!("{bad} failed")).color(Color32::from_rgb(230, 150, 60))).on_hover_text("points that did not converge even after the automatic recovery steps (they are listed in the log and excluded from the results)");
                            }
                        }
                        ui.label(RichText::new(stage).small());
                    });
                }
            });
        } else {
            egui::ScrollArea::vertical().id_salt("log").stick_to_bottom(true).auto_shrink([false, false]).show(ui, |ui| {
                for l in self.log.iter().filter(|l| l.level >= self.log_min) {
                    let col = match l.level {
                        LogLevel::Error => Color32::from_rgb(230, 80, 70),
                        LogLevel::Warn => Color32::from_rgb(230, 160, 60),
                        LogLevel::Debug => ui.visuals().weak_text_color(),
                        LogLevel::Info => ui.visuals().text_color(),
                    };
                    let prefix = if l.job.is_empty() { String::new() } else { format!("[{}] ", l.job) };
                    ui.label(RichText::new(format!("{prefix}{}", l.text)).color(col).monospace());
                }
            });
        }
    }
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // background IO
        while let Ok(m) = self.io_rx.try_recv() {
            self.io_busy = false;
            match m {
                IoMsg::Loaded(Ok((p, path))) => {
                    self.snapshot();
                    self.project = p;
                    self.path = Some(path.clone());
                    self.dirty = false;
                    self.sel_foil = self.project.foils.first().map(|f| f.name.clone());
                    self.sel_plane = 0;
                    self.details.clear();
                    self.pui = Default::default();
                    self.notify(format!("opened {}", path.display()), false);
                }
                IoMsg::Loaded(Err(e)) => self.notify(e, true),
                IoMsg::Saved(Ok(p)) => {
                    if p.extension().map(|e| e == "foiley").unwrap_or(false) && !p.to_string_lossy().contains("autosave") {
                        self.path = Some(p.clone());
                        self.dirty = false;
                        self.notify(format!("saved {}", p.display()), false);
                    }
                }
                IoMsg::Saved(Err(e)) => self.notify(format!("save failed: {e}"), true),
            }
        }
        self.process_jobs();
        // crash safety: autosave every 3 minutes when modified
        if self.dirty && self.autosave_at.elapsed() > Duration::from_secs(180) && !self.io_busy {
            self.autosave_at = Instant::now();
            let p = self.project.clone();
            let tx = self.io_tx.clone();
            let path = std::env::temp_dir().join("foiley_autosave.foiley");
            std::thread::spawn(move || {
                let _ = tx.send(IoMsg::Saved(p.save(&path).map(|_| path)));
            });
        }
        // keyboard shortcuts
        let (s, o, z, y) = ctx.input(|i| (i.modifiers.command && i.key_pressed(egui::Key::S), i.modifiers.command && i.key_pressed(egui::Key::O), i.modifiers.command && i.key_pressed(egui::Key::Z), i.modifiers.command && i.key_pressed(egui::Key::Y)));
        if s {
            self.save(false);
        }
        if o {
            self.open_dialog();
        }
        if z && !ctx.wants_keyboard_input() {
            self.undo();
        }
        if y && !ctx.wants_keyboard_input() {
            self.redo();
        }
        let title = format!("Foiley Studio - {}{}", self.path.as_ref().map(|p| p.file_name().unwrap_or_default().to_string_lossy().to_string()).unwrap_or_else(|| self.project.name.clone()), if self.dirty { " *" } else { "" });
        ctx.send_viewport_cmd(egui::ViewportCommand::Title(title));

        egui::TopBottomPanel::top("menu").show(ctx, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("New project").clicked() {
                        self.snapshot();
                        self.project = Project::starter();
                        self.path = None;
                        self.sel_foil = self.project.foils.first().map(|f| f.name.clone());
                        self.pui = Default::default();
                        ui.close();
                    }
                    if ui.button("Open...        Ctrl+O").clicked() {
                        ui.close();
                        self.open_dialog();
                    }
                    if ui.button("Save           Ctrl+S").clicked() {
                        ui.close();
                        self.save(false);
                    }
                    if ui.button("Save as...").clicked() {
                        ui.close();
                        self.save(true);
                    }
                    ui.separator();
                    if ui.button("Import foils / planes / polars...").on_hover_text(".dat foils (Selig or Lednicer), xflr5/flow5 plane .xml, XFoil polar .txt, other projects").clicked() {
                        ui.close();
                        self.import_dialog();
                    }
                    let recovery = std::env::temp_dir().join("foiley_autosave.foiley");
                    if recovery.exists() && ui.button("Recover autosave").clicked() {
                        ui.close();
                        self.load_async(recovery);
                    }
                });
                ui.menu_button("Edit", |ui| {
                    if ui.add_enabled(!self.undo.is_empty(), egui::Button::new("Undo   Ctrl+Z")).clicked() {
                        self.undo();
                        ui.close();
                    }
                    if ui.add_enabled(!self.redo.is_empty(), egui::Button::new("Redo   Ctrl+Y")).clicked() {
                        self.redo();
                        ui.close();
                    }
                });
                ui.menu_button("View", |ui| {
                    if ui.checkbox(&mut self.dark, "Dark theme").changed() {
                        ctx.set_visuals(if self.dark { egui::Visuals::dark() } else { egui::Visuals::light() });
                    }
                    ui.label("Zoom: Ctrl +/- ");
                });
                ui.menu_button("Help", |ui| {
                    ui.label("Foiley Studio");
                    widgets::hint(ui, "2D: XFoil 6.97 engine (Rust port).\n3D: LLT, VLM1/VLM2, constant source/doublet panels.\nAll analyses run in the background: watch the Jobs panel.\nHover any field for an explanation.");
                });
                ui.separator();
                for (m, name, help) in Module::ALL {
                    if ui.selectable_label(self.module == m, RichText::new(name).strong()).on_hover_text(help).clicked() {
                        self.module = m;
                    }
                }
            });
        });
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                let r = self.jobs.running();
                if r > 0 {
                    ui.spinner();
                    ui.label(format!("{r} job(s) running"));
                } else {
                    ui.label("ready");
                }
                if self.io_busy {
                    ui.spinner();
                    ui.label("reading / writing file...");
                }
                if let Some((m, t, err)) = &self.toast {
                    if t.elapsed() < Duration::from_secs(8) {
                        ui.separator();
                        ui.label(RichText::new(m).color(if *err { Color32::from_rgb(230, 80, 70) } else { ui.visuals().text_color() }));
                    }
                }
            });
        });
        egui::TopBottomPanel::bottom("jobs").resizable(true).default_height(150.0).min_height(60.0).show(ctx, |ui| self.jobs_panel(ui));
        egui::SidePanel::left("explorer").resizable(true).default_width(200.0).show(ctx, |ui| self.explorer(ui));
        egui::CentralPanel::default().show(ctx, |ui| match self.module {
            Module::FoilDesign => self.foil_design(ui),
            Module::FoilAnalysis => self.foil_analysis(ui),
            Module::PlaneDesign => self.plane_design(ui),
            Module::PlaneAnalysis => self.plane_analysis(ui),
            Module::Stability => self.stability(ui),
        });
        if self.jobs.running() > 0 || self.io_busy {
            ctx.request_repaint_after(Duration::from_millis(100));
        }
    }

    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        storage.set_string("dark", if self.dark { "1".into() } else { "0".into() });
    }

    fn on_exit(&mut self, _gl: Option<&eframe::glow::Context>) {
        self.jobs.cancel_all();
    }
}

fn main() -> eframe::Result {
    let opts = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([1400.0, 900.0]).with_min_inner_size([900.0, 600.0]).with_title("Foiley Studio"),
        ..Default::default()
    };
    eframe::run_native("Foiley Studio", opts, Box::new(|cc| Ok(Box::new(App::new(cc)))))
}
