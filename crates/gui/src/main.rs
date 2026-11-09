//! Foiley Studio - graphical interface.
//!
//! The user interface thread never runs a solver: every analysis (and
//! meshing, saving, loading) is a background job.  The bottom panel shows
//! every job with a live progress bar, stage, point counts, elapsed time,
//! ETA and a Cancel button, so the window stays responsive at all times.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod body_list;
mod foil_ui;
mod plane_ui;
mod view3d;
mod widgets;

use eframe::egui::{self, Color32, RichText};
use foiley_app::{JobId, JobKind, JobManager, JobOutput, JobStatus, Project};
use foiley_core::progress::LogLevel;
use foiley_foil::analysis::PointDetail;
use foiley_plane::{Plane, PlanePolarDef};
use std::collections::{HashMap, VecDeque};
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
    /// apply an imported CAD body to the plane of this name
    CadImport { plane: String },
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct ProjectRevision {
    project: u64,
    revision: u64,
}

pub struct PendingJob {
    id: JobId,
    generation: u64,
    action: Pending,
}

enum IoMsg {
    Loaded { project: u64, result: Result<(Project, PathBuf), String> },
    Saved { version: ProjectRevision, autosave: bool, result: Result<PathBuf, String> },
}

enum IoTask {
    Load { owner: u64, path: PathBuf },
    Save { version: ProjectRevision, project: Project, path: PathBuf, autosave: bool },
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
    pub pending: Vec<PendingJob>,
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
    version: ProjectRevision,
    job_generation: u64,
    queued_io: VecDeque<IoTask>,
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
        let project = Project::starter();
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
        let mut app = Self::from_project(project, jobs, dark);
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
                            app.fui.visible.insert(q.key());
                        }
                        for q in &app.project.plane_polars {
                            app.pui.visible.insert((q.plane.clone(), q.name.clone()));
                            app.pui.sel_polar = Some((q.plane.clone(), q.name.clone()));
                        }
                    }
                    Err(e) => app.notify(e, true),
                }
            } else {
                match import_file(&mut app.project, &f) {
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

    fn from_project(project: Project, jobs: JobManager, dark: bool) -> Self {
        let (io_tx, io_rx) = mpsc::channel();
        let sel_foil = project.foils.first().map(|f| f.name.clone());
        App {
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
            version: ProjectRevision::default(),
            job_generation: 0,
            queued_io: VecDeque::new(),
            autosave_at: Instant::now(),
            details: HashMap::new(),
            point_detail: None,
            bottom_tab: 0,
        }
    }

    /// Record the project state before a modification (undo).
    pub fn snapshot(&mut self) {
        self.undo.push(self.project.clone());
        if self.undo.len() > 40 {
            self.undo.remove(0);
        }
        self.redo.clear();
        self.version.revision += 1;
        self.dirty = true;
    }
    fn undo(&mut self) {
        if let Some(p) = self.undo.pop() {
            self.redo.push(std::mem::replace(&mut self.project, p));
            self.version.revision += 1;
            self.dirty = true;
            self.notify("undone", false);
        }
    }
    fn redo(&mut self) {
        if let Some(p) = self.redo.pop() {
            self.undo.push(std::mem::replace(&mut self.project, p));
            self.version.revision += 1;
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
        self.track_job(j.id, then);
        self.bottom_tab = 0;
    }

    fn track_job(&mut self, id: JobId, action: Pending) {
        self.pending.push(PendingJob { id, generation: self.job_generation, action });
    }

    fn current_job(&self, id: JobId) -> bool {
        self.pending.iter().any(|p| p.id == id && p.generation == self.job_generation)
    }

    fn invalidate_project_jobs(&mut self) {
        self.job_generation += 1;
        for pending in &self.pending {
            if let Some(job) = self.jobs.get(pending.id) {
                job.cancel();
            }
        }
        // Keep their ownership records until completion so cancelled partial
        // results are collected and discarded, never applied to the successor.
        self.pui.mesh_job = None;
        self.fui.opt_job = None;
        self.details.clear();
        self.point_detail = None;
    }

    fn replace_project(&mut self, project: Project, path: Option<PathBuf>, dirty: bool) {
        self.invalidate_project_jobs();
        self.version.project += 1;
        self.version.revision = 0;
        self.project = project;
        self.path = path;
        self.dirty = dirty;
        self.undo.clear();
        self.redo.clear();
        self.sel_foil = self.project.foils.first().map(|f| f.name.clone());
        self.sel_plane = 0;
        self.fui = Default::default();
        for p in &self.project.foil_polars {
            self.fui.visible.insert(p.key());
        }
        self.pui = Default::default();
        for p in &self.project.plane_polars {
            let key = (p.plane.clone(), p.name.clone());
            self.pui.visible.insert(key.clone());
            self.pui.sel_polar = Some(key);
        }
    }

    fn open_dialog(&mut self) {
        if let Some(p) = rfd::FileDialog::new().add_filter("Foiley project", &["foiley", "json"]).pick_file() {
            self.load_async(p);
        }
    }

    fn load_async(&mut self, p: PathBuf) {
        self.queued_io.push_back(IoTask::Load { owner: self.version.project, path: p });
        self.start_next_io();
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
        self.save_to(path, false);
    }

    fn save_to(&mut self, path: PathBuf, autosave: bool) {
        self.queued_io.push_back(IoTask::Save { version: self.version, project: self.project.clone(), path, autosave });
        self.start_next_io();
    }

    fn start_next_io(&mut self) {
        // Requests retain their captured project/revision while writes run
        // serially. Ctrl+S and Open remain effective during an autosave.
        if self.io_busy { return; }
        while let Some(task) = self.queued_io.pop_front() {
            if matches!(&task, IoTask::Load { owner, .. } if *owner != self.version.project) { continue; }
            self.io_busy = true;
            let tx = self.io_tx.clone();
            std::thread::spawn(move || {
                let message = match task {
                    IoTask::Load { owner, path } => {
                        let result = Project::load(&path).map(|project| (project, path));
                        IoMsg::Loaded { project: owner, result }
                    }
                    IoTask::Save { version, project, path, autosave } => {
                        let result = project.save(&path).map(|_| path);
                        IoMsg::Saved { version, autosave, result }
                    }
                };
                let _ = tx.send(message);
            });
            break;
        }
    }

    fn process_io(&mut self, message: IoMsg) {
        self.io_busy = false;
        match message {
            IoMsg::Loaded { project, result } if project == self.version.project => match result {
                Ok((project, path)) => {
                    self.replace_project(project, Some(path.clone()), false);
                    self.notify(format!("opened {}", path.display()), false);
                }
                Err(e) => self.notify(e, true),
            },
            IoMsg::Saved { version, autosave, result } if version.project == self.version.project => match result {
                Ok(path) if !autosave => {
                    self.path = Some(path.clone());
                    if version == self.version { self.dirty = false; }
                    self.notify(format!("saved {}", path.display()), false);
                }
                Err(e) => self.notify(format!("save failed: {e}"), true),
                _ => {}
            },
            _ => {}
        }
        self.start_next_io();
    }

    fn import_dialog(&mut self) {
        if let Some(files) = rfd::FileDialog::new()
            .add_filter("All supported", &["dat", "txt", "xml", "pol", "foiley", "json", "xfl", "fl5", "plr", "avl"])
            .add_filter("xflr5 / flow5 projects and polars", &["xfl", "fl5", "plr"])
            .add_filter("AVL geometry", &["avl"])
            .add_filter("Airfoil coordinates", &["dat", "txt"])
            .add_filter("xflr5 / flow5 plane", &["xml"])
            .add_filter("XFoil polar", &["txt", "pol"])
            .pick_files()
        {
            self.snapshot();
            for f in files {
                match import_file(&mut self.project, &f) {
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
        for (i, pending) in self.pending.iter().enumerate() {
            if let Some(j) = self.jobs.get(pending.id) {
                if j.status().is_finished() {
                    done.push(i);
                }
            }
        }
        for i in done.into_iter().rev() {
            let PendingJob { id, generation, action: what } = self.pending.remove(i);
            let Some(job) = self.jobs.get(id) else { continue };
            let st = job.status();
            let out = job.take_output();
            if generation != self.job_generation { continue; }
            match (st, out) {
                (JobStatus::Done, Some(out)) | (JobStatus::Cancelled, Some(out)) => self.apply_output(generation, &job.title, what, out),
                (JobStatus::Failed(e), _) => self.notify(format!("{} failed: {e}", job.title), true),
                (JobStatus::Cancelled, None) => self.notify(format!("{} cancelled", job.title), false),
                _ => {}
            }
        }
    }

    fn apply_output(&mut self, generation: u64, title: &str, what: Pending, out: JobOutput) {
        if generation != self.job_generation { return; }
        match (what, out) {
            (Pending::FoilPolar, JobOutput::FoilPolar { polar, details }) => {
                let msg = format!("{}: {} points{}", polar.name, polar.points.len(), if polar.failed.is_empty() { String::new() } else { format!(", {} not converged", polar.failed.len()) });
                self.snapshot();
                let key = self.project.add_foil_polar(polar);
                if !details.is_empty() {
                    let stored = self.details.entry(key.clone()).or_default();
                    for detail in details {
                        stored.retain(|d| d.point != detail.point);
                        stored.push(detail);
                    }
                }
                self.fui.visible.insert(key);
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
                    let key = self.project.add_foil_polar(p);
                    self.fui.visible.insert(key);
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
            (Pending::CadImport { plane }, JobOutput::CadImport(r)) => match self.project.planes.iter().position(|p| p.name == plane) {
                Some(i) => {
                    let mut edited = self.project.planes[i].clone();
                    match r.apply(&mut edited) {
                        Ok(msg) => {
                            self.snapshot();
                            self.project.planes[i] = edited;
                            self.notify(format!("{title}: {msg}"), false);
                        }
                        Err(e) => self.notify(format!("{title}: {e}"), true),
                    }
                }
                None => self.notify(format!("{title}: plane '{plane}' no longer exists"), true),
            },
            (Pending::Mesh(key), JobOutput::Mesh(m)) if self.pui.mesh_job.map(|(_, pending_key)| pending_key) == Some(key) => {
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
        let mut name = self.project.name.clone();
        if ui.add(egui::TextEdit::singleline(&mut name).hint_text("project name")).changed() {
            // A rename is not worth an undo snapshot per keystroke.
            self.version.revision += 1;
            self.dirty = true;
            self.project.name = name;
        }
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
        // Background completions carry the exact project/revision captured
        // when the operation started.
        while let Ok(message) = self.io_rx.try_recv() {
            self.process_io(message);
        }
        self.process_jobs();
        // crash safety: autosave every 3 minutes when modified
        if self.dirty && self.autosave_at.elapsed() > Duration::from_secs(180) && !self.io_busy {
            self.autosave_at = Instant::now();
            self.save_to(std::env::temp_dir().join("foiley_autosave.foiley"), true);
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
                        self.replace_project(Project::starter(), None, false);
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

/// Import a file into the project.  xflr5 / flow5 projects, .plr polars and
/// AVL geometry return their full conversion report as the message.
fn import_file(p: &mut Project, f: &std::path::Path) -> Result<String, String> {
    if foiley_app::legacy_io::is_interchange(f) {
        foiley_app::legacy_io::import_into(p, f).map(|r| r.text().trim_end().to_string())
    } else {
        p.import(f)
    }
}

#[cfg(test)]
mod review_regressions {
    use super::*;
    use foiley_foil::{AnalysisSettings, FoilPoint, FoilPolar};

    fn app() -> App {
        App::from_project(Project::starter(), JobManager::new(1), true)
    }

    #[test]
    fn save_completion_tracks_exact_project_revision_and_autosave_kind() {
        let mut app = app();
        app.snapshot();
        let saved = app.version;
        app.snapshot();
        app.project.name = "edit made during save".into();
        // A legitimate manual filename containing "autosave" still saves.
        let path = PathBuf::from("my_autosave_research.json");
        app.process_io(IoMsg::Saved { version: saved, autosave: false, result: Ok(path.clone()) });
        assert_eq!(app.path, Some(path.clone()));
        assert!(app.dirty, "the saved older revision must not clear new edits");
        app.process_io(IoMsg::Saved { version: app.version, autosave: true, result: Ok(PathBuf::from("recovery.foiley")) });
        assert!(app.dirty, "autosave is recovery, not a user save");
        assert_eq!(app.path, Some(path.clone()));
        app.process_io(IoMsg::Saved { version: app.version, autosave: false, result: Ok(path) });
        assert!(!app.dirty);
        app.snapshot();
        app.process_io(IoMsg::Saved { version: app.version, autosave: false, result: Err("disk full".into()) });
        assert!(app.dirty);
    }

    #[test]
    fn replacement_ignores_old_save_and_old_load_completion() {
        let mut app = app();
        app.snapshot();
        let old = app.version;
        let successor_path = PathBuf::from("successor.foiley");
        app.replace_project(Project::new("successor"), Some(successor_path.clone()), true);
        app.process_io(IoMsg::Saved { version: old, autosave: false, result: Ok(PathBuf::from("old.foiley")) });
        assert!(app.dirty);
        assert_eq!(app.path, Some(successor_path.clone()));
        app.process_io(IoMsg::Loaded { project: old.project, result: Ok((Project::new("old load"), PathBuf::from("old.foiley"))) });
        assert_eq!(app.project.name, "successor");
        assert_eq!(app.path, Some(successor_path));
    }

    #[test]
    fn overlapping_io_is_serialized_and_new_edits_stay_dirty_on_disk() {
        let mut app = app();
        let path = std::env::temp_dir().join(format!("foiley_save_revision_{}.foiley", std::process::id()));
        let second_path = path.with_extension("queued.foiley");
        app.snapshot();
        app.project.name = "written revision".into();
        app.save_to(path.clone(), false);
        app.snapshot();
        app.project.name = "newer revision".into();
        app.save_to(second_path.clone(), false);
        assert!(app.io_busy);
        assert_eq!(app.queued_io.len(), 1);
        let completion = app.io_rx.recv_timeout(Duration::from_secs(10)).unwrap();
        app.process_io(completion);
        assert!(app.dirty);
        assert_eq!(Project::load(&path).unwrap().name, "written revision");
        assert!(app.io_busy, "the queued explicit save must start after the first completion");
        let completion = app.io_rx.recv_timeout(Duration::from_secs(10)).unwrap();
        app.process_io(completion);
        assert!(!app.dirty);
        assert_eq!(Project::load(&second_path).unwrap().name, "newer revision");
        std::fs::remove_file(path).unwrap();
        std::fs::remove_file(second_path).unwrap();
    }

    #[test]
    fn undo_during_save_as_preserves_document_identity_and_destination() {
        let mut app = app();
        app.snapshot();
        app.project.name = "written".into();
        let saved = app.version;
        app.undo();
        assert_eq!(saved.project, app.version.project);
        assert_ne!(saved.revision, app.version.revision);
        let path = PathBuf::from("save_as.foiley");
        app.process_io(IoMsg::Saved { version: saved, autosave: false, result: Ok(path.clone()) });
        assert_eq!(app.path, Some(path));
        assert!(app.dirty);
    }

    #[test]
    fn queued_manual_save_writes_its_snapshot_after_new_project() {
        let mut app = app();
        let dir = std::env::temp_dir().join(format!("foiley_queued_save_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let autosave = dir.join("recovery.foiley");
        let manual = dir.join("explicit.foiley");
        app.snapshot();
        app.project.name = "original snapshot".into();
        app.save_to(autosave.clone(), true);
        app.save_to(manual.clone(), false);
        app.replace_project(Project::new("new document"), None, true);
        for _ in 0..2 {
            let completion = app.io_rx.recv_timeout(Duration::from_secs(10)).unwrap();
            app.process_io(completion);
        }
        assert_eq!(Project::load(&manual).unwrap().name, "original snapshot");
        assert_eq!(app.project.name, "new document");
        assert!(app.path.is_none());
        assert!(app.dirty);
        assert!(!app.io_busy);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn open_requested_during_autosave_runs_after_the_write() {
        let mut app = app();
        let dir = std::env::temp_dir().join(format!("foiley_queued_open_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let incoming = dir.join("incoming.foiley");
        Project::new("opened after autosave").save(&incoming).unwrap();
        app.snapshot();
        app.save_to(dir.join("recovery.foiley"), true);
        app.load_async(incoming.clone());
        assert_eq!(app.queued_io.len(), 1);
        for _ in 0..2 {
            let completion = app.io_rx.recv_timeout(Duration::from_secs(10)).unwrap();
            app.process_io(completion);
        }
        assert_eq!(app.project.name, "opened after autosave");
        assert_eq!(app.path, Some(incoming));
        assert!(!app.dirty);
        assert!(!app.io_busy);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn replaced_project_rejects_old_results_and_chained_jobs() {
        let mut app = app();
        let old = app.job_generation;
        let plane = app.project.planes[0].clone();
        app.track_job(123, Pending::AutoPolarsThen(Box::new((plane.clone(), PlanePolarDef::default()))));
        assert_eq!(app.pending[0].generation, old);
        assert!(app.current_job(123));
        app.replace_project(Project::new("successor"), None, false);
        assert!(!app.current_job(123), "partial previews must belong to the current project too");
        let polar = FoilPolar { foil: "old foil".into(), ..Default::default() };
        app.apply_output(old, "old 2D", Pending::FoilPolar, JobOutput::FoilPolar { polar: polar.clone(), details: vec![PointDetail::default()] });
        app.apply_output(old, "old batch", Pending::Batch, JobOutput::Batch(vec![polar.clone()]));
        app.apply_output(old, "old chain", Pending::AutoPolarsThen(Box::new((plane, PlanePolarDef::default()))), JobOutput::Batch(vec![polar]));
        app.apply_output(old, "old point", Pending::FoilPoint, JobOutput::FoilPoint(Some(PointDetail::default())));
        app.apply_output(old, "old CAD", Pending::CadImport { plane: "successor".into() }, JobOutput::CadImport(foiley_app::CadImportResult { body: 0, name: "stale".into(), triangles: vec![], notes: vec![] }));
        assert!(app.project.foil_polars.is_empty());
        assert!(app.jobs.jobs().is_empty(), "an old automatic batch must not launch a 3D successor job");
        assert!(app.details.is_empty());
        assert!(app.point_detail.is_none());
        assert!(!app.dirty);
    }

    #[test]
    fn completed_pending_job_is_collected_without_applying_to_loaded_project() {
        let mut app = app();
        let plane = app.project.planes[0].clone();
        app.submit("automatic batch".into(), JobKind::Batch { foils: vec![], cases: vec![], settings: AnalysisSettings::default(), sweep: foiley_foil::Sweep::AlphaList(vec![]) }, Pending::AutoPolarsThen(Box::new((plane, PlanePolarDef::default()))));
        let id = app.pending[0].id;
        app.jobs.wait(id, |_| {});
        assert!(app.jobs.get(id).unwrap().has_output());
        app.process_io(IoMsg::Loaded { project: app.version.project, result: Ok((Project::new("loaded"), PathBuf::from("loaded.foiley"))) });
        app.process_jobs();
        assert!(app.pending.is_empty());
        assert!(!app.jobs.get(id).unwrap().has_output());
        assert_eq!(app.jobs.jobs().len(), 1, "no chained 3D job may be created");
        assert_eq!(app.project.name, "loaded");
        assert!(app.project.foil_polars.is_empty());
        assert!(!app.dirty);
    }

    #[test]
    fn near_re_details_and_visibility_are_independent() {
        let mut app = app();
        let foil = app.project.foils[0].clone();
        for re in [300_000.0, 300_100.0] {
            let settings = AnalysisSettings { re, ..Default::default() };
            let point = FoilPoint { alpha: 2.0, re, cl: re / 1e6, ..Default::default() };
            let polar = FoilPolar { name: settings.polar_name(), foil: foil.name.clone(), foil_hash: Some(foil.geom_hash()), re,
                settings: Some(settings), points: vec![point], ..Default::default() };
            app.apply_output(app.job_generation, "polar", Pending::FoilPolar, JobOutput::FoilPolar { polar, details: vec![PointDetail { point, ..Default::default() }] });
        }
        assert_eq!(app.project.foil_polars.len(), 2);
        assert_eq!(app.fui.visible.len(), 2);
        assert_eq!(app.details.len(), 2);
        let a = app.project.foil_polars[0].clone();
        let b = app.project.foil_polars[1].clone();
        app.show_point(&foil.name, &a.identity().hex(), a.points[0]);
        assert_eq!(app.point_detail.as_ref().unwrap().point.re, 300_000.0);
        app.show_point(&foil.name, &b.identity().hex(), b.points[0]);
        assert_eq!(app.point_detail.as_ref().unwrap().point.re, 300_100.0);
        app.fui.visible.remove(&a.key());
        assert!(app.fui.visible.contains(&b.key()));
    }
}

