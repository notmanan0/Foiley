//! Foil design and 2D analysis modules.

use crate::widgets::{self, check, heading, hint, int, num, section, series_color};
use crate::{App, Pending};
use eframe::egui::{self, Color32, RichText};
use egui_plot::{Legend, Line, Plot, PlotPoints, Points};
use foiley_app::JobKind;
use foiley_foil::analysis::{PointDetail, PolarType};
use foiley_foil::batch::{re_list_geometric, BatchCase};
use foiley_foil::foil::Paneling;
use foiley_foil::{AnalysisSettings, Foil, FoilPoint, Sweep};
use foiley_foil::inverse::{modify, InverseMode, ModKind, Modification, SpeedDist};
use foiley_foil::optim::{Algorithm, Goal, OptimResult, OptimSettings};
use foiley_app::JobId;
use std::collections::HashSet;

pub struct FoilUi {
    pub visible: HashSet<(String, String)>,
    pub settings: AnalysisSettings,
    pub sweep_mode: usize,
    pub range: [[f64; 3]; 4],
    pub point_value: f64,
    pub point_is_alpha: bool,
    pub batch_re: [f64; 2],
    pub batch_count: usize,
    pub batch_mach: String,
    pub batch_ncrit: String,
    pub batch_foils: HashSet<String>,
    pub plots: [(usize, usize); 4],
    pub naca: String,
    pub flap: [f64; 3],
    pub te_gap: [f64; 2],
    pub le: [f64; 2],
    pub tc: [f64; 2],
    pub high: [f64; 2],
    pub repanel: Paneling,
    pub mix_with: String,
    pub mix_frac: f64,
    pub preview: Option<Foil>,
    pub new_name: String,
    pub bl_var: usize,
    pub cp_inviscid: bool,
    pub show_failed: bool,
    pub inv_mode: InverseMode,
    pub inv_alpha: f64,
    pub qdist: Option<SpeedDist>,
    pub inv_target: Option<(Vec<f64>, Option<(usize, usize)>)>,
    pub inv_upper: bool,
    pub inv_x: [f64; 2],
    pub inv_kind: usize,
    pub inv_val: [f64; 2],
    pub inv_iter: usize,
    pub opt: OptimSettings,
    pub opt_job: Option<JobId>,
    pub opt_result: Option<OptimResult>,
}

impl Default for FoilUi {
    fn default() -> Self {
        FoilUi {
            visible: HashSet::new(),
            settings: AnalysisSettings::default(),
            sweep_mode: 0,
            range: [[-5.0, 15.0, 0.5], [0.0, 1.2, 0.1], [1e5, 1e6, 1e5], [-10.0, 20.0, 2.0]],
            point_value: 4.0,
            point_is_alpha: true,
            batch_re: [5e4, 1e6],
            batch_count: 6,
            batch_mach: "0".into(),
            batch_ncrit: "9".into(),
            batch_foils: HashSet::new(),
            plots: [(0, 1), (2, 1), (0, 5), (0, 8)],
            naca: "4412".into(),
            flap: [0.75, 0.5, 10.0],
            te_gap: [0.0, 0.8],
            le: [1.0, 0.1],
            tc: [12.0, 2.0],
            high: [0.3, 0.4],
            repanel: Paneling::default(),
            mix_with: String::new(),
            mix_frac: 0.5,
            preview: None,
            new_name: String::new(),
            bl_var: 0,
            cp_inviscid: true,
            show_failed: true,
            inv_mode: InverseMode::Mixed,
            inv_alpha: 2.0,
            qdist: None,
            inv_target: None,
            inv_upper: true,
            inv_x: [0.3, 0.7],
            inv_kind: 0,
            inv_val: [0.97, 0.02],
            inv_iter: 10,
            opt: OptimSettings::default(),
            opt_job: None,
            opt_result: None,
        }
    }
}

pub const VARS: [&str; 13] = ["alpha (deg)", "CL", "CD", "CDp", "CDf", "Cm", "Xtr top", "Xtr bottom", "CL/CD", "CL^1.5/CD", "Re", "Cp min", "flap (deg)"];

pub fn var(p: &FoilPoint, i: usize) -> f64 {
    match i {
        0 => p.alpha,
        1 => p.cl,
        2 => p.cd,
        3 => p.cdp,
        4 => p.cdf,
        5 => p.cm,
        6 => p.xtr_top,
        7 => p.xtr_bot,
        8 => p.ld(),
        9 => {
            if p.cd > 0.0 {
                p.cl.max(0.0).powf(1.5) / p.cd
            } else {
                0.0
            }
        }
        10 => p.re,
        11 => p.cpmin,
        _ => p.flap,
    }
}

fn foil_lines(f: &Foil) -> Vec<[f64; 2]> {
    f.x.iter().zip(&f.y).map(|(x, y)| [*x, *y]).collect()
}

impl App {
    fn selected_foil(&self) -> Option<Foil> {
        self.sel_foil.as_ref().and_then(|n| self.project.foil(n)).cloned()
    }

    fn foil_picker(&mut self, ui: &mut egui::Ui) {
        let names: Vec<String> = self.project.foils.iter().map(|f| f.name.clone()).collect();
        egui::ComboBox::from_id_salt("foilpick").width(220.0).selected_text(self.sel_foil.clone().unwrap_or_else(|| "(no foil)".into())).show_ui(ui, |ui| {
            for n in names {
                if ui.selectable_label(self.sel_foil.as_deref() == Some(n.as_str()), &n).clicked() {
                    self.sel_foil = Some(n);
                    self.fui.preview = None;
                }
            }
        });
    }

    pub fn foil_design(&mut self, ui: &mut egui::Ui) {
        egui::SidePanel::left("foil_design_left").resizable(true).default_width(340.0).show_inside(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                heading(ui, "Foil");
                self.foil_picker(ui);
                let foil = self.selected_foil();
                if let Some(f) = &foil {
                    let p = f.props();
                    egui::Grid::new("props").num_columns(2).striped(true).show(ui, |ui| {
                        let mut row = |a: &str, b: String| {
                            ui.label(a);
                            ui.label(b);
                            ui.end_row();
                        };
                        row("points", p.points.to_string());
                        row("max thickness", format!("{:.2}% at x/c {:.3}", p.thickness * 100.0, p.x_thickness));
                        row("max camber", format!("{:.2}% at x/c {:.3}", p.camber * 100.0, p.x_camber));
                        row("LE radius", format!("{:.4} c", p.le_radius));
                        row("TE gap / angle", format!("{:.4} c / {:.1} deg", p.te_gap, p.te_angle));
                        row("area", format!("{:.4} c^2", p.area));
                        if let Some((x, y)) = f.hinge {
                            row("control-surface hinge", format!("x/c {x:.2}, y/t {y:.2}"));
                        }
                    });
                    if !f.note.is_empty() {
                        hint(ui, &f.note);
                    }
                }
                ui.separator();
                heading(ui, "Create");
                ui.horizontal(|ui| {
                    ui.label("NACA");
                    ui.add(egui::TextEdit::singleline(&mut self.fui.naca).desired_width(70.0)).on_hover_text("4-digit (e.g. 2412: 2% camber at 40%, 12% thick) or 5-digit (230xx)");
                    if ui.button("Add").clicked() {
                        match Foil::naca(&self.fui.naca) {
                            Ok(mut f) => {
                                self.snapshot();
                                f.name = self.project.unique_foil_name(&f.name);
                                self.sel_foil = Some(f.name.clone());
                                self.project.foils.push(f);
                            }
                            Err(e) => self.notify(e.to_string(), true),
                        }
                    }
                });
                if ui.button("Load .dat files...").on_hover_text("Selig or Lednicer coordinate files; clockwise files are reordered automatically").clicked() {
                    if let Some(files) = rfd::FileDialog::new().add_filter("Airfoil", &["dat", "txt"]).pick_files() {
                        self.snapshot();
                        for p in files {
                            match Foil::load(&p) {
                                Ok(mut f) => {
                                    f.name = self.project.unique_foil_name(&f.name);
                                    self.sel_foil = Some(f.name.clone());
                                    self.project.foils.push(f);
                                }
                                Err(e) => self.notify(format!("{}: {e}", p.display()), true),
                            }
                        }
                    }
                }
                let Some(foil) = foil else { return };
                ui.separator();
                heading(ui, "Modify (preview, then apply)");
                hint(ui, "Each operation shows a preview (orange). Apply replaces the foil (its polars are cleared because they no longer match); 'Save as new' keeps the original.");
                let mut preview: Option<Result<Foil, String>> = None;
                let e = |r: Result<Foil, foiley_foil::FoilError>| r.map_err(|e| e.to_string());
                section(ui, "Trailing-edge flap", true, |ui| {
                    egui::Grid::new("flap").num_columns(2).show(ui, |ui| {
                        num(ui, "hinge x/c", &mut self.fui.flap[0], 0.01, 0.05..=0.98, "chordwise hinge position");
                        num(ui, "hinge y/t", &mut self.fui.flap[1], 0.01, 0.0..=1.0, "hinge height: 0 = lower surface, 1 = upper surface");
                        num(ui, "deflection (deg)", &mut self.fui.flap[2], 0.5, -45.0..=45.0, "positive = trailing edge down");
                    });
                    ui.horizontal(|ui| {
                        if ui.button("Preview").clicked() {
                            preview = Some(e(foil.with_flap(self.fui.flap[0], self.fui.flap[1], self.fui.flap[2])));
                        }
                        if ui.button("Mark as control surface").on_hover_text("Use this foil's hinge for 3D control surfaces: wings using it get deflectable flaps").clicked() {
                            let mut f = foil.clone();
                            f.hinge = Some((self.fui.flap[0], self.fui.flap[1]));
                            preview = Some(Ok(f));
                        }
                    });
                });
                section(ui, "Thickness and camber", false, |ui| {
                    egui::Grid::new("tc").num_columns(2).show(ui, |ui| {
                        num(ui, "max thickness (%c)", &mut self.fui.tc[0], 0.1, 0.5..=40.0, "new maximum thickness");
                        num(ui, "max camber (%c)", &mut self.fui.tc[1], 0.1, -10.0..=15.0, "new maximum camber (cannot be created on a symmetric foil - use interpolation)");
                        num(ui, "x of max thickness", &mut self.fui.high[0], 0.01, 0.05..=0.8, "move the point of maximum thickness");
                        num(ui, "x of max camber", &mut self.fui.high[1], 0.01, 0.05..=0.9, "move the point of maximum camber");
                    });
                    ui.horizontal(|ui| {
                        if ui.button("Set t & c").clicked() {
                            let p = foil.props();
                            let c = if p.camber.abs() > 1e-6 { Some(self.fui.tc[1] / 100.0) } else { None };
                            preview = Some(e(foil.with_tc(Some(self.fui.tc[0] / 100.0), c)));
                        }
                        if ui.button("Move max points").clicked() {
                            preview = Some(e(foil.with_high_points(self.fui.high[0], self.fui.high[1])));
                        }
                    });
                });
                section(ui, "Leading / trailing edge", false, |ui| {
                    egui::Grid::new("lete").num_columns(2).show(ui, |ui| {
                        num(ui, "LE radius factor", &mut self.fui.le[0], 0.01, 0.1..=5.0, "new/old leading-edge radius");
                        num(ui, "LE blend x/c", &mut self.fui.le[1], 0.01, 0.01..=0.5, "distance over which the change is blended");
                        num(ui, "TE gap (c)", &mut self.fui.te_gap[0], 0.001, 0.0..=0.05, "trailing-edge thickness");
                        num(ui, "TE blend x/c", &mut self.fui.te_gap[1], 0.01, 0.1..=1.0, "blending distance from the TE");
                    });
                    ui.horizontal(|ui| {
                        if ui.button("Set LE radius").clicked() {
                            preview = Some(e(foil.with_le_radius(self.fui.le[0], self.fui.le[1])));
                        }
                        if ui.button("Set TE gap").clicked() {
                            preview = Some(e(foil.with_te_gap(self.fui.te_gap[0], self.fui.te_gap[1])));
                        }
                    });
                });
                section(ui, "Paneling and normalisation", false, |ui| {
                    egui::Grid::new("pane").num_columns(2).show(ui, |ui| {
                        int(ui, "panel nodes", &mut self.fui.repanel.panels, 20..=350, "XFoil PANGEN node count (160 default)");
                        num(ui, "LE bunching", &mut self.fui.repanel.le_bunching, 0.05, 0.0..=3.0, "curvature-based bunching strength");
                        num(ui, "TE/LE density", &mut self.fui.repanel.te_bunching, 0.01, 0.01..=1.0, "trailing-edge panel density relative to the LE");
                    });
                    ui.horizontal(|ui| {
                        if ui.button("Repanel").clicked() {
                            preview = Some(e(foil.repaneled(&self.fui.repanel)));
                        }
                        if ui.button("Normalise").on_hover_text("de-rotate, unit chord, LE at origin").clicked() {
                            preview = Some(Ok(foil.normalized()));
                        }
                    });
                });
                section(ui, "Interpolate with another foil", false, |ui| {
                    let names: Vec<String> = self.project.foils.iter().map(|f| f.name.clone()).filter(|n| *n != foil.name).collect();
                    egui::ComboBox::from_id_salt("mixwith").selected_text(if self.fui.mix_with.is_empty() { "(choose)" } else { &self.fui.mix_with }).show_ui(ui, |ui| {
                        for n in names {
                            ui.selectable_value(&mut self.fui.mix_with, n.clone(), n);
                        }
                    });
                    ui.add(egui::Slider::new(&mut self.fui.mix_frac, 0.0..=1.0).text("fraction of the other foil"));
                    if ui.button("Preview mix").clicked() {
                        if let Some(o) = self.project.foil(&self.fui.mix_with) {
                            preview = Some(e(Foil::interpolate(&foil, o, self.fui.mix_frac, &foil.name)));
                        }
                    }
                });
                section(ui, "Inverse design (QDES / MDES)", false, |ui| {
                    hint(ui, "Change the target surface speed and solve for the foil that produces it (XFoil's mixed and full inverse methods).");
                    ui.horizontal(|ui| {
                        ui.selectable_value(&mut self.fui.inv_mode, InverseMode::Mixed, "mixed (one segment)");
                        ui.selectable_value(&mut self.fui.inv_mode, InverseMode::Full, "full (whole surface)");
                    });
                    egui::Grid::new("inv").num_columns(2).show(ui, |ui| {
                        num(ui, "alpha (deg)", &mut self.fui.inv_alpha, 0.1, -20.0..=20.0, "design angle of attack (inviscid speed distribution)");
                    });
                    if ui.button("Get speed distribution").clicked() {
                        self.fui.inv_target = None;
                        let k = JobKind::Inverse { foil: foil.clone(), alpha: self.fui.inv_alpha, mode: self.fui.inv_mode, paneling: self.fui.repanel, target: None, iterations: 0 };
                        self.submit(format!("{} speed distribution", foil.name), k, Pending::Inverse);
                    }
                    if let Some(d) = self.fui.qdist.clone() {
                        ui.horizontal(|ui| {
                            ui.selectable_value(&mut self.fui.inv_upper, true, "upper surface");
                            ui.selectable_value(&mut self.fui.inv_upper, false, "lower surface");
                        });
                        egui::Grid::new("invmod").num_columns(2).show(ui, |ui| {
                            num(ui, "segment from x/c", &mut self.fui.inv_x[0], 0.01, 0.0..=1.0, "start of the modified segment");
                            num(ui, "segment to x/c", &mut self.fui.inv_x[1], 0.01, 0.0..=1.0, "end of the modified segment");
                            ui.label("change");
                            widgets::pick(ui, "invkind", &mut self.fui.inv_kind, &["scale speed", "add bump", "smooth"]);
                            ui.end_row();
                            if self.fui.inv_kind < 2 {
                                let (lbl, sp) = if self.fui.inv_kind == 0 { ("factor", 0.002) } else { ("bump (Q/Vinf)", 0.002) };
                                num(ui, lbl, &mut self.fui.inv_val[self.fui.inv_kind], sp, -2.0..=3.0, "scale: 0.97 slows the flow by 3%; bump: adds a smooth speed bump");
                            }
                            int(ui, "iterations", &mut self.fui.inv_iter, 1..=100, "mixed-inverse Newton iterations");
                        });
                        let kind = match self.fui.inv_kind {
                            0 => ModKind::Scale(self.fui.inv_val[0]),
                            1 => ModKind::Bump(self.fui.inv_val[1]),
                            _ => ModKind::Smooth,
                        };
                        let m = Modification { x_from: self.fui.inv_x[0], x_to: self.fui.inv_x[1], upper: self.fui.inv_upper, kind };
                        let (t, seg) = modify(&d, &m);
                        self.fui.inv_target = Some((t.clone(), seg));
                        if seg.is_none() {
                            ui.label(RichText::new("the segment contains fewer than 2 points").color(Color32::from_rgb(230, 80, 70)));
                        } else if ui.button("Solve for this target").clicked() {
                            let k = JobKind::Inverse { foil: foil.clone(), alpha: self.fui.inv_alpha, mode: self.fui.inv_mode, paneling: self.fui.repanel, target: Some((t, seg)), iterations: self.fui.inv_iter };
                            self.submit(format!("{} inverse design", foil.name), k, Pending::Inverse);
                        }
                    }
                });
                section(ui, "Optimiser (MOPSO / genetic)", false, |ui| {
                    hint(ui, "Hicks-Henne bumps on both surfaces (and optionally a flap) are optimised with XFoil at one operating point. Every evaluation runs in the background.");
                    let o = &mut self.fui.opt;
                    ui.horizontal(|ui| {
                        ui.selectable_value(&mut o.algorithm, Algorithm::Mopso, "particle swarm (MOPSO)");
                        ui.selectable_value(&mut o.algorithm, Algorithm::Genetic, "genetic");
                    });
                    egui::Grid::new("opt").num_columns(2).show(ui, |ui| {
                        num(ui, "alpha (deg)", &mut o.alpha, 0.1, -10.0..=20.0, "operating point angle of attack");
                        let sp = o.analysis.re * 0.01;
                        num(ui, "Reynolds", &mut o.analysis.re, sp, 1e4..=1e8, "");
                        num(ui, "Ncrit", &mut o.analysis.ncrit, 0.1, 0.5..=15.0, "");
                        int(ui, "population", &mut o.population, 4..=100, "candidates per iteration (7-15 recommended)");
                        int(ui, "iterations", &mut o.iterations, 1..=500, "");
                        int(ui, "bumps per side", &mut o.bumps_per_side, 1..=10, "Hicks-Henne functions on each surface");
                        num(ui, "max bump (c)", &mut o.max_amplitude, 0.001, 0.0005..=0.05, "maximum bump amplitude, fraction of chord");
                        let mut fl = o.flap.is_some();
                        if check(ui, "optimise a TE flap angle", &mut fl, "adds the flap deflection as a design variable") {
                            o.flap = if fl { Some((0.0, 10.0)) } else { None };
                        }
                        if let Some((a, b)) = &mut o.flap {
                            num(ui, "  flap min (deg)", a, 0.5, -30.0..=30.0, "");
                            num(ui, "  flap max (deg)", b, 0.5, -30.0..=30.0, "");
                        }
                    });
                    ui.label(RichText::new("objectives").strong());
                    for (k, spec) in o.objectives.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.checkbox(&mut spec.enabled, spec.objective.label());
                            let mut gi = match spec.goal {
                                Goal::Minimise => 0,
                                Goal::Maximise => 1,
                                Goal::Target { .. } => 2,
                            };
                            if widgets::pick(ui, &format!("goal{k}"), &mut gi, &["minimise", "maximise", "target"]) {
                                spec.goal = match gi {
                                    0 => Goal::Minimise,
                                    1 => Goal::Maximise,
                                    _ => Goal::Target { value: 0.0, tolerance: 1e-3 },
                                };
                            }
                            if let Goal::Target { value, tolerance } = &mut spec.goal {
                                ui.add(egui::DragValue::new(value).speed(0.001).prefix("= "));
                                ui.add(egui::DragValue::new(tolerance).speed(0.0001).prefix("+/- "));
                            }
                        });
                    }
                    let running = self.fui.opt_job.and_then(|id| self.jobs.get(id)).map(|j| !j.status().is_finished()).unwrap_or(false);
                    ui.add_enabled_ui(!running, |ui| {
                        if ui.button(RichText::new("Run optimiser").strong()).clicked() {
                            let k = JobKind::Optimise { foil: foil.clone(), settings: self.fui.opt.clone() };
                            let j = self.jobs.submit(format!("optimise {}", foil.name), k);
                            self.fui.opt_job = Some(j.id);
                            self.fui.opt_result = None;
                            self.pending.push((j.id, Pending::Optimise));
                        }
                    });
                    let live = self.fui.opt_job.and_then(|id| self.jobs.get(id)).and_then(|j| j.partial.lock().unwrap().optim.clone()).or_else(|| self.fui.opt_result.clone());
                    if let Some(r) = live {
                        if let Some(b) = &r.best {
                            ui.label(format!("best: CL {:.4}  CD {:.5}  CL/CD {:.1}  Cm {:.4}  ({} evaluations, {} Pareto)", b.cl, b.cd, b.cl / b.cd.max(1e-12), b.cm, r.evaluations, r.pareto.len()));
                            if ui.button("Preview best design").clicked() {
                                let mut f = b.foil.clone();
                                f.note = format!("optimised: CL {:.3} CD {:.5}", b.cl, b.cd);
                                preview = Some(Ok(f));
                            }
                        }
                    }
                });
                if let Some(r) = preview {
                    match r {
                        Ok(p) => {
                            self.fui.new_name = self.project.unique_foil_name(&format!("{} mod", foil.name));
                            self.fui.preview = Some(p);
                        }
                        Err(err) => self.notify(err, true),
                    }
                }
                if let Some(pv) = self.fui.preview.clone() {
                    ui.separator();
                    ui.label(RichText::new(format!("Preview: {}", pv.note)).color(Color32::from_rgb(255, 160, 60)));
                    ui.horizontal(|ui| {
                        if ui.button("Replace foil").clicked() {
                            self.snapshot();
                            let mut f = pv.clone();
                            f.name = foil.name.clone();
                            self.project.upsert_foil(f);
                            self.fui.preview = None;
                        }
                        ui.add(egui::TextEdit::singleline(&mut self.fui.new_name).desired_width(110.0));
                        if ui.button("Save as new").clicked() {
                            self.snapshot();
                            let mut f = pv.clone();
                            f.name = self.project.unique_foil_name(&self.fui.new_name);
                            self.sel_foil = Some(f.name.clone());
                            self.project.foils.push(f);
                            self.fui.preview = None;
                        }
                        if ui.button("Discard").clicked() {
                            self.fui.preview = None;
                        }
                    });
                }
                ui.separator();
                ui.horizontal(|ui| {
                    if ui.button("Export .dat...").clicked() {
                        if let Some(p) = rfd::FileDialog::new().set_file_name(format!("{}.dat", foil.name)).save_file() {
                            if let Err(e) = foil.save(&p) {
                                self.notify(e.to_string(), true);
                            }
                        }
                    }
                    let mut name = foil.name.clone();
                    if ui.add(egui::TextEdit::singleline(&mut name).desired_width(120.0)).on_hover_text("rename (updates polars and wings using it)").lost_focus() && name != foil.name && !name.trim().is_empty() && self.project.foil(&name).is_none() {
                        self.snapshot();
                        let old = foil.name.clone();
                        if let Some(f) = self.project.foils.iter_mut().find(|f| f.name == old) {
                            f.name = name.clone();
                        }
                        for p in &mut self.project.foil_polars {
                            if p.foil == old {
                                p.foil = name.clone();
                            }
                        }
                        for pl in &mut self.project.planes {
                            for w in &mut pl.wings {
                                for s in &mut w.sections {
                                    if s.foil_left == old {
                                        s.foil_left = name.clone();
                                    }
                                    if s.foil_right == old {
                                        s.foil_right = name.clone();
                                    }
                                }
                            }
                        }
                        self.sel_foil = Some(name);
                    }
                });
            });
        });
        egui::CentralPanel::default().show_inside(ui, |ui| {
            let foil = self.selected_foil();
            let pv = self.fui.preview.clone();
            let live = self.fui.opt_job.and_then(|id| self.jobs.get(id)).and_then(|j| j.partial.lock().unwrap().optim.clone()).or_else(|| self.fui.opt_result.clone());
            let extra = self.fui.qdist.is_some() || live.is_some();
            let h = ui.available_height();
            Plot::new("foilshape").data_aspect(1.0).height(if extra { h * 0.55 } else { h }).legend(Legend::default()).show(ui, |pu| {
                if let Some(best) = live.as_ref().and_then(|r| r.best.as_ref()) {
                    pu.line(Line::new("optimiser best", PlotPoints::from(foil_lines(&best.foil))).color(Color32::from_rgb(90, 220, 120)).width(1.5));
                }
                if let Some(f) = &foil {
                    pu.line(Line::new(f.name.clone(), PlotPoints::from(foil_lines(f))).color(series_color(0)).width(2.0));
                    let (xs, c, _) = f.camber_thickness(60);
                    pu.line(Line::new("camber line", PlotPoints::from(xs.iter().zip(&c).map(|(x, y)| [*x, *y]).collect::<Vec<_>>())).color(Color32::GRAY).style(egui_plot::LineStyle::dashed_dense()));
                    pu.points(Points::new("nodes", PlotPoints::from(foil_lines(f))).radius(1.5).color(series_color(0)));
                    if let Some((x, yr)) = f.hinge {
                        let (u, l) = f.surface_y(x);
                        pu.points(Points::new("hinge", PlotPoints::from(vec![[x, l + yr * (u - l)]])).radius(5.0).color(Color32::YELLOW));
                    }
                }
                if let Some(p) = &pv {
                    pu.line(Line::new("preview", PlotPoints::from(foil_lines(p))).color(Color32::from_rgb(255, 160, 60)).width(2.0));
                }
            });
            if extra {
                ui.horizontal(|ui| {
                    let w = ui.available_width();
                    let hh = ui.available_height();
                    if let Some(d) = &self.fui.qdist {
                        ui.vertical(|ui| {
                            ui.set_width(if live.is_some() { w * 0.5 } else { w });
                            ui.label(RichText::new(format!("Surface speed Q/Vinf at alpha {:.2} (CL {:.3})", d.alpha, d.cl)).strong());
                            Plot::new("qdist").height(hh - 24.0).legend(Legend::default()).x_axis_label("x/c").show(ui, |pu| {
                                let side = |upper: bool, q: &[f64]| -> Vec<[f64; 2]> { (0..d.q.len()).filter(|&i| d.is_upper(i) == upper).map(|i| [d.x[i], q[i].abs()]).collect() };
                                pu.line(Line::new("upper", PlotPoints::from(side(true, &d.q))).color(series_color(0)));
                                pu.line(Line::new("lower", PlotPoints::from(side(false, &d.q))).color(series_color(1)));
                                if let Some((t, _)) = &self.fui.inv_target {
                                    if t.len() == d.q.len() {
                                        pu.line(Line::new("target", PlotPoints::from(side(self.fui.inv_upper, t))).color(Color32::from_rgb(255, 160, 60)).style(egui_plot::LineStyle::dashed_dense()));
                                    }
                                }
                            });
                        });
                    }
                    if let Some(r) = &live {
                        ui.vertical(|ui| {
                            ui.label(RichText::new("Optimiser: best score per iteration").strong());
                            Plot::new("opthist").height(hh - 24.0).x_axis_label("iteration").show(ui, |pu| {
                                let pts: Vec<[f64; 2]> = r.history.iter().map(|(i, v)| [*i as f64, *v]).collect();
                                pu.line(Line::new("best", PlotPoints::from(pts.clone())).color(series_color(2)));
                                pu.points(Points::new("best", PlotPoints::from(pts)).radius(3.0).color(series_color(2)));
                            });
                        });
                    }
                });
            }
        });
    }

    pub fn foil_analysis(&mut self, ui: &mut egui::Ui) {
        egui::SidePanel::left("foil_an_left").resizable(true).default_width(330.0).show_inside(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                heading(ui, "Foil");
                self.foil_picker(ui);
                ui.separator();
                heading(ui, "Analysis settings");
                let s = &mut self.fui.settings;
                let mut pt = PolarType::ALL.iter().position(|t| *t == s.polar_type).unwrap_or(0);
                let labels: Vec<&str> = PolarType::ALL.iter().map(|t| t.label()).collect();
                if widgets::pick(ui, "ptype", &mut pt, &labels) {
                    s.polar_type = PolarType::ALL[pt];
                    self.fui.sweep_mode = match s.polar_type {
                        PolarType::FixedAlpha => 2,
                        PolarType::FlapSweep => 3,
                        _ if self.fui.sweep_mode >= 2 => 0,
                        _ => self.fui.sweep_mode,
                    };
                }
                egui::Grid::new("xset").num_columns(2).show(ui, |ui| {
                    let re_label = match s.polar_type {
                        PolarType::FixedLift => "Re*sqrt(CL)",
                        PolarType::RubberChord => "Re*CL",
                        _ => "Reynolds",
                    };
                    if s.polar_type != PolarType::FixedAlpha {
                        let sp = s.re * 0.01;
                        num(ui, re_label, &mut s.re, sp, 1e3..=1e9, "chord Reynolds number V c / nu (for type 2/3 the product held constant)");
                    }
                    num(ui, "Mach", &mut s.mach, 0.01, 0.0..=0.9, "freestream Mach number (Karman-Tsien compressibility correction)");
                    num(ui, "Ncrit", &mut s.ncrit, 0.1, 0.1..=20.0, "e^N transition criterion: 9 = average wind tunnel, 11-12 = quiet free flight, gliders, 4-7 = turbulent or rough conditions");
                    num(ui, "trip top x/c", &mut s.xtr_top, 0.01, 0.0..=1.0, "forced transition on the upper surface (1 = free transition)");
                    num(ui, "trip bottom x/c", &mut s.xtr_bot, 0.01, 0.0..=1.0, "forced transition on the lower surface (1 = free transition)");
                    if matches!(s.polar_type, PolarType::FixedAlpha | PolarType::FlapSweep) {
                        num(ui, "alpha (deg)", &mut s.alpha_fixed, 0.1, -30.0..=30.0, "fixed angle of attack");
                    }
                    if s.polar_type == PolarType::FlapSweep {
                        num(ui, "hinge x/c", &mut s.flap_hinge.0, 0.01, 0.05..=0.98, "flap hinge position");
                        num(ui, "hinge y/t", &mut s.flap_hinge.1, 0.01, 0.0..=1.0, "flap hinge height (0 lower .. 1 upper surface)");
                    }
                    check(ui, "viscous (boundary layer)", &mut s.viscous, "unchecked = inviscid panel solution only");
                });
                ui.separator();
                heading(ui, "Sweep");
                let modes = ["alpha", "CL", "Re (type 4)", "flap (type 6)"];
                ui.horizontal(|ui| {
                    for (i, m) in modes.iter().enumerate() {
                        let allowed = match self.fui.settings.polar_type {
                            PolarType::FixedAlpha => i == 2,
                            PolarType::FlapSweep => i == 3,
                            _ => i < 2,
                        };
                        ui.add_enabled_ui(allowed, |ui| {
                            ui.selectable_value(&mut self.fui.sweep_mode, i, *m);
                        });
                    }
                });
                let m = self.fui.sweep_mode;
                egui::Grid::new("sweep").num_columns(2).show(ui, |ui| {
                    let r = &mut self.fui.range[m];
                    let sp = if m == 2 { 1e4 } else { 0.1 };
                    num(ui, "from", &mut r[0], sp, -1e9..=1e9, "first value");
                    num(ui, "to", &mut r[1], sp, -1e9..=1e9, "last value");
                    num(ui, "step", &mut r[2], sp * 0.1, 1e-6..=1e9, "increment");
                });
                let r = self.fui.range[m];
                let npts = ((r[1] - r[0]).abs() / r[2].max(1e-12)).floor() as usize + 1;
                hint(ui, &format!("{npts} operating points"));
                section(ui, "Solver options", false, |ui| {
                    let s = &mut self.fui.settings;
                    egui::Grid::new("solver").num_columns(2).show(ui, |ui| {
                        int(ui, "iterations", &mut s.max_iter, 10..=1000, "Newton iteration limit per point (XFoil ITER)");
                        int(ui, "panel nodes", &mut s.paneling.panels, 40..=350, "airfoil panel nodes (XFoil PANE)");
                        num(ui, "VACC", &mut s.vaccel, 0.001, 0.0..=0.1, "viscous acceleration parameter (XFoil VACC, default 0.01)");
                        check(ui, "recover unconverged points", &mut s.recover, "re-initialise the boundary layer and march in sub-steps from a converged neighbour when a point fails");
                        check(ui, "sweep outwards from zero lift", &mut s.from_zero, "every point warm-starts from an attached-flow neighbour (more robust near stall)");
                        check(ui, "parallel chunks", &mut s.parallel, "split long sweeps over several cores");
                        check(ui, "keep Cp / BL for every point", &mut s.keep_details, "store distributions so clicking a polar point shows them instantly");
                        num(ui, "discard CD below", &mut s.min_cd, 0.0001, 0.0..=0.01, "drop converged points with an implausibly low drag (0 = keep all)");
                    });
                });
                ui.separator();
                let foil = self.selected_foil();
                ui.add_enabled_ui(foil.is_some(), |ui| {
                    if ui.add(egui::Button::new(RichText::new("Run polar").strong()).min_size(egui::vec2(140.0, 28.0))).clicked() {
                        let f = foil.clone().unwrap();
                        let r = self.fui.range[m];
                        let sweep = match m {
                            0 => Sweep::Alpha { from: r[0], to: r[1], step: r[2] },
                            1 => Sweep::Cl { from: r[0], to: r[1], step: r[2] },
                            2 => Sweep::Re { from: r[0], to: r[1], step: r[2] },
                            _ => Sweep::Flap { from: r[0], to: r[1], step: r[2] },
                        };
                        let st = self.fui.settings.clone();
                        if let Err(e) = st.validate() {
                            self.notify(e, true);
                        } else {
                            let t = format!("{} {}", f.name, st.polar_name());
                            self.submit(t, JobKind::FoilPolar { foil: f, settings: st, sweep }, Pending::FoilPolar);
                        }
                    }
                    ui.horizontal(|ui| {
                        ui.selectable_value(&mut self.fui.point_is_alpha, true, "alpha");
                        ui.selectable_value(&mut self.fui.point_is_alpha, false, "CL");
                        ui.add(egui::DragValue::new(&mut self.fui.point_value).speed(0.1));
                        if ui.button("Run point").on_hover_text("single operating point with Cp and boundary-layer plots").clicked() {
                            let f = foil.clone().unwrap();
                            let t = format!("{} {} = {}", f.name, if self.fui.point_is_alpha { "alpha" } else { "CL" }, self.fui.point_value);
                            self.submit(t, JobKind::FoilPoint { foil: f, settings: self.fui.settings.clone(), value: self.fui.point_value, is_alpha: self.fui.point_is_alpha }, Pending::FoilPoint);
                        }
                    });
                });
                ui.separator();
                section(ui, "Batch polar mesh", false, |ui| {
                    hint(ui, "Polars for several foils and Reynolds numbers at once (used by 3D viscous analyses). Uses the alpha sweep and solver options above.");
                    egui::Grid::new("batch").num_columns(2).show(ui, |ui| {
                        num(ui, "Re from", &mut self.fui.batch_re[0], 1e3, 1e3..=1e9, "lowest Reynolds number");
                        num(ui, "Re to", &mut self.fui.batch_re[1], 1e4, 1e3..=1e9, "highest Reynolds number");
                        int(ui, "count (geometric)", &mut self.fui.batch_count, 1..=usize::MAX, "number of Reynolds numbers, geometrically spaced (rounded to 1000; duplicates are dropped)");
                        ui.label("Mach list");
                        ui.text_edit_singleline(&mut self.fui.batch_mach);
                        ui.end_row();
                        ui.label("Ncrit list");
                        ui.text_edit_singleline(&mut self.fui.batch_ncrit);
                        ui.end_row();
                    });
                    ui.label("foils:");
                    let names: Vec<String> = self.project.foils.iter().map(|f| f.name.clone()).collect();
                    ui.horizontal_wrapped(|ui| {
                        for n in names {
                            let mut on = self.fui.batch_foils.contains(&n);
                            if ui.checkbox(&mut on, &n).changed() {
                                if on {
                                    self.fui.batch_foils.insert(n);
                                } else {
                                    self.fui.batch_foils.remove(&n);
                                }
                            }
                        }
                    });
                    let parse = |s: &str| -> Vec<f64> { s.split(',').filter_map(|t| t.trim().parse().ok()).collect() };
                    let res = re_list_geometric(self.fui.batch_re[0], self.fui.batch_re[1], self.fui.batch_count);
                    let (machs, ncs) = (parse(&self.fui.batch_mach), parse(&self.fui.batch_ncrit));
                    let foils: Vec<Foil> = self.project.foils.iter().filter(|f| self.fui.batch_foils.contains(&f.name)).cloned().collect();
                    let n = foils.len() * res.len() * machs.len().max(1) * ncs.len().max(1);
                    if ui.add_enabled(n > 0, egui::Button::new(format!("Run batch ({n} polars)"))).clicked() {
                        let mut cases = vec![];
                        for &re in &res {
                            for &m in if machs.is_empty() { &[0.0][..] } else { &machs[..] } {
                                for &nc in if ncs.is_empty() { &[9.0][..] } else { &ncs[..] } {
                                    cases.push(BatchCase { re, mach: m, ncrit: nc });
                                }
                            }
                        }
                        let r = self.fui.range[0];
                        let st = AnalysisSettings { polar_type: PolarType::FixedSpeed, ..self.fui.settings.clone() };
                        self.submit(format!("batch: {n} polars"), JobKind::Batch { foils, cases, settings: st, sweep: Sweep::Alpha { from: r[0], to: r[1], step: r[2] } }, Pending::Batch);
                    }
                });
            });
        });
        // polar list
        egui::SidePanel::right("polar_list").resizable(true).default_width(250.0).show_inside(ui, |ui| {
            heading(ui, "Polars");
            ui.horizontal(|ui| {
                if ui.small_button("show all").clicked() {
                    for p in &self.project.foil_polars {
                        self.fui.visible.insert((p.foil.clone(), p.name.clone()));
                    }
                }
                if ui.small_button("hide all").clicked() {
                    self.fui.visible.clear();
                }
                if ui.small_button("only this foil").clicked() {
                    self.fui.visible.clear();
                    for p in &self.project.foil_polars {
                        if Some(&p.foil) == self.sel_foil.as_ref() {
                            self.fui.visible.insert((p.foil.clone(), p.name.clone()));
                        }
                    }
                }
            });
            ui.checkbox(&mut self.fui.show_failed, "mark unconverged points");
            egui::ScrollArea::vertical().show(ui, |ui| {
                let mut del = None;
                let mut vis_idx = 0;
                for (i, p) in self.project.foil_polars.iter().enumerate() {
                    let key = (p.foil.clone(), p.name.clone());
                    let mut on = self.fui.visible.contains(&key);
                    let swatch = if on {
                        vis_idx += 1;
                        Some(series_color(vis_idx - 1))
                    } else {
                        None
                    };
                    ui.horizontal(|ui| {
                        let (r, painter) = ui.allocate_painter(egui::vec2(10.0, 10.0), egui::Sense::hover());
                        if let Some(c) = swatch {
                            painter.rect_filled(r.rect, 2.0, c);
                        }
                        if ui.checkbox(&mut on, "").changed() {
                            if on {
                                self.fui.visible.insert(key.clone());
                            } else {
                                self.fui.visible.remove(&key);
                            }
                        }
                        let r = ui.label(RichText::new(format!("{}\n{}", p.foil, p.name)).small());
                        r.on_hover_text(format!("{} points, {} not converged\n{}", p.points.len(), p.failed.len(), if p.failed.is_empty() { String::new() } else { format!("not converged at {:?}", p.failed) }));
                        ui.menu_button("...", |ui| {
                            if ui.button("Export CSV").clicked() {
                                if let Some(path) = rfd::FileDialog::new().set_file_name(format!("{}_{}.csv", p.foil, p.name)).save_file() {
                                    let _ = std::fs::write(path, p.to_csv());
                                }
                                ui.close();
                            }
                            if ui.button("Export XFoil text").clicked() {
                                if let Some(path) = rfd::FileDialog::new().set_file_name(format!("{}_{}.txt", p.foil, p.name)).save_file() {
                                    let _ = std::fs::write(path, p.to_xfoil_text());
                                }
                                ui.close();
                            }
                            if ui.button("Delete").clicked() {
                                del = Some(i);
                                ui.close();
                            }
                        });
                    });
                }
                if let Some(i) = del {
                    self.snapshot();
                    self.project.foil_polars.remove(i);
                }
            });
        });
        egui::CentralPanel::default().show_inside(ui, |ui| {
            let has_detail = self.point_detail.is_some();
            let h = ui.available_height();
            let ph = if has_detail { h * 0.58 } else { h };
            let w = ui.available_width();
            let polars: Vec<_> = self.project.foil_polars.iter().filter(|p| self.fui.visible.contains(&(p.foil.clone(), p.name.clone()))).cloned().collect();
            let mut clicked: Option<(String, String, f64)> = None;
            egui::Grid::new("polarplots").num_columns(2).spacing([6.0, 6.0]).show(ui, |ui| {
                for k in 0..4 {
                    ui.vertical(|ui| {
                        ui.set_width(w / 2.0 - 8.0);
                        ui.horizontal(|ui| {
                            let (mut x, mut y) = self.fui.plots[k];
                            widgets::pick(ui, &format!("py{k}"), &mut y, &VARS);
                            ui.label("vs");
                            widgets::pick(ui, &format!("px{k}"), &mut x, &VARS);
                            self.fui.plots[k] = (x, y);
                        });
                        let (xv, yv) = self.fui.plots[k];
                        let mut plot = Plot::new(format!("polar{k}")).height(ph / 2.0 - 40.0);
                        if polars.len() <= 5 {
                            plot = plot.legend(Legend::default().position(egui_plot::Corner::RightBottom));
                        }
                        let r = plot.x_axis_label(VARS[xv]).y_axis_label(VARS[yv]).show(ui, |pu| {
                            for (i, p) in polars.iter().enumerate() {
                                let pts: Vec<[f64; 2]> = p.points.iter().map(|q| [var(q, xv), var(q, yv)]).collect();
                                let name = format!("{} {}", p.foil, p.name);
                                pu.line(Line::new(name.clone(), PlotPoints::from(pts.clone())).color(series_color(i)));
                                pu.points(Points::new(name.clone(), PlotPoints::from(pts)).radius(2.0).color(series_color(i)));
                                if self.fui.show_failed && xv == 0 && !p.failed.is_empty() && p.polar_type <= 3 {
                                    // unconverged alphas as crosses on the x axis value line
                                    let lo = p.points.iter().map(|q| var(q, yv)).fold(f64::MAX, f64::min);
                                    let fpts: Vec<[f64; 2]> = p.failed.iter().map(|a| [*a, lo]).collect();
                                    pu.points(Points::new(format!("{name} not converged"), PlotPoints::from(fpts)).shape(egui_plot::MarkerShape::Cross).radius(5.0).color(Color32::from_rgb(230, 80, 70)));
                                }
                            }
                            if pu.response().clicked() {
                                if let Some(c) = pu.pointer_coordinate() {
                                    let b = pu.plot_bounds();
                                    let (sx, sy) = (b.width().max(1e-12), b.height().max(1e-12));
                                    let mut best = (f64::MAX, None);
                                    for p in &polars {
                                        for q in &p.points {
                                            let d = ((var(q, xv) - c.x) / sx).powi(2) + ((var(q, yv) - c.y) / sy).powi(2);
                                            if d < best.0 {
                                                best = (d, Some((p.foil.clone(), p.name.clone(), q.alpha)));
                                            }
                                        }
                                    }
                                    if best.0 < 0.002 {
                                        clicked = best.1;
                                    }
                                }
                            }
                        });
                        r.response.on_hover_text("click a point to show its Cp and boundary layer");
                    });
                    if k % 2 == 1 {
                        ui.end_row();
                    }
                }
            });
            if let Some((foil, pol, alpha)) = clicked {
                self.show_point(&foil, &pol, alpha);
            }
            if let Some(d) = self.point_detail.clone() {
                ui.separator();
                self.detail_plots(ui, &d, w);
            }
        });
    }

    /// Show the distributions of a polar point (stored, or recomputed in the background).
    fn show_point(&mut self, foil: &str, polar: &str, alpha: f64) {
        if let Some(ds) = self.details.get(&(foil.to_string(), polar.to_string())) {
            if let Some(d) = ds.iter().min_by(|a, b| (a.point.alpha - alpha).abs().total_cmp(&(b.point.alpha - alpha).abs())) {
                self.point_detail = Some(d.clone());
                return;
            }
        }
        let Some(p) = self.project.foil_polars.iter().find(|p| p.foil == foil && p.name == polar).cloned() else { return };
        let Some(f) = self.project.foil(foil).cloned() else { return };
        let st = AnalysisSettings {
            polar_type: match p.polar_type {
                2 => PolarType::FixedLift,
                3 => PolarType::RubberChord,
                _ => PolarType::FixedSpeed,
            },
            re: p.re,
            mach: p.mach,
            ncrit: p.ncrit,
            xtr_top: p.xtr_top,
            xtr_bot: p.xtr_bot,
            viscous: p.viscous,
            ..self.fui.settings.clone()
        };
        self.submit(format!("{foil} alpha = {alpha:.2} ({polar})"), JobKind::FoilPoint { foil: f, settings: st, value: alpha, is_alpha: true }, Pending::FoilPoint);
    }

    fn detail_plots(&mut self, ui: &mut egui::Ui, d: &PointDetail, w: f32) {
        let q = &d.point;
        ui.horizontal(|ui| {
            ui.label(RichText::new(format!("alpha {:.2}  CL {:.4}  CD {:.5}  Cm {:.4}  L/D {:.1}  Xtr {:.3}/{:.3}  Re {:.0}", q.alpha, q.cl, q.cd, q.cm, q.ld(), q.xtr_top, q.xtr_bot, q.re)).strong());
            ui.checkbox(&mut self.fui.cp_inviscid, "inviscid Cp");
            ui.label("BL variable:");
            widgets::pick(ui, "blvar", &mut self.fui.bl_var, &["Ue/Vinf", "shape factor H", "Cf", "displacement thickness", "momentum thickness", "N / sqrt(Ctau)"]);
            if ui.small_button("close").clicked() {
                self.point_detail = None;
            }
        });
        let h = ui.available_height() - 4.0;
        ui.horizontal(|ui| {
            ui.vertical(|ui| {
                ui.set_width(w / 2.0 - 6.0);
                Plot::new("cp").height(h).legend(Legend::default()).y_axis_label("-Cp").x_axis_label("x/c").show(ui, |pu| {
                    let pts: Vec<[f64; 2]> = d.x.iter().zip(&d.cp_visc).map(|(x, c)| [*x, -c]).collect();
                    pu.line(Line::new("-Cp viscous", PlotPoints::from(pts)).color(series_color(0)).width(2.0));
                    if self.fui.cp_inviscid {
                        let pts: Vec<[f64; 2]> = d.x.iter().zip(&d.cp_inv).map(|(x, c)| [*x, -c]).collect();
                        pu.line(Line::new("-Cp inviscid", PlotPoints::from(pts)).color(Color32::GRAY).style(egui_plot::LineStyle::dashed_loose()));
                    }
                });
            });
            ui.vertical(|ui| {
                ui.set_width(w / 2.0 - 6.0);
                let v = self.fui.bl_var;
                Plot::new("bl").height(h).legend(Legend::default()).x_axis_label("x/c").show(ui, |pu| {
                    for (i, (name, b)) in [("top", &d.top), ("bottom", &d.bot), ("wake", &d.wake)].iter().enumerate() {
                        let ys = match v {
                            0 => &b.ue,
                            1 => &b.h,
                            2 => &b.cf,
                            3 => &b.dstar,
                            4 => &b.theta,
                            _ => &b.n_or_ctau,
                        };
                        let pts: Vec<[f64; 2]> = b.x.iter().zip(ys).map(|(x, y)| [*x, *y]).collect();
                        pu.line(Line::new(*name, PlotPoints::from(pts)).color(series_color(i)));
                    }
                });
            });
        });
    }
}
