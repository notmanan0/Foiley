//! Plane design, 3D analysis and stability modules.

use crate::view3d::{self, Camera, ColorField, Overlay};
use crate::widgets::{self, check, heading, hint, int, num, num_u, section, series_color};
use crate::{App, Pending};
use eframe::egui::{self, Color32, RichText};
use egui_plot::{Legend, Line, MarkerShape, Plot, PlotPoints, Points};
use foiley_app::{auto, JobId, JobKind};
use foiley_core::atmosphere::{isa, Fluid};
use foiley_core::{Spacing, Vec3};
use foiley_plane::analysis::{ControlRange, ControlTarget, ExtraDrag, Method, PlaneOpp, PolarKind, RefDims, Setting};
use foiley_plane::geom::{Fuse, FuseFrame, PointMass, Wing, WingKind, WingSection};
use foiley_plane::mesh::{Mesh, MeshOptions};
use foiley_plane::stability::{modal_initial_state, time_response};
use foiley_plane::{Plane, PlanePolar, PlanePolarDef};
use std::collections::HashSet;
use std::hash::{Hash, Hasher};

pub struct PlaneUi {
    pub sel_wing: usize,
    pub design_tab: usize,
    pub def: PlanePolarDef,
    pub altitude: f64,
    pub water: bool,
    pub auto_polars: usize,
    pub visible: HashSet<(String, String)>,
    pub sel_polar: Option<(String, String)>,
    pub sel_point: usize,
    pub plots: [(usize, usize); 4],
    pub span_var: usize,
    pub cam: Camera,
    pub mesh: Option<(u64, Mesh)>,
    pub mesh_job: Option<JobId>,
    pub view_thick: bool,
    pub color_cp: bool,
    pub edges: bool,
    pub show_3d: bool,
    pub stab_long: bool,
    pub stab_mode: usize,
    pub resp_time: f64,
    pub resp_var: usize,
    /// 0 = free modal response, 1 step, 2 pulse, 3 doublet, 4 sine
    pub resp_kind: usize,
    pub resp_amp: f64,
    pub resp_arg: f64,
}

impl Default for PlaneUi {
    fn default() -> Self {
        PlaneUi {
            sel_wing: 0,
            design_tab: 0,
            def: PlanePolarDef::default(),
            altitude: 0.0,
            water: false,
            auto_polars: 6,
            visible: HashSet::new(),
            sel_polar: None,
            sel_point: 0,
            plots: [(0, 4), (5, 4), (0, 9), (0, 13)],
            span_var: 0,
            cam: Camera::default(),
            mesh: None,
            mesh_job: None,
            view_thick: false,
            color_cp: true,
            edges: true,
            show_3d: true,
            stab_long: true,
            stab_mode: 0,
            resp_time: 10.0,
            resp_var: 0,
            resp_kind: 0,
            resp_amp: 0.2,
            resp_arg: 0.5,
        }
    }
}

pub const VARS3: [&str; 20] = ["alpha (deg)", "beta (deg)", "V (m/s)", "control", "CL", "CD", "CDi", "CDv", "CY", "Cm", "Cl (roll)", "Cn (yaw)", "CL/CD", "CL^1.5/CD", "Vz sink (m/s)", "efficiency e", "lift (N)", "drag (N)", "root bending (N.m)", "Xcp (m)"];

/// Value of plot variable `i`; NaN when the point does not have it (LLT
/// lateral coefficients), and such points are left out of the plots.
fn var3(o: &PlaneOpp, i: usize) -> f64 {
    match i {
        0 => o.alpha,
        1 => o.beta,
        2 => o.velocity,
        3 => o.ctrl,
        4 => o.cl,
        5 => o.cd,
        6 => o.cdi,
        7 => o.cdv,
        8 => o.cy.unwrap_or(f64::NAN),
        9 => o.cm,
        10 => o.cr.unwrap_or(f64::NAN),
        11 => o.cn.unwrap_or(f64::NAN),
        12 => o.ld(),
        13 => {
            if o.cd > 0.0 {
                o.cl.max(0.0).powf(1.5) / o.cd
            } else {
                0.0
            }
        }
        14 => o.vz(),
        15 => o.efficiency,
        16 => o.lift,
        17 => o.drag,
        18 => o.root_bending,
        _ => o.x_cp,
    }
}

const SPANVARS: [&str; 11] = ["local Cl", "Cl c / MAC", "induced Cd", "viscous Cd", "induced angle (deg)", "Reynolds", "transition top x/c", "transition bottom x/c", "lift per span (N/m)", "bending moment (N.m)", "effective alpha (deg)"];

fn plane_key(p: &Plane, o: &MeshOptions) -> u64 {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    serde_json::to_string(p).unwrap_or_default().hash(&mut h);
    o.thick.hash(&mut h);
    o.include_fuse.hash(&mut h);
    let mut f: Vec<_> = o.flaps.iter().map(|(k, v)| (*k, v.to_bits())).collect();
    f.sort();
    f.hash(&mut h);
    h.finish()
}

impl App {
    /// Ask for a mesh of `plane` (built in the background).
    fn ensure_mesh(&mut self, plane: &Plane, opts: MeshOptions) {
        let key = plane_key(plane, &opts);
        if self.pui.mesh.as_ref().map(|m| m.0) == Some(key) {
            return;
        }
        if let Some(id) = self.pui.mesh_job {
            if self.jobs.get(id).map(|j| !j.status().is_finished()).unwrap_or(false) {
                return;
            }
            self.pui.mesh_job = None;
        }
        let db = self.project.foil_db();
        let missing: Vec<String> = plane.foil_names().into_iter().filter(|n| !db.has_foil(n)).collect();
        let opts = MeshOptions { allow_missing_foils: true, ..opts };
        if !missing.is_empty() {
            self.toast = Some((format!("foils not loaded: {} (drawn as flat plates)", missing.join(", ")), std::time::Instant::now(), true));
        }
        let j = self.jobs.submit(format!("mesh {}", plane.name), JobKind::Mesh { plane: plane.clone(), db, options: opts });
        self.pui.mesh_job = Some(j.id);
        self.pending.push((j.id, Pending::Mesh(key)));
    }

    pub fn plane_design(&mut self, ui: &mut egui::Ui) {
        if self.project.planes.is_empty() {
            ui.label("No plane in the project. Use '+ New plane' in the explorer or import an xflr5/flow5 .xml file.");
            return;
        }
        self.sel_plane = self.sel_plane.min(self.project.planes.len() - 1);
        let foil_names: Vec<String> = self.project.foils.iter().map(|f| f.name.clone()).collect();
        let before = self.project.planes[self.sel_plane].clone();
        let mut edited = before.clone();
        let db = self.project.foil_db();
        egui::SidePanel::left("plane_left").resizable(true).default_width(560.0).show_inside(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                let p = &mut edited;
                ui.horizontal(|ui| {
                    ui.label("Plane");
                    ui.text_edit_singleline(&mut p.name);
                });
                ui.horizontal(|ui| {
                    ui.selectable_value(&mut self.pui.design_tab, 0, "Wings");
                    ui.selectable_value(&mut self.pui.design_tab, 1, "Fuselage");
                    ui.selectable_value(&mut self.pui.design_tab, 2, "Masses");
                    ui.selectable_value(&mut self.pui.design_tab, 3, "Summary");
                });
                ui.separator();
                match self.pui.design_tab {
                    0 => wings_editor(ui, p, &mut self.pui.sel_wing, &foil_names),
                    1 => fuse_editor(ui, p),
                    2 => masses_editor(ui, p),
                    _ => summary(ui, p, &db),
                }
            });
        });
        if edited != before {
            if let Err(e) = edited.validate() {
                self.toast = Some((e, std::time::Instant::now(), true));
            }
            // keep plane polars attached when renaming
            if edited.name != before.name {
                for q in &mut self.project.plane_polars {
                    if q.plane == before.name {
                        q.plane = edited.name.clone();
                    }
                }
            }
            self.snapshot();
            self.project.planes[self.sel_plane] = edited.clone();
        }
        egui::CentralPanel::default().show_inside(ui, |ui| {
            ui.horizontal(|ui| {
                ui.checkbox(&mut self.pui.view_thick, "thick surfaces (panel-method mesh)");
                ui.checkbox(&mut self.pui.edges, "panel edges");
                let plane_now = self.project.planes[self.sel_plane].clone();
                ui.menu_button("Export", |ui| {
                    let db = self.project.foil_db();
                    if ui.button("AVL geometry (.avl + foil files)").clicked() {
                        if let Some(path) = rfd::FileDialog::new().set_file_name(format!("{}.avl", plane_now.name)).save_file() {
                            let mut ok = std::fs::write(&path, foiley_plane::export::avl(&plane_now, &db, 0.0)).is_ok();
                            let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
                            for (n, t) in foiley_plane::export::avl_foil_files(&plane_now, &self.project.foils) {
                                ok &= std::fs::write(dir.join(n), t).is_ok();
                            }
                            self.notify(if ok { format!("wrote {}", path.display()) } else { "export failed".into() }, !ok);
                        }
                        ui.close();
                    }
                    if ui.button("AVL mass file (.mass)").clicked() {
                        if let Some(path) = rfd::FileDialog::new().set_file_name(format!("{}.mass", plane_now.name)).save_file() {
                            let _ = std::fs::write(&path, foiley_plane::export::avl_mass(&plane_now, &db));
                        }
                        ui.close();
                    }
                    if ui.button("Panel mesh (.stl)").clicked() {
                        if let Some(path) = rfd::FileDialog::new().set_file_name(format!("{}.stl", plane_now.name)).save_file() {
                            match foiley_plane::export::stl_mesh(&plane_now, &db, self.pui.view_thick) {
                                Ok(b) => {
                                    let _ = std::fs::write(&path, b);
                                }
                                Err(e) => self.notify(e.to_string(), true),
                            }
                        }
                        ui.close();
                    }
                    if ui.button("xflr5 plane XML").clicked() {
                        if let Some(path) = rfd::FileDialog::new().set_file_name(format!("{}.xml", plane_now.name)).save_file() {
                            let _ = std::fs::write(&path, foiley_plane::xml::write_plane(&plane_now));
                        }
                        ui.close();
                    }
                });
                hint(ui, "orange = control surfaces; the mesh is rebuilt in the background after each edit");
            });
            let plane = self.project.planes[self.sel_plane].clone();
            let opts = MeshOptions { thick: self.pui.view_thick, include_fuse: true, fuse_in_thin: true, ..Default::default() };
            self.ensure_mesh(&plane, opts.clone());
            let key = plane_key(&plane, &MeshOptions { allow_missing_foils: true, ..opts });
            let _ = key;
            let mesh = self.pui.mesh.as_ref().map(|m| &m.1);
            let inertia = plane.inertia(&self.project.foil_db());
            let mut ov = Overlay { lines: vec![], points: vec![] };
            if inertia.mass > 0.0 {
                ov.points.push((inertia.cog, Color32::YELLOW, format!("CoG {:.3} kg", inertia.mass)));
            }
            for pm in &plane.point_masses {
                ov.points.push((pm.position, Color32::from_rgb(255, 120, 200), pm.tag.clone()));
            }
            view3d::show(ui, &mut self.pui.cam, mesh, None, &ov, self.pui.edges);
        });
    }

    fn selected_plane_polar(&self) -> Option<PlanePolar> {
        let (pl, name) = self.pui.sel_polar.as_ref()?;
        self.project.plane_polars.iter().find(|p| &p.plane == pl && &p.name == name).cloned()
    }

    pub fn plane_analysis(&mut self, ui: &mut egui::Ui) {
        if self.project.planes.is_empty() {
            ui.label("No plane in the project.");
            return;
        }
        self.sel_plane = self.sel_plane.min(self.project.planes.len() - 1);
        let plane = self.project.planes[self.sel_plane].clone();
        egui::SidePanel::left("an3d_left").resizable(true).default_width(340.0).show_inside(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                heading(ui, &format!("Analysis of '{}'", plane.name));
                self.polar_def_editor(ui, &plane);
            });
        });
        egui::SidePanel::right("an3d_right").resizable(true).default_width(240.0).show_inside(ui, |ui| {
            heading(ui, "3D polars");
            egui::ScrollArea::vertical().show(ui, |ui| {
                let mut del = None;
                let mut vis_idx = 0;
                for (i, p) in self.project.plane_polars.iter().enumerate() {
                    let key = (p.plane.clone(), p.name.clone());
                    let swatch = if self.pui.visible.contains(&key) {
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
                        let mut on = self.pui.visible.contains(&key);
                        if ui.checkbox(&mut on, "").changed() {
                            if on {
                                self.pui.visible.insert(key.clone());
                            } else {
                                self.pui.visible.remove(&key);
                            }
                        }
                        let sel = self.pui.sel_polar.as_ref() == Some(&key);
                        if ui.selectable_label(sel, RichText::new(format!("{}\n{}", p.plane, p.name)).small()).on_hover_text(format!("{} points, {} failed, {} unknowns", p.points.len(), p.failed.len(), p.unknowns)).clicked() {
                            self.pui.sel_polar = Some(key.clone());
                            self.pui.sel_point = 0;
                            self.pui.visible.insert(key.clone());
                            if let Some(d) = &p.def {
                                self.pui.def = d.clone();
                            }
                        }
                        ui.menu_button("...", |ui| {
                            if ui.button("Export CSV").clicked() {
                                if let Some(path) = rfd::FileDialog::new().set_file_name(format!("{}_{}.csv", p.plane, p.name)).save_file() {
                                    let _ = std::fs::write(path, p.to_csv());
                                }
                                ui.close();
                            }
                            if ui.button("Export JSON (all results)").clicked() {
                                if let Some(path) = rfd::FileDialog::new().set_file_name(format!("{}_{}.json", p.plane, p.name)).save_file() {
                                    let _ = std::fs::write(path, serde_json::to_string_pretty(p).unwrap_or_default());
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
                    self.project.plane_polars.remove(i);
                }
            });
        });
        egui::CentralPanel::default().show_inside(ui, |ui| self.results_3d(ui, &plane));
    }

    fn polar_def_editor(&mut self, ui: &mut egui::Ui, plane: &Plane) {
        let d = &mut self.pui.def;
        let mut k = PolarKind::ALL.iter().position(|x| *x == d.kind).unwrap_or(0);
        let kl: Vec<&str> = PolarKind::ALL.iter().map(|k| k.label()).collect();
        let k_on: Vec<bool> = PolarKind::ALL.iter().map(|k| d.method.supports(*k, Setting::PolarType)).collect();
        ui.label("Polar type");
        if widgets::pick_enabled(ui, "kind3", &mut k, &kl, &k_on) {
            d.kind = PolarKind::ALL[k];
            d.clear_unsupported();
            d.range = match d.kind {
                PolarKind::FixedAlpha => foiley_plane::analysis::Range { from: 8.0, to: 30.0, step: 2.0 },
                PolarKind::Sideslip => foiley_plane::analysis::Range { from: -10.0, to: 10.0, step: 2.0 },
                PolarKind::Control | PolarKind::Stability => foiley_plane::analysis::Range { from: 0.0, to: 1.0, step: 0.25 },
                _ => foiley_plane::analysis::Range { from: -4.0, to: 10.0, step: 1.0 },
            };
        }
        let mut m = Method::ALL.iter().position(|x| *x == d.method).unwrap_or(2);
        let ml: Vec<&str> = Method::ALL.iter().map(|m| m.label()).collect();
        let m_on: Vec<bool> = Method::ALL.iter().map(|m| m.supports(d.kind, Setting::PolarType)).collect();
        ui.label("Method");
        if widgets::pick_enabled(ui, "meth3", &mut m, &ml, &m_on) {
            d.method = Method::ALL[m];
            d.clear_unsupported();
        }
        let on = |s: Setting| d.method.supports(d.kind, s);
        let (on_beta, on_ground, on_fdrag, on_deriv, on_fwake, on_flaps) = (on(Setting::Sideslip), on(Setting::GroundEffect), on(Setting::FuseDrag), on(Setting::Derivatives), on(Setting::FreeWake), on(Setting::Flaps));
        hint(
            ui,
            match d.method {
                Method::Llt => "Fast and nonlinear (captures stall from the 2D polars); main wing only, no sweep/dihedral effects.",
                Method::Vlm1 | Method::Vlm2 => "Thin-surface vortex lattice: whole plane, fuselage ignored. Good default for performance and stability.",
                Method::Panel => "Thick surfaces with source/doublet panels; includes the fuselage. Slower; use fewer panels.",
                Method::TriPanel => "Like 3D panels, on triangles (as flow5). STL/STEP fuselages need a triangle method: this one or linear doublets.",
                Method::TriLinear => "Linear doublets on triangles with Galerkin testing and a linear trailing-edge wake. Three unknowns per triangle; use a coarser mesh initially.",
            },
        );
        egui::Grid::new("def3").num_columns(2).show(ui, |ui| {
            match d.kind {
                PolarKind::FixedSpeed | PolarKind::Sideslip | PolarKind::List | PolarKind::Control => {
                    num_u(ui, "speed", &mut d.velocity, 0.1, 0.1..=400.0, " m/s", "freestream speed");
                }
                _ => {}
            }
            if matches!(d.kind, PolarKind::FixedAlpha | PolarKind::Sideslip | PolarKind::Control) {
                num_u(ui, "alpha", &mut d.alpha, 0.1, -30.0..=30.0, " deg", "angle of attack");
            }
            if matches!(d.kind, PolarKind::FixedSpeed | PolarKind::FixedLift | PolarKind::Glide | PolarKind::FixedAlpha | PolarKind::Stability) {
                widgets::num_u_on(ui, on_beta, "sideslip", &mut d.beta, 0.1, -30.0..=30.0, " deg", "sideslip angle (+ = wind from the right)");
            }
            if d.kind != PolarKind::List {
                ui.label(RichText::new(format!("sweep: {}", d.kind.swept())).strong());
                ui.end_row();
                let sp = if d.kind == PolarKind::FixedAlpha { 0.5 } else { 0.1 };
                num(ui, "  from", &mut d.range.from, sp, -1e4..=1e4, "first value of the swept variable");
                num(ui, "  to", &mut d.range.to, sp, -1e4..=1e4, "last value");
                num(ui, "  step", &mut d.range.step, sp * 0.1, 1e-4..=1e4, "increment");
            }
        });
        if d.kind == PolarKind::List {
            ui.label("points (alpha, beta, V):");
            let mut del = None;
            for (i, pt) in d.list.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.add(egui::DragValue::new(&mut pt.0).speed(0.1).prefix("a "));
                    ui.add(egui::DragValue::new(&mut pt.1).speed(0.1).prefix("b "));
                    ui.add(egui::DragValue::new(&mut pt.2).speed(0.1).prefix("V "));
                    if ui.small_button("x").clicked() {
                        del = Some(i);
                    }
                });
            }
            if let Some(i) = del {
                d.list.remove(i);
            }
            if ui.button("+ point").clicked() {
                d.list.push((2.0, 0.0, d.velocity));
            }
        }
        if d.kind == PolarKind::Control {
            section(ui, "Control ranges (value = min + t (max - min))", true, |ui| {
                ranges_editor(ui, &mut d.ranges, plane);
            });
        }
        if d.kind == PolarKind::Stability {
            section(ui, "Control gains (deflection per unit of control parameter)", true, |ui| {
                hint(ui, "Each point is trimmed (Cm = 0, lift = weight); flaps deflect by gain x control.");
                gains_editor(ui, &mut d.gains, plane);
            });
        }
        section(ui, "Aerodynamics", true, |ui| {
            egui::Grid::new("aero3").num_columns(2).show(ui, |ui| {
                check(ui, "viscous drag from 2D polars", &mut d.viscous, "adds profile drag interpolated from the foils' polars at each strip's Re and CL");
                if d.viscous || d.method == Method::Llt {
                    int(ui, "auto polar mesh (Re count)", &mut self.pui.auto_polars, 0..=20, "missing 2D polars are computed automatically over the Re range of this analysis (0 = never)");
                    check(ui, "interpolate by CL (xflr5)", &mut d.viscous_from_cl, "unchecked: interpolate by the section's effective angle of attack");
                    check(ui, "strict 2D polar envelope", &mut d.strict_envelope, "fail a point when any strip's Re or CL falls outside its 2D polars (instead of clamping with a warning)");
                }
                widgets::check_on(ui, on_fdrag, "fuselage friction drag", &mut d.fuse_drag, "turbulent skin friction x form factor on the fuselage wetted area");
                if d.method.is_panel() {
                    check(ui, "include fuselage panels", &mut d.include_fuse, "model the fuselage surface in the panel method");
                    num(ui, "wake length (x MAC)", &mut d.wake_length, 1.0, 1.0..=200.0, "length of the flat wake behind the trailing edges");
                    int(ui, "wake panels", &mut d.wake_panels, 1..=50, "panels per wake column");
                }
                if d.method == Method::Llt {
                    int(ui, "LLT stations", &mut d.llt.stations, 8..=200, "spanwise stations");
                    int(ui, "LLT max iterations", &mut d.llt.max_iter, 10..=5000, "");
                }
                let mut ge = d.ground_height.is_some();
                if widgets::check_on(ui, on_ground, "ground effect", &mut ge, "mirror image below the plane") {
                    d.ground_height = if ge { Some(0.5) } else { None };
                }
                if let Some(h) = &mut d.ground_height {
                    num_u(ui, "  ground below origin", h, 0.01, 0.0..=100.0, " m", "height of the plane origin above the ground");
                }
                widgets::check_on(ui, on_deriv, "stability derivatives at every point", &mut d.derivatives, "computes derivatives and modes for each point (needs mass); always on for T7");
                check(ui, "keep panel Cp for 3D colour maps", &mut d.keep_cp, "uses memory for large meshes");
                if matches!(d.method, Method::Vlm1 | Method::Vlm2) && on_fwake {
                    int(ui, "free-wake iterations", &mut d.free_wake, 0..=usize::MAX,"0 = rigid wake along x (fast, linear). >0 relaxes the trailing vortices along the local flow at every point (slower; shown in the 3D view)");
                }
            });
        });
        section(ui, "Fluid and reference", false, |ui| {
            egui::Grid::new("fluid3").num_columns(2).show(ui, |ui| {
                let mut ch = num_u(ui, "altitude (ISA)", &mut self.pui.altitude, 10.0, -500.0..=20000.0, " m", "sets density and viscosity from the standard atmosphere");
                ch |= check(ui, "water (hydrofoil)", &mut self.pui.water, "fresh water at 15 C");
                if ch {
                    d.fluid = if self.pui.water {
                        Fluid::water()
                    } else {
                        let a = isa(self.pui.altitude, 0.0);
                        Fluid { rho: a.rho, mu: a.mu, a: a.a }
                    };
                }
                num_u(ui, "density", &mut d.fluid.rho, 0.001, 0.01..=2000.0, " kg/m3", "");
                let mut nu = d.fluid.nu() * 1e5;
                if num_u(ui, "kinematic viscosity", &mut nu, 0.01, 0.01..=100.0, " e-5 m2/s", "") {
                    d.fluid.mu = nu * 1e-5 * d.fluid.rho;
                }
                let mut rk = match d.reference {
                    RefDims::Projected => 0,
                    RefDims::Planform => 1,
                    RefDims::Custom { .. } => 2,
                };
                ui.label("reference");
                if widgets::pick(ui, "ref3", &mut rk, &["main wing, projected", "main wing, planform", "custom"]) {
                    let (s, b, c) = plane.reference(true);
                    d.reference = match rk {
                        0 => RefDims::Projected,
                        1 => RefDims::Planform,
                        _ => RefDims::Custom { area: s, span: b, chord: c },
                    };
                }
                ui.end_row();
                if let RefDims::Custom { area, span, chord } = &mut d.reference {
                    num_u(ui, "  area", area, 0.001, 1e-6..=1e4, " m2", "");
                    num_u(ui, "  span", span, 0.01, 1e-6..=1e4, " m", "");
                    num_u(ui, "  chord", chord, 0.001, 1e-6..=1e4, " m", "");
                }
            });
        });
        section(ui, "Mass override", false, |ui| {
            let mut ov = d.inertia.is_some();
            if ui.checkbox(&mut ov, "override the plane's mass properties").changed() {
                d.inertia = if ov { Some(plane.inertia(&self.project.foil_db())) } else { None };
            }
            if let Some(i) = &mut d.inertia {
                egui::Grid::new("inert").num_columns(2).show(ui, |ui| {
                    num_u(ui, "mass", &mut i.mass, 0.01, 0.0..=1e6, " kg", "");
                    num_u(ui, "CoG x", &mut i.cog.x, 0.001, -1e3..=1e3, " m", "");
                    num_u(ui, "CoG z", &mut i.cog.z, 0.001, -1e3..=1e3, " m", "");
                    num_u(ui, "Ixx", &mut i.ixx, 0.001, 0.0..=1e6, " kg m2", "");
                    num_u(ui, "Iyy", &mut i.iyy, 0.001, 0.0..=1e6, " kg m2", "");
                    num_u(ui, "Izz", &mut i.izz, 0.001, 0.0..=1e6, " kg m2", "");
                    num_u(ui, "Ixz", &mut i.ixz, 0.001, -1e6..=1e6, " kg m2", "");
                });
            }
        });
        section(ui, "Flap presets and extra drag", false, |ui| {
            if !on_flaps {
                hint(ui, "This method ignores flap deflections.");
            }
            hint(ui, "Preset deflections (all polar types). Flaps exist where both section foils are marked as control surfaces (Foil design > Mark as control surface).");
            let flaps = flap_targets(plane);
            let mut del = None;
            for (i, ((w, f), deg)) in d.flaps.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.label(format!("{} flap {}", plane.wings.get(*w).map(|w| w.name.as_str()).unwrap_or("?"), *f + 1));
                    ui.add(egui::DragValue::new(deg).speed(0.2).suffix(" deg"));
                    if ui.small_button("x").clicked() {
                        del = Some(i);
                    }
                });
            }
            if let Some(i) = del {
                d.flaps.remove(i);
            }
            ui.add_enabled_ui(on_flaps, |ui| ui.menu_button("+ flap preset", |ui| {
                for (w, f, name) in &flaps {
                    if ui.button(name).clicked() {
                        d.flaps.push(((*w, *f), 5.0));
                        ui.close();
                    }
                }
                if flaps.is_empty() {
                    ui.label("no control surfaces on this plane");
                }
            }));
            ui.separator();
            let mut del = None;
            for (i, e) in d.extra_drag.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.add(egui::TextEdit::singleline(&mut e.name).desired_width(80.0));
                    ui.add(egui::DragValue::new(&mut e.area).speed(0.001).prefix("A ").suffix(" m2"));
                    ui.add(egui::DragValue::new(&mut e.cd).speed(0.001).prefix("CD "));
                    if ui.small_button("x").clicked() {
                        del = Some(i);
                    }
                });
            }
            if let Some(i) = del {
                d.extra_drag.remove(i);
            }
            if ui.button("+ extra drag").on_hover_text("e.g. landing gear: frontal area and drag coefficient").clicked() {
                d.extra_drag.push(ExtraDrag { name: "gear".into(), area: 0.001, cd: 1.0 });
            }
        });
        ui.separator();
        ui.horizontal(|ui| {
            ui.label("name");
            let mut auto_name = d.name.is_empty();
            let hint_name = d.auto_name();
            ui.add(egui::TextEdit::singleline(&mut d.name).hint_text(hint_name).desired_width(150.0));
            if ui.checkbox(&mut auto_name, "auto").changed() && auto_name {
                d.name.clear();
            }
        });
        if ui.add(egui::Button::new(RichText::new("Run analysis").strong()).min_size(egui::vec2(160.0, 30.0))).clicked() {
            self.run_3d(plane);
        }
    }

    pub fn run_3d_public(&mut self, plane: &Plane) {
        self.run_3d(plane)
    }

    fn run_3d(&mut self, plane: &Plane) {
        let mut def = self.pui.def.clone();
        if def.name.is_empty() {
            def.name = def.auto_name();
        }
        if let Err(e) = def.validate() {
            self.notify(e, true);
            return;
        }
        let db = self.project.foil_db();
        let missing: Vec<String> = plane.foil_names().into_iter().filter(|n| !db.has_foil(n)).collect();
        if !missing.is_empty() {
            self.notify(format!("load these foils first: {}", missing.join(", ")), true);
            return;
        }
        let title = format!("{} - {}", plane.name, def.name);
        if (def.viscous || def.method == Method::Llt) && self.pui.auto_polars > 0 {
            match auto::plan(plane, &self.project.foils, &db, &def, self.pui.auto_polars) {
                Ok(plan) if !plan.is_empty() => {
                    self.notify(format!("computing {} first", plan.describe()), false);
                    self.submit(format!("2D polars for {}", plane.name), JobKind::AutoPolars { batches: plan.batches }, Pending::AutoPolarsThen(Box::new((plane.clone(), def))));
                    return;
                }
                Ok(_) => {}
                Err(e) => {
                    self.notify(e, true);
                    return;
                }
            }
        }
        self.submit(title, JobKind::PlanePolar { plane: plane.clone(), db, def }, Pending::PlanePolar);
    }

    fn results_3d(&mut self, ui: &mut egui::Ui, plane: &Plane) {
        // polars to plot, including live partial results of running jobs
        let mut series: Vec<(String, Vec<PlaneOpp>, bool)> = self.project.plane_polars.iter().filter(|p| self.pui.visible.contains(&(p.plane.clone(), p.name.clone()))).map(|p| (format!("{} {}", p.plane, p.name), p.points.clone(), false)).collect();
        for j in self.jobs.jobs() {
            if j.kind == "3D polar" && !j.status().is_finished() {
                let pts = j.partial.lock().unwrap().plane_points.clone();
                if !pts.is_empty() {
                    let mut pts = pts;
                    pts.sort_by(|a, b| (a.alpha + a.velocity + a.ctrl + a.beta).total_cmp(&(b.alpha + b.velocity + b.ctrl + b.beta)));
                    series.push((format!("{} (running)", j.title), pts, true));
                }
            }
        }
        let sel = self.selected_plane_polar();
        let w = ui.available_width();
        let h = ui.available_height();
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.pui.show_3d, "3D view");
            if let Some(p) = &sel {
                if !p.points.is_empty() {
                    let n = p.points.len();
                    self.pui.sel_point = self.pui.sel_point.min(n - 1);
                    ui.label("point:");
                    ui.add(egui::Slider::new(&mut self.pui.sel_point, 0..=n - 1).show_value(false));
                    let o = &p.points[self.pui.sel_point];
                    ui.label(RichText::new(format!("alpha {:.2}  V {:.2}  CL {:.4}  CD {:.5}  Cm {:.4}  L/D {:.1}", o.alpha, o.velocity, o.cl, o.cd, o.cm, o.ld())).strong());
                    for wn in o.warnings.iter().take(2) {
                        ui.label(RichText::new("!").color(Color32::from_rgb(230, 160, 60))).on_hover_text(wn);
                    }
                    if let Some(e) = &o.stability_error {
                        ui.label(RichText::new("no stability").color(Color32::from_rgb(230, 160, 60))).on_hover_text(e);
                    }
                }
            }
        });
        let top_h = if self.pui.show_3d { h * 0.5 } else { h * 0.62 };
        let mut clicked: Option<usize> = None;
        egui::Grid::new("plots3").num_columns(2).show(ui, |ui| {
            for k in 0..4 {
                ui.vertical(|ui| {
                    ui.set_width(w / 2.0 - 8.0);
                    ui.horizontal(|ui| {
                        let (mut x, mut y) = self.pui.plots[k];
                        widgets::pick(ui, &format!("p3y{k}"), &mut y, &VARS3);
                        ui.label("vs");
                        widgets::pick(ui, &format!("p3x{k}"), &mut x, &VARS3);
                        self.pui.plots[k] = (x, y);
                    });
                    let (xv, yv) = self.pui.plots[k];
                    let mut plot = Plot::new(format!("p3{k}")).height(top_h / 2.0 - 34.0);
                    if series.len() <= 5 {
                        plot = plot.legend(Legend::default().position(egui_plot::Corner::RightBottom));
                    }
                    plot.x_axis_label(VARS3[xv]).y_axis_label(VARS3[yv]).show(ui, |pu| {
                        for (i, (name, pts, running)) in series.iter().enumerate() {
                            let xy: Vec<[f64; 2]> = pts.iter().map(|o| [var3(o, xv), var3(o, yv)]).filter(|p| p[0].is_finite() && p[1].is_finite()).collect();
                            let col = series_color(i);
                            let mut l = Line::new(name.clone(), PlotPoints::from(xy.clone())).color(col);
                            if *running {
                                l = l.style(egui_plot::LineStyle::dashed_dense());
                            }
                            pu.line(l);
                            pu.points(Points::new(name.clone(), PlotPoints::from(xy)).radius(2.5).color(col));
                        }
                        if let Some(p) = &sel {
                            if let Some(o) = p.points.get(self.pui.sel_point) {
                                pu.points(Points::new("selected", PlotPoints::from(vec![[var3(o, xv), var3(o, yv)]])).radius(6.0).shape(MarkerShape::Circle).filled(false).color(Color32::YELLOW));
                            }
                            if pu.response().clicked() {
                                if let Some(c) = pu.pointer_coordinate() {
                                    let b = pu.plot_bounds();
                                    let best = p.points.iter().enumerate().map(|(i, o)| (i, ((var3(o, xv) - c.x) / b.width()).powi(2) + ((var3(o, yv) - c.y) / b.height()).powi(2))).filter(|d| d.1.is_finite()).min_by(|a, b| a.1.total_cmp(&b.1));
                                    if let Some((i, d)) = best {
                                        if d < 0.003 {
                                            clicked = Some(i);
                                        }
                                    }
                                }
                            }
                        }
                    });
                });
                if k % 2 == 1 {
                    ui.end_row();
                }
            }
        });
        if let Some(i) = clicked {
            self.pui.sel_point = i;
        }
        ui.separator();
        let rest = ui.available_height();
        ui.horizontal(|ui| {
            ui.vertical(|ui| {
                ui.set_width(if self.pui.show_3d { w * 0.45 } else { w });
                ui.horizontal(|ui| {
                    ui.label("span distribution:");
                    widgets::pick(ui, "spanvar", &mut self.pui.span_var, &SPANVARS);
                });
                Plot::new("span").height(rest - 30.0).legend(Legend::default()).x_axis_label("span position (m)").show(ui, |pu| {
                    if let Some(p) = &sel {
                        if let Some(o) = p.points.get(self.pui.sel_point) {
                            let mac = p.ref_dims.2.max(1e-9);
                            let nw = o.strips.iter().map(|s| s.wing).max().map(|m| m + 1).unwrap_or(0);
                            for wi in 0..nw {
                                let pts: Vec<[f64; 2]> = o
                                    .strips
                                    .iter()
                                    .filter(|s| s.wing == wi)
                                    .map(|s| {
                                        let v = match self.pui.span_var {
                                            0 => s.cl,
                                            1 => s.cl * s.chord / mac,
                                            2 => s.cdi,
                                            3 => s.cdv,
                                            4 => s.ai,
                                            5 => s.re,
                                            6 => s.xtr_top,
                                            7 => s.xtr_bot,
                                            8 => s.lift_per_span,
                                            9 => s.bending,
                                            _ => s.alpha_eff,
                                        };
                                        [s.span_pos, v]
                                    })
                                    .collect();
                                let name = plane.wings.get(wi).map(|w| w.name.clone()).unwrap_or_else(|| format!("wing {wi}"));
                                pu.line(Line::new(name.clone(), PlotPoints::from(pts.clone())).color(series_color(wi)));
                                pu.points(Points::new(name, PlotPoints::from(pts)).radius(1.5).color(series_color(wi)));
                            }
                        }
                    }
                });
            });
            if self.pui.show_3d {
                ui.vertical(|ui| {
                    ui.horizontal(|ui| {
                        ui.checkbox(&mut self.pui.color_cp, "colour by Cp");
                        ui.checkbox(&mut self.pui.edges, "edges");
                    });
                    let def = sel.as_ref().and_then(|p| p.def.clone()).unwrap_or_else(|| self.pui.def.clone());
                    let opts = foiley_plane::analysis::mesh_options(&def, &Default::default(), &Default::default());
                    self.ensure_mesh(plane, opts);
                    let mesh = self.pui.mesh.as_ref().map(|m| &m.1);
                    let mut ov = Overlay { lines: vec![], points: vec![] };
                    let mut cpv: Vec<f32> = vec![];
                    if let (Some(p), Some(m)) = (&sel, mesh) {
                        if let Some(o) = p.points.get(self.pui.sel_point) {
                            ov.points.push((o.cog, Color32::YELLOW, "CoG".into()));
                            for l in &o.wake {
                                for w in l.windows(2) {
                                    let v = |a: [f32; 3]| Vec3::new(a[0] as f64, a[1] as f64, a[2] as f64);
                                    ov.lines.push((v(w[0]), v(w[1]), Color32::from_rgb(170, 200, 255)));
                                }
                            }
                            if o.cp.len() == m.panels.len() {
                                cpv = o.cp.clone();
                            }
                            // lift vectors per strip
                            let lmax = o.strips.iter().map(|s| s.lift_per_span.abs()).fold(0.0, f64::max).max(1e-9);
                            let span = plane.main_wing().map(|w| w.planform_span()).unwrap_or(1.0);
                            for (s, r) in m.strips.iter().zip(&o.strips) {
                                let len = 0.15 * span * r.lift_per_span / lmax;
                                ov.lines.push((s.c4, s.c4 + s.normal * len, Color32::from_rgb(120, 230, 120)));
                            }
                        }
                    }
                    let field = if self.pui.color_cp && !cpv.is_empty() {
                        let mut v: Vec<f32> = cpv.iter().cloned().filter(|x| x.is_finite()).collect();
                        v.sort_by(|a, b| a.total_cmp(b));
                        let lo = v.get(v.len() / 50).cloned().unwrap_or(-1.0);
                        let hi = v.get(v.len() - 1 - v.len() / 50).cloned().unwrap_or(1.0);
                        Some(ColorField { values: &cpv, min: lo, max: hi, label: if def.method.is_panel() { "Cp".into() } else { "dCp (loading)".into() } })
                    } else {
                        None
                    };
                    view3d::show(ui, &mut self.pui.cam, mesh, field.as_ref(), &ov, self.pui.edges);
                });
            }
        });
    }

    pub fn stability(&mut self, ui: &mut egui::Ui) {
        let stab_polars: Vec<PlanePolar> = self.project.plane_polars.iter().filter(|p| p.points.iter().any(|o| o.stability.is_some() || o.stability_error.is_some())).cloned().collect();
        if stab_polars.is_empty() {
            ui.vertical_centered(|ui| {
                ui.add_space(40.0);
                ui.label(RichText::new("No stability results yet").size(18.0));
                hint(ui, "Run a 'T7 - stability' polar in 3D analysis (or enable 'stability derivatives at every point'). The plane needs a mass: set wing/fuselage masses or point masses in Plane design > Masses.");
                if ui.button("Set up a T7 stability analysis").clicked() {
                    self.pui.def.kind = PolarKind::Stability;
                    self.pui.def.range = foiley_plane::analysis::Range { from: 0.0, to: 0.0, step: 1.0 };
                    self.module = crate::Module::PlaneAnalysis;
                }
            });
            return;
        }
        let key = self.pui.sel_polar.clone().filter(|k| stab_polars.iter().any(|p| (p.plane.clone(), p.name.clone()) == *k)).unwrap_or_else(|| (stab_polars[0].plane.clone(), stab_polars[0].name.clone()));
        let pol = stab_polars.iter().find(|p| (p.plane.clone(), p.name.clone()) == key).unwrap().clone();
        let pts: Vec<&PlaneOpp> = pol.points.iter().filter(|o| o.stability.is_some()).collect();
        let missing: Vec<&PlaneOpp> = pol.points.iter().filter(|o| o.stability.is_none() && o.stability_error.is_some()).collect();
        let missing_text = missing.iter().map(|o| format!("alpha {:.2}  V {:.2}  ctrl {:.3}: {}", o.alpha, o.velocity, o.ctrl, o.stability_error.as_deref().unwrap_or(""))).collect::<Vec<_>>().join("\n");
        if pts.is_empty() {
            ui.vertical_centered(|ui| {
                ui.add_space(40.0);
                ui.label(RichText::new(format!("Stability failed at every point of {} / {}", key.0, key.1)).size(18.0).color(Color32::from_rgb(230, 160, 60)));
                ui.label(missing_text);
            });
            return;
        }
        self.pui.sel_point = self.pui.sel_point.min(pts.len() - 1);
        egui::SidePanel::left("stab_left").resizable(true).default_width(430.0).show_inside(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label("polar");
                    egui::ComboBox::from_id_salt("stabpol").selected_text(format!("{} / {}", key.0, key.1)).show_ui(ui, |ui| {
                        for p in &stab_polars {
                            let k = (p.plane.clone(), p.name.clone());
                            if ui.selectable_label(k == key, format!("{} / {}", p.plane, p.name)).clicked() {
                                self.pui.sel_polar = Some(k);
                                self.pui.sel_point = 0;
                            }
                        }
                    });
                });
                if pts.len() > 1 {
                    ui.add(egui::Slider::new(&mut self.pui.sel_point, 0..=pts.len() - 1).text("point"));
                }
                if !missing.is_empty() {
                    ui.label(RichText::new(format!("{} point(s) without stability", missing.len())).color(Color32::from_rgb(230, 160, 60))).on_hover_text(&missing_text);
                }
                let o = pts[self.pui.sel_point];
                let s = o.stability.as_ref().unwrap();
                heading(ui, "Trim");
                egui::Grid::new("trim").num_columns(2).striped(true).show(ui, |ui| {
                    let mut row = |a: &str, b: String, h: &str| {
                        ui.label(a).on_hover_text(h);
                        ui.label(b);
                        ui.end_row();
                    };
                    row("control", format!("{:.3}", o.ctrl), "");
                    row("trim alpha", format!("{:.3} deg", s.alpha), "angle of attack with zero pitching moment");
                    row("trim speed", format!("{:.3} m/s", s.u0), "speed where lift equals weight");
                    row("CL / CD", format!("{:.4} / {:.5}", o.cl, o.cd), "");
                    row("neutral point", format!("{:.4} m", s.x_np), "CoG position giving neutral static stability");
                    row("static margin", format!("{:.1} % MAC", s.static_margin * 100.0), "(x_NP - x_CoG) / MAC; positive is stable, 5-15% typical");
                    row("CoG", format!("({:.4}, {:.4}, {:.4}) m", o.cog.x, o.cog.y, o.cog.z), "");
                    row("mass", format!("{:.3} kg", o.mass), "");
                });
                heading(ui, "Modes");
                egui::Grid::new("modes").num_columns(5).striped(true).show(ui, |ui| {
                    for h in ["mode", "eigenvalue", "freq (Hz)", "damping", "t1/2 / t2 (s)"] {
                        ui.label(RichText::new(h).strong());
                    }
                    ui.end_row();
                    for m in &s.modes {
                        ui.label(&m.name);
                        ui.label(format!("{:+.4} {:+.4}i", m.eigenvalue.re, m.eigenvalue.im));
                        ui.label(format!("{:.4}", m.freq_hz));
                        ui.label(RichText::new(format!("{:+.3}", m.damping)).color(if m.stable { Color32::from_rgb(80, 190, 90) } else { Color32::from_rgb(230, 80, 70) }));
                        ui.label(format!("{} {:.2}", if m.stable { "half" } else { "double" }, m.t_half_or_double));
                        ui.end_row();
                    }
                });
                heading(ui, "Non-dimensional derivatives (per rad)");
                let n = &s.nondim;
                egui::Grid::new("nd").num_columns(6).striped(true).show(ui, |ui| {
                    let items = [("CXu", n.cxu), ("CZu", n.czu), ("Cmu", n.cmu), ("CXa", n.cxa), ("CZa", n.cza), ("Cma", n.cma), ("CXq", n.cxq), ("CZq", n.czq), ("Cmq", n.cmq), ("CYb", n.cyb), ("CYp", n.cyp), ("CYr", n.cyr), ("Clb", n.clb), ("Clp", n.clp), ("Clr", n.clr), ("Cnb", n.cnb), ("Cnp", n.cnp), ("Cnr", n.cnr)];
                    for (i, (k, v)) in items.iter().enumerate() {
                        ui.label(RichText::new(*k).strong());
                        ui.label(format!("{v:+.4}"));
                        if i % 3 == 2 {
                            ui.end_row();
                        }
                    }
                });
                section(ui, "State matrices", false, |ui| {
                    for (name, a, vars) in [("longitudinal  [u, w, q, theta]", &s.a_long, ["u", "w", "q", "th"]), ("lateral  [v, p, r, phi]", &s.a_lat, ["v", "p", "r", "ph"])] {
                        ui.label(name);
                        egui::Grid::new(name).num_columns(5).show(ui, |ui| {
                            for (r, row) in a.iter().enumerate() {
                                ui.label(vars[r]);
                                for v in row {
                                    ui.monospace(format!("{v:+10.4}"));
                                }
                                ui.end_row();
                            }
                        });
                    }
                });
            });
        });
        egui::CentralPanel::default().show_inside(ui, |ui| {
            let h = ui.available_height();
            let o = pts[self.pui.sel_point];
            let s = o.stability.as_ref().unwrap();
            ui.label(RichText::new("Root locus (all points of the polar)").strong());
            Plot::new("rootlocus").height(h * 0.45).legend(Legend::default()).x_axis_label("real (1/s)").y_axis_label("imaginary (rad/s)").show(ui, |pu| {
                let mut by_mode: std::collections::BTreeMap<String, Vec<[f64; 2]>> = Default::default();
                for p in &pts {
                    for m in &p.stability.as_ref().unwrap().modes {
                        let e = by_mode.entry(m.name.clone()).or_default();
                        e.push([m.eigenvalue.re, m.eigenvalue.im]);
                        e.push([m.eigenvalue.re, -m.eigenvalue.im]);
                    }
                }
                for (i, (name, v)) in by_mode.into_iter().enumerate() {
                    pu.points(Points::new(name, PlotPoints::from(v)).radius(4.0).shape(MarkerShape::Cross).color(series_color(i)));
                }
                let cur: Vec<[f64; 2]> = s.modes.iter().flat_map(|m| [[m.eigenvalue.re, m.eigenvalue.im], [m.eigenvalue.re, -m.eigenvalue.im]]).collect();
                pu.points(Points::new("selected point", PlotPoints::from(cur)).radius(7.0).filled(false).color(Color32::YELLOW));
                pu.vline(egui_plot::VLine::new("stability limit", 0.0).color(Color32::GRAY));
            });
            ui.horizontal(|ui| {
                ui.label(RichText::new("Time response").strong());
                ui.selectable_value(&mut self.pui.stab_long, true, "longitudinal");
                ui.selectable_value(&mut self.pui.stab_long, false, "lateral");
                let modes: Vec<&str> = s.modes.iter().filter(|m| m.longitudinal == self.pui.stab_long).map(|m| m.name.as_str()).collect();
                ui.label("excite mode:");
                let mut idx = self.pui.stab_mode.min(modes.len().saturating_sub(1));
                widgets::pick(ui, "stabmode", &mut idx, &modes);
                self.pui.stab_mode = idx;
                ui.add(egui::DragValue::new(&mut self.pui.resp_time).range(0.5..=300.0).suffix(" s").prefix("duration "));
            });
            if s.control.is_some() {
                ui.horizontal(|ui| {
                    ui.label("control input:");
                    widgets::pick(ui, "respkind", &mut self.pui.resp_kind, &["none (free modal response)", "step", "pulse", "doublet", "sine"]);
                    if self.pui.resp_kind > 0 {
                        ui.add(egui::DragValue::new(&mut self.pui.resp_amp).speed(0.01).prefix("amplitude "));
                        if self.pui.resp_kind > 1 {
                            let lbl = if self.pui.resp_kind == 4 { "frequency (Hz) " } else { "duration (s) " };
                            ui.add(egui::DragValue::new(&mut self.pui.resp_arg).speed(0.01).range(0.01..=50.0).prefix(lbl));
                        }
                        hint(ui, "amplitude in units of the T7 control parameter (deflection = gain x amplitude)");
                    }
                });
            } else {
                hint(ui, "Forced responses need control gains in the T7 analysis (control derivatives).");
            }
            if self.pui.resp_kind > 0 && s.control.is_some() {
                use foiley_plane::stability::{forced_response, ControlInput};
                let input = match self.pui.resp_kind {
                    1 => ControlInput::Step,
                    2 => ControlInput::Pulse { duration: self.pui.resp_arg },
                    3 => ControlInput::Doublet { duration: self.pui.resp_arg },
                    _ => ControlInput::Sine { hz: self.pui.resp_arg },
                };
                let (a, b) = if self.pui.stab_long { (&s.a_long, &s.b_long) } else { (&s.a_lat, &s.b_lat) };
                let dt = (self.pui.resp_time / 2000.0).max(1e-4);
                let r = forced_response(a, b, input, self.pui.resp_amp, self.pui.resp_time, dt);
                let names = if self.pui.stab_long { ["u (m/s)", "w (m/s)", "q (deg/s)", "theta (deg)"] } else { ["v (m/s)", "p (deg/s)", "r (deg/s)", "phi (deg)"] };
                Plot::new("resp_forced").height(ui.available_height() - 4.0).legend(Legend::default()).x_axis_label("time (s)").show(ui, |pu| {
                    for k in 0..4 {
                        let sc = if k >= 2 || (!self.pui.stab_long && k == 1) { 180.0 / std::f64::consts::PI } else { 1.0 };
                        let pts: Vec<[f64; 2]> = r.iter().map(|(t, x)| [*t, x[k] * sc]).collect();
                        pu.line(Line::new(names[k], PlotPoints::from(pts)).color(series_color(k)));
                    }
                    let inp: Vec<[f64; 2]> = r.iter().map(|(t, _)| [*t, self.pui.resp_amp * input.value(*t)]).collect();
                    pu.line(Line::new("control input", PlotPoints::from(inp)).color(Color32::GRAY).style(egui_plot::LineStyle::dashed_dense()));
                });
                return;
            }
            let modes: Vec<_> = s.modes.iter().filter(|m| m.longitudinal == self.pui.stab_long).collect();
            if let Some(m) = modes.get(self.pui.stab_mode) {
                let a = if self.pui.stab_long { &s.a_long } else { &s.a_lat };
                let x0 = modal_initial_state(m, 1.0);
                let dt = (self.pui.resp_time / 2000.0).max(1e-4);
                let r = time_response(a, x0, self.pui.resp_time, dt);
                let names = if self.pui.stab_long { ["u (m/s)", "w (m/s)", "q (deg/s)", "theta (deg)"] } else { ["v (m/s)", "p (deg/s)", "r (deg/s)", "phi (deg)"] };
                Plot::new("resp").height(ui.available_height() - 4.0).legend(Legend::default()).x_axis_label("time (s)").show(ui, |pu| {
                    for k in 0..4 {
                        let sc = if k >= 2 || (!self.pui.stab_long && k == 1) { 180.0 / std::f64::consts::PI } else { 1.0 };
                        let pts: Vec<[f64; 2]> = r.iter().map(|(t, x)| [*t, x[k] * sc]).collect();
                        pu.line(Line::new(names[k], PlotPoints::from(pts)).color(series_color(k)));
                    }
                });
            }
        });
    }
}

fn flap_targets(p: &Plane) -> Vec<(usize, usize, String)> {
    // flaps: surfaces whose both foils carry a hinge; numbered left tip -> right tip per wing
    let mut out = vec![];
    for (wi, w) in p.wings.iter().enumerate() {
        let n = w.sections.len().saturating_sub(1);
        let sides = if w.two_sided { 2 } else { 1 };
        for k in 0..n * sides {
            out.push((wi, k, format!("{} flap {}", w.name, k + 1)));
        }
    }
    out
}

fn ranges_editor(ui: &mut egui::Ui, ranges: &mut Vec<ControlRange>, plane: &Plane) {
    let mut del = None;
    for (i, r) in ranges.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            ui.label(r.target.label(plane));
            ui.add(egui::DragValue::new(&mut r.min).speed(0.1).prefix("min "));
            ui.add(egui::DragValue::new(&mut r.max).speed(0.1).prefix("max "));
            if ui.small_button("x").clicked() {
                del = Some(i);
            }
        });
    }
    if let Some(i) = del {
        ranges.remove(i);
    }
    ui.menu_button("+ range", |ui| {
        let mut opts = vec![(ControlTarget::Velocity, 10.0, 20.0), (ControlTarget::Alpha, 0.0, 6.0), (ControlTarget::Beta, -5.0, 5.0), (ControlTarget::CogX, 0.1, 0.2)];
        for (w, f, _) in flap_targets(plane) {
            opts.push((ControlTarget::Flap { wing: w, index: f }, -10.0, 10.0));
        }
        for (wi, _) in plane.wings.iter().enumerate() {
            opts.push((ControlTarget::Tilt { wing: wi }, -3.0, 3.0));
        }
        for (t, a, b) in opts {
            if ui.button(t.label(plane)).clicked() {
                ranges.push(ControlRange { target: t, min: a, max: b });
                ui.close();
            }
        }
    });
}

fn gains_editor(ui: &mut egui::Ui, gains: &mut Vec<(ControlTarget, f64)>, plane: &Plane) {
    let mut del = None;
    for (i, (t, g)) in gains.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            ui.label(t.label(plane));
            ui.add(egui::DragValue::new(g).speed(0.1).suffix(" deg/unit"));
            if ui.small_button("x").clicked() {
                del = Some(i);
            }
        });
    }
    if let Some(i) = del {
        gains.remove(i);
    }
    ui.menu_button("+ gain", |ui| {
        for (w, f, name) in flap_targets(plane) {
            if ui.button(name).clicked() {
                gains.push((ControlTarget::Flap { wing: w, index: f }, 5.0));
                ui.close();
            }
        }
        for (wi, w) in plane.wings.iter().enumerate() {
            if ui.button(format!("{} tilt", w.name)).clicked() {
                gains.push((ControlTarget::Tilt { wing: wi }, 2.0));
                ui.close();
            }
        }
    });
}

fn spacing_pick(ui: &mut egui::Ui, id: String, s: &mut Spacing) {
    egui::ComboBox::from_id_salt(id).width(70.0).selected_text(format!("{s:?}")).show_ui(ui, |ui| {
        for v in Spacing::ALL {
            ui.selectable_value(s, v, v.label());
        }
    });
}

fn wings_editor(ui: &mut egui::Ui, p: &mut Plane, sel: &mut usize, foils: &[String]) {
    ui.horizontal_wrapped(|ui| {
        for (i, w) in p.wings.iter().enumerate() {
            ui.selectable_value(sel, i, &w.name).on_hover_text(w.kind.label());
        }
    });
    ui.horizontal(|ui| {
        let f0 = foils.first().cloned().unwrap_or_default();
        let f1 = foils.get(1).cloned().unwrap_or_else(|| f0.clone());
        ui.menu_button("+ wing", |ui| {
            if ui.button("main wing").clicked() {
                p.wings.push(Wing::default_main(&f0));
                ui.close();
            }
            if ui.button("elevator").clicked() {
                p.wings.push(Wing::default_elevator(&f1));
                ui.close();
            }
            if ui.button("fin").clicked() {
                p.wings.push(Wing::default_fin(&f1));
                ui.close();
            }
        });
        if !p.wings.is_empty() && ui.button("duplicate").clicked() {
            let mut w = p.wings[*sel].clone();
            w.name = format!("{} copy", w.name);
            p.wings.push(w);
        }
        if p.wings.len() > 1 && ui.button("delete").clicked() {
            p.wings.remove(*sel);
        }
    });
    *sel = (*sel).min(p.wings.len().saturating_sub(1));
    let Some(w) = p.wings.get_mut(*sel) else { return };
    egui::Grid::new("wingprops").num_columns(4).show(ui, |ui| {
        ui.label("name");
        ui.text_edit_singleline(&mut w.name);
        ui.label("type");
        egui::ComboBox::from_id_salt("wkind").selected_text(w.kind.label()).show_ui(ui, |ui| {
            for k in WingKind::ALL {
                ui.selectable_value(&mut w.kind, k, k.label());
            }
        });
        ui.end_row();
        ui.label("position x, y, z (m)").on_hover_text("root leading-edge reference point in plane axes (x aft, y right, z up)");
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut w.position.x).speed(0.005));
            ui.add(egui::DragValue::new(&mut w.position.y).speed(0.005));
            ui.add(egui::DragValue::new(&mut w.position.z).speed(0.005));
        });
        ui.label("tilt (deg)").on_hover_text("incidence about y; + = nose up");
        ui.add(egui::DragValue::new(&mut w.ry).speed(0.1));
        ui.end_row();
        ui.label("rotation about x (deg)").on_hover_text("90 turns a one-sided wing into a vertical fin");
        ui.add(egui::DragValue::new(&mut w.rx).speed(1.0));
        ui.label("structural mass (kg)").on_hover_text("distributed over the wing volume for inertia");
        ui.add(egui::DragValue::new(&mut w.mass).speed(0.01).range(0.0..=1e6));
        ui.end_row();
        ui.checkbox(&mut w.two_sided, "two-sided").on_hover_text("build both halves (unchecked: single surface, e.g. fin)");
        ui.checkbox(&mut w.symmetric, "symmetric foils").on_hover_text("left foils follow the right foils");
        ui.checkbox(&mut w.close_inner, "close root").on_hover_text("close the root side (panel method)");
        ui.end_row();
    });
    ui.label(RichText::new(format!("span {:.3} m | area {:.4} m2 | MAC {:.4} m | AR {:.2} | taper {:.2} | sweep {:.1} deg", w.planform_span(), w.planform_area(), w.mac(), w.aspect_ratio(), w.taper(), w.sweep())).small());
    ui.separator();
    ui.label(RichText::new("Sections (root to tip)").strong());
    hint(ui, "dihedral, panel counts and distributions apply to the panel outboard of a section");
    let mut action: Option<(usize, u8)> = None;
    egui::ScrollArea::horizontal().show(ui, |ui| {
        egui::Grid::new("sections").striped(true).num_columns(12).show(ui, |ui| {
            for h in ["", "y (m)", "chord (m)", "offset (m)", "dihedral", "twist", "foil (right)", "foil (left)", "nx", "x-dist", "ny", "y-dist"] {
                ui.label(RichText::new(h).strong());
            }
            ui.end_row();
            let n = w.sections.len();
            let sym = w.symmetric;
            for (i, s) in w.sections.iter_mut().enumerate() {
                ui.menu_button(format!("{}", i + 1), |ui| {
                    if ui.button("insert after (split)").clicked() {
                        action = Some((i, 0));
                        ui.close();
                    }
                    if n > 2 && ui.button("delete").clicked() {
                        action = Some((i, 1));
                        ui.close();
                    }
                });
                ui.add(egui::DragValue::new(&mut s.y).speed(0.005).range(0.0..=1e4));
                ui.add(egui::DragValue::new(&mut s.chord).speed(0.002).range(0.001..=1e3));
                ui.add(egui::DragValue::new(&mut s.offset).speed(0.002));
                ui.add(egui::DragValue::new(&mut s.dihedral).speed(0.1).suffix("°"));
                ui.add(egui::DragValue::new(&mut s.twist).speed(0.1).suffix("°"));
                for (side, f) in [("r", &mut s.foil_right), ("l", &mut s.foil_left)] {
                    if side == "l" && sym {
                        ui.label("=");
                        continue;
                    }
                    egui::ComboBox::from_id_salt(format!("foil{side}{i}")).width(95.0).selected_text(f.clone()).show_ui(ui, |ui| {
                        for name in foils {
                            ui.selectable_value(f, name.clone(), name);
                        }
                    });
                }
                if sym {
                    s.foil_left = s.foil_right.clone();
                }
                ui.add(egui::DragValue::new(&mut s.nx).range(1..=100));
                spacing_pick(ui, format!("xd{i}"), &mut s.xdist);
                ui.add(egui::DragValue::new(&mut s.ny).range(1..=200));
                spacing_pick(ui, format!("yd{i}"), &mut s.ydist);
                ui.end_row();
            }
        });
    });
    match action {
        Some((i, 0)) => {
            if i + 1 < w.sections.len() {
                w.split_section(i);
            } else {
                let mut s = w.sections[i].clone();
                s.y += 0.2 * s.chord.max(0.05);
                w.sections.push(s);
            }
        }
        Some((i, 1)) => {
            w.sections.remove(i);
        }
        _ => {}
    }
    if ui.button("+ tip section").clicked() {
        let last = w.sections.last().cloned().unwrap_or_else(|| WingSection::new(0.0, 0.2, 0.0, ""));
        let mut s = last.clone();
        s.y += 0.25 * last.chord.max(0.1);
        s.chord *= 0.8;
        w.sections.push(s);
    }
}

fn fuse_editor(ui: &mut egui::Ui, p: &mut Plane) {
    if p.fuses.is_empty() {
        hint(ui, "No fuselage. VLM ignores fuselages; the 3D panel method models them.");
        if ui.button("+ add a default pod").clicked() {
            p.fuses.push(Fuse::default_pod());
        }
        return;
    }
    let mut remove = false;
    let f = &mut p.fuses[0];
    ui.label(RichText::new(format!("type: {}", f.kind_label())).strong());
    ui.horizontal(|ui| {
        static UNITS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let mut u = UNITS.load(std::sync::atomic::Ordering::Relaxed);
        if widgets::pick(ui, "fuseunits", &mut u, &["units: auto (STEP) / m", "mm", "cm", "inch", "foot"]) {
            UNITS.store(u, std::sync::atomic::Ordering::Relaxed);
        }
        if ui.button("Import STL / STEP...").on_hover_text("replace this fuselage by a CAD surface; use the triangle panel method to analyse it").clicked() {
            if let Some(path) = rfd::FileDialog::new().add_filter("CAD surface", &["stl", "step", "stp"]).pick_file() {
                let ext = path.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default();
                let res = if ext == "stl" {
                    std::fs::read(&path).map_err(|e| e.to_string()).and_then(|b| foiley_plane::cad::read_stl(&b))
                } else {
                    std::fs::read_to_string(&path).map_err(|e| e.to_string()).and_then(|t| foiley_plane::cad::read_step_target(&t, 3000).map(|r| r.0))
                };
                if let Ok(mut tris) = res {
                    let s = [1.0, 1e-3, 1e-2, 0.0254, 0.3048][u];
                    if s != 1.0 {
                        for t in &mut tris {
                            for v in t.iter_mut() {
                                *v = *v * s;
                            }
                        }
                    }
                    let mass = f.mass;
                    let name = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
                    *f = Fuse::from_triangles(&name, tris, mass);
                }
            }
        }
        if f.triangles.is_empty() && f.nurbs.is_none() && ui.button("Convert to NURBS").on_hover_text("edit control points like xflr5 / flow5").clicked() {
            f.nurbs = f.to_nurbs();
        }
        if f.nurbs.is_some() && ui.button("Back to elliptic frames").clicked() {
            f.nurbs = None;
        }
    });
    if !f.triangles.is_empty() {
        ui.label(format!("{} triangles, length {:.3} m, wetted area {:.4} m2", f.triangles.len(), f.length(), f.wetted_area()));
        egui::Grid::new("fusetri").num_columns(2).show(ui, |ui| {
            ui.label("position x, y, z (m)");
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut f.position.x).speed(0.005));
                ui.add(egui::DragValue::new(&mut f.position.y).speed(0.005));
                ui.add(egui::DragValue::new(&mut f.position.z).speed(0.005));
            });
            ui.end_row();
            num_u(ui, "structural mass", &mut f.mass, 0.01, 0.0..=1e6, " kg", "distributed over the skin");
        });
        if ui.button("remove fuselage").clicked() {
            p.fuses.clear();
        }
        return;
    }
    if let Some(n) = &mut f.nurbs {
        egui::Grid::new("nurbsdeg").num_columns(2).show(ui, |ui| {
            int(ui, "degree along x", &mut n.degree_u, 1..=5, "");
            int(ui, "degree around", &mut n.degree_v, 1..=5, "");
            int(ui, "panels along x", &mut f.nx, 2..=200, "");
            int(ui, "panels around (per side)", &mut f.nh, 2..=60, "");
            num_u(ui, "structural mass", &mut f.mass, 0.01, 0.0..=1e6, " kg", "distributed over the skin");
        });
        ui.label(RichText::new("NURBS frames (right half, points from top to bottom)").strong());
        let np = n.frames.first().map(|fr| fr.points.len()).unwrap_or(0);
        let mut del = None;
        for (i, fr) in n.frames.iter_mut().enumerate() {
            egui::CollapsingHeader::new(format!("frame {} at x = {:.3}", i + 1, fr.x)).id_salt(format!("nf{i}")).show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label("x");
                    if ui.add(egui::DragValue::new(&mut fr.x).speed(0.002)).changed() {
                        let x = fr.x;
                        fr.points.iter_mut().for_each(|pt| pt.x = x);
                    }
                    if ui.small_button("delete frame").clicked() {
                        del = Some(i);
                    }
                });
                for (k, pt) in fr.points.iter_mut().enumerate() {
                    ui.horizontal(|ui| {
                        ui.label(format!("{}", k + 1));
                        ui.add(egui::DragValue::new(&mut pt.y).speed(0.001).prefix("y "));
                        ui.add(egui::DragValue::new(&mut pt.z).speed(0.001).prefix("z "));
                    });
                }
            });
        }
        if let Some(i) = del {
            if n.frames.len() > 2 {
                n.frames.remove(i);
            }
        }
        if ui.button("+ frame").clicked() {
            if let Some(last) = n.frames.last().cloned() {
                let mut fr = last;
                fr.x += 0.1;
                let x = fr.x;
                fr.points.iter_mut().for_each(|pt| pt.x = x);
                n.frames.push(fr);
            }
        }
        let _ = np;
        if ui.button("remove fuselage").clicked() {
            remove = true;
        }
        if remove {
            p.fuses.clear();
        }
        return;
    }
    egui::Grid::new("fuse").num_columns(2).show(ui, |ui| {
        ui.label("name");
        ui.text_edit_singleline(&mut f.name);
        ui.end_row();
        ui.label("position x, y, z (m)");
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut f.position.x).speed(0.005));
            ui.add(egui::DragValue::new(&mut f.position.y).speed(0.005));
            ui.add(egui::DragValue::new(&mut f.position.z).speed(0.005));
        });
        ui.end_row();
        num(ui, "section exponent", &mut f.exponent, 0.1, 1.0..=10.0, "2 = elliptic sections, larger = more rectangular");
        int(ui, "panels along x", &mut f.nx, 2..=200, "");
        int(ui, "panels around (per side)", &mut f.nh, 2..=60, "");
        num_u(ui, "structural mass", &mut f.mass, 0.01, 0.0..=1e6, " kg", "distributed over the skin");
    });
    ui.label(RichText::new("Frames").strong());
    let mut del = None;
    egui::Grid::new("frames").striped(true).num_columns(6).show(ui, |ui| {
        for h in ["", "x (m)", "z centre (m)", "half width (m)", "half height (m)", ""] {
            ui.label(RichText::new(h).strong());
        }
        ui.end_row();
        for (i, fr) in f.frames.iter_mut().enumerate() {
            ui.label(format!("{}", i + 1));
            ui.add(egui::DragValue::new(&mut fr.x).speed(0.002));
            ui.add(egui::DragValue::new(&mut fr.z).speed(0.001));
            ui.add(egui::DragValue::new(&mut fr.half_width).speed(0.001).range(0.0..=100.0));
            ui.add(egui::DragValue::new(&mut fr.half_height).speed(0.001).range(0.0..=100.0));
            if ui.small_button("x").clicked() {
                del = Some(i);
            }
            ui.end_row();
        }
    });
    if let Some(i) = del {
        if f.frames.len() > 2 {
            f.frames.remove(i);
        }
    }
    ui.horizontal(|ui| {
        if ui.button("+ frame").clicked() {
            let last = *f.frames.last().unwrap();
            f.frames.push(FuseFrame { x: last.x + 0.1, ..last });
        }
        if ui.button("remove fuselage").clicked() {
            remove = true;
        }
    });
    if f.validate().is_err() {
        ui.label(RichText::new("frame x positions must increase").color(Color32::from_rgb(230, 80, 70)));
    }
    if remove {
        p.fuses.clear();
    }
}

fn masses_editor(ui: &mut egui::Ui, p: &mut Plane) {
    hint(ui, "Structural masses are set on each wing / fuselage. Add batteries, servos, ballast... here (plane axes).");
    let mut del = None;
    egui::Grid::new("pm").striped(true).num_columns(6).show(ui, |ui| {
        for h in ["tag", "mass (kg)", "x (m)", "y (m)", "z (m)", ""] {
            ui.label(RichText::new(h).strong());
        }
        ui.end_row();
        for (i, m) in p.point_masses.iter_mut().enumerate() {
            ui.add(egui::TextEdit::singleline(&mut m.tag).desired_width(110.0));
            ui.add(egui::DragValue::new(&mut m.mass).speed(0.005).range(0.0..=1e6));
            ui.add(egui::DragValue::new(&mut m.position.x).speed(0.002));
            ui.add(egui::DragValue::new(&mut m.position.y).speed(0.002));
            ui.add(egui::DragValue::new(&mut m.position.z).speed(0.002));
            if ui.small_button("x").clicked() {
                del = Some(i);
            }
            ui.end_row();
        }
    });
    if let Some(i) = del {
        p.point_masses.remove(i);
    }
    if ui.button("+ point mass").clicked() {
        p.point_masses.push(PointMass { tag: "mass".into(), mass: 0.1, position: Vec3::ZERO });
    }
}

fn summary(ui: &mut egui::Ui, p: &Plane, db: &foiley_plane::FoilDb) {
    let i = p.inertia(db);
    egui::Grid::new("summ").striped(true).num_columns(2).show(ui, |ui| {
        let mut row = |a: &str, b: String| {
            ui.label(a);
            ui.label(b);
            ui.end_row();
        };
        if let Some(w) = p.main_wing() {
            row("reference wing", w.name.clone());
            row("span (planform / projected)", format!("{:.4} / {:.4} m", w.planform_span(), w.projected_span()));
            row("area (planform / projected)", format!("{:.5} / {:.5} m2", w.planform_area(), w.projected_area()));
            row("MAC", format!("{:.4} m", w.mac()));
            row("aspect ratio", format!("{:.3}", w.aspect_ratio()));
            let (wl, _) = (i.mass / w.planform_area().max(1e-9), 0);
            row("wing loading", format!("{:.2} kg/m2  ({:.1} g/dm2)", wl, wl * 10.0));
        }
        let (h, v) = p.tail_volumes();
        row("tail volume (horizontal / vertical)", format!("{h:.3} / {v:.4}"));
        row("total mass", format!("{:.4} kg", i.mass));
        row("CoG", format!("({:.4}, {:.4}, {:.4}) m", i.cog.x, i.cog.y, i.cog.z));
        row("Ixx / Iyy / Izz", format!("{:.5} / {:.5} / {:.5} kg m2", i.ixx, i.iyy, i.izz));
        row("Ixz", format!("{:.5} kg m2", i.ixz));
        row("foils used", p.foil_names().join(", "));
    });
}
