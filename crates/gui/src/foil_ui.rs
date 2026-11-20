//! Foil design and 2D analysis modules.

use crate::widgets::{self, check, heading, hint, int, int_keep, num, num_keep, num_u, section, series_color};
use crate::{App, Pending};
use eframe::egui::{self, Color32, RichText};
use egui_plot::{Legend, Line, Plot, PlotPoints, Points};
use foiley_app::JobKind;
use foiley_foil::analysis::{PointDetail, PolarType};
use foiley_foil::batch::{re_list_geometric, BatchCase};
use foiley_foil::foil::{node_hash, NodeInput, Paneling};
use foiley_foil::{AnalysisSettings, Foil, FoilPoint, FoilPolar, Sweep};
use foiley_foil::inverse::{modify, InverseMode, ModKind, Modification, SpeedDist};
use foiley_foil::gen::cst::{Cst, CstFit};
use foiley_foil::gen::kt::KarmanTrefftz;
use foiley_foil::gen::parsec::Parsec;
use foiley_foil::gen::plate::{Plate, PlateCamber, PlateEdge};
use foiley_foil::optim::{Algorithm, Candidate, DragGoal, Goal, OperatingPoint, OptimResult, OptimSettings, Parameterization, SparConstraint};
use foiley_foil::Recipe;
use foiley_app::JobId;
use std::collections::{BTreeSet, HashMap, HashSet};

struct PolarLabel {
    key: (String, String),
    name: String,
}

/// Compute identities once per frame and reuse them across selectors and all
/// four plots. Only settings collisions need a short hash in their label.
fn polar_labels(polars: &[FoilPolar]) -> Vec<PolarLabel> {
    let mut labels: Vec<_> = polars.iter().map(|p| PolarLabel { key: p.key(), name: p.display_name() }).collect();
    let mut counts = HashMap::new();
    for label in &labels { *counts.entry((label.key.0.clone(), label.name.clone())).or_insert(0) += 1; }
    for label in &mut labels {
        if counts[&(label.key.0.clone(), label.name.clone())] > 1 {
            label.name.push_str(&format!(" / {}", &label.key.1[..12]));
        }
    }
    labels
}

/// What an inverse-design result was computed from: the foil (name and
/// node hash) and the solver settings.  A speed distribution is shown and
/// solved only while all of them still hold.
#[derive(Clone, Debug, PartialEq)]
pub struct InvKey {
    pub foil: String,
    pub hash: u64,
    pub mode: InverseMode,
    pub alpha: f64,
    pub paneling: Paneling,
}

pub const GEN_FAMILIES: [&str; 5] = ["NACA", "CST (Kulfan)", "PARSEC", "Karman-Trefftz", "cambered plate"];

const NACA_HELP: &str = "4-digit 2412 \u{2022} 5-digit 23012, reflexed 23112, any design CL 43012 \u{2022} 4-digit modified 0012-64 \u{2022} 16-series 16-212 \u{2022} 6-series 63-206, 64(2)-415, 65-415 a=0.5 \u{2022} 6A-series 64A010";

/// The Create section's generator: one parameter set per family, and the
/// last generated preview (regenerated only when a parameter changes).
pub struct GenUi {
    pub family: usize,
    pub closed_te: bool,
    pub cst: Cst,
    pub fit_order: usize,
    pub fit: Option<Result<(String, CstFit), String>>,
    pub parsec: Parsec,
    pub kt: KarmanTrefftz,
    pub kt_alpha: f64,
    pub plate: Plate,
    pub show: bool,
    cache: Option<(Recipe, Result<Foil, String>)>,
    #[allow(clippy::type_complexity)]
    kt_exact: Option<((KarmanTrefftz, f64), Result<(f64, f64), String>)>,
}

impl Default for GenUi {
    fn default() -> Self {
        GenUi { family: 0, closed_te: false, cst: Cst::default(), fit_order: 8, fit: None, parsec: Parsec::default(), kt: KarmanTrefftz::default(), kt_alpha: 4.0, plate: Plate::default(), show: true, cache: None, kt_exact: None }
    }
}

impl GenUi {
    /// The recipe the editor currently describes (`naca` is the NACA text field).
    pub fn recipe(&self, naca: &str) -> Recipe {
        match self.family {
            1 => Recipe::Cst(self.cst.clone()),
            2 => Recipe::Parsec(self.parsec),
            3 => Recipe::KarmanTrefftz(self.kt),
            4 => Recipe::Plate(self.plate),
            _ => {
                let code = naca.trim();
                let closed = self.closed_te && !code.to_ascii_lowercase().contains("closed");
                Recipe::Naca { designation: if closed { format!("{code} closed TE") } else { code.to_string() } }
            }
        }
    }

    /// Put a foil's recipe into the editor, if it is valid (a saved project
    /// could hold, say, CST weight lists of different lengths).
    pub fn load(&mut self, r: &Recipe, naca: &mut String) -> Result<(), String> {
        r.validate().map_err(|e| format!("cannot edit this recipe: {e}"))?;
        match r {
            Recipe::Naca { designation } => {
                self.family = 0;
                let lower = designation.to_ascii_lowercase();
                self.closed_te = lower.ends_with("closed te");
                *naca = if self.closed_te { designation[..designation.len() - "closed TE".len()].trim().to_string() } else { designation.clone() };
            }
            Recipe::Cst(c) => {
                self.family = 1;
                self.cst = c.clone();
            }
            Recipe::Parsec(p) => {
                self.family = 2;
                self.parsec = *p;
            }
            Recipe::KarmanTrefftz(k) => {
                self.family = 3;
                self.kt = *k;
            }
            Recipe::Plate(p) => {
                self.family = 4;
                self.plate = *p;
            }
        }
        Ok(())
    }

    /// The generated foil for the current parameters (cached).
    pub fn preview(&mut self, naca: &str) -> Result<Foil, String> {
        let r = self.recipe(naca);
        if self.cache.as_ref().map(|c| &c.0) != Some(&r) {
            let name = match &r {
                Recipe::Naca { .. } => None,
                _ => Some(["NACA", "CST", "PARSEC", "Karman-Trefftz", "Plate"][self.family.min(4)]),
            };
            let out = Foil::from_recipe(&r, name).map_err(|e| e.to_string());
            self.cache = Some((r, out));
        }
        self.cache.as_ref().unwrap().1.clone()
    }

    /// The cached preview, if it is valid and shown.
    pub fn shown(&self) -> Option<&Foil> {
        self.cache.as_ref().filter(|_| self.show).and_then(|c| c.1.as_ref().ok())
    }
}

/// UI state for "Use points from 3D analysis" inside the optimiser.
#[derive(Clone, Debug, Default)]
pub struct From3dPicker {
    pub force_open: Option<bool>,
    pub polar: Option<(String, String)>,
    pub wing: usize,
    pub span_pos: Option<f64>,
    pub selected_points: BTreeSet<usize>,
}

/// Which 3D plane, polar and wing station an operating-point table came from.
#[derive(Clone, Debug, PartialEq)]
pub struct SectionOptimOrigin {
    pub plane: String,
    pub polar_name: String,
    pub polar_def: Option<foiley_plane::PlanePolarDef>,
    pub wing: usize,
    pub wing_name: String,
    pub span_pos: f64,
    pub foil: String,
    /// Resolved 2D analysis settings of `foil` in the source 3D polar (Mach,
    /// Ncrit, trips, wake controls) so the optimiser and comparison rerun keep
    /// identical transition settings.
    pub source_analysis: Option<AnalysisSettings>,
}

/// Frozen metadata captured when an optimiser job is submitted, keeping the
/// result, spar constraint, base foil name and 3D origin bound to that run.
#[derive(Clone, Debug)]
pub struct OptimRunMeta {
    pub base_foil: String,
    pub spar: Option<SparConstraint>,
    pub origin: Option<SectionOptimOrigin>,
    pub all_target_goals: bool,
}

impl OptimRunMeta {
    pub fn follow_foil_rename(&mut self, old: &str, new: &str) {
        if self.base_foil == old {
            self.base_foil = new.to_string();
        }
        if let Some(o) = &mut self.origin {
            if o.foil == old {
                o.foil = new.to_string();
            }
        }
    }

    pub fn follow_plane_rename(&mut self, old: &str, new: &str) {
        if let Some(o) = &mut self.origin {
            if o.plane == old {
                o.plane = new.to_string();
            }
        }
    }
}

type Apply3dRequest = ((String, String), usize, f64, BTreeSet<usize>);

/// Safe default spar minimum thickness (`0.5%` ..= `cur_t`) for any valid foil,
/// including thin sections below 4% thickness.
fn default_spar_min(cur_t: f64) -> f64 {
    if cur_t.is_finite() && cur_t > 0.005 {
        let stepped = ((cur_t * 200.0).floor() / 200.0) - 0.005;
        let floor = (cur_t * 0.85).max(0.005);
        stepped.max(floor).clamp(0.005, cur_t)
    } else {
        0.005
    }
}

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
    pub bl_arc_length: bool,
    pub cp_inviscid: bool,
    pub show_failed: bool,
    pub inv_mode: InverseMode,
    pub inv_alpha: f64,
    pub qdist: Option<SpeedDist>,
    pub inv_target: Option<foiley_app::InverseTarget>,
    pub inv_upper: bool,
    pub inv_x: [f64; 2],
    pub inv_kind: usize,
    pub inv_val: [f64; 2],
    pub inv_iter: usize,
    /// What the speed distribution was computed from.
    pub inv_src: Option<InvKey>,
    /// Name and node hash of the foil the preview was made from.
    pub preview_for: Option<(String, u64)>,
    pub opt: OptimSettings,
    pub opt_multipoint: bool,
    pub opt_points_saved: Vec<OperatingPoint>,
    pub opt_force_open: Option<bool>,
    pub opt_scroll_to: bool,
    pub opt_sel_pareto: Option<usize>,
    pub opt_save_name: String,
    pub opt_from_3d: From3dPicker,
    pub opt_3d_origin: Option<SectionOptimOrigin>,
    pub opt_job: Option<(JobId, OptimRunMeta)>,
    pub opt_result: Option<(OptimResult, OptimRunMeta)>,
    pub gen: GenUi,
}

impl FoilUi {
    pub fn follow_foil_rename(&mut self, old: &str, new: &str) {
        if let Some(o) = &mut self.opt_3d_origin {
            if o.foil == old {
                o.foil = new.to_string();
            }
        }
        if let Some((_, m)) = &mut self.opt_job {
            m.follow_foil_rename(old, new);
        }
        if let Some((_, m)) = &mut self.opt_result {
            m.follow_foil_rename(old, new);
        }
    }

    pub fn follow_plane_rename(&mut self, old: &str, new: &str) {
        if let Some((pl, _)) = &mut self.opt_from_3d.polar {
            if pl == old {
                *pl = new.to_string();
            }
        }
        if let Some(o) = &mut self.opt_3d_origin {
            if o.plane == old {
                o.plane = new.to_string();
            }
        }
        if let Some((_, m)) = &mut self.opt_job {
            m.follow_plane_rename(old, new);
        }
        if let Some((_, m)) = &mut self.opt_result {
            m.follow_plane_rename(old, new);
        }
    }
}

impl Default for FoilUi {
    fn default() -> Self {
        let opt = OptimSettings {
            parameterization: Parameterization::Cst,
            cst_order: 4,
            cst_max_delta: 0.05,
            operating_points: OptimSettings::default_operating_points(),
            spar: None,
            population: 12,
            iterations: 25,
            ..OptimSettings::default()
        };
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
            bl_arc_length: false,
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
            inv_src: None,
            preview_for: None,
            opt,
            opt_multipoint: true,
            opt_points_saved: Vec::new(),
            opt_force_open: None,
            opt_scroll_to: false,
            opt_sel_pareto: None,
            opt_save_name: String::new(),
            opt_from_3d: From3dPicker::default(),
            opt_3d_origin: None,
            opt_job: None,
            opt_result: None,
            gen: GenUi::default(),
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

    /// The Create section: pick a family, set its parameters with a live
    /// preview in the foil plot, then add it or replace the selected foil.
    fn generator_ui(&mut self, ui: &mut egui::Ui, selected: Option<&Foil>) {
        let fui = &mut self.fui;
        ui.horizontal(|ui| {
            ui.label("family");
            widgets::pick(ui, "genfam", &mut fui.gen.family, &GEN_FAMILIES);
            ui.checkbox(&mut fui.gen.show, "preview").on_hover_text("draw the new foil (dashed) in the plot");
        });
        let g = &mut fui.gen;
        match g.family {
            1 => {
                hint(ui, "y = x^N1 (1-x)^N2 sum A_i B_i(x) + x dz on each surface. A0 sets the LE radius (r = A0^2/2), the last weight the TE wedge.");
                let mut order = g.cst.order();
                egui::Grid::new("cstp").num_columns(2).show(ui, |ui| {
                    if int_keep(ui, "order", &mut order, 0..=30, "Bernstein polynomial order (order + 1 weights per surface)") {
                        for w in [&mut g.cst.upper, &mut g.cst.lower] {
                            let last = *w.last().unwrap_or(&0.1);
                            w.resize(order + 1, last);
                        }
                    }
                    num_keep(ui, "TE y upper", &mut g.cst.te_upper, 0.0002, -10.0..=10.0, "trailing-edge height of the upper surface");
                    num_keep(ui, "TE y lower", &mut g.cst.te_lower, 0.0002, -10.0..=10.0, "trailing-edge height of the lower surface");
                    num_keep(ui, "N1", &mut g.cst.n1, 0.01, 0.0..=1.0, "class exponent at the nose (0.5 = round; must be > 0)");
                    num_keep(ui, "N2", &mut g.cst.n2, 0.01, 0.0..=2.0, "class exponent at the tail (1 = sharp; must be > 0)");
                });
                egui::Grid::new("cstw").num_columns(3).striped(true).show(ui, |ui| {
                    ui.label("");
                    ui.label("upper");
                    ui.label("lower");
                    ui.end_row();
                    for i in 0..g.cst.upper.len().min(g.cst.lower.len()) {
                        ui.label(format!("A{i}"));
                        ui.add(egui::DragValue::new(&mut g.cst.upper[i]).speed(0.001).fixed_decimals(4));
                        ui.add(egui::DragValue::new(&mut g.cst.lower[i]).speed(0.001).fixed_decimals(4));
                        ui.end_row();
                    }
                });
                ui.horizontal(|ui| {
                    ui.add_enabled_ui(selected.is_some(), |ui| {
                        if ui.button("Fit to selected foil").on_hover_text("least-squares CST weights of this order for the selected foil").on_disabled_hover_text("select a foil first").clicked() {
                            if let Some(f) = selected {
                                match Cst::fit(f, g.fit_order) {
                                    Ok((c, fit)) => {
                                        g.cst = c;
                                        g.fit = Some(Ok((f.name.clone(), fit)));
                                    }
                                    Err(e) => g.fit = Some(Err(format!("{}: {e}", f.name))),
                                }
                            }
                        }
                    });
                    ui.label("order");
                    ui.add(egui::DragValue::new(&mut g.fit_order).range(0..=30).clamp_existing_to_range(false));
                });
                match &g.fit {
                    Some(Ok((n, fit))) => hint(ui, &format!("fit of {n}: rms {:.1e}, max {:.1e} (chord fractions)", fit.rms, fit.max)),
                    Some(Err(e)) => {
                        ui.label(RichText::new(format!("fit failed: {e}")).color(Color32::from_rgb(230, 80, 70)));
                    }
                    None => {}
                }
            }
            2 => {
                hint(ui, "z = sum a_n x^(n-1/2) per surface, from the LE radius, the crests and the trailing edge.");
                let p = &mut g.parsec;
                egui::Grid::new("parsec").num_columns(2).show(ui, |ui| {
                    num_keep(ui, "LE radius", &mut p.r_le, 0.0005, 0.0..=0.2, "leading-edge radius / chord");
                    num_keep(ui, "upper crest x", &mut p.x_up, 0.005, 0.01..=0.99, "");
                    num_keep(ui, "upper crest z", &mut p.z_up, 0.001, -1.0..=1.0, "");
                    num_keep(ui, "upper crest z''", &mut p.zxx_up, 0.01, -100.0..=100.0, "curvature at the crest (negative)");
                    num_keep(ui, "lower crest x", &mut p.x_lo, 0.005, 0.01..=0.99, "");
                    num_keep(ui, "lower crest z", &mut p.z_lo, 0.001, -1.0..=1.0, "");
                    num_keep(ui, "lower crest z''", &mut p.zxx_lo, 0.01, -100.0..=100.0, "curvature at the crest (positive)");
                    num_keep(ui, "TE height", &mut p.z_te, 0.001, -1.0..=1.0, "height of the trailing-edge midpoint");
                    num_keep(ui, "TE thickness", &mut p.dz_te, 0.0005, 0.0..=0.5, "");
                    num_keep(ui, "TE direction (deg)", &mut p.alpha_te, 0.1, -44.0..=44.0, "bisector angle below the chord line");
                    num_keep(ui, "TE wedge (deg)", &mut p.beta_te, 0.1, 0.0..=89.0, "included angle");
                });
            }
            3 => {
                hint(ui, "Conformal map of a circle through zeta = 1 centred at (-x, y); TE angle 0 gives the Joukowsky section. Its inviscid flow is known exactly: a reference for the panel solvers.");
                let k = &mut g.kt;
                egui::Grid::new("kt").num_columns(2).show(ui, |ui| {
                    num_keep(ui, "centre x offset", &mut k.x_center, 0.002, 0.001..=0.5, "leftward offset of the circle centre: thickness");
                    num_keep(ui, "centre y offset", &mut k.y_center, 0.002, -0.5..=0.5, "upward offset: camber");
                    num_keep(ui, "TE angle (deg)", &mut k.te_angle, 0.5, 0.0..=60.0, "");
                    num_keep(ui, "alpha (deg)", &mut g.kt_alpha, 0.1, -15.0..=15.0, "angle for the exact solution below");
                });
                // the exact solution costs about half a millisecond: solve it
                // once per parameter set, not on every repaint
                let key = (*k, g.kt_alpha);
                if g.kt_exact.as_ref().map(|c| c.0) != Some(key) {
                    g.kt_exact = Some((key, k.exact(g.kt_alpha).map(|s| (s.cl, s.cm)).map_err(|e| e.to_string())));
                }
                if let Some((_, Ok((cl, cm)))) = &g.kt_exact {
                    hint(ui, &format!("exact inviscid CL {cl:.4}, Cm(c/4) {cm:.4} at {:.1} deg", g.kt_alpha));
                }
            }
            4 => {
                let p = &mut g.plate;
                let mut t = p.thickness * 100.0;
                let mut c = p.camber * 100.0;
                let mut kind = usize::from(p.kind == PlateCamber::Parabolic);
                egui::Grid::new("plate").num_columns(2).show(ui, |ui| {
                    if num_keep(ui, "thickness (%c)", &mut t, 0.1, 0.2..=15.0, "constant plate thickness; the edges are semicircles") {
                        p.thickness = t / 100.0;
                    }
                    if num_keep(ui, "camber (%c)", &mut c, 0.1, -15.0..=15.0, "") {
                        p.camber = c / 100.0;
                    }
                    ui.label("camber line");
                    if widgets::pick(ui, "platekind", &mut kind, &["circular arc", "parabolic"]) {
                        p.kind = if kind == 1 { PlateCamber::Parabolic } else { PlateCamber::CircularArc };
                    }
                    ui.end_row();
                    if p.kind == PlateCamber::Parabolic {
                        num_keep(ui, "x of max camber", &mut p.x_camber, 0.01, 0.1..=0.9, "");
                    }
                    let mut te = usize::from(p.trailing_edge == PlateEdge::Round);
                    ui.label("trailing edge").on_hover_text("square: blunt, cut square to the plate (XFoil can analyse it); round: a semicircle, geometry only (XFoil's viscous solution does not converge on it)");
                    if widgets::pick(ui, "platete", &mut te, &["square", "round"]) {
                        p.trailing_edge = if te == 1 { PlateEdge::Round } else { PlateEdge::Square };
                    }
                    ui.end_row();
                });
            }
            _ => {
                ui.horizontal(|ui| {
                    ui.label("NACA");
                    ui.add(egui::TextEdit::singleline(&mut fui.naca).desired_width(150.0)).on_hover_text(NACA_HELP);
                });
                ui.checkbox(&mut g.closed_te, "closed trailing edge").on_hover_text("4- and 5-digit sections: x^4 coefficient -0.1036 instead of -0.1015");
                hint(ui, NACA_HELP);
            }
        }
        let res = fui.gen.preview(&fui.naca);
        let mut add = None;
        let mut replace = None;
        match &res {
            Ok(f) => {
                let p = f.props();
                hint(ui, &format!("{}: t/c {:.2}% at {:.3}, camber {:.2}% at {:.3}, LE radius {:.4}", f.name, p.thickness * 100.0, p.x_thickness, p.camber * 100.0, p.x_camber, p.le_radius));
                ui.horizontal(|ui| {
                    if ui.button("Add").on_hover_text("add as a new foil").clicked() {
                        add = Some(f.clone());
                    }
                    if let Some(sel) = selected {
                        if ui.button("Replace selected").on_hover_text(format!("give '{}' this shape (keeps its name and hinge; its polars are cleared)", sel.name)).clicked() {
                            replace = Some((sel.clone(), f.clone()));
                        }
                    }
                });
            }
            Err(e) => {
                ui.label(RichText::new(e).color(Color32::from_rgb(230, 80, 70)));
            }
        }
        if let Some(r) = selected.and_then(|f| f.recipe.as_ref()) {
            if ui.button(format!("Edit parameters of '{}'", selected.unwrap().name)).on_hover_text(format!("load its {} recipe into the editor", r.family_label())).clicked() {
                let r = r.clone();
                if let Err(e) = self.fui.gen.load(&r, &mut self.fui.naca) {
                    self.notify(e, true);
                }
            }
        }
        if let Some(mut f) = add {
            self.snapshot_foil_geometry();
            f.name = self.project.unique_foil_name(&f.name);
            self.sel_foil = Some(f.name.clone());
            self.project.foils.push(f);
        }
        if let Some((sel, mut f)) = replace {
            self.snapshot_foil_geometry();
            f.name = sel.name.clone();
            f.hinge = sel.hinge;
            self.project.upsert_foil(f);
            self.fui.preview = None;
            self.fui.preview_for = None;
            self.fui.qdist = None;
            self.fui.inv_target = None;
            self.fui.inv_src = None;
        }
    }

    pub fn foil_design(&mut self, ui: &mut egui::Ui) {
        egui::SidePanel::left("foil_design_left").resizable(true).default_width(390.0).show_inside(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                heading(ui, "Foil");
                self.foil_picker(ui);
                let foil = self.selected_foil();
                // a preview and a speed distribution belong to the foil and
                // geometry they were made from (and a distribution to its
                // solver settings): drop them as soon as any of that changes
                // (another foil selected or added, an edit, undo, redo),
                // before any panel or plot shows or applies them
                let geom = foil.as_ref().map(|f| (f.name.clone(), node_hash(&f.x, &f.y)));
                if self.fui.preview.is_some() && self.fui.preview_for != geom {
                    self.fui.preview = None;
                    self.fui.preview_for = None;
                }
                let key = geom.as_ref().map(|(n, h)| InvKey { foil: n.clone(), hash: *h, mode: self.fui.inv_mode, alpha: self.fui.inv_alpha, paneling: self.fui.repanel });
                if self.fui.inv_src.is_some() && self.fui.inv_src != key {
                    self.fui.qdist = None;
                    self.fui.inv_target = None;
                    self.fui.inv_src = None;
                }
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
                self.generator_ui(ui, foil.as_ref());
                if ui.button("Load .dat files...").on_hover_text("Selig or Lednicer coordinate files; clockwise files are reordered automatically").clicked() {
                    if let Some(files) = crate::file_dialog(ui.ctx()).and_then(|dialog| dialog.add_filter("Airfoil", &["dat", "txt"]).pick_files()) {
                        self.snapshot_foil_geometry();
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
                    let key = InvKey { foil: foil.name.clone(), hash: node_hash(&foil.x, &foil.y), mode: self.fui.inv_mode, alpha: self.fui.inv_alpha, paneling: self.fui.repanel };
                    if ui.button("Get speed distribution").clicked() {
                        self.fui.inv_target = None;
                        let k = JobKind::Inverse { foil: foil.clone(), alpha: self.fui.inv_alpha, mode: self.fui.inv_mode, paneling: self.fui.repanel, target: None, iterations: 0 };
                        self.submit(format!("{} speed distribution", foil.name), k, Pending::Inverse(key.clone()));
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
                            self.submit(format!("{} inverse design", foil.name), k, Pending::Inverse(key.clone()));
                        }
                    }
                });
                let mut save_and_rerun: Option<(Candidate, OptimRunMeta)> = None;
                let mut save_new_cand: Option<(Candidate, OptimRunMeta)> = None;
                let force_open = self.fui.opt_force_open.take();
                let opt_header = egui::CollapsingHeader::new(RichText::new("Optimiser (MOPSO / genetic)").strong())
                    .default_open(false)
                    .open(force_open)
                    .show(ui, |ui| {
                        self.optimiser_ui(ui, &foil, &mut preview, &mut save_new_cand, &mut save_and_rerun);
                    });
                if self.fui.opt_scroll_to {
                    opt_header.header_response.scroll_to_me(Some(egui::Align::TOP));
                    self.fui.opt_scroll_to = false;
                }
                if let Some((cand, meta)) = save_new_cand {
                    let name = self.save_optim_candidate(&cand, &meta.base_foil);
                    self.sel_foil = Some(name.clone());
                    self.fui.preview = None;
                    self.notify(format!("saved optimised foil '{name}'"), false);
                }
                if let Some((cand, meta)) = save_and_rerun {
                    self.save_and_rerun_3d_study(&cand, &meta);
                }
                if let Some(r) = preview {
                    match r {
                        Ok(p) => {
                            self.fui.new_name = self.project.unique_foil_name(&format!("{} mod", foil.name));
                            self.fui.preview = Some(p);
                            self.fui.preview_for = Some((foil.name.clone(), node_hash(&foil.x, &foil.y)));
                        }
                        Err(err) => self.notify(err, true),
                    }
                }
                if let Some(pv) = self.fui.preview.clone() {
                    ui.separator();
                    ui.label(RichText::new(format!("Preview: {}", pv.note)).color(Color32::from_rgb(255, 160, 60)));
                    ui.horizontal(|ui| {
                        if ui.button("Replace foil").clicked() {
                            self.snapshot_foil_geometry();
                            let mut f = pv.clone();
                            f.name = foil.name.clone();
                            self.project.upsert_foil(f);
                            self.fui.preview = None;
                        }
                        ui.add(egui::TextEdit::singleline(&mut self.fui.new_name).desired_width(110.0));
                        if ui.button("Save as new").clicked() {
                            self.snapshot_foil_geometry();
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
                        if let Some(p) = crate::file_dialog(ui.ctx()).and_then(|dialog| dialog.set_file_name(format!("{}.dat", foil.name)).save_file()) {
                            if let Err(e) = foil.save(&p) {
                                self.notify(e.to_string(), true);
                            }
                        }
                    }
                    let mut name = foil.name.clone();
                    if ui.add(egui::TextEdit::singleline(&mut name).desired_width(120.0)).on_hover_text("rename (updates polars and wings using it)").lost_focus() && name != foil.name && !name.trim().is_empty() && self.project.foil(&name).is_none() {
                        self.snapshot_foil_geometry();
                        self.rename_foil(&foil.name, &name);
                        self.sel_foil = Some(name);
                    }
                });
            });
        });
        egui::CentralPanel::default().show_inside(ui, |ui| {
            let foil = self.selected_foil();
            let pv = self.fui.preview.clone();
            let new_foil = self.fui.gen.shown().cloned();
            let live = self.live_optim();
            let sel_cand = live.as_ref().and_then(|(r, _)| selected_optim_candidate(r, self.fui.opt_sel_pareto)).cloned();
            let active_spar = live.as_ref().and_then(|(_, m)| m.spar).or(self.fui.opt.spar);
            let extra = self.fui.qdist.is_some() || live.is_some();
            let h = ui.available_height();
            let plot_h = if extra { (h * 0.50).max(140.0) } else { h };
            Plot::new("foilshape").data_aspect(1.0).height(plot_h).legend(Legend::default()).show(ui, |pu| {
                if let Some(cand) = &sel_cand {
                    pu.line(Line::new("optimiser candidate", PlotPoints::from(foil_lines(&cand.foil))).color(Color32::from_rgb(90, 220, 120)).width(1.8_f32));
                }
                if let Some(f) = &foil {
                    pu.line(Line::new(f.name.clone(), PlotPoints::from(foil_lines(f))).color(series_color(0)).width(2.0_f32));
                    let (xs, c, _) = f.camber_thickness(60);
                    pu.line(Line::new("camber line", PlotPoints::from(xs.iter().zip(&c).map(|(x, y)| [*x, *y]).collect::<Vec<_>>())).color(Color32::GRAY).style(egui_plot::LineStyle::dashed_dense()));
                    pu.points(Points::new("nodes", PlotPoints::from(foil_lines(f))).radius(1.5_f32).color(series_color(0)));
                    if let Some((x, yr)) = f.hinge {
                        let (u, l) = f.surface_y(x);
                        pu.points(Points::new("hinge", PlotPoints::from(vec![[x, l + yr * (u - l)]])).radius(5.0_f32).color(Color32::YELLOW));
                    }
                    if let Some(sp) = active_spar {
                        let ref_foil = sel_cand.as_ref().map(|c| &c.foil).unwrap_or(f);
                        let (u, l) = ref_foil.surface_y(sp.x);
                        let mid = 0.5 * (u + l);
                        let half = 0.5 * sp.min_thickness;
                        let spar_pts = vec![[sp.x, mid - half], [sp.x, mid + half]];
                        pu.line(Line::new(format!("spar min {:.1}%", sp.min_thickness * 100.0), PlotPoints::from(spar_pts.clone())).color(Color32::from_rgb(255, 210, 70)).width(3.0_f32));
                        pu.points(Points::new("spar limits", PlotPoints::from(spar_pts)).radius(3.5_f32).color(Color32::from_rgb(255, 210, 70)));
                    }
                }
                if let Some(p) = &pv {
                    pu.line(Line::new("preview", PlotPoints::from(foil_lines(p))).color(Color32::from_rgb(255, 160, 60)).width(2.0_f32));
                }
                if let Some(g) = &new_foil {
                    pu.line(Line::new(format!("new: {}", g.name), PlotPoints::from(foil_lines(g))).color(Color32::from_rgb(230, 110, 230)).width(1.5_f32).style(egui_plot::LineStyle::dashed_loose()));
                }
            });
            if extra {
                let w = ui.available_width();
                let hh = ui.available_height().max(130.0);
                ui.allocate_ui_with_layout(egui::vec2(w, hh), egui::Layout::left_to_right(egui::Align::TOP), |ui| {
                    if let Some(d) = &self.fui.qdist {
                        let qw = if live.is_some() { w * 0.40 } else { w };
                        ui.allocate_ui_with_layout(egui::vec2(qw, hh), egui::Layout::top_down(egui::Align::LEFT), |ui| {
                            ui.label(RichText::new(format!("Surface speed Q/Vinf at alpha {:.2} (CL {:.3})", d.alpha, d.cl)).strong());
                            Plot::new("qdist").height((hh - 24.0).max(80.0)).legend(Legend::default()).x_axis_label("x/c").show(ui, |pu| {
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
                    if let Some((r, meta)) = &live {
                        let avail_w = if self.fui.qdist.is_some() { (w * 0.58).max(200.0) } else { w };
                        self.optim_results_panel(ui, r, meta, avail_w, hh);
                    }
                });
            }
        });
    }

    pub fn foil_analysis(&mut self, ui: &mut egui::Ui) {
        let mut labels = polar_labels(&self.project.foil_polars);
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
                    let mut separate = s.ncrit_bot.is_some();
                    if ui.checkbox(&mut separate, "separate bottom Ncrit").changed() {
                        s.ncrit_bot = separate.then_some(s.ncrit);
                    }
                    if let Some(n) = &mut s.ncrit_bot {
                        num(ui, "Ncrit bottom", n, 0.1, 0.1..=20.0, "e^N transition criterion on the lower surface (XFoil 6.99 VPAR NB)");
                    }
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
                        let mut keep = s.paneling == NodeInput::Preserve;
                        if check(ui, "keep input nodes", &mut keep, "solve on the foil's own coordinates (XFoil LOAD + PCOP) instead of repaneling them (PANE)") {
                            s.paneling = if keep { NodeInput::Preserve } else { NodeInput::default() };
                        }
                        if let NodeInput::Pangen(p) = &mut s.paneling {
                            int(ui, "panel nodes", &mut p.panels, 40..=350, "airfoil panel nodes (XFoil PANE)");
                        }
                        num(ui, "wake length / chord", &mut s.wake_length, 0.1, 1.0..=20.0, "wake extent behind the trailing edge (XFoil 6.99 VPAR WAKE, default 1); beyond about 4 chords the wake array limit spreads the same number of points further apart");
                        num(ui, "wake dissipation", &mut s.wake_dissipation, 0.01, 0.01..=10.0, "wake inverse dissipation-length factor (XFoil 6.99 VPAR LAG, default 0.9)");
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
                    for label in &labels {
                        self.fui.visible.insert(label.key.clone());
                    }
                }
                if ui.small_button("hide all").clicked() {
                    self.fui.visible.clear();
                }
                if ui.small_button("only this foil").clicked() {
                    self.fui.visible.clear();
                    for (p, label) in self.project.foil_polars.iter().zip(&labels) {
                        if Some(&p.foil) == self.sel_foil.as_ref() {
                            self.fui.visible.insert(label.key.clone());
                        }
                    }
                }
            });
            ui.checkbox(&mut self.fui.show_failed, "mark unconverged points");
            egui::ScrollArea::vertical().show(ui, |ui| {
                let mut del = None;
                let mut vis_idx = 0;
                for (i, p) in self.project.foil_polars.iter().enumerate() {
                    let key = labels[i].key.clone();
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
                        let r = ui.label(RichText::new(format!("{}\n{}", p.foil, labels[i].name)).small());
                        r.on_hover_text(format!("{} points, {} not converged\n{}", p.points.len(), p.failed.len(), if p.failed.is_empty() { String::new() } else { format!("not converged at {:?}", p.failed) }));
                        ui.menu_button("...", |ui| {
                            if ui.button("Export CSV").clicked() {
                                if let Some(path) = crate::file_dialog(ui.ctx()).and_then(|dialog| dialog.set_file_name(format!("{}.csv", p.export_stem())).save_file()) {
                                    let _ = std::fs::write(path, p.to_csv());
                                }
                                ui.close();
                            }
                            if ui.button("Export XFoil text").clicked() {
                                if let Some(path) = crate::file_dialog(ui.ctx()).and_then(|dialog| dialog.set_file_name(format!("{}.txt", p.export_stem())).save_file()) {
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
                    self.snapshot_foil_polars();
                    self.project.foil_polars.remove(i);
                    labels.remove(i);
                }
            });
        });
        egui::CentralPanel::default().show_inside(ui, |ui| {
            let has_detail = self.point_detail.is_some();
            let h = ui.available_height();
            let ph = if has_detail { h * 0.58 } else { h };
            let w = ui.available_width();
            let polars: Vec<_> = self.project.foil_polars.iter().zip(&labels).filter(|(_, label)| self.fui.visible.contains(&label.key)).map(|(p, label)| (p.clone(), label)).collect();
            let mut clicked: Option<(String, String, FoilPoint)> = None;
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
                            for (i, (p, label)) in polars.iter().enumerate() {
                                let pts: Vec<[f64; 2]> = p.points.iter().map(|q| [var(q, xv), var(q, yv)]).collect();
                                let name = format!("{} {}", p.foil, label.name);
                                pu.line(Line::new(name.clone(), PlotPoints::from(pts.clone())).color(series_color(i)));
                                pu.points(Points::new(name.clone(), PlotPoints::from(pts)).radius(2.0_f32).color(series_color(i)));
                                if self.fui.show_failed && xv == 0 && !p.failed.is_empty() && p.polar_type <= 3 {
                                    // unconverged alphas as crosses on the x axis value line
                                    let lo = p.points.iter().map(|q| var(q, yv)).fold(f64::MAX, f64::min);
                                    let fpts: Vec<[f64; 2]> = p.failed.iter().map(|a| [*a, lo]).collect();
                                    pu.points(Points::new(format!("{name} not converged"), PlotPoints::from(fpts)).shape(egui_plot::MarkerShape::Cross).radius(5.0_f32).color(Color32::from_rgb(230, 80, 70)));
                                }
                            }
                            if pu.response().clicked() {
                                if let Some(c) = pu.pointer_coordinate() {
                                    let b = pu.plot_bounds();
                                    let (sx, sy) = (b.width().max(1e-12), b.height().max(1e-12));
                                    let mut best = (f64::MAX, None);
                                    for (p, label) in &polars {
                                        for q in &p.points {
                                            let d = ((var(q, xv) - c.x) / sx).powi(2) + ((var(q, yv) - c.y) / sy).powi(2);
                                            if d < best.0 {
                                                best = (d, Some((label.key.0.clone(), label.key.1.clone(), *q)));
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
            if let Some((foil, pol, point)) = clicked {
                self.show_point(&foil, &pol, point);
            }
            if let Some(d) = self.point_detail.clone() {
                ui.separator();
                self.detail_plots(ui, &d, w);
            }
        });
    }

    /// Show the distributions of a polar point (stored, or recomputed in the background).
    pub(crate) fn show_point(&mut self, foil: &str, identity: &str, point: FoilPoint) {
        if let Some(ds) = self.details.get(&(foil.to_string(), identity.to_string())) {
            if let Some(d) = ds.iter().find(|d| d.point == point) {
                self.point_detail = Some(d.clone());
                return;
            }
        }
        let Some(p) = self.project.foil_polars.iter().find(|p| p.key() == (foil.to_string(), identity.to_string())).cloned() else { return };
        let Some(mut f) = self.project.foil(foil).cloned() else { return };
        let mut st = p.settings.clone().unwrap_or_else(|| AnalysisSettings {
            re: p.re, mach: p.mach, ncrit: p.ncrit, ncrit_bot: p.ncrit_bot, xtr_top: p.xtr_top, xtr_bot: p.xtr_bot, viscous: p.viscous,
            ..Default::default()
        });
        // Recompute this operating point, including its exact Re/Mach and
        // deflected geometry, rather than the current editor's settings.
        let flap = if p.polar_type == 6 { point.flap } else { p.flap_deg };
        if flap != 0.0 {
            let hinge = if p.polar_type == 6 { st.flap_hinge } else { f.hinge.unwrap_or(st.flap_hinge) };
            match f.with_flap(hinge.0, hinge.1, flap) {
                Ok(deflected) => f = deflected,
                Err(e) => { self.notify(e.to_string(), true); return; }
            }
        }
        st.polar_type = PolarType::FixedSpeed;
        st.re = point.re;
        st.mach = point.mach;
        self.submit(format!("{foil} alpha = {:.2} ({})", point.alpha, p.display_name()), JobKind::FoilPoint { foil: f, settings: st, value: point.alpha, is_alpha: true }, Pending::FoilPoint);
    }

    fn detail_plots(&mut self, ui: &mut egui::Ui, d: &PointDetail, w: f32) {
        let q = &d.point;
        ui.horizontal(|ui| {
            ui.label(RichText::new(format!("alpha {:.2}  CL {:.4}  CD {:.5}  Cm {:.4}  L/D {:.1}  Xtr {:.3}/{:.3}  Re {:.0}", q.alpha, q.cl, q.cd, q.cm, q.ld(), q.xtr_top, q.xtr_bot, q.re)).strong());
            ui.checkbox(&mut self.fui.cp_inviscid, "inviscid Cp");
            ui.checkbox(&mut self.fui.bl_arc_length, "BL versus arc length");
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
                    pu.line(Line::new("-Cp viscous", PlotPoints::from(pts)).color(series_color(0)).width(2.0_f32));
                    if self.fui.cp_inviscid {
                        let pts: Vec<[f64; 2]> = d.x.iter().zip(&d.cp_inv).map(|(x, c)| [*x, -c]).collect();
                        pu.line(Line::new("-Cp inviscid", PlotPoints::from(pts)).color(Color32::GRAY).style(egui_plot::LineStyle::dashed_loose()));
                    }
                });
            });
            ui.vertical(|ui| {
                ui.set_width(w / 2.0 - 6.0);
                let v = self.fui.bl_var;
                Plot::new("bl").height(h).legend(Legend::default()).x_axis_label(if self.fui.bl_arc_length { "s/c" } else { "x/c" }).show(ui, |pu| {
                    for (i, (name, b)) in [("top", &d.top), ("bottom", &d.bot), ("wake", &d.wake)].iter().enumerate() {
                        let ys = match v {
                            0 => &b.ue,
                            1 => &b.h,
                            2 => &b.cf,
                            3 => &b.dstar,
                            4 => &b.theta,
                            _ => &b.n_or_ctau,
                        };
                        let coordinate = if self.fui.bl_arc_length { &b.s } else { &b.x };
                        let pts: Vec<[f64; 2]> = coordinate.iter().zip(ys).map(|(x, y)| [*x, *y]).collect();
                        pu.line(Line::new(*name, PlotPoints::from(pts)).color(series_color(i)));
                    }
                });
            });
        });
    }

    /// Currently displayed optimiser result (live running partial or finished
    /// result) paired with the frozen metadata of the job that produced it.
    fn live_optim(&self) -> Option<(OptimResult, OptimRunMeta)> {
        if let Some((id, meta)) = &self.fui.opt_job {
            if let Some(r) = self.jobs.get(*id).and_then(|j| j.partial.lock().unwrap().optim.clone()) {
                return Some((r, meta.clone()));
            }
        }
        self.fui.opt_result.clone()
    }

    /// Resolve the 2D analysis settings (Mach, Ncrit, bottom Ncrit, trips, wake
    /// controls) of the polar family actually used for `foil_name` in `polar`.
    fn resolve_source_foil_analysis(&mut self, polar: &foiley_plane::PlanePolar, foil_name: &str) -> Option<AnalysisSettings> {
        if foil_name.is_empty() {
            return None;
        }
        let sel = polar.def.as_ref().map(|d| d.family).unwrap_or_default();
        let db = self.foil_db();
        let base = AnalysisSettings::default();
        if let Ok(Some(id)) = db.resolve(foil_name, &sel) {
            let fam = db.family(id);
            let k = fam.key;
            let mut st = fam.polars.first().and_then(|p| p.settings.clone()).unwrap_or(base);
            st.mach = k.mach;
            st.ncrit = k.ncrit;
            st.ncrit_bot = k.ncrit_bot;
            st.xtr_top = k.xtr_top;
            st.xtr_bot = k.xtr_bot;
            return Some(st);
        }
        // Fallback: if the polar definition specifies explicit family settings, carry them.
        if sel != foiley_plane::FamilySelector::default() {
            return Some(AnalysisSettings {
                mach: sel.mach.unwrap_or(base.mach),
                ncrit: sel.ncrit.unwrap_or(base.ncrit),
                ncrit_bot: sel.ncrit_bot,
                xtr_top: sel.xtr_top.unwrap_or(base.xtr_top),
                xtr_bot: sel.xtr_bot.unwrap_or(base.xtr_bot),
                ..base
            });
        }
        None
    }

    /// Jump from 3D analysis to the 2D airfoil optimiser, pre-loading slow flight,
    /// cruise and fast glide from `(wing, span_pos)` on `polar_key`.
    pub(crate) fn open_section_optimiser(&mut self, polar_key: (String, String), wing: usize, span_pos: f64) {
        let Some(polar) = self.project.plane_polars.iter().find(|p| p.plane == polar_key.0 && p.name == polar_key.1).cloned() else {
            return;
        };
        let plane = self.project.planes.iter().find(|pl| pl.name == polar.plane);
        let stations = polar.wing_stations(plane, wing);
        let station = stations
            .iter()
            .min_by(|a, b| (a.span_pos - span_pos).abs().total_cmp(&(b.span_pos - span_pos).abs()))
            .cloned();
        let actual_y = station.as_ref().map(|s| s.span_pos).unwrap_or(span_pos);
        let regimes = polar.pick_flight_regimes(wing, actual_y);
        let selected: BTreeSet<usize> = regimes.iter().map(|r| r.point_index).collect();

        self.fui.opt_from_3d.polar = Some(polar_key.clone());
        self.fui.opt_from_3d.wing = wing;
        self.fui.opt_from_3d.span_pos = Some(actual_y);
        self.fui.opt_from_3d.selected_points = selected.clone();
        self.fui.opt_from_3d.force_open = Some(true);
        self.fui.opt_force_open = Some(true);
        self.fui.opt_scroll_to = true;
        self.module = crate::Module::FoilDesign;

        if selected.is_empty() {
            self.notify("selected 3D polar has no valid operating points at this station", true);
            return;
        }
        self.fui.opt_multipoint = true;
        self.apply_3d_station_points(&polar_key, wing, actual_y, &selected);
    }

    /// Populate `self.fui.opt.operating_points` from the selected points of a 3D polar
    /// at `(wing, span_pos)`, aligning the selected foil and 2D flow settings.
    fn apply_3d_station_points(&mut self, polar_key: &(String, String), wing: usize, span_pos: f64, selected_points: &BTreeSet<usize>) {
        let Some(polar) = self.project.plane_polars.iter().find(|p| p.plane == polar_key.0 && p.name == polar_key.1).cloned() else {
            return;
        };
        let plane = self.project.planes.iter().find(|pl| pl.name == polar.plane);
        let stations = polar.wing_stations(plane, wing);
        let station = stations
            .iter()
            .min_by(|a, b| (a.span_pos - span_pos).abs().total_cmp(&(b.span_pos - span_pos).abs()))
            .cloned();
        let actual_y = station.as_ref().map(|s| s.span_pos).unwrap_or(span_pos);
        let wing_name = station.as_ref().map(|s| s.wing_name.clone()).unwrap_or_else(|| format!("wing {wing}"));
        let foil_name = station.as_ref().and_then(|s| s.foil.clone()).unwrap_or_default();

        // Keep exact regime ordering, labels and default goals when the selected
        // points match `pick_flight_regimes` (guaranteed agreement with CLI).
        let regimes = polar.pick_flight_regimes(wing, actual_y);
        let regime_set: BTreeSet<usize> = regimes.iter().map(|r| r.point_index).collect();
        let rows: Vec<OperatingPoint> = if !regimes.is_empty() && *selected_points == regime_set {
            regimes.iter().map(|sp| sp.to_default_operating_point()).collect()
        } else {
            let regime_map: HashMap<usize, (String, DragGoal)> = regimes.into_iter().map(|r| (r.point_index, (r.label, r.default_goal))).collect();
            let mut pts: Vec<foiley_plane::StationPoint> = selected_points
                .iter()
                .filter_map(|&idx| {
                    let fallback = polar.points.get(idx).map(|o| format!("a={:.1} V={:.1}", o.alpha, o.velocity)).unwrap_or_default();
                    let (label, goal) = regime_map.get(&idx).cloned().unwrap_or((fallback, DragGoal::Minimise));
                    let mut sp = polar.station_point(idx, wing, actual_y, label)?;
                    sp.default_goal = goal;
                    Some(sp)
                })
                .collect();
            if pts.is_empty() {
                self.notify("pick at least one 3D operating point first", true);
                return;
            }
            // Order consistently with the polar's flight envelope (ascending speed
            // when speed varies, otherwise descending whole-aircraft CL).
            let (v_min, v_max) = pts.iter().fold((f64::INFINITY, 0.0f64), |(lo, hi), p| (lo.min(p.velocity), hi.max(p.velocity)));
            let v_varies = (v_max - v_min) / v_max.max(1e-9) > 0.05;
            if v_varies {
                pts.sort_by(|a, b| a.velocity.total_cmp(&b.velocity).then_with(|| b.plane_cl.total_cmp(&a.plane_cl)));
            } else {
                pts.sort_by(|a, b| b.plane_cl.total_cmp(&a.plane_cl));
            }
            let any_cap = pts.iter().any(|p| p.default_goal.is_cap());
            let n = pts.len();
            pts.iter()
                .enumerate()
                .map(|(i, sp)| {
                    let goal = if any_cap {
                        sp.default_goal
                    } else if n >= 2 && i == 0 {
                        DragGoal::AtMostOriginal
                    } else {
                        DragGoal::Minimise
                    };
                    sp.to_operating_point(goal, 1.0)
                })
                .collect()
        };
        let n = rows.len();

        // Detach any previous optimiser job/result so an earlier candidate can
        // never be installed onto this newly loaded wing station.
        if let Some((old_id, _)) = self.fui.opt_job.take() {
            if let Some(j) = self.jobs.get(old_id) {
                if !j.status().is_finished() {
                    j.cancel();
                }
            }
        }
        self.fui.opt.operating_points = rows;
        self.fui.opt_multipoint = true;
        self.fui.opt_sel_pareto = None;
        self.fui.opt_result = None;

        let source_analysis = self.resolve_source_foil_analysis(&polar, &foil_name);
        if !foil_name.is_empty() && self.project.foil(&foil_name).is_some() {
            self.sel_foil = Some(foil_name.clone());
            self.fui.preview = None;
            self.fui.opt_save_name = self.project.unique_foil_name(&format!("{foil_name} opt"));
            // Inherit resolved Mach/Ncrit/trips/wake controls so 2D and 3D agree.
            // With no resolvable source, reset the controls the optimiser panel does not
            // show so an earlier station's trips/wake cannot silently carry over.
            let d = AnalysisSettings::default();
            let src = source_analysis.as_ref().unwrap_or(&d);
            if source_analysis.is_some() {
                self.fui.opt.analysis.mach = src.mach;
                self.fui.opt.analysis.ncrit = src.ncrit;
            }
            self.fui.opt.analysis.ncrit_bot = src.ncrit_bot;
            self.fui.opt.analysis.xtr_top = src.xtr_top;
            self.fui.opt.analysis.xtr_bot = src.xtr_bot;
            self.fui.opt.analysis.wake_length = src.wake_length;
            self.fui.opt.analysis.wake_dissipation = src.wake_dissipation;
            // Set a practical default spar thickness constraint without panicking on thin foils.
            if let Some(f) = self.project.foil(&foil_name) {
                let spar_x = self.fui.opt.spar.map(|sp| sp.x).unwrap_or(0.30);
                let cur_t = f.thickness_at(spar_x);
                if cur_t > 0.005 {
                    let default_min = default_spar_min(cur_t);
                    match &mut self.fui.opt.spar {
                        None => self.fui.opt.spar = Some(SparConstraint { x: spar_x, min_thickness: default_min }),
                        Some(sp) if sp.min_thickness > cur_t + 1e-4 => sp.min_thickness = default_min,
                        _ => {}
                    }
                }
            }
        }
        self.fui.opt_3d_origin = Some(SectionOptimOrigin {
            plane: polar.plane.clone(),
            polar_name: polar.name.clone(),
            polar_def: polar.def.clone(),
            wing,
            wing_name: wing_name.clone(),
            span_pos: actual_y,
            foil: foil_name,
            source_analysis,
        });
        self.notify(format!("loaded {} operating points from {} ({wing_name} y = {:.3} m)", n, polar.name, actual_y), false);
    }

    fn save_optim_candidate(&mut self, cand: &Candidate, base_name: &str) -> String {
        self.snapshot_foil_geometry();
        let desired = if self.fui.opt_save_name.trim().is_empty() {
            format!("{base_name} opt")
        } else {
            self.fui.opt_save_name.trim().to_string()
        };
        let name = self.project.unique_foil_name(&desired);
        let mut f = cand.foil.clone();
        f.name = name.clone();
        if !cand.op_points.is_empty() {
            let summary: Vec<String> = cand.op_points.iter().map(|e| format!("{} CD {:.5}", e.label, e.cd)).collect();
            f.note = format!("optimised from {base_name}: {}", summary.join(", "));
        } else {
            f.note = format!("optimised from {base_name}: CL {:.3} CD {:.5}", cand.cl, cand.cd);
        }
        self.project.foils.push(f);
        self.fui.opt_save_name = self.project.unique_foil_name(&format!("{base_name} opt"));
        name
    }

    /// Check whether `origin` still points to a valid source plane, wing,
    /// matching wing section foil, and 3D polar definition.
    fn validate_rerun_origin(&self, origin: &SectionOptimOrigin, base_name: &str) -> Result<(foiley_plane::Plane, foiley_plane::PlanePolarDef, usize, String), String> {
        let source_plane = self
            .project
            .planes
            .iter()
            .find(|p| p.name == origin.plane)
            .cloned()
            .ok_or_else(|| format!("source plane '{}' is no longer in the project", origin.plane))?;
        let source_polar = self.project.plane_polars.iter().find(|p| p.plane == origin.plane && p.name == origin.polar_name);
        let def = source_polar
            .and_then(|p| p.def.clone())
            .or_else(|| origin.polar_def.clone())
            .ok_or_else(|| format!("source 3D polar '{}' definition is no longer available", origin.polar_name))?;
        let target_wing = if def.method == foiley_plane::Method::Llt && origin.wing == 0 {
            source_plane.wings.iter().position(|w| w.kind == foiley_plane::WingKind::Main).unwrap_or(0)
        } else {
            origin.wing
        };
        let wing = source_plane
            .wings
            .get(target_wing)
            .ok_or_else(|| format!("wing {} no longer exists on '{}'", target_wing, source_plane.name))?;
        let target_foil = if origin.foil.is_empty() { base_name.to_string() } else { origin.foil.clone() };
        let has_match = wing.sections.iter().any(|sec| sec.foil_right == target_foil || sec.foil_left == target_foil);
        if !has_match {
            return Err(format!("wing '{}' on '{}' no longer uses foil '{}'", wing.name, source_plane.name, target_foil));
        }
        Ok((source_plane, def, target_wing, target_foil))
    }

    /// Save the candidate as a new foil, duplicate the source 3D plane using that
    /// foil on the optimised wing, keep the original 3D polar visible, and rerun
    /// the same 3D analysis with matching transition settings so whole-aircraft
    /// performance can be compared directly.
    fn save_and_rerun_3d_study(&mut self, cand: &Candidate, meta: &OptimRunMeta) {
        if !cand.constraints_met() {
            self.notify("candidate violates a drag cap or spar constraint and cannot be saved for rerun", true);
            return;
        }
        let Some(origin) = meta.origin.clone() else {
            self.notify("no source 3D study linked to this run; save the foil and assign it in Plane design", true);
            return;
        };
        let (source_plane, def, target_wing, target_foil) = match self.validate_rerun_origin(&origin, &meta.base_foil) {
            Ok(v) => v,
            Err(e) => {
                self.notify(e, true);
                return;
            }
        };

        let new_foil = self.save_optim_candidate(cand, &meta.base_foil);
        let mut new_plane = source_plane.clone();
        new_plane.name = self.project.unique_plane_name(&format!("{} ({new_foil})", source_plane.name));

        let mut swapped = 0usize;
        if let Some(w) = new_plane.wings.get_mut(target_wing) {
            for sec in &mut w.sections {
                if sec.foil_right == target_foil {
                    sec.foil_right = new_foil.clone();
                    swapped += 1;
                }
                if sec.foil_left == target_foil {
                    sec.foil_left = new_foil.clone();
                    swapped += 1;
                }
            }
        }
        if swapped == 0 {
            self.notify(format!("no wing section matched '{target_foil}'"), true);
            return;
        }
        self.project.planes.push(new_plane.clone());
        self.sel_plane = self.project.planes.len() - 1;
        // Keep the original 3D polar visible for direct comparison on every plot.
        if self.project.plane_polars.iter().any(|p| p.plane == origin.plane && p.name == origin.polar_name) {
            self.pui.visible.insert((origin.plane.clone(), origin.polar_name.clone()));
        }
        let mut foil_defaults = HashMap::new();
        if let Some(src_an) = origin.source_analysis.clone() {
            foil_defaults.insert(new_foil.clone(), src_an);
        }
        let auto_polars = self.pui.auto_polars.max(6);
        let def_name = if def.name.is_empty() { def.auto_name() } else { def.name.clone() };
        self.module = crate::Module::PlaneAnalysis;
        self.notify(format!("saved '{new_foil}' -> rerunning {def_name} on '{}'", new_plane.name), false);
        self.run_3d_with_def(&new_plane, def, auto_polars, &foil_defaults);
    }

    fn optimiser_ui(
        &mut self,
        ui: &mut egui::Ui,
        foil: &Foil,
        preview: &mut Option<Result<Foil, String>>,
        save_new_cand: &mut Option<(Candidate, OptimRunMeta)>,
        save_and_rerun: &mut Option<(Candidate, OptimRunMeta)>,
    ) {
        hint(ui, "Optimise the foil shape with XFoil across one or several operating points (Re, target CL) with drag goals and a spar thickness floor. Every evaluation runs in the background.");

        // Mode toggle: Operating-point table (multi-point) vs Single alpha
        let prev_mp = self.fui.opt_multipoint;
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.fui.opt_multipoint, true, "Operating points");
            ui.selectable_value(&mut self.fui.opt_multipoint, false, "Single alpha");
        });
        if self.fui.opt_multipoint != prev_mp {
            if self.fui.opt_multipoint {
                self.fui.opt.operating_points = if self.fui.opt_points_saved.is_empty() {
                    OptimSettings::default_operating_points()
                } else {
                    std::mem::take(&mut self.fui.opt_points_saved)
                };
            } else {
                self.fui.opt_points_saved = std::mem::take(&mut self.fui.opt.operating_points);
            }
        }

        // "Use points from 3D analysis" panel
        if self.fui.opt_multipoint {
            self.picker_3d_ui(ui);
        }

        // Operating-point table (multi-point mode)
        if self.fui.opt_multipoint {
            ui.add_space(4.0);
            ui.label(RichText::new("Operating points").strong());
            let mut remove_row: Option<usize> = None;
            egui::ScrollArea::horizontal().id_salt("opt_table_scroll").show(ui, |ui| {
                egui::Grid::new("opt_pts").striped(true).spacing([4.0, 3.0]).show(ui, |ui| {
                    for h in ["", "Point", "Re", "CL", "Drag goal", "Wt", ""] {
                        ui.label(RichText::new(h).small().strong());
                    }
                    ui.end_row();
                    let n = self.fui.opt.operating_points.len();
                    for (i, op) in self.fui.opt.operating_points.iter_mut().enumerate() {
                        ui.checkbox(&mut op.enabled, "");
                        ui.add(egui::TextEdit::singleline(&mut op.label).desired_width(68.0).hint_text(format!("OP {}", i + 1)));
                        let re_sp = (op.re * 0.01).max(100.0);
                        ui.add(egui::DragValue::new(&mut op.re).speed(re_sp).range(1.001e3..=1e8).max_decimals(0));
                        ui.add(egui::DragValue::new(&mut op.cl).speed(0.01).range(-3.0..=3.5).fixed_decimals(3));
                        ui.horizontal(|ui| {
                            let mut gi = match op.goal {
                                DragGoal::Minimise => 0,
                                DragGoal::AtMostOriginal => 1,
                                DragGoal::AtMost { .. } => 2,
                                DragGoal::Target { .. } => 3,
                            };
                            if widgets::pick(ui, &format!("dgoal{i}"), &mut gi, &["minimise", "<= original", "<= max CD", "target CD"]) {
                                op.goal = match gi {
                                    0 => DragGoal::Minimise,
                                    1 => DragGoal::AtMostOriginal,
                                    2 => DragGoal::AtMost { cd: 0.012 },
                                    _ => DragGoal::Target { cd: 0.008 },
                                };
                            }
                            match &mut op.goal {
                                DragGoal::AtMost { cd } | DragGoal::Target { cd } => {
                                    ui.add(egui::DragValue::new(cd).speed(0.0001).range(0.0005..=0.5).fixed_decimals(5));
                                }
                                _ => {}
                            }
                        });
                        ui.add(egui::DragValue::new(&mut op.weight).speed(0.1).range(0.1..=20.0).fixed_decimals(1));
                        if n > 1 && ui.small_button("x").on_hover_text("remove row").clicked() {
                            remove_row = Some(i);
                        } else {
                            ui.label("");
                        }
                        ui.end_row();
                    }
                });
            });
            if let Some(i) = remove_row {
                self.fui.opt.operating_points.remove(i);
            }
            ui.horizontal(|ui| {
                if ui.small_button("+ Add point").clicked() {
                    let next = self.fui.opt.operating_points.last().cloned().unwrap_or_default();
                    self.fui.opt.operating_points.push(OperatingPoint { label: String::new(), ..next });
                }
                if ui.small_button("Glider defaults").on_hover_text("slow flight (<= orig), cruise (min), fast glide (min)").clicked() {
                    self.fui.opt.operating_points = OptimSettings::default_operating_points();
                    self.fui.opt_3d_origin = None;
                }
            });
        }

        // Minimum thickness at spar
        ui.add_space(4.0);
        let mut spar_on = self.fui.opt.spar.is_some();
        if ui.checkbox(&mut spar_on, RichText::new("Minimum thickness at spar").strong()).on_hover_text("reject candidates thinner than this at the spar chord station before calling XFoil").changed() {
            if spar_on {
                let cur = foil.thickness_at(0.30);
                let min_t = default_spar_min(cur);
                self.fui.opt.spar = Some(SparConstraint { x: 0.30, min_thickness: min_t });
            } else {
                self.fui.opt.spar = None;
            }
        }
        if let Some(sp) = &mut self.fui.opt.spar {
            let cur_t = foil.thickness_at(sp.x);
            egui::Grid::new("opt_spar").num_columns(2).show(ui, |ui| {
                num(ui, "  spar x/c", &mut sp.x, 0.01, 0.05..=0.80, "chordwise spar location");
                let mut min_pct = sp.min_thickness * 100.0;
                if num_u(ui, "  min thickness", &mut min_pct, 0.1, 0.5..=35.0, " %c", "minimum airfoil thickness t/c at spar x/c") {
                    sp.min_thickness = min_pct / 100.0;
                }
            });
            ui.horizontal(|ui| {
                hint(ui, &format!("{} thickness at x/c {:.2}: {:.2}%", foil.name, sp.x, cur_t * 100.0));
                if cur_t > 0.005 && ui.small_button("match").on_hover_text("set minimum thickness to the current foil's thickness at this spar station").clicked() {
                    sp.min_thickness = ((cur_t * 1000.0).floor() / 1000.0).max(0.005);
                }
            });
            if sp.min_thickness > cur_t + 1e-4 {
                ui.label(RichText::new("current foil is thinner than the spar floor (click 'match' or lower min thickness)").small().color(Color32::from_rgb(235, 165, 70)));
            }
        }

        ui.separator();
        let o = &mut self.fui.opt;
        ui.horizontal(|ui| {
            ui.selectable_value(&mut o.algorithm, Algorithm::Mopso, "particle swarm (MOPSO)");
            ui.selectable_value(&mut o.algorithm, Algorithm::Genetic, "genetic");
        });
        ui.horizontal(|ui| {
            ui.label("shape variables");
            ui.selectable_value(&mut o.parameterization, Parameterization::Cst, "CST weights");
            ui.selectable_value(&mut o.parameterization, Parameterization::HicksHenne, "Hicks-Henne bumps");
        });
        egui::Grid::new("opt").num_columns(2).show(ui, |ui| {
            if !self.fui.opt_multipoint {
                num(ui, "alpha (deg)", &mut o.alpha, 0.1, -10.0..=20.0, "operating point angle of attack");
                let sp = o.analysis.re * 0.01;
                num(ui, "Reynolds", &mut o.analysis.re, sp, 1e4..=1e8, "");
            }
            num(ui, "Ncrit", &mut o.analysis.ncrit, 0.1, 0.5..=15.0, "");
            num(ui, "Mach", &mut o.analysis.mach, 0.01, 0.0..=0.8, "");
            int(ui, "population", &mut o.population, 4..=100, "candidates per iteration (8-16 recommended)");
            int(ui, "iterations", &mut o.iterations, 1..=500, "");
            if o.parameterization == Parameterization::Cst {
                int(ui, "CST order", &mut o.cst_order, 1..=12, "order + 1 weight changes on each surface");
                num(ui, "max weight change", &mut o.cst_max_delta, 0.002, 0.002..=0.3, "largest change of a CST weight (0.05 moves a 12% section by roughly 0.5% of chord)");
            } else {
                int(ui, "bumps per side", &mut o.bumps_per_side, 1..=10, "Hicks-Henne functions on each surface");
                num(ui, "max bump (c)", &mut o.max_amplitude, 0.001, 0.0005..=0.05, "maximum bump amplitude, fraction of chord");
            }
            let mut fl = o.flap.is_some();
            if check(ui, "optimise a TE flap angle", &mut fl, "adds the flap deflection as a design variable") {
                o.flap = if fl { Some((0.0, 10.0)) } else { None };
            }
            if let Some((a, b)) = &mut o.flap {
                num(ui, "  flap min (deg)", a, 0.5, -30.0..=30.0, "");
                num(ui, "  flap max (deg)", b, 0.5, -30.0..=30.0, "");
            }
        });
        // Trips, bottom Ncrit and wake are inherited from a 3D station's polar but have
        // no control here; show them whenever they differ from the defaults.
        let d = AnalysisSettings::default();
        let a = &o.analysis;
        let hidden_changed = a.ncrit_bot.is_some()
            || a.xtr_top != d.xtr_top
            || a.xtr_bot != d.xtr_bot
            || a.wake_length != d.wake_length
            || a.wake_dissipation != d.wake_dissipation;
        if hidden_changed {
            let nb = a.ncrit_bot.map(|n| format!("Ncrit bottom {n:.1}, ")).unwrap_or_default();
            let text = format!("inherited from 3D polar: {nb}trips x/c {:.2}/{:.2}, wake {:.2}/{:.2}", a.xtr_top, a.xtr_bot, a.wake_length, a.wake_dissipation);
            ui.horizontal(|ui| {
                ui.label(RichText::new(text).weak());
                if ui.small_button("reset").on_hover_text("restore free transition and default wake").clicked() {
                    o.analysis.ncrit_bot = d.ncrit_bot;
                    o.analysis.xtr_top = d.xtr_top;
                    o.analysis.xtr_bot = d.xtr_bot;
                    o.analysis.wake_length = d.wake_length;
                    o.analysis.wake_dissipation = d.wake_dissipation;
                }
            });
        }

        if !self.fui.opt_multipoint {
            ui.label(RichText::new("objectives").strong());
            for (k, spec) in self.fui.opt.objectives.iter_mut().enumerate() {
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
        }

        let running = self.fui.opt_job.as_ref().and_then(|(id, _)| self.jobs.get(*id)).map(|j| !j.status().is_finished()).unwrap_or(false);
        ui.add_enabled_ui(!running, |ui| {
            if ui.button(RichText::new("Run optimiser").strong()).clicked() {
                let st = self.fui.opt.clone();
                if let Err(e) = st.validate().and_then(|_| st.analysis.validate()) {
                    self.notify(e, true);
                } else {
                    if self.fui.opt_save_name.trim().is_empty() {
                        self.fui.opt_save_name = self.project.unique_foil_name(&format!("{} opt", foil.name));
                    }
                    let enabled_ops = st.enabled_operating_points();
                    let obj_ops: Vec<_> = enabled_ops.iter().filter(|o| !o.goal.is_cap()).collect();
                    let all_target_goals = !obj_ops.is_empty() && obj_ops.iter().all(|o| matches!(o.goal, DragGoal::Target { .. }));
                    // Only attach the 3D origin if this run actually optimises the linked station's foil.
                    let origin = self
                        .fui
                        .opt_3d_origin
                        .clone()
                        .filter(|orig| orig.foil.is_empty() || orig.foil == foil.name);
                    let meta = OptimRunMeta {
                        base_foil: foil.name.clone(),
                        spar: st.spar,
                        origin,
                        all_target_goals,
                    };
                    let k = JobKind::Optimise { foil: foil.clone(), settings: st };
                    let j = self.jobs.submit(format!("optimise {}", foil.name), k);
                    self.fui.opt_job = Some((j.id, meta.clone()));
                    self.fui.opt_result = None;
                    self.fui.opt_sel_pareto = None;
                    self.track_job(j.id, Pending::Optimise { id: j.id, meta: Box::new(meta) });
                }
            }
        });

        let live = self.live_optim();
        if let Some((r, meta)) = live {
            if let Some(cand) = selected_optim_candidate(&r, self.fui.opt_sel_pareto).cloned() {
                ui.separator();
                if cand.op_points.is_empty() {
                    ui.label(format!("best: CL {:.4}  CD {:.5}  CL/CD {:.1}  Cm {:.4}  ({} evals, {} Pareto)", cand.cl, cand.cd, cand.cl / cand.cd.max(1e-12), cand.cm, r.evaluations, r.pareto.len()));
                } else {
                    ui.label(format!("{} evals, {} Pareto (from {}) — see comparison table", r.evaluations, r.pareto.len(), meta.base_foil));
                }
                if self.fui.opt_save_name.trim().is_empty() {
                    self.fui.opt_save_name = self.project.unique_foil_name(&format!("{} opt", meta.base_foil));
                }
                ui.horizontal_wrapped(|ui| {
                    ui.add(egui::TextEdit::singleline(&mut self.fui.opt_save_name).desired_width(130.0).hint_text("new foil name"));
                    if ui.button("Save as new foil").on_hover_text("add the candidate airfoil to the project").clicked() {
                        *save_new_cand = Some((cand.clone(), meta.clone()));
                    }
                    if ui.button("Preview").on_hover_text("show the candidate as the orange preview").clicked() {
                        let mut f = cand.foil.clone();
                        f.note = format!("optimised candidate (CL {:.3} CD {:.5})", cand.cl, cand.cd);
                        *preview = Some(Ok(f));
                    }
                });
                if let Some(origin) = &meta.origin {
                    match self.validate_rerun_origin(origin, &meta.base_foil) {
                        Ok(_) => {
                            if ui
                                .button(RichText::new("Save foil & rerun 3D study").strong())
                                .on_hover_text(format!("save as a new foil, duplicate '{}' with it on '{}', and rerun '{}' alongside the original polar", origin.plane, origin.wing_name, origin.polar_name))
                                .clicked()
                            {
                                *save_and_rerun = Some((cand, meta));
                            }
                        }
                        Err(reason) => {
                            ui.add_enabled(false, egui::Button::new("Save foil & rerun 3D study")).on_disabled_hover_text(reason);
                        }
                    }
                }
            } else if !running {
                ui.separator();
                ui.label(RichText::new("No candidate satisfied all drag caps and spar constraints.").small().color(Color32::from_rgb(235, 90, 75)));
            }
        }
    }

    /// "Use points from 3D analysis" sub-panel inside the optimiser.
    fn picker_3d_ui(&mut self, ui: &mut egui::Ui) {
        let force_open = self.fui.opt_from_3d.force_open.take();
        let polars: Vec<(String, String)> = self
            .project
            .plane_polars
            .iter()
            .filter(|p| p.points.iter().any(|o| !o.strips.is_empty()))
            .map(|p| (p.plane.clone(), p.name.clone()))
            .collect();
        let mut go_3d = false;
        let mut apply_req: Option<Apply3dRequest> = None;

        egui::CollapsingHeader::new(RichText::new("Use points from 3D analysis").strong())
            .default_open(false)
            .open(force_open)
            .show(ui, |ui| {
                if polars.is_empty() {
                    hint(ui, "Run a 3D analysis first (e.g. a T2 fixed-lift or T3 glide polar on the glider) to extract local Re and CL from a wing station.");
                    if ui.button("Go to 3D analysis").clicked() {
                        go_3d = true;
                    }
                    return;
                }
                // Keep selected polar valid
                if self.fui.opt_from_3d.polar.as_ref().is_none_or(|k| !polars.contains(k)) {
                    let prefer = self.pui.sel_polar.clone().filter(|k| polars.contains(k)).unwrap_or_else(|| polars[0].clone());
                    self.fui.opt_from_3d.polar = Some(prefer);
                    self.fui.opt_from_3d.span_pos = None;
                    self.fui.opt_from_3d.selected_points.clear();
                }
                let polar_key = self.fui.opt_from_3d.polar.clone().unwrap();
                let mut polar_changed = false;
                ui.horizontal(|ui| {
                    ui.label("3D polar");
                    egui::ComboBox::from_id_salt("opt3d_pol")
                        .width(210.0)
                        .selected_text(format!("{} / {}", polar_key.0, polar_key.1))
                        .show_ui(ui, |ui| {
                            for k in &polars {
                                if ui.selectable_label(k == &polar_key, format!("{} / {}", k.0, k.1)).clicked() && k != &polar_key {
                                    self.fui.opt_from_3d.polar = Some(k.clone());
                                    polar_changed = true;
                                }
                            }
                        });
                });
                if polar_changed {
                    self.fui.opt_from_3d.span_pos = None;
                    self.fui.opt_from_3d.selected_points.clear();
                }
                let polar_key = self.fui.opt_from_3d.polar.clone().unwrap();
                let Some(polar) = self.project.plane_polars.iter().find(|p| p.plane == polar_key.0 && p.name == polar_key.1) else {
                    return;
                };
                let plane = self.project.planes.iter().find(|pl| pl.name == polar.plane);
                let wings = polar.strip_wings();
                if wings.is_empty() {
                    hint(ui, "This polar has no saved strip distributions.");
                    return;
                }
                if !wings.contains(&self.fui.opt_from_3d.wing) {
                    self.fui.opt_from_3d.wing = wings[0];
                    self.fui.opt_from_3d.span_pos = None;
                }
                if wings.len() > 1 {
                    let cur_w = self.fui.opt_from_3d.wing;
                    let cur_name = polar.wing_name_for_strip(plane, cur_w);
                    ui.horizontal(|ui| {
                        ui.label("wing");
                        egui::ComboBox::from_id_salt("opt3d_wing").selected_text(cur_name).show_ui(ui, |ui| {
                            for &wi in &wings {
                                let wname = polar.wing_name_for_strip(plane, wi);
                                if ui.selectable_value(&mut self.fui.opt_from_3d.wing, wi, wname).changed() {
                                    self.fui.opt_from_3d.span_pos = None;
                                }
                            }
                        });
                    });
                }
                let wing = self.fui.opt_from_3d.wing;
                let stations = polar.wing_stations(plane, wing);
                if stations.is_empty() {
                    return;
                }
                let sel_st_idx = match self.fui.opt_from_3d.span_pos {
                    Some(y) => stations.iter().enumerate().min_by(|a, b| (a.1.span_pos - y).abs().total_cmp(&(b.1.span_pos - y).abs())).map(|(i, _)| i).unwrap_or(0),
                    None => {
                        let idx = stations
                            .iter()
                            .enumerate()
                            .filter(|(_, s)| s.span_pos >= 0.0)
                            .min_by(|a, b| (a.1.span_frac - 0.45).abs().total_cmp(&(b.1.span_frac - 0.45).abs()))
                            .map(|(i, _)| i)
                            .unwrap_or(stations.len() / 2);
                        self.fui.opt_from_3d.span_pos = Some(stations[idx].span_pos);
                        idx
                    }
                };
                let st = &stations[sel_st_idx];
                let actual_y = st.span_pos;
                ui.horizontal(|ui| {
                    ui.label("station");
                    let cur_txt = format!("y {:+.3} m ({:.0}%, c {:.3} m)", st.span_pos, st.span_frac * 100.0, st.chord);
                    egui::ComboBox::from_id_salt("opt3d_st").width(210.0).selected_text(cur_txt).show_ui(ui, |ui| {
                        for s in &stations {
                            let foil_s = s.foil.as_deref().unwrap_or("");
                            let txt = format!("y {:+.3} m ({:.0}%, c {:.3} m) {foil_s}", s.span_pos, s.span_frac * 100.0, s.chord);
                            if ui.selectable_label((s.span_pos - actual_y).abs() < 1e-6, txt).clicked() {
                                self.fui.opt_from_3d.span_pos = Some(s.span_pos);
                            }
                        }
                    });
                });
                if let Some(f) = &st.foil {
                    if let Some((a, b, t)) = &st.foil_blend {
                        hint(ui, &format!("station foil: {:.0}% {a} / {:.0}% {b} (dominant: {f}; rerun replaces sections using {f})", (1.0 - t) * 100.0, t * 100.0));
                    } else {
                        hint(ui, &format!("station foil: {f}"));
                    }
                }
                // Default point selection if empty: slow / cruise / fast
                if self.fui.opt_from_3d.selected_points.is_empty() {
                    for r in polar.pick_flight_regimes(wing, actual_y) {
                        self.fui.opt_from_3d.selected_points.insert(r.point_index);
                    }
                }
                let regime_tags: HashMap<usize, String> = polar
                    .pick_flight_regimes(wing, actual_y)
                    .into_iter()
                    .map(|r| (r.point_index, r.label))
                    .collect();

                ui.horizontal(|ui| {
                    if ui.small_button("Pick slow / cruise / fast").clicked() {
                        self.fui.opt_from_3d.selected_points = regime_tags.keys().copied().collect();
                    }
                    if ui.small_button("Clear").clicked() {
                        self.fui.opt_from_3d.selected_points.clear();
                    }
                });

                egui::ScrollArea::vertical().id_salt("opt3d_pts").max_height(135.0).show(ui, |ui| {
                    egui::Grid::new("opt3d_grid").striped(true).spacing([6.0, 2.0]).show(ui, |ui| {
                        for h in ["", "alpha", "V", "L/D", "local Cl", "local Re", ""] {
                            ui.label(RichText::new(h).small().strong());
                        }
                        ui.end_row();
                        for (idx, opp) in polar.points.iter().enumerate() {
                            let Some(sr) = polar.strip_at(idx, wing, actual_y) else { continue };
                            let mut on = self.fui.opt_from_3d.selected_points.contains(&idx);
                            if ui.checkbox(&mut on, "").changed() {
                                if on {
                                    self.fui.opt_from_3d.selected_points.insert(idx);
                                } else {
                                    self.fui.opt_from_3d.selected_points.remove(&idx);
                                }
                            }
                            ui.label(RichText::new(format!("{:.1}°", opp.alpha)).small());
                            ui.label(RichText::new(format!("{:.1}", opp.velocity)).small());
                            ui.label(RichText::new(format!("{:.1}", opp.ld())).small());
                            ui.label(RichText::new(format!("{:.3}", sr.cl)).small());
                            ui.label(RichText::new(format!("{:.0}k", sr.re / 1e3)).small());
                            let tag_txt = regime_tags.get(&idx).map(|s| s.as_str()).unwrap_or("");
                            if sr.out_of_polar {
                                let env = foiley_plane::foildb::Envelope::from_u8(sr.envelope).label();
                                let combined = if tag_txt.is_empty() { "clamped".into() } else { format!("{tag_txt} (!)") };
                                ui.label(RichText::new(combined).small().color(Color32::from_rgb(235, 165, 70)))
                                    .on_hover_text(format!("2D polar lookup at this station was {env}"));
                            } else if !tag_txt.is_empty() {
                                ui.label(RichText::new(tag_txt).small().color(Color32::from_rgb(120, 200, 255)));
                            } else {
                                ui.label("");
                            }
                            ui.end_row();
                        }
                    });
                });

                let n_sel = self.fui.opt_from_3d.selected_points.len();
                if ui.add_enabled(n_sel > 0, egui::Button::new(format!("Use {n_sel} selected point(s)"))).clicked() {
                    apply_req = Some((polar_key, wing, actual_y, self.fui.opt_from_3d.selected_points.clone()));
                }
            });

        if go_3d {
            self.module = crate::Module::PlaneAnalysis;
        }
        if let Some((pk, wing, y, pts)) = apply_req {
            self.apply_3d_station_points(&pk, wing, y, &pts);
        }
    }

    /// Bottom central panel for optimiser results: original-vs-candidate drag
    /// comparison table (every operating-point row + spar thickness) alongside
    /// the iteration convergence history plot.
    fn optim_results_panel(&mut self, ui: &mut egui::Ui, r: &OptimResult, meta: &OptimRunMeta, w: f32, hh: f32) {
        let is_multi = r.baseline.as_ref().is_some_and(|b| !b.op_points.is_empty())
            || r.best.as_ref().is_some_and(|b| !b.op_points.is_empty());

        if !is_multi {
            ui.allocate_ui_with_layout(egui::vec2(w, hh), egui::Layout::top_down(egui::Align::LEFT), |ui| {
                ui.label(RichText::new("Optimiser: best score per iteration").strong());
                Plot::new("opthist").height((hh - 24.0).max(80.0)).x_axis_label("iteration").show(ui, |pu| {
                    let pts: Vec<[f64; 2]> = r.history.iter().map(|(i, v)| [*i as f64, *v]).collect();
                    pu.line(Line::new("best", PlotPoints::from(pts.clone())).color(series_color(2)));
                    pu.points(Points::new("best", PlotPoints::from(pts)).radius(3.0_f32).color(series_color(2)));
                });
            });
            return;
        }

        let cand = selected_optim_candidate(r, self.fui.opt_sel_pareto);
        let table_w = (w * 0.64).clamp(260.0, (w - 140.0).max(260.0));
        let plot_w = (w - table_w - 10.0).max(120.0);

        ui.allocate_ui_with_layout(egui::vec2(table_w, hh), egui::Layout::top_down(egui::Align::LEFT), |ui| {
            ui.horizontal(|ui| {
                ui.label(RichText::new("Original vs candidate drag").strong());
                if r.pareto.len() > 1 {
                    let cur_label = match self.fui.opt_sel_pareto {
                        None => "best compromise".into(),
                        Some(i) => format!("Pareto #{}", i + 1),
                    };
                    egui::ComboBox::from_id_salt("opt_pareto_pick").width(120.0).selected_text(cur_label).show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.fui.opt_sel_pareto, None, "best compromise");
                        for i in 0..r.pareto.len() {
                            ui.selectable_value(&mut self.fui.opt_sel_pareto, Some(i), format!("Pareto #{}", i + 1));
                        }
                    });
                }
                match cand {
                    Some(c) if c.constraints_met() => {
                        let base_ops = r.baseline.as_ref().map(|b| b.op_points.as_slice()).unwrap_or(&[]);
                        let all_improved_or_met = c.op_points.iter().enumerate().all(|(i, ev)| {
                            let bo = base_ops.get(i);
                            ev.meets_goal(bo).unwrap_or(true)
                                && match ev.goal {
                                    DragGoal::Minimise => bo.filter(|o| o.converged).is_none_or(|o| ev.cd <= o.cd + 1e-7),
                                    _ => true,
                                }
                        });
                        let badge = if all_improved_or_met { "all goals met" } else { "constraints met" };
                        ui.label(RichText::new(badge).small().color(Color32::from_rgb(80, 195, 100)));
                    }
                    Some(_) => {
                        ui.label(RichText::new("drag cap violated").small().color(Color32::from_rgb(235, 130, 60)));
                    }
                    None => {
                        ui.label(RichText::new("no feasible design").small().color(Color32::from_rgb(235, 90, 75)));
                    }
                }
            });
            let scroll_h = ui.available_height().max(80.0);
            egui::ScrollArea::both().id_salt("opt_cmp_scroll").max_height(scroll_h).show(ui, |ui| {
                let base_ops = r.baseline.as_ref().map(|b| b.op_points.as_slice()).unwrap_or(&[]);
                let cand_ops = cand.map(|c| c.op_points.as_slice()).unwrap_or(&[]);
                let nrows = base_ops.len().max(cand_ops.len());
                egui::Grid::new("opt_cmp_grid").striped(true).spacing([10.0, 4.0]).show(ui, |ui| {
                    for h in ["Point", "Re", "CL", "Goal", "Original CD", "Candidate CD", "Delta", "alpha"] {
                        ui.label(RichText::new(h).small().strong());
                    }
                    ui.end_row();
                    for i in 0..nrows {
                        let bo = base_ops.get(i);
                        let co = cand_ops.get(i);
                        let lbl = co.or(bo).map(|e| e.label.as_str()).unwrap_or("-");
                        let re = co.or(bo).map(|e| e.re).unwrap_or(0.0);
                        let cl = co.or(bo).map(|e| e.target_cl).unwrap_or(0.0);
                        let goal = co.or(bo).map(|e| e.goal).unwrap_or(DragGoal::Minimise);
                        ui.label(lbl);
                        ui.monospace(format!("{re:.0}"));
                        ui.monospace(format!("{cl:.3}"));
                        let goal_txt = match goal {
                            DragGoal::Minimise => "minimise".into(),
                            DragGoal::AtMostOriginal => "<= orig".into(),
                            DragGoal::AtMost { cd } => format!("<= {cd:.5}"),
                            DragGoal::Target { cd } => format!("= {cd:.5}"),
                        };
                        ui.label(RichText::new(goal_txt).small());
                        match bo {
                            Some(b) if b.converged => ui.monospace(format!("{:.5}", b.cd)),
                            _ => ui.label("-"),
                        };
                        match co {
                            Some(c) if c.converged => ui.monospace(format!("{:.5}", c.cd)),
                            _ => ui.label("-"),
                        };
                        match (bo, co) {
                            (Some(b), Some(c)) if b.converged && c.converged && b.cd > 1e-9 => {
                                let d = c.cd - b.cd;
                                let pct = d / b.cd * 100.0;
                                let met = c.meets_goal(Some(b)).unwrap_or(true);
                                let col = if !met {
                                    Color32::from_rgb(235, 90, 75)
                                } else if d < -1e-6 {
                                    Color32::from_rgb(80, 195, 100)
                                } else if d > 1e-6 {
                                    Color32::from_rgb(235, 165, 70)
                                } else {
                                    ui.visuals().text_color()
                                };
                                ui.label(RichText::new(format!("{d:+.5} ({pct:+.1}%)")).monospace().color(col));
                            }
                            _ => {
                                ui.label("-");
                            }
                        }
                        match co {
                            Some(c) if c.converged => ui.monospace(format!("{:.2}°", c.alpha)),
                            _ => ui.label("-"),
                        };
                        ui.end_row();
                    }
                });
                if let Some(sp) = meta.spar {
                    let bt = r.baseline.as_ref().and_then(|b| b.spar_thickness);
                    let ct = cand.and_then(|c| c.spar_thickness);
                    match (bt, ct) {
                        (Some(b), Some(c)) => {
                            let ok = c + 1e-6 >= sp.min_thickness;
                            let col = if ok { Color32::from_rgb(80, 195, 100) } else { Color32::from_rgb(235, 90, 75) };
                            ui.add_space(4.0);
                            ui.label(
                                RichText::new(format!(
                                    "Spar thickness at x/c {:.2}: original {:.2}%  ->  candidate {:.2}%  (min {:.2}% {})",
                                    sp.x,
                                    b * 100.0,
                                    c * 100.0,
                                    sp.min_thickness * 100.0,
                                    if ok { "OK" } else { "TOO THIN" }
                                ))
                                .small()
                                .color(col),
                            );
                        }
                        (Some(b), None) => {
                            ui.add_space(4.0);
                            ui.label(
                                RichText::new(format!(
                                    "Spar thickness at x/c {:.2}: original {:.2}% (min {:.2}%)",
                                    sp.x,
                                    b * 100.0,
                                    sp.min_thickness * 100.0,
                                ))
                                .small(),
                            );
                        }
                        _ => {}
                    }
                }
            });
        });

        ui.allocate_ui_with_layout(egui::vec2(plot_w, hh), egui::Layout::top_down(egui::Align::LEFT), |ui| {
            let title = if meta.all_target_goals { "Target CD error (normalised)" } else { "Weighted CD / original" };
            ui.label(RichText::new(title).strong());
            let plot_h = ui.available_height().max(80.0);
            Plot::new("opthist").height(plot_h).x_axis_label("iteration").show(ui, |pu| {
                let pts: Vec<[f64; 2]> = r.history.iter().filter(|(_, v)| v.is_finite()).map(|(i, v)| [*i as f64, *v]).collect();
                if !meta.all_target_goals {
                    pu.hline(egui_plot::HLine::new("original", 1.0).color(Color32::GRAY).style(egui_plot::LineStyle::dashed_dense()));
                }
                pu.line(Line::new("best", PlotPoints::from(pts.clone())).color(series_color(2)));
                pu.points(Points::new("best", PlotPoints::from(pts)).radius(3.0_f32).color(series_color(2)));
            });
        });
    }
}

fn selected_optim_candidate(r: &OptimResult, sel: Option<usize>) -> Option<&Candidate> {
    sel.and_then(|i| r.pareto.get(i)).or(r.best.as_ref()).filter(|c| c.constraints_met())
}

#[cfg(test)]
mod polar_label_regressions {
    use super::*;
    use foiley_app::{JobManager, Project};
    use foiley_plane::analysis::{Method, PlaneOpp, PlanePolar, PlanePolarDef, PolarKind, StripResult};

    #[test]
    fn labels_show_exact_re_and_only_disambiguate_settings_collisions() {
        let a = AnalysisSettings { re: 300_000.0, ..Default::default() };
        let b = AnalysisSettings { re: 300_100.0, ..a.clone() };
        let c = AnalysisSettings { cm_ref: (0.5, 0.0), ..a.clone() };
        let polars: Vec<_> = [a, b, c].into_iter().map(|settings| FoilPolar {
            name: settings.polar_name(), foil: "NACA 0012".into(), re: settings.re, settings: Some(settings), ..Default::default()
        }).collect();
        let labels = polar_labels(&polars);
        assert!(labels[0].name.contains(" / "));
        assert!(!labels[1].name.contains(" / "));
        assert!(labels[2].name.contains(" / "));
        assert!(labels[0].name.contains("Re=300000"));
        assert!(labels[1].name.contains("Re=300100"));
        assert_ne!(labels[0].name, labels[2].name);
        assert_ne!(labels[0].key, labels[2].key);
    }

    #[test]
    fn section_optimiser_loads_3d_points_and_reruns_study_on_plane_copy() {
        let mut project = Project::starter();
        let plane = project.planes[0].clone();
        let make_opp = |alpha: f64, v: f64, cl: f64, cd: f64, local_cl: f64| PlaneOpp {
            alpha,
            velocity: v,
            cl,
            cd,
            strips: vec![
                StripResult { wing: 0, span_pos: 0.72, chord: 0.247, cl: local_cl, re: v * 0.247 / 1.5e-5, cdv: 0.009, alpha_eff: alpha, ..Default::default() },
            ],
            ..Default::default()
        };
        let polar = PlanePolar {
            name: "T2-Vlm2".into(),
            plane: plane.name.clone(),
            def: Some(PlanePolarDef { kind: PolarKind::FixedLift, method: Method::Vlm2, viscous: false, ..Default::default() }),
            points: vec![
                make_opp(0.0, 14.5, 0.16, 0.013, 0.18),
                make_opp(4.0, 8.2, 0.52, 0.019, 0.58),
                make_opp(8.0, 6.4, 0.88, 0.038, 0.96),
            ],
            ..Default::default()
        };
        project.add_plane_polar(polar);

        let mut app = App::from_project(project, JobManager::new(1), true);
        app.open_section_optimiser((plane.name.clone(), "T2-Vlm2".into()), 0, 0.72);

        assert_eq!(app.module, crate::Module::FoilDesign);
        assert_eq!(app.sel_foil.as_deref(), Some("NACA 2412"));
        assert!(app.fui.opt_multipoint);
        assert_eq!(app.fui.opt.operating_points.len(), 3);
        assert_eq!(app.fui.opt.operating_points[0].label, "Slow flight");
        assert_eq!(app.fui.opt.operating_points[0].goal, DragGoal::AtMostOriginal);
        assert_eq!(app.fui.opt.operating_points[1].label, "Cruise");
        assert_eq!(app.fui.opt.operating_points[1].goal, DragGoal::Minimise);
        assert_eq!(app.fui.opt.operating_points[2].label, "Fast glide");
        assert!(app.fui.opt.spar.is_some());
        assert!(app.fui.opt_3d_origin.is_some());

        // Saving and rerunning creates a new foil + plane copy without mutating the original plane
        let cand = Candidate {
            genes: vec![],
            foil: app.project.foil("NACA 2412").unwrap().clone(),
            cl: 0.58,
            cd: 0.008,
            cm: -0.05,
            cpmin: -1.1,
            objectives: vec![0.008],
            scores: vec![0.008],
            feasible: true,
            op_points: vec![],
            spar_thickness: Some(0.118),
            violation: 0.0,
        };
        let meta = OptimRunMeta {
            base_foil: "NACA 2412".into(),
            spar: app.fui.opt.spar,
            origin: app.fui.opt_3d_origin.clone(),
            all_target_goals: false,
        };
        app.save_and_rerun_3d_study(&cand, &meta);
        assert_eq!(app.module, crate::Module::PlaneAnalysis);
        assert!(app.project.foil("NACA 2412 opt").is_some());
        assert_eq!(app.project.planes.len(), 2);
        assert_eq!(app.project.planes[0].wings[0].sections[0].foil_right, "NACA 2412");
        assert_eq!(app.project.planes[1].wings[0].sections[0].foil_right, "NACA 2412 opt");
        assert!(app.pui.visible.contains(&(plane.name.clone(), "T2-Vlm2".into())));
        assert!(!app.pending.is_empty(), "a 3D rerun job should have been queued");
        app.jobs.cancel_all();
    }

    #[test]
    fn thin_foil_spar_default_and_stale_candidate_and_transition_settings() {
        // Thin foil (NACA 0003, ~3% t/c) must never panic when computing default spar min
        let thin = Foil::naca("0003").unwrap();
        let t = thin.thickness_at(0.30);
        assert!(t > 0.02 && t < 0.04);
        let d = default_spar_min(t);
        assert!(d >= 0.005 && d <= t, "default_spar_min({t}) = {d}");

        // Setup project with Ncrit=5 polar family AND Ncrit=9 polar family on NACA 2412
        let mut project = Project::starter();
        let mut fp9 = FoilPolar { foil: "NACA 2412".into(), name: "T1-Re200k-N9".into(), re: 2.0e5, ncrit: 9.0, viscous: true, polar_type: 1, ..Default::default() };
        fp9.points.push(FoilPoint { alpha: 0.0, cl: 0.2, cd: 0.01, ..Default::default() });
        fp9.points.push(FoilPoint { alpha: 4.0, cl: 0.6, cd: 0.012, ..Default::default() });
        let mut fp5 = FoilPolar { foil: "NACA 2412".into(), name: "T1-Re200k-N5".into(), re: 2.0e5, ncrit: 5.0, viscous: true, polar_type: 1, ..Default::default() };
        fp5.points.push(FoilPoint { alpha: 0.0, cl: 0.2, cd: 0.011, ..Default::default() });
        fp5.points.push(FoilPoint { alpha: 4.0, cl: 0.6, cd: 0.014, ..Default::default() });
        project.add_foil_polar(fp9);
        project.add_foil_polar(fp5);

        let plane = project.planes[0].clone();
        let make_opp = |alpha: f64, v: f64, cl: f64, cd: f64, local_cl: f64| PlaneOpp {
            alpha,
            velocity: v,
            cl,
            cd,
            strips: vec![
                StripResult { wing: 0, span_pos: 0.72, chord: 0.247, cl: local_cl, re: v * 0.247 / 1.5e-5, cdv: 0.009, alpha_eff: alpha, ..Default::default() },
                StripResult { wing: 1, span_pos: 0.15, chord: 0.11, cl: 0.05, re: v * 0.11 / 1.5e-5, cdv: 0.011, alpha_eff: alpha, ..Default::default() },
            ],
            ..Default::default()
        };
        let mut def = PlanePolarDef { kind: PolarKind::FixedLift, method: Method::Vlm2, viscous: true, ..Default::default() };
        def.family.ncrit = Some(5.0);
        let polar = PlanePolar {
            name: "T2-N5".into(),
            plane: plane.name.clone(),
            def: Some(def),
            points: vec![
                make_opp(0.0, 14.5, 0.16, 0.013, 0.18),
                make_opp(4.0, 8.2, 0.52, 0.019, 0.58),
                make_opp(8.0, 6.4, 0.88, 0.038, 0.96),
            ],
            ..Default::default()
        };
        project.add_plane_polar(polar);

        let mut app = App::from_project(project, JobManager::new(1), true);
        app.open_section_optimiser((plane.name.clone(), "T2-N5".into()), 0, 0.72);
        // Must pick Ncrit=5 from the resolved family, NOT Ncrit=9 from the first polar!
        assert!((app.fui.opt.analysis.ncrit - 5.0).abs() < 1e-9, "got ncrit {}", app.fui.opt.analysis.ncrit);
        assert_eq!(app.fui.opt_3d_origin.as_ref().unwrap().source_analysis.as_ref().unwrap().ncrit, 5.0);

        // Simulate a finished main-wing job stored in opt_job + opt_result, then switch to elevator (wing 1)
        let dummy_job = app.jobs.submit("dummy", JobKind::FoilPoint { foil: Foil::naca("2412").unwrap(), settings: AnalysisSettings::default(), value: 2.0, is_alpha: true });
        let main_meta = OptimRunMeta {
            base_foil: "NACA 2412".into(),
            spar: app.fui.opt.spar,
            origin: app.fui.opt_3d_origin.clone(),
            all_target_goals: false,
        };
        app.fui.opt_job = Some((dummy_job.id, main_meta.clone()));
        app.fui.opt_result = Some((OptimResult::default(), main_meta));
        let elev_pts: BTreeSet<usize> = [0, 1, 2].into_iter().collect();
        app.apply_3d_station_points(&(plane.name.clone(), "T2-N5".into()), 1, 0.15, &elev_pts);
        // Applying elevator points MUST detach the earlier main-wing job and result!
        assert!(app.fui.opt_job.is_none());
        assert!(app.fui.opt_result.is_none());
        assert!(app.live_optim().is_none());
        assert_eq!(app.sel_foil.as_deref(), Some("NACA 0009"));

        // Foil rename updates opt_3d_origin and validation rejects stale origin after foil removal
        app.rename_foil("NACA 0009", "Tail Foil");
        assert_eq!(app.fui.opt_3d_origin.as_ref().unwrap().foil, "Tail Foil");
        let origin_now = app.fui.opt_3d_origin.clone().unwrap();
        assert!(app.validate_rerun_origin(&origin_now, "Tail Foil").is_ok());
        // Removing all matching sections from the elevator makes rerun validation fail cleanly
        for sec in &mut app.project.planes[0].wings[1].sections {
            sec.foil_right = "NACA 2412".into();
            sec.foil_left = "NACA 2412".into();
        }
        assert!(app.validate_rerun_origin(&origin_now, "Tail Foil").is_err());
        app.jobs.cancel_all();
    }
}

