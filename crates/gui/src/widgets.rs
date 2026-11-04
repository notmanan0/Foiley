//! Small UI helpers: labelled numeric fields with help text, colour maps,
//! plot-variable pickers.

use eframe::egui::{self, Color32, RichText, Ui};

/// A labelled numeric field in a two-column grid, with hover help.
pub fn num(ui: &mut Ui, label: &str, v: &mut f64, speed: f64, range: std::ops::RangeInclusive<f64>, help: &str) -> bool {
    let l = ui.label(label);
    if !help.is_empty() {
        l.on_hover_text(help);
    }
    let r = ui.add(egui::DragValue::new(v).speed(speed).range(range).max_decimals(6));
    let r = if help.is_empty() { r } else { r.on_hover_text(help) };
    ui.end_row();
    r.changed()
}

/// Same with a unit suffix.
pub fn num_u(ui: &mut Ui, label: &str, v: &mut f64, speed: f64, range: std::ops::RangeInclusive<f64>, unit: &str, help: &str) -> bool {
    let l = ui.label(label);
    if !help.is_empty() {
        l.on_hover_text(help);
    }
    let r = ui.add(egui::DragValue::new(v).speed(speed).range(range).suffix(unit).max_decimals(6));
    let r = if help.is_empty() { r } else { r.on_hover_text(help) };
    ui.end_row();
    r.changed()
}

pub fn int(ui: &mut Ui, label: &str, v: &mut usize, range: std::ops::RangeInclusive<usize>, help: &str) -> bool {
    let l = ui.label(label);
    if !help.is_empty() {
        l.on_hover_text(help);
    }
    let r = ui.add(egui::DragValue::new(v).range(range));
    let r = if help.is_empty() { r } else { r.on_hover_text(help) };
    ui.end_row();
    r.changed()
}

/// `check`, greyed out when `on` is false.
pub fn check_on(ui: &mut Ui, on: bool, label: &str, v: &mut bool, help: &str) -> bool {
    ui.label("");
    let r = ui.add_enabled(on, egui::Checkbox::new(v, label)).on_hover_text(help).on_disabled_hover_text("not supported by this method / polar type");
    ui.end_row();
    r.changed()
}

/// `num_u`, greyed out when `on` is false.
#[allow(clippy::too_many_arguments)]
pub fn num_u_on(ui: &mut Ui, on: bool, label: &str, v: &mut f64, speed: f64, range: std::ops::RangeInclusive<f64>, unit: &str, help: &str) -> bool {
    let l = ui.add_enabled(on, egui::Label::new(label));
    if !help.is_empty() {
        l.on_hover_text(help);
    }
    let r = ui.add_enabled(on, egui::DragValue::new(v).speed(speed).range(range).suffix(unit).max_decimals(6)).on_disabled_hover_text("not supported by this method / polar type");
    ui.end_row();
    r.changed()
}

/// `pick` with some entries greyed out (`enabled[i] == false`).
pub fn pick_enabled(ui: &mut Ui, id: &str, current: &mut usize, labels: &[&str], enabled: &[bool]) -> bool {
    let mut changed = false;
    egui::ComboBox::from_id_salt(id).selected_text(labels.get(*current).copied().unwrap_or("?")).show_ui(ui, |ui| {
        for (i, l) in labels.iter().enumerate() {
            let on = enabled.get(i).copied().unwrap_or(true);
            if ui.add_enabled(on, egui::Button::selectable(*current == i, *l)).clicked() && *current != i {
                *current = i;
                changed = true;
            }
        }
    });
    changed
}

pub fn check(ui: &mut Ui, label: &str, v: &mut bool, help: &str) -> bool {
    ui.label("");
    let r = ui.checkbox(v, label).on_hover_text(help);
    ui.end_row();
    r.changed()
}

/// A short explanation in a muted colour.
pub fn hint(ui: &mut Ui, text: &str) {
    ui.label(RichText::new(text).small().color(ui.visuals().weak_text_color()));
}

pub fn heading(ui: &mut Ui, text: &str) {
    ui.add_space(4.0);
    ui.label(RichText::new(text).strong().size(15.0));
}

/// Collapsible section with a title.
pub fn section(ui: &mut Ui, title: &str, open: bool, add: impl FnOnce(&mut Ui)) {
    egui::CollapsingHeader::new(RichText::new(title).strong()).default_open(open).show(ui, add);
}

/// Perceptually ordered colour map (blue - cyan - green - yellow - red), t in 0..1.
pub fn colormap(t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    let stops: [(f32, [f32; 3]); 5] = [(0.0, [0.10, 0.20, 0.80]), (0.25, [0.10, 0.70, 0.90]), (0.5, [0.20, 0.80, 0.30]), (0.75, [0.95, 0.85, 0.15]), (1.0, [0.90, 0.15, 0.10])];
    for w in stops.windows(2) {
        let (a, ca) = w[0];
        let (b, cb) = w[1];
        if t <= b {
            let u = (t - a) / (b - a);
            let c: Vec<f32> = (0..3).map(|i| ca[i] + u * (cb[i] - ca[i])).collect();
            return Color32::from_rgb((c[0] * 255.0) as u8, (c[1] * 255.0) as u8, (c[2] * 255.0) as u8);
        }
    }
    Color32::RED
}

/// Distinct colours for overlaid curves.
pub fn series_color(i: usize) -> Color32 {
    const C: [Color32; 10] = [
        Color32::from_rgb(31, 119, 180),
        Color32::from_rgb(255, 127, 14),
        Color32::from_rgb(44, 160, 44),
        Color32::from_rgb(214, 39, 40),
        Color32::from_rgb(148, 103, 189),
        Color32::from_rgb(140, 86, 75),
        Color32::from_rgb(227, 119, 194),
        Color32::from_rgb(127, 127, 127),
        Color32::from_rgb(188, 189, 34),
        Color32::from_rgb(23, 190, 207),
    ];
    C[i % C.len()]
}

/// Combo box over a list of labels; returns true when changed.
pub fn pick(ui: &mut Ui, id: &str, current: &mut usize, labels: &[&str]) -> bool {
    let mut changed = false;
    egui::ComboBox::from_id_salt(id).selected_text(labels.get(*current).copied().unwrap_or("?")).show_ui(ui, |ui| {
        for (i, l) in labels.iter().enumerate() {
            if ui.selectable_value(current, i, *l).changed() {
                changed = true;
            }
        }
    });
    changed
}

pub fn fmt_secs(d: std::time::Duration) -> String {
    let s = d.as_secs_f64();
    if s < 60.0 {
        format!("{s:.1} s")
    } else {
        format!("{}:{:02} min", (s / 60.0) as u64, s as u64 % 60)
    }
}
