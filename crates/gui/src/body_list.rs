//! The plane editor's body (fuselage) list: select, add, duplicate, remove,
//! and put a CAD surface into one body.  Plain state on `Plane::fuses` and
//! a selected index, with no egui, so `crates/gui/tests/f004b.rs` drives it
//! without rendering.

use foiley_plane::geom::{Fuse, Plane};

/// Lateral spacing of a new body from the previous one, m (a new default pod
/// placed on top of an existing body would interpenetrate it).
pub const NEW_BODY_SPACING: f64 = 0.3;

#[derive(Clone, Debug, PartialEq)]
pub enum BodyAction {
    Select(usize),
    /// a default pod, beside the last body
    AddPod,
    /// a copy of the selected body, beside the last body
    Duplicate,
    /// remove the selected body
    Remove,
}

/// The selected body index, clamped to the list (None when there is no body).
pub fn selected(p: &Plane, sel: &mut usize) -> Option<usize> {
    if p.fuses.is_empty() {
        *sel = 0;
        return None;
    }
    *sel = (*sel).min(p.fuses.len() - 1);
    Some(*sel)
}

/// The selected body, for editing.
pub fn selected_mut<'a>(p: &'a mut Plane, sel: &mut usize) -> Option<&'a mut Fuse> {
    let i = selected(p, sel)?;
    p.fuses.get_mut(i)
}

/// A body name not used yet ("Fuselage", "Fuselage 2", ...).
pub fn unique_name(p: &Plane, base: &str) -> String {
    let taken = |n: &str| p.fuses.iter().any(|f| f.name == n);
    if !taken(base) {
        return base.to_string();
    }
    (2..).map(|k| format!("{base} {k}")).find(|n| !taken(n)).unwrap()
}

/// Label of body `i` in the list.
pub fn label(p: &Plane, i: usize) -> String {
    p.fuses.get(i).map(|f| format!("{}. {}", i + 1, f.name)).unwrap_or_default()
}

fn beside_last(p: &Plane) -> f64 {
    p.fuses.iter().map(|f| f.position.y).fold(f64::NEG_INFINITY, f64::max).max(-NEW_BODY_SPACING) + NEW_BODY_SPACING
}

pub fn apply(p: &mut Plane, sel: &mut usize, a: BodyAction) {
    match a {
        BodyAction::Select(i) => {
            *sel = i;
            selected(p, sel);
        }
        BodyAction::AddPod => {
            let mut f = Fuse::default_pod();
            if !p.fuses.is_empty() {
                f.position.y = beside_last(p);
            }
            f.name = unique_name(p, &f.name);
            p.fuses.push(f);
            *sel = p.fuses.len() - 1;
        }
        BodyAction::Duplicate => {
            if let Some(i) = selected(p, sel) {
                let mut f = p.fuses[i].clone();
                f.position.y = beside_last(p);
                f.name = unique_name(p, &format!("{} copy", f.name));
                p.fuses.push(f);
                *sel = p.fuses.len() - 1;
            }
        }
        BodyAction::Remove => {
            if let Some(i) = selected(p, sel) {
                p.fuses.remove(i);
                selected(p, sel);
            }
        }
    }
}
