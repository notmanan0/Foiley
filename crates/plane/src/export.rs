//! Exchange: AVL geometry (.avl) export and import, AVL mass (.mass) files,
//! STL meshes, and the conversion report every import returns.

use crate::analysis::ControlTarget;
use crate::foildb::FoilDb;
use crate::geom::{Plane, Wing, WingKind, WingSection};
use crate::mesh::{build, FlapInfo, MeshOptions};
use crate::PlaneError;
use foiley_core::{Spacing, Vec3, DEG};
use foiley_foil::Foil;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

fn naca_digits(name: &str) -> Option<String> {
    let l = name.to_lowercase().replace([' ', '-', '_'], "");
    let d = l.strip_prefix("naca")?;
    (d.len() == 4 && d.chars().all(|c| c.is_ascii_digit())).then(|| d.to_string())
}

/// Shortest text that reads back to the same f64 (no exponent, no "-0").
fn num(v: f64) -> String {
    if v == 0.0 {
        "0".into()
    } else {
        format!("{v}")
    }
}

/// AVL spacing parameter (0 equal, 1 cosine, 2 sine = bunched at the first
/// section, -2 = bunched at the last).  `reversed`: the sections are listed
/// tip to root (explicit left halves).
fn avl_space(s: Spacing, reversed: bool) -> f64 {
    match s {
        Spacing::Uniform => 0.0,
        Spacing::Cosine | Spacing::Tanh => 1.0,
        Spacing::Sine => {
            if reversed {
                2.0
            } else {
                -2.0
            }
        }
        Spacing::MinusSine => {
            if reversed {
                -2.0
            } else {
                2.0
            }
        }
    }
}

fn spacing_from_avl(v: f64) -> (Spacing, bool) {
    let r = v.round().clamp(-3.0, 3.0);
    let s = match r as i64 {
        1 | -1 => Spacing::Cosine,
        2 => Spacing::MinusSine,
        -2 => Spacing::Sine,
        _ => Spacing::Uniform,
    };
    (s, (v - r).abs() < 1e-9)
}

fn sanitize(s: &str) -> String {
    let t: String = s.chars().map(|c| if c.is_ascii_alphanumeric() { c } else { '_' }).collect();
    let t = t.trim_matches('_').to_string();
    if t.is_empty() {
        "wing".into()
    } else {
        t
    }
}

/// One AVL `CONTROL` line on the sections bounding wing surface `surf`.
struct AvlControl {
    surf: usize,
    name: String,
    gain: f64,
    xhinge: f64,
    sgn_dup: f64,
}

fn flap_gain(gains: &[(ControlTarget, f64)], wing: usize, index: usize) -> f64 {
    gains.iter().find(|(t, _)| *t == ControlTarget::Flap { wing, index }).map(|g| g.1).unwrap_or(1.0)
}

fn control_word(kind: WingKind, aileron: bool) -> &'static str {
    if aileron {
        return "aileron";
    }
    match kind {
        WingKind::Elevator => "elevator",
        WingKind::Fin => "rudder",
        _ => "flap",
    }
}

fn unique(used: &mut HashSet<String>, name: String, iw: usize) -> String {
    let mut n = name.clone();
    let mut k = 0;
    while used.contains(&n) {
        k += 1;
        n = if k == 1 { format!("{name}_w{}", iw + 1) } else { format!("{name}_w{}_{k}", iw + 1) };
    }
    used.insert(n.clone());
    n
}

/// How one wing is written: one surface (mirrored with `YDUPLICATE 0.0`, or
/// one-sided) or an explicit right and left surface.
struct WingPlan {
    mirrored: bool,
    right: Vec<AvlControl>,
    left: Vec<AvlControl>,
}

fn partner<'a>(f: &FlapInfo, other: &[&'a FlapInfo]) -> Option<&'a FlapInfo> {
    other.iter().copied().find(|o| o.surf == f.surf)
}

fn plan_wing(iw: usize, w: &Wing, flaps: &[&FlapInfo], gains: &[(ControlTarget, f64)], used: &mut HashSet<String>) -> WingPlan {
    let side = |left: bool| {
        let mut v: Vec<&FlapInfo> = flaps.iter().copied().filter(|f| f.left == left).collect();
        v.sort_by_key(|f| f.surf);
        v
    };
    let (rf, lf) = (side(false), side(true));
    let g = |f: &FlapInfo| flap_gain(gains, iw, f.index);
    let antisym = |a: &FlapInfo, b: Option<&FlapInfo>| b.map(|b| g(a) != 0.0 && g(b) == -g(a)).unwrap_or(false);
    let n = w.sections.len();
    let foils_mirror = (0..n).all(|i| w.foil(i, true) == w.foil(i, false));
    let controls_mirror = rf.len() == lf.len()
        && rf.iter().all(|r| match partner(r, &lf) {
            Some(l) => g(l).abs() == g(r).abs() && (l.xhinge - r.xhinge).abs() < 1e-12,
            None => false,
        });
    let mirrored = w.two_sided && w.position.y.abs() < 1e-12 && w.rx.abs() < 1e-12 && foils_mirror && controls_mirror;
    let base = sanitize(&w.name);
    let mut plan = WingPlan { mirrored, right: vec![], left: vec![] };
    for (k, r) in rf.iter().enumerate() {
        let l = partner(r, &lf);
        let aileron = antisym(r, l);
        let word = control_word(w.kind, aileron);
        let suffix = if mirrored || !w.two_sided { "" } else { "_R" };
        let sgn = if mirrored && aileron { -1.0 } else { 1.0 };
        let name = unique(used, format!("{base}_{word}{}{suffix}", k + 1), iw);
        plan.right.push(AvlControl { surf: r.surf, name, gain: g(r), xhinge: r.xhinge, sgn_dup: sgn });
    }
    if !mirrored {
        for (k, l) in lf.iter().enumerate() {
            let aileron = antisym(l, partner(l, &rf));
            let name = unique(used, format!("{base}_{}{}_L", control_word(w.kind, aileron), k + 1), iw);
            plan.left.push(AvlControl { surf: l.surf, name, gain: g(l), xhinge: l.xhinge, sgn_dup: 1.0 });
        }
    }
    plan
}

fn write_surface(o: &mut String, title: &str, w: &Wing, left: bool, ydup: bool, ctrls: &[AvlControl]) {
    let nx = w.sections.first().map(|s| s.nx).unwrap_or(8);
    let cspace = avl_space(w.sections.first().map(|s| s.xdist).unwrap_or_default(), false);
    o.push_str(&format!("#{}\nSURFACE\n{title}\n#Nchord Cspace\n{nx} {}\n", "=".repeat(60), num(cspace)));
    if ydup {
        o.push_str("YDUPLICATE\n0.0\n");
    }
    if w.kind == WingKind::Fin && !w.two_sided {
        o.push_str("#vertical surface\n");
    }
    let fr = w.frames(left);
    let n = w.sections.len();
    let order: Vec<usize> = if left { (0..n).rev().collect() } else { (0..n).collect() };
    for &k in &order {
        let sec = &w.sections[k];
        let le = fr[k].le();
        // strips from this section to the next one listed
        let seg = if left { &w.sections[k.saturating_sub(1)] } else { sec };
        o.push_str(&format!(
            "SECTION\n#Xle Yle Zle Chord Ainc Nspan Sspace\n{} {} {} {} {} {} {}\n",
            num(le.x),
            num(le.y),
            num(le.z),
            num(sec.chord),
            num(sec.twist + w.ry),
            seg.ny,
            num(avl_space(seg.ydist, left))
        ));
        let foil = w.foil(k, left);
        match naca_digits(foil) {
            Some(d) => o.push_str(&format!("NACA\n{d}\n")),
            None => o.push_str(&format!("AFILE\n{}.dat\n", foil.replace(' ', "_"))),
        }
        for c in ctrls.iter().filter(|c| c.surf == k || c.surf + 1 == k) {
            o.push_str(&format!("CONTROL\n#name gain Xhinge XYZhvec SgnDup\n{} {} {} 0 0 0 {}\n", c.name, num(c.gain), num(c.xhinge), num(c.sgn_dup)));
        }
    }
}

/// AVL geometry file, every flap with unit gain.  See [`avl_with_gains`].
pub fn avl(plane: &Plane, db: &FoilDb, mach: f64) -> String {
    avl_with_gains(plane, db, mach, &[])
}

/// AVL geometry file.  Foils that are not NACA 4-digit are referenced as
/// `AFILE <name>.dat`; write those files next to the .avl (see
/// [`avl_foil_files`]).
///
/// * A wing is written as one surface with `YDUPLICATE 0.0` only when it is
///   two-sided, centred, unrotated about x, and mirror-identical in foils,
///   twist and controls.  Otherwise each side is its own `SURFACE` (the left
///   one listed tip to root), so offset two-sided wings always get an
///   explicit left half.
/// * Every flap (a wing surface whose two foils have hinges) is one
///   `CONTROL`, with a unique name.  `gains` is the definition's mixing (T7
///   gains, degrees per unit control); flaps it does not list get gain 1.
///   A mirrored pair with opposite gains is an aileron (`SgnDup -1`); other
///   pairs (flaps, elevators) get `SgnDup +1`.
/// * Numbers are written in full (shortest round-trip form).
pub fn avl_with_gains(plane: &Plane, db: &FoilDb, mach: f64, gains: &[(ControlTarget, f64)]) -> String {
    let (s, b, c) = plane.reference(true);
    let i = plane.inertia(db);
    let mut o = String::new();
    o.push_str(&format!("{}\n#Mach\n{mach:.3}\n#IYsym IZsym Zsym\n0 0 0.0\n#Sref Cref Bref\n{s:.6} {c:.6} {b:.6}\n#Xref Yref Zref\n{:.6} {:.6} {:.6}\n#CDp\n0.0\n", plane.name, i.cog.x, i.cog.y, i.cog.z));
    let flaps = build(plane, db, &MeshOptions { allow_missing_foils: true, ..Default::default() }).map(|m| m.flaps).unwrap_or_default();
    let mut used = HashSet::new();
    for (iw, w) in plane.wings.iter().enumerate() {
        let wf: Vec<&FlapInfo> = flaps.iter().filter(|f| f.wing == iw).collect();
        let plan = plan_wing(iw, w, &wf, gains, &mut used);
        let title = format!("{} ({})", w.name, w.kind.label());
        if plan.mirrored || !w.two_sided {
            write_surface(&mut o, &title, w, false, plan.mirrored, &plan.right);
        } else {
            write_surface(&mut o, &format!("{title} right"), w, false, false, &plan.right);
            write_surface(&mut o, &format!("{title} left"), w, true, false, &plan.left);
        }
    }
    o
}

/// Foil files referenced by `AFILE` entries: (file name, .dat text).
pub fn avl_foil_files(plane: &Plane, foils: &[foiley_foil::Foil]) -> Vec<(String, String)> {
    let mut out = Vec::new();
    for n in plane.foil_names() {
        if naca_digits(&n).is_some() {
            continue;
        }
        if let Some(f) = foils.iter().find(|f| f.name == n) {
            out.push((format!("{}.dat", n.replace(' ', "_")), f.to_dat()));
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Conversion report
// ---------------------------------------------------------------------------

/// One item of a conversion report and why it was approximated or skipped.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct ReportEntry {
    pub item: String,
    pub reason: String,
}

/// What an import (AVL, xflr5 `.xfl`, flow5 `.fl5`) mapped, approximated
/// and skipped, with reasons.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct ConversionReport {
    /// the file or format the report is about
    pub source: String,
    pub mapped: Vec<String>,
    pub approximated: Vec<ReportEntry>,
    pub skipped: Vec<ReportEntry>,
}

impl ConversionReport {
    pub fn new(source: impl Into<String>) -> Self {
        ConversionReport { source: source.into(), ..Default::default() }
    }
    pub fn map(&mut self, item: impl Into<String>) {
        self.mapped.push(item.into());
    }
    pub fn approx(&mut self, item: impl Into<String>, reason: impl Into<String>) {
        self.approximated.push(ReportEntry { item: item.into(), reason: reason.into() });
    }
    pub fn skip(&mut self, item: impl Into<String>, reason: impl Into<String>) {
        self.skipped.push(ReportEntry { item: item.into(), reason: reason.into() });
    }
    /// Append another report's entries.
    pub fn merge(&mut self, o: ConversionReport) {
        self.mapped.extend(o.mapped);
        self.approximated.extend(o.approximated);
        self.skipped.extend(o.skipped);
    }
    /// Nothing approximated or skipped.
    pub fn is_lossless(&self) -> bool {
        self.approximated.is_empty() && self.skipped.is_empty()
    }
    /// Plain-text report for the CLI and the GUI.
    pub fn text(&self) -> String {
        let mut o = format!("Conversion report: {}\n  mapped: {}, approximated: {}, skipped: {}\n", self.source, self.mapped.len(), self.approximated.len(), self.skipped.len());
        for m in &self.mapped {
            o.push_str(&format!("  mapped       {m}\n"));
        }
        for e in &self.approximated {
            o.push_str(&format!("  approximated {}: {}\n", e.item, e.reason));
        }
        for e in &self.skipped {
            o.push_str(&format!("  skipped      {}: {}\n", e.item, e.reason));
        }
        o
    }
}

// ---------------------------------------------------------------------------
// AVL import
// ---------------------------------------------------------------------------

/// An imported AVL geometry.
#[derive(Clone, Debug)]
pub struct AvlImport {
    pub plane: Plane,
    /// NACA foils the file names, and copies of `AFILE` foils (matched by
    /// name against the foils passed to [`import_avl`]) that needed a hinge,
    /// hinged where a `CONTROL` needs it; add them to the project.
    pub foils: Vec<Foil>,
    /// controls as flap mixing: (AVL control name, flap, gain in degrees per unit)
    pub mixing: Vec<(String, ControlTarget, f64)>,
    pub report: ConversionReport,
}

#[derive(Clone, Debug)]
struct ACtrl {
    name: String,
    gain: f64,
    xhinge: f64,
    sgn_dup: f64,
}

#[derive(Clone, Debug)]
struct ASec {
    le: Vec3,
    chord: f64,
    ainc: f64,
    nspan: Option<usize>,
    sspace: Option<f64>,
    /// foil name, NACA digits when `naca`
    foil: String,
    naca: bool,
    ctrls: Vec<ACtrl>,
}

#[derive(Clone, Debug)]
struct ASurf {
    title: String,
    nchord: usize,
    cspace: f64,
    nspan: Option<usize>,
    sspace: Option<f64>,
    ydup: Option<f64>,
    scale: Vec3,
    translate: Vec3,
    angle: f64,
    secs: Vec<ASec>,
}

fn nums(line: &str) -> Vec<f64> {
    line.split(|c: char| c.is_whitespace() || c == ',').filter(|t| !t.is_empty()).map_while(|t| t.replace(['D', 'd'], "e").parse::<f64>().ok()).collect()
}

fn strip_comment(l: &str) -> &str {
    let cut = l.find(['!', '#']).unwrap_or(l.len());
    l[..cut].trim()
}

const KEYWORDS: [&str; 19] = ["SURF", "BODY", "COMP", "INDE", "YDUP", "SCAL", "TRAN", "ANGL", "NOWA", "NOAL", "NOLO", "SECT", "NACA", "AIRF", "AFIL", "BFIL", "CLAF", "CDCL", "CONT"];

fn keyword(l: &str) -> Option<String> {
    let w = l.split_whitespace().next()?.to_uppercase();
    let k: String = w.chars().take(4).collect();
    (KEYWORDS.contains(&k.as_str()) || k == "DESI").then_some(k)
}

/// Read an AVL geometry file.  `foils` are the project's foils: `AFILE`
/// references are matched against their names (spaces written as `_`).
pub fn import_avl(text: &str, foils: &[Foil]) -> Result<AvlImport, String> {
    let mut rep = ConversionReport::new("AVL geometry");
    // logical lines: (comment-stripped, raw)
    let lines: Vec<(&str, &str)> = text.lines().map(|l| (strip_comment(l), l.trim())).filter(|(s, _)| !s.is_empty()).collect();
    let mut it = 0usize;
    let title = lines.first().map(|l| l.1.to_string()).ok_or("empty AVL file")?;
    it += 1;
    let mut header: Vec<Vec<f64>> = Vec::new();
    while it < lines.len() && keyword(lines[it].0).is_none() {
        header.push(nums(lines[it].0));
        it += 1;
    }
    if header.len() < 4 {
        return Err("AVL header incomplete: expected Mach, IYsym IZsym Zsym, Sref Cref Bref and Xref Yref Zref lines".into());
    }
    if header[0].first().map(|m| *m != 0.0).unwrap_or(false) {
        rep.skip(format!("Mach {}", header[0][0]), "the Mach number belongs to the analysis definition");
    }
    if header[1].first().map(|v| *v != 0.0).unwrap_or(false) {
        rep.skip("IYsym (image about y = Ysym)", "the geometry is imported as written; make wings two-sided to model the image");
    }
    if header[1].get(1).map(|v| *v != 0.0).unwrap_or(false) {
        rep.skip("IZsym (ground or ceiling image)", "set the ground height in the analysis definition");
    }
    rep.skip("Sref, Cref, Bref", "Foiley derives the reference dimensions from the main wing (or the definition's custom ones)");
    rep.skip("Xref, Yref, Zref", "the moment reference is the plane's centre of gravity");
    if header.get(4).map(|c| c.first().map(|v| *v != 0.0).unwrap_or(false)).unwrap_or(false) {
        rep.skip("CDp", "profile drag comes from the 2D polars; add an extra drag to the definition if needed");
    }

    let mut surfs: Vec<ASurf> = Vec::new();
    let mut in_body = false;
    let next = |it: &mut usize| -> Result<(&str, &str), String> {
        *it += 1;
        lines.get(*it).copied().ok_or_else(|| "AVL file ends inside a keyword block".to_string())
    };
    while it < lines.len() {
        let (l, _) = lines[it];
        let Some(kw) = keyword(l) else {
            rep.skip(format!("line '{l}'"), "not an AVL keyword");
            it += 1;
            continue;
        };
        match kw.as_str() {
            "SURF" => {
                in_body = false;
                let name = next(&mut it)?.1.to_string();
                let v = nums(next(&mut it)?.0);
                if v.len() < 2 {
                    return Err(format!("surface '{name}': expected Nchord Cspace"));
                }
                surfs.push(ASurf {
                    title: name,
                    nchord: v[0].max(1.0) as usize,
                    cspace: v[1],
                    nspan: v.get(2).map(|n| n.max(1.0) as usize),
                    sspace: v.get(3).copied(),
                    ydup: None,
                    scale: Vec3::new(1.0, 1.0, 1.0),
                    translate: Vec3::ZERO,
                    angle: 0.0,
                    secs: vec![],
                });
            }
            "BODY" => {
                in_body = true;
                let name = next(&mut it)?.1.to_string();
                next(&mut it)?;
                rep.skip(format!("BODY '{name}'"), "AVL slender bodies are not imported; add a fuselage and import its shape");
            }
            "COMP" | "INDE" => {
                next(&mut it)?;
            }
            "YDUP" | "SCAL" | "TRAN" | "ANGL" => {
                let v = nums(next(&mut it)?.0);
                let s = if in_body { None } else { surfs.last_mut() };
                if let Some(s) = s {
                    match kw.as_str() {
                        "YDUP" => s.ydup = v.first().copied(),
                        "SCAL" if v.len() >= 3 => s.scale = Vec3::new(v[0], v[1], v[2]),
                        "TRAN" if v.len() >= 3 => s.translate = Vec3::new(v[0], v[1], v[2]),
                        "ANGL" => s.angle = v.first().copied().unwrap_or(0.0),
                        _ => return Err(format!("{kw}: expected 3 values")),
                    }
                }
            }
            "NOWA" | "NOAL" | "NOLO" => {
                let s = surfs.last().map(|s| s.title.clone()).unwrap_or_default();
                rep.skip(format!("{} on '{s}'", l.split_whitespace().next().unwrap_or("")), "an AVL solver flag with no Foiley equivalent");
            }
            "SECT" => {
                let v = nums(next(&mut it)?.0);
                let s = surfs.last_mut().ok_or("SECTION before any SURFACE")?;
                if v.len() < 5 {
                    return Err(format!("surface '{}': SECTION needs Xle Yle Zle Chord Ainc", s.title));
                }
                s.secs.push(ASec { le: Vec3::new(v[0], v[1], v[2]), chord: v[3], ainc: v[4], nspan: v.get(5).map(|n| n.max(0.0) as usize), sspace: v.get(6).copied(), foil: String::new(), naca: false, ctrls: vec![] });
            }
            "NACA" => {
                if nums(l.split_once(char::is_whitespace).map(|x| x.1).unwrap_or("")).len() >= 2 {
                    rep.approx("NACA X1 X2 range", "the whole foil is used");
                }
                let d = next(&mut it)?.0.split_whitespace().next().unwrap_or("").to_string();
                if let Some(sec) = surfs.last_mut().and_then(|s| s.secs.last_mut()) {
                    sec.foil = d;
                    sec.naca = true;
                }
            }
            "AFIL" => {
                let f = next(&mut it)?.1.trim().trim_matches('"').to_string();
                let stem = std::path::Path::new(&f).file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or(f.clone());
                let hit = foils.iter().find(|x| x.name == stem || x.name.replace(' ', "_") == stem).or_else(|| foils.iter().find(|x| x.name.replace(' ', "_").eq_ignore_ascii_case(&stem)));
                let name = match hit {
                    Some(x) => x.name.clone(),
                    None => {
                        rep.approx(format!("AFILE {f}"), format!("the .dat file is not read here: load it as foil '{stem}'"));
                        stem
                    }
                };
                if let Some(sec) = surfs.last_mut().and_then(|s| s.secs.last_mut()) {
                    sec.foil = name;
                }
            }
            "AIRF" => {
                let s = surfs.last().map(|s| s.title.clone()).unwrap_or_default();
                while it + 1 < lines.len() && keyword(lines[it + 1].0).is_none() {
                    it += 1;
                }
                rep.skip(format!("inline AIRFOIL coordinates on '{s}'"), "load the foil as a .dat file and assign it to the section");
                if let Some(sec) = surfs.last_mut().and_then(|s| s.secs.last_mut()) {
                    sec.foil = format!("{s} airfoil");
                }
            }
            "BFIL" => {
                next(&mut it)?;
            }
            "CLAF" | "CDCL" | "DESI" => {
                next(&mut it)?;
                rep.skip(kw.clone(), "AVL section polar / design data; Foiley uses the 2D polars of the section foils");
            }
            "CONT" => {
                let raw = next(&mut it)?.0;
                let mut tok = raw.split_whitespace();
                let name = tok.next().unwrap_or("control").to_string();
                let v = nums(&tok.collect::<Vec<_>>().join(" "));
                if v.len() < 2 {
                    return Err(format!("CONTROL '{name}': expected gain and Xhinge"));
                }
                if v.len() >= 5 && v[2..5].iter().any(|x| *x != 0.0) {
                    rep.approx(format!("CONTROL '{name}' hinge vector"), "Foiley hinges run along the surface's hinge line");
                }
                let c = ACtrl { name, gain: v[0], xhinge: v[1], sgn_dup: v.get(5).copied().unwrap_or(1.0) };
                if let Some(sec) = surfs.last_mut().and_then(|s| s.secs.last_mut()) {
                    sec.ctrls.push(c);
                }
            }
            _ => {}
        }
        it += 1;
    }

    // surfaces -> wings
    let mut plane = Plane { name: title, description: "imported from AVL".into(), wings: vec![], fuses: vec![], point_masses: vec![] };
    let mut mixing = Vec::new();
    let mut naca_hinge: HashMap<String, Vec<Option<f64>>> = HashMap::new();
    let mut afile_hinge: HashMap<String, Vec<Option<f64>>> = HashMap::new();
    let mut i = 0;
    while i < surfs.len() {
        let mut s = surfs[i].clone();
        apply_transforms(&mut s);
        if s.secs.len() < 2 {
            rep.skip(format!("surface '{}'", s.title), "fewer than 2 sections");
            i += 1;
            continue;
        }
        let pair = if s.ydup.is_none() && i + 1 < surfs.len() && surfs[i + 1].ydup.is_none() {
            let mut l = surfs[i + 1].clone();
            apply_transforms(&mut l);
            mirror_plane(&s.secs, &l.secs).map(|ym| (l, ym))
        } else {
            None
        };
        let iw = plane.wings.len();
        let (wing, ctrl_sides) = match pair {
            Some((l, ym)) => {
                let lsecs: Vec<ASec> = l.secs.iter().rev().cloned().collect();
                let title = strip_side(&s.title);
                let w = make_wing(&title, &s, &s.secs, Some(&lsecs), Some(ym), &mut rep);
                rep.map(format!("surfaces '{}' and '{}' -> two-sided wing '{}'", s.title, l.title, w.name));
                i += 2;
                (w, vec![(true, lsecs, false), (false, s.secs.clone(), false)])
            }
            None => {
                let two = s.ydup.is_some();
                let w = make_wing(&s.title, &s, &s.secs, None, s.ydup, &mut rep);
                rep.map(format!("surface '{}' -> {} wing '{}'", s.title, if two { "two-sided" } else { "one-sided" }, w.name));
                i += 1;
                let mut sides = vec![(false, s.secs.clone(), false)];
                if two {
                    sides.insert(0, (true, s.secs.clone(), true));
                }
                (w, sides)
            }
        };
        // flaps: wing surfaces whose two sections share a control; numbered
        // like the mesh (left side tip to root, then right side root to tip)
        let mut index = 0usize;
        for (left, secs, dup) in &ctrl_sides {
            let order: Vec<usize> = if *left { (0..secs.len() - 1).rev().collect() } else { (0..secs.len() - 1).collect() };
            for k in order {
                if let Some(c) = secs[k].ctrls.iter().find(|c| secs[k + 1].ctrls.iter().any(|d| d.name == c.name)) {
                    let gain = if *dup { c.gain * c.sgn_dup } else { c.gain };
                    mixing.push((c.name.clone(), ControlTarget::Flap { wing: iw, index }, gain));
                    rep.map(format!("CONTROL '{}' -> {} flap {} ({} side, gain {})", c.name, wing.name, index + 1, if *left { "left" } else { "right" }, num(gain)));
                    index += 1;
                }
            }
        }
        // the YDUPLICATE image repeats the same sections: count them once
        for (_, secs, _) in ctrl_sides.iter().filter(|c| !c.2) {
            for sec in secs {
                let h = sec.ctrls.first().map(|c| c.xhinge);
                if sec.naca {
                    naca_hinge.entry(sec.foil.clone()).or_default().push(h);
                } else {
                    afile_hinge.entry(sec.foil.clone()).or_default().push(h);
                }
            }
        }
        plane.wings.push(wing);
    }
    // NACA foils, hinged where their sections carry a control
    let mut names: Vec<&String> = naca_hinge.keys().collect();
    names.sort();
    let mut out_foils = Vec::new();
    for d in names {
        let uses = &naca_hinge[d];
        let mut f = Foil::naca(d).map_err(|e| format!("NACA {d}: {e}"))?;
        let hinges: Vec<f64> = uses.iter().flatten().copied().collect();
        if let Some(&h) = hinges.first() {
            f.hinge = Some((h, 0.5));
            if hinges.iter().any(|x| (x - h).abs() > 1e-9) || hinges.len() != uses.len() {
                rep.approx(format!("foil {}", f.name), format!("one hinge at x/c = {} is used wherever this foil is (the file uses it with different or no controls)", num(h)));
            }
        }
        out_foils.push(f);
    }
    // AFILE foils carrying a control need a hinge at its x/c
    let mut names: Vec<&String> = afile_hinge.keys().collect();
    names.sort();
    for name in names {
        let uses = &afile_hinge[name];
        let hinges: Vec<f64> = uses.iter().flatten().copied().collect();
        let Some(&h) = hinges.first() else { continue };
        let mixed = hinges.iter().any(|x| (x - h).abs() > 1e-9) || hinges.len() != uses.len();
        match foils.iter().find(|f| &f.name == name) {
            Some(f) if f.hinge.map(|g| (g.0 - h).abs() < 1e-6).unwrap_or(false) && !mixed => {}
            Some(f) => {
                let mut f = f.clone();
                f.hinge = Some((h, 0.5));
                rep.approx(format!("foil {}", f.name), format!("hinge set at x/c = {} for its CONTROL (Foiley flaps come from foil hinges){}", num(h), if mixed { "; the file uses this foil with different or no controls too" } else { "" }));
                out_foils.push(f);
            }
            None => rep.approx(format!("foil {name}"), format!("its sections carry a CONTROL: give the foil a hinge at x/c = {} once it is loaded", num(h))),
        }
    }
    if plane.wings.is_empty() {
        return Err("the AVL file has no wing surface".into());
    }
    Ok(AvlImport { plane, foils: out_foils, mixing, report: rep })
}

fn apply_transforms(s: &mut ASurf) {
    let (sc, tr, an) = (s.scale, s.translate, s.angle);
    for sec in &mut s.secs {
        sec.le = Vec3::new(sec.le.x * sc.x + tr.x, sec.le.y * sc.y + tr.y, sec.le.z * sc.z + tr.z);
        sec.chord *= sc.x;
        sec.ainc += an;
    }
    s.scale = Vec3::new(1.0, 1.0, 1.0);
    s.translate = Vec3::ZERO;
    s.angle = 0.0;
}

/// `l` (listed tip to root) is the mirror image of `r` (root to tip) about a
/// plane y = ym: returns ym.
fn mirror_plane(r: &[ASec], l: &[ASec]) -> Option<f64> {
    let n = r.len();
    if n < 2 || l.len() != n {
        return None;
    }
    let ym = 0.5 * (r[0].le.y + l[n - 1].le.y);
    let scale = r.iter().chain(l).map(|s| s.le.norm() + s.chord).fold(1.0, f64::max);
    let tol = 1e-9 * scale;
    let ok = (0..n).all(|k| {
        let (a, b) = (&r[k], &l[n - 1 - k]);
        (a.le.x - b.le.x).abs() < tol && ((2.0 * ym - a.le.y) - b.le.y).abs() < tol && (a.le.z - b.le.z).abs() < tol && (a.chord - b.chord).abs() < tol && (a.ainc - b.ainc).abs() < 1e-9
    });
    (ok && r[n - 1].le.y > ym + tol).then_some(ym)
}

fn strip_side(t: &str) -> String {
    t.strip_suffix(" right").or_else(|| t.strip_suffix(" left")).unwrap_or(t).to_string()
}

/// Wing name and kind from a surface title ("name (Main wing)" as written by
/// the exporter, or a plain AVL name).
fn name_kind(title: &str, vertical: bool, two_sided: bool) -> (String, WingKind) {
    for k in WingKind::ALL {
        if let Some(n) = title.strip_suffix(&format!(" ({})", k.label())) {
            return (n.to_string(), k);
        }
    }
    let l = title.to_lowercase();
    let kind = if vertical {
        WingKind::Fin
    } else if ["elev", "stab", "tail", "htp"].iter().any(|k| l.contains(k)) {
        WingKind::Elevator
    } else if two_sided {
        WingKind::Main
    } else {
        WingKind::Other
    };
    (title.to_string(), kind)
}

fn foil_name(s: &ASec) -> String {
    if s.naca {
        format!("NACA {}", s.foil)
    } else {
        s.foil.clone()
    }
}

/// A Foiley wing through the AVL sections.  AVL places each section at its
/// (rotated) leading edge and turns it about the LE by Ainc; Foiley twists
/// about the quarter chord, so the untwisted leading edges are found by a
/// fixed-point iteration through `Wing::frames`.
fn make_wing(title: &str, s: &ASurf, right: &[ASec], left: Option<&[ASec]>, mirror_y: Option<f64>, rep: &mut ConversionReport) -> Wing {
    let n = right.len();
    let two_sided = mirror_y.is_some();
    let (dy, dz) = (right[n - 1].le.y - right[0].le.y, right[n - 1].le.z - right[0].le.z);
    let (name, kind) = name_kind(title, dz.abs() > dy.abs() && !two_sided, two_sided);
    let (xdist, exact) = spacing_from_avl(s.cspace);
    if !exact {
        rep.approx(format!("Cspace {} of '{name}'", num(s.cspace)), "rounded to the nearest AVL spacing type");
    }
    // spanwise panels: per section, or the surface's Nspan shared by span length
    let len = |k: usize| ((right[k + 1].le.y - right[k].le.y).powi(2) + (right[k + 1].le.z - right[k].le.z).powi(2)).sqrt();
    let total: f64 = (0..n - 1).map(len).sum::<f64>().max(1e-12);
    let mut sections = Vec::with_capacity(n);
    for (k, a) in right.iter().enumerate() {
        let mut sec = WingSection::new(0.0, a.chord, 0.0, &foil_name(a));
        if let Some(l) = left {
            sec.foil_left = foil_name(&l[k]);
        }
        sec.twist = a.ainc;
        sec.nx = s.nchord.clamp(1, 200);
        sec.xdist = xdist;
        let (ny, sp) = match (s.nspan, a.nspan) {
            (Some(ns), _) => (((ns as f64) * if k + 1 < n { len(k) / total } else { 1.0 }).round().max(1.0) as usize, s.sspace.unwrap_or(1.0)),
            (None, Some(ny)) => (ny.max(1), a.sspace.unwrap_or(1.0)),
            (None, None) => (8, 1.0),
        };
        sec.ny = ny.clamp(1, 500);
        sec.ydist = spacing_from_avl(sp).0;
        sections.push(sec);
    }
    let symmetric = sections.iter().all(|x| x.foil_left == x.foil_right);
    let mut w = Wing::new(&name, kind, sections);
    w.two_sided = two_sided;
    w.symmetric = symmetric;
    w.close_inner = !two_sided;
    let target: Vec<Vec3> = right.iter().map(|a| a.le).collect();
    let scale = target.iter().map(|p| p.norm()).fold(1.0, f64::max);
    let mut u = target.clone();
    for _ in 0..200 {
        place_sections(&mut w, &u, mirror_y);
        let fr = w.frames(false);
        let mut err: f64 = 0.0;
        for k in 0..n {
            let d = target[k] - fr[k].le();
            u[k] = u[k] + d;
            err = err.max(d.norm());
        }
        if err <= 1e-15 * scale {
            break;
        }
    }
    place_sections(&mut w, &u, mirror_y);
    w
}

/// Wing position, developed span stations, dihedrals and offsets through the
/// untwisted leading edges `u` (plane axes, rx = ry = 0).
fn place_sections(w: &mut Wing, u: &[Vec3], mirror_y: Option<f64>) {
    let y0 = mirror_y.unwrap_or(u[0].y);
    w.position = Vec3::new(0.0, y0, u[0].z);
    let mut y = u[0].y - y0;
    let n = u.len();
    for k in 0..n {
        if k > 0 {
            let (dy, dz) = (u[k].y - u[k - 1].y, u[k].z - u[k - 1].z);
            let l = (dy * dy + dz * dz).sqrt();
            w.sections[k - 1].dihedral = if l > 0.0 { dz.atan2(dy) / DEG } else { 0.0 };
            y += l;
        }
        w.sections[k].y = y;
        w.sections[k].offset = u[k].x;
    }
    if n >= 2 {
        w.sections[n - 1].dihedral = w.sections[n - 2].dihedral;
    }
}

/// AVL mass file: one lumped line per point mass plus distributed masses
/// as their equivalent lumped inertia.
pub fn avl_mass(plane: &Plane, db: &FoilDb) -> String {
    let i = plane.inertia(db);
    let mut o = String::from("# AVL mass file written by Foiley\nLunit = 1.0 m\nMunit = 1.0 kg\nTunit = 1.0 s\ng = 9.81\nrho = 1.225\n#  mass   x   y   z   Ixx   Iyy   Izz  Ixy  Ixz  Iyz\n");
    o.push_str(&format!(
        "{:.6} {:.6} {:.6} {:.6} {:.6e} {:.6e} {:.6e} {:.6e} {:.6e} {:.6e}  ! whole aircraft (about its CoG)\n",
        i.mass, i.cog.x, i.cog.y, i.cog.z, i.ixx, i.iyy, i.izz, i.ixy, i.ixz, i.iyz
    ));
    o
}

/// The panel mesh (thick surfaces + fuselage) as binary STL.
pub fn stl_mesh(plane: &Plane, db: &FoilDb, thick: bool) -> Result<Vec<u8>, PlaneError> {
    let m = build(plane, db, &MeshOptions { thick, include_fuse: true, allow_missing_foils: true, ..Default::default() })?;
    let mut tris = Vec::with_capacity(m.panels.len() * 2);
    for p in &m.panels {
        let [a, b, c, d] = p.nodes;
        tris.push([a, b, c]);
        if !p.is_tri() {
            tris.push([a, c, d]);
        }
    }
    if !thick {
        // thin surfaces: also show the fuselage
        for f in &plane.fuses {
            tris.extend(f.world_triangles(f.nx.max(8), f.nh.max(6)));
        }
    }
    Ok(crate::cad::write_stl(&tris, &plane.name))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn avl_has_surfaces() {
        let mut db = FoilDb::new();
        db.add_foil(&foiley_foil::Foil::naca("2412").unwrap());
        db.add_foil(&foiley_foil::Foil::naca("0009").unwrap());
        let p = Plane::default_plane("NACA 2412", "NACA 0009");
        let t = avl(&p, &db, 0.0);
        assert_eq!(t.matches("SURFACE").count(), 3);
        assert!(t.contains("NACA\n2412") && t.contains("YDUPLICATE"));
        assert!(avl_mass(&p, &db).contains("whole aircraft"));
        let stl = stl_mesh(&p, &db, true).unwrap();
        assert!(crate::cad::read_stl(&stl).unwrap().len() > 500);
    }
}
