//! Import (and export) of plane definitions in the xflr5 (`<explane>`) and
//! flow5 (`<xflplane>`) XML formats.

use crate::geom::{Fuse, FuseFrame, NurbsBody, NurbsFrame, Plane, PointMass, Wing, WingKind, WingSection};
use crate::PlaneError;
use foiley_core::{Spacing, Vec3};
use quick_xml::events::Event;
use quick_xml::Reader;

#[derive(Debug, Default, Clone)]
struct Node {
    name: String,
    text: String,
    children: Vec<Node>,
}

impl Node {
    fn child(&self, n: &str) -> Option<&Node> {
        self.children.iter().find(|c| c.name == n)
    }
    fn all<'a>(&'a self, n: &'a str) -> impl Iterator<Item = &'a Node> + 'a {
        self.children.iter().filter(move |c| c.name == n)
    }
    fn txt(&self, n: &str) -> Option<&str> {
        self.child(n).map(|c| c.text.trim())
    }
    fn f(&self, n: &str) -> Option<f64> {
        self.txt(n).and_then(|t| t.parse().ok())
    }
    fn b(&self, n: &str) -> Option<bool> {
        self.txt(n).map(|t| t.eq_ignore_ascii_case("true") || t == "1")
    }
    fn v3(&self, n: &str) -> Option<Vec3> {
        let t = self.txt(n)?;
        let v: Vec<f64> = t.split(|c: char| c == ',' || c.is_whitespace()).filter(|s| !s.is_empty()).filter_map(|s| s.parse().ok()).collect();
        if v.len() >= 3 {
            Some(Vec3::new(v[0], v[1], v[2]))
        } else {
            None
        }
    }
}

fn unescape(s: &str) -> String {
    s.replace("&lt;", "<").replace("&gt;", ">").replace("&quot;", "\"").replace("&apos;", "'").replace("&amp;", "&")
}

fn parse(text: &str) -> Result<Node, PlaneError> {
    let mut r = Reader::from_str(text);
    r.config_mut().trim_text(true);
    let mut stack: Vec<Node> = vec![Node { name: "#root".into(), ..Default::default() }];
    loop {
        match r.read_event() {
            Ok(Event::Start(e)) => stack.push(Node { name: e.name().as_ref().to_lowercase(), ..Default::default() }),
            Ok(Event::Empty(e)) => {
                let n = Node { name: e.name().as_ref().to_lowercase(), ..Default::default() };
                stack.last_mut().unwrap().children.push(n);
            }
            Ok(Event::Text(t)) => {
                let s = unescape(&t.into_inner());
                if let Some(top) = stack.last_mut() {
                    top.text.push_str(&s);
                }
            }
            Ok(Event::CData(t)) => {
                if let Some(top) = stack.last_mut() {
                    top.text.push_str(&t.into_inner());
                }
            }
            Ok(Event::GeneralRef(g)) => {
                let s = match &*g.into_inner() {
                    "amp" => "&",
                    "lt" => "<",
                    "gt" => ">",
                    "quot" => "\"",
                    "apos" => "'",
                    _ => "",
                };
                if let Some(top) = stack.last_mut() {
                    top.text.push_str(s);
                }
            }
            Ok(Event::End(_)) => {
                let n = stack.pop().unwrap();
                match stack.last_mut() {
                    Some(p) => p.children.push(n),
                    None => return Err(PlaneError::Input("unbalanced XML".into())),
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => return Err(PlaneError::Input(format!("XML error at byte {}: {e}", r.error_position()))),
            _ => {}
        }
    }
    stack.pop().ok_or_else(|| PlaneError::Input("empty XML".into()))
}

fn dist(s: &str) -> Spacing {
    match s.trim().to_uppercase().as_str() {
        "COSINE" => Spacing::Cosine,
        "SINE" | "EXP" => Spacing::MinusSine,
        "INVERSESINE" | "INV_SINE" | "INV_EXP" => Spacing::Sine,
        "TANH" | "INV_SINH" => Spacing::Tanh,
        _ => Spacing::Uniform,
    }
}

fn dist_name(s: Spacing) -> &'static str {
    match s {
        Spacing::Uniform => "UNIFORM",
        Spacing::Cosine => "COSINE",
        Spacing::Sine => "INVERSESINE",
        Spacing::MinusSine => "SINE",
        Spacing::Tanh => "TANH",
    }
}

fn point_masses(n: &Node, lf: f64, mf: f64) -> Vec<PointMass> {
    n.all("point_mass")
        .map(|p| PointMass { tag: p.txt("tag").unwrap_or("").to_string(), mass: p.f("mass").unwrap_or(0.0) * mf, position: p.v3("coordinates").unwrap_or_default() * lf })
        .collect()
}

/// Read a plane from xflr5 / flow5 XML text.
pub fn read_plane(text: &str) -> Result<Plane, PlaneError> {
    let root = parse(text)?;
    let top = root.children.iter().find(|c| c.name == "explane" || c.name == "xflplane").ok_or_else(|| PlaneError::Input("not an xflr5 <explane> or flow5 <xflplane> file".into()))?;
    let flow5 = top.name == "xflplane";
    let (mut lf, mut mf) = (1.0, 1.0);
    if let Some(u) = top.child("units") {
        if let Some(v) = u.f("length_unit_to_meter") {
            lf = v;
        }
        if let Some(v) = u.f("meter_to_length_unit") {
            lf = 1.0 / v;
        }
        if let Some(v) = u.f("mass_unit_to_kg") {
            mf = v;
        }
        if let Some(v) = u.f("kg_to_mass_unit") {
            mf = 1.0 / v;
        }
    }
    let p = top.child("plane").ok_or_else(|| PlaneError::Input("missing <Plane>".into()))?;
    let mut plane = Plane { name: p.txt("name").unwrap_or("Imported plane").to_string(), description: p.txt("description").unwrap_or("").to_string(), wings: vec![], fuses: vec![], point_masses: vec![] };
    if let Some(i) = p.child("inertia") {
        plane.point_masses = point_masses(i, lf, mf);
    }
    for w in p.all("wing") {
        let kind = match w.txt("type").unwrap_or("").to_uppercase().as_str() {
            "MAINWING" => WingKind::Main,
            "ELEVATOR" => WingKind::Elevator,
            "FIN" => WingKind::Fin,
            _ => WingKind::Other,
        };
        let mut secs = Vec::new();
        if let Some(ss) = w.child("sections") {
            for s in ss.all("section") {
                let fl = s.txt("left_side_foilname").or(s.txt("left_side_foil_file")).unwrap_or("").trim_end_matches(".dat").to_string();
                let fr = s.txt("right_side_foilname").or(s.txt("right_side_foil_file")).unwrap_or("").trim_end_matches(".dat").to_string();
                secs.push(WingSection {
                    y: s.f("y_position").unwrap_or(0.0) * lf,
                    chord: s.f("chord").unwrap_or(0.1) * lf,
                    offset: s.f("xoffset").unwrap_or(0.0) * lf,
                    dihedral: s.f("dihedral").unwrap_or(0.0),
                    twist: s.f("twist").unwrap_or(0.0),
                    foil_left: fl.clone(),
                    foil_right: if fr.is_empty() { fl } else { fr },
                    nx: s.f("x_number_of_panels").unwrap_or(7.0) as usize,
                    ny: s.f("y_number_of_panels").unwrap_or(7.0) as usize,
                    xdist: dist(s.txt("x_panel_distribution").unwrap_or("COSINE")),
                    ydist: dist(s.txt("y_panel_distribution").unwrap_or("UNIFORM")),
                });
            }
        }
        let mut wing = Wing::new(w.txt("name").unwrap_or("wing"), kind, secs);
        wing.position = w.v3("position").unwrap_or_default() * lf;
        wing.symmetric = w.b("symetric").or(w.b("symmetric")).unwrap_or(true);
        if flow5 {
            wing.two_sided = w.b("two_sided").unwrap_or(true);
            wing.close_inner = w.b("closed_inner_side").unwrap_or(false);
            let rx = w.f("rx_angle").unwrap_or(0.0);
            // flow5 builds one-sided wings on the left half; we build them on the right
            wing.rx = if wing.two_sided { rx } else { -rx };
            wing.ry = w.f("ry_angle").unwrap_or(0.0);
        } else {
            wing.ry = w.f("tilt_angle").unwrap_or(0.0);
            let is_fin = w.b("isfin").unwrap_or(kind == WingKind::Fin);
            if is_fin {
                wing.rx = 90.0;
                wing.two_sided = w.b("issymfin").unwrap_or(false);
                wing.close_inner = true;
            }
        }
        if let Some(i) = w.child("inertia") {
            wing.mass = i.f("volume_mass").or(i.f("mass")).unwrap_or(0.0) * mf;
            wing.point_masses = point_masses(i, lf, mf);
        }
        let double_fin = !flow5 && w.b("isdoublefin").unwrap_or(false);
        if double_fin {
            let mut a = wing.clone();
            let mut b = wing.clone();
            let y = wing.position.y.abs();
            a.position.y = -y;
            b.position.y = y;
            a.name = format!("{} (left)", wing.name);
            b.name = format!("{} (right)", wing.name);
            plane.wings.push(a);
            plane.wings.push(b);
        } else {
            plane.wings.push(wing);
        }
    }
    for b in p.all("body") {
        let frames_parent = b.child("nurbs").unwrap_or(b);
        let body_pos = b.v3("position").unwrap_or_default() * lf;
        let mut nframes: Vec<NurbsFrame> = Vec::new();
        for f in frames_parent.all("frame") {
            let pos = f.v3("position").unwrap_or_default() * lf;
            let pts: Vec<Vec3> = f
                .all("point")
                .filter_map(|n| {
                    let v: Vec<f64> = n.text.split(|c: char| c == ',' || c.is_whitespace()).filter(|s| !s.is_empty()).filter_map(|s| s.parse().ok()).collect();
                    (v.len() >= 3).then(|| Vec3::new(v[0], v[1], v[2]) * lf)
                })
                .collect();
            if pts.len() < 2 {
                continue;
            }
            // points carry their own x in flow5; xflr5 points may omit it: use the frame x
            let fx = if pts.iter().all(|p| p.x == 0.0) { pos.x } else { pts[0].x };
            let pts = pts.into_iter().map(|p| Vec3::new(fx, p.y.abs(), p.z)).collect();
            nframes.push(NurbsFrame { x: fx, points: pts });
        }
        nframes.sort_by(|a, b| a.x.total_cmp(&b.x));
        if nframes.len() >= 2 {
            // equal point counts are required for the tensor-product surface
            let np = nframes.iter().map(|f| f.points.len()).min().unwrap_or(0);
            for f in &mut nframes {
                f.points.truncate(np);
            }
            let frames: Vec<FuseFrame> = nframes
                .iter()
                .map(|f| {
                    let hw = f.points.iter().map(|p| p.y.abs()).fold(0.0, f64::max);
                    let zmax = f.points.iter().map(|p| p.z).fold(f64::MIN, f64::max);
                    let zmin = f.points.iter().map(|p| p.z).fold(f64::MAX, f64::min);
                    FuseFrame { x: f.x, z: 0.5 * (zmax + zmin), half_width: hw, half_height: 0.5 * (zmax - zmin) }
                })
                .collect();
            let flat = b.txt("type").map(|t| t.to_uppercase().contains("FLAT")).unwrap_or(false);
            let du = b.f("x_degree").or(frames_parent.f("u_degree")).unwrap_or(3.0) as usize;
            let dv = b.f("hoop_degree").or(frames_parent.f("v_degree")).unwrap_or(3.0) as usize;
            let nurbs = NurbsBody { frames: nframes, degree_u: if flat { 1 } else { du }, degree_v: if flat { 1 } else { dv } };
            plane.fuses.push(Fuse {
                name: b.txt("name").unwrap_or("Fuselage").to_string(),
                position: body_pos,
                frames,
                exponent: 2.0,
                nx: b.f("x_panels").unwrap_or(19.0) as usize,
                nh: (b.f("hoop_panels").unwrap_or(11.0) as usize).max(4),
                mass: b.child("inertia").and_then(|i| i.f("volume_mass")).unwrap_or(0.0) * mf,
                nurbs: Some(nurbs),
                triangles: vec![],
            });
        }
    }
    plane.validate().map_err(PlaneError::Geometry)?;
    Ok(plane)
}

/// Write a plane in the xflr5 `<explane>` format (SI units).
pub fn write_plane(p: &Plane) -> String {
    let mut s = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!DOCTYPE explane>\n<explane version=\"1.0\">\n    <Units>\n        <length_unit_to_meter>1</length_unit_to_meter>\n        <mass_unit_to_kg>1</mass_unit_to_kg>\n    </Units>\n    <Plane>\n");
    let esc = |t: &str| t.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
    s.push_str(&format!("        <Name>{}</Name>\n        <Description>{}</Description>\n        <Inertia>\n", esc(&p.name), esc(&p.description)));
    for m in &p.point_masses {
        s.push_str(&format!("            <Point_Mass><Tag>{}</Tag><Mass>{}</Mass><coordinates>{},{},{}</coordinates></Point_Mass>\n", esc(&m.tag), m.mass, m.position.x, m.position.y, m.position.z));
    }
    s.push_str("        </Inertia>\n");
    s.push_str(&format!("        <has_body>{}</has_body>\n", !p.fuses.is_empty()));
    for w in &p.wings {
        let ty = match w.kind {
            WingKind::Main => "MAINWING",
            WingKind::Elevator => "ELEVATOR",
            WingKind::Fin => "FIN",
            WingKind::Other => "SECONDWING",
        };
        let is_fin = w.rx.abs() > 45.0;
        s.push_str(&format!(
            "        <wing>\n            <Name>{}</Name>\n            <Type>{ty}</Type>\n            <Position>{},{},{}</Position>\n            <Tilt_angle>{}</Tilt_angle>\n            <Symetric>{}</Symetric>\n            <isFin>{is_fin}</isFin>\n            <isDoubleFin>false</isDoubleFin>\n            <isSymFin>{}</isSymFin>\n            <Inertia>\n                <Volume_Mass>{}</Volume_Mass>\n",
            esc(&w.name), w.position.x, w.position.y, w.position.z, w.ry, w.symmetric, is_fin && w.two_sided, w.mass
        ));
        for m in &w.point_masses {
            s.push_str(&format!("                <Point_Mass><Tag>{}</Tag><Mass>{}</Mass><coordinates>{},{},{}</coordinates></Point_Mass>\n", esc(&m.tag), m.mass, m.position.x, m.position.y, m.position.z));
        }
        s.push_str("            </Inertia>\n            <Sections>\n");
        for sec in &w.sections {
            s.push_str(&format!(
                "                <Section>\n                    <y_position>{}</y_position>\n                    <Chord>{}</Chord>\n                    <xOffset>{}</xOffset>\n                    <Dihedral>{}</Dihedral>\n                    <Twist>{}</Twist>\n                    <x_number_of_panels>{}</x_number_of_panels>\n                    <x_panel_distribution>{}</x_panel_distribution>\n                    <y_number_of_panels>{}</y_number_of_panels>\n                    <y_panel_distribution>{}</y_panel_distribution>\n                    <Left_Side_FoilName>{}</Left_Side_FoilName>\n                    <Right_Side_FoilName>{}</Right_Side_FoilName>\n                </Section>\n",
                sec.y, sec.chord, sec.offset, sec.dihedral, sec.twist, sec.nx, dist_name(sec.xdist), sec.ny, dist_name(sec.ydist), esc(&sec.foil_left), esc(&sec.foil_right)
            ));
        }
        s.push_str("            </Sections>\n        </wing>\n");
    }
    for f in &p.fuses {
        let Some(n) = f.nurbs.clone().or_else(|| f.to_nurbs()) else { continue };
        s.push_str(&format!(
            "        <body>\n            <Name>{}</Name>\n            <Type>NURBS</Type>\n            <Position>{},{},{}</Position>\n            <x_degree>{}</x_degree>\n            <hoop_degree>{}</hoop_degree>\n            <x_panels>{}</x_panels>\n            <hoop_panels>{}</hoop_panels>\n            <Inertia>\n                <Volume_Mass>{}</Volume_Mass>\n            </Inertia>\n",
            esc(&f.name), f.position.x, f.position.y, f.position.z, n.degree_u, n.degree_v, f.nx, f.nh, f.mass
        ));
        for fr in &n.frames {
            s.push_str(&format!("            <frame>\n                <Position>{},0,0</Position>\n", fr.x));
            for pt in &fr.points {
                s.push_str(&format!("                <point>{},{},{}</point>\n", pt.x, pt.y, pt.z));
            }
            s.push_str("            </frame>\n");
        }
        s.push_str("        </body>\n");
    }
    s.push_str("    </Plane>\n</explane>\n");
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn roundtrip() {
        let p = Plane::default_plane("NACA 2412", "NACA 0009");
        let x = write_plane(&p);
        let q = read_plane(&x).unwrap();
        assert_eq!(q.wings.len(), 3);
        assert_eq!(q.wings[0].sections.len(), 2);
        assert!((q.wings[0].sections[1].chord - 0.19).abs() < 1e-12);
        assert!(q.wings[2].rx > 45.0 && !q.wings[2].two_sided);
        assert!((q.point_masses[0].mass - 0.6).abs() < 1e-12);
    }
}
