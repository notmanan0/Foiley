//! Exports: AVL geometry (.avl) and mass (.mass) files, STL meshes.

use crate::foildb::FoilDb;
use crate::geom::{Plane, WingKind};
use crate::mesh::{build, MeshOptions};
use crate::PlaneError;

fn naca_digits(name: &str) -> Option<String> {
    let l = name.to_lowercase().replace([' ', '-', '_'], "");
    let d = l.strip_prefix("naca")?;
    (d.len() == 4 && d.chars().all(|c| c.is_ascii_digit())).then(|| d.to_string())
}

/// AVL geometry file.  Foils that are not NACA 4-digit are referenced as
/// `AFILE <name>.dat`; write those files next to the .avl (see
/// [`avl_foil_files`]).  Control surfaces become AVL CONTROL entries.
pub fn avl(plane: &Plane, db: &FoilDb, mach: f64) -> String {
    let (s, b, c) = plane.reference(true);
    let i = plane.inertia(db);
    let mut o = String::new();
    o.push_str(&format!("{}\n#Mach\n{mach:.3}\n#IYsym IZsym Zsym\n0 0 0.0\n#Sref Cref Bref\n{s:.6} {c:.6} {b:.6}\n#Xref Yref Zref\n{:.6} {:.6} {:.6}\n#CDp\n0.0\n", plane.name, i.cog.x, i.cog.y, i.cog.z));
    for w in &plane.wings {
        let nx = w.sections.first().map(|s| s.nx).unwrap_or(8);
        let ny: usize = w.sections.iter().take(w.sections.len() - 1).map(|s| s.ny).sum();
        o.push_str(&format!("#{}\nSURFACE\n{} ({})\n#Nchord Cspace Nspan Sspace\n{} 1.0 {} -2.0\n", "=".repeat(60), w.name, w.kind.label(), nx, ny.max(1)));
        if w.two_sided && w.position.y.abs() < 1e-9 {
            o.push_str("YDUPLICATE\n0.0\n");
        }
        if w.kind == WingKind::Fin && !w.two_sided {
            o.push_str("#vertical surface\n");
        }
        let fr = w.frames(false);
        for (k, sec) in w.sections.iter().enumerate() {
            let f = &fr[k];
            let le = f.le();
            let inc = sec.twist + w.ry;
            o.push_str(&format!("SECTION\n#Xle Yle Zle Chord Ainc\n{:.6} {:.6} {:.6} {:.6} {:.3}\n", le.x, le.y, le.z, sec.chord, inc));
            match naca_digits(&sec.foil_right) {
                Some(d) => o.push_str(&format!("NACA\n{d}\n")),
                None => o.push_str(&format!("AFILE\n{}.dat\n", sec.foil_right.replace(' ', "_"))),
            }
            if let Some(sh) = db.shape(&sec.foil_right) {
                if let Some((xh, _)) = sh.hinge {
                    let name = match w.kind {
                        WingKind::Elevator => "elevator",
                        WingKind::Fin => "rudder",
                        _ => "flap",
                    };
                    let sgn = if w.kind == WingKind::Main { 1.0 } else { 1.0 };
                    o.push_str(&format!("CONTROL\n#name gain Xhinge XYZhvec SgnDup\n{name} 1.0 {xh:.4} 0.0 0.0 0.0 {sgn:.1}\n"));
                }
            }
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
