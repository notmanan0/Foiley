//! Reproduce the code-review findings tracked in docs/REMEDIATION_PLAN.md.
//! This diagnostic prints current behavior; it does not assert correctness.
use foiley_core::{Progress, Vec3};
use foiley_foil::{Foil, FoilPoint, FoilPolar};
use foiley_plane::analysis::*;
use foiley_plane::geom::*;
use foiley_plane::FoilDb;

fn run(p: &Plane, db: &FoilDb, d: &PlanePolarDef) -> PlaneOpp {
    match run_polar(p, db, d, &Progress::none(), &|_| {}) {
        Ok(mut pol) if !pol.points.is_empty() => pol.points.remove(0),
        Ok(pol) => {
            println!("  (no point: {:?})", pol.failed);
            PlaneOpp::default()
        }
        Err(e) => {
            println!("  (rejected: {e})");
            PlaneOpp::default()
        }
    }
}

fn polar(foil: &str, name: &str, re: f64, ncrit: f64, cd: f64) -> FoilPolar {
    FoilPolar {
        foil: foil.into(),
        name: name.into(),
        re,
        ncrit,
        polar_type: 1,
        viscous: true,
        points: (-20..=20)
            .map(|a| FoilPoint {
                alpha: a as f64,
                cl: 0.1 * a as f64,
                cd,
                re,
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    }
}

fn main() {
    let f = Foil::naca("0012").unwrap();
    let mut db = FoilDb::new();
    db.add_foil(&f);
    let mut root = WingSection::new(0.0, 1.0, 0.0, &f.name);
    root.nx = 3;
    root.ny = 8;
    let tip = WingSection::new(4.0, 1.0, 0.0, &f.name);
    let p = Plane {
        name: "review rectangle".into(),
        description: String::new(),
        wings: vec![Wing::new("wing", WingKind::Main, vec![root, tip])],
        fuses: vec![],
        point_masses: vec![],
    };
    let inertia = Inertia {
        mass: 10.0,
        cog: Vec3::new(0.25, 0.0, 0.0),
        ixx: 10.0,
        iyy: 10.0,
        izz: 10.0,
        ..Default::default()
    };
    let d = PlanePolarDef {
        viscous: false,
        free_wake: 4,
        wake_length: 10.0,
        inertia: Some(inertia),
        alpha: 8.0,
        range: Range {
            from: 8.0,
            to: 8.0,
            step: 1.0,
        },
        ..Default::default()
    };
    let fixed = run(&p, &db, &d);
    let t2 = run(
        &p,
        &db,
        &PlanePolarDef {
            kind: PolarKind::FixedLift,
            ..d.clone()
        },
    );
    println!(
        "T2 relaxed wake: lift/weight={:.9}",
        t2.lift / (inertia.mass * G)
    );
    let t3 = run(
        &p,
        &db,
        &PlanePolarDef {
            kind: PolarKind::Glide,
            ..d.clone()
        },
    );
    let t3rigid = run(
        &p,
        &db,
        &PlanePolarDef {
            kind: PolarKind::Glide,
            free_wake: 0,
            ..d.clone()
        },
    );
    println!(
        "T3 relaxed vs rigid: delta_CL={:.12}, retained_wake_lines={}",
        t3.cl - t3rigid.cl,
        t3.wake.len()
    );
    let t6 = run(
        &p,
        &db,
        &PlanePolarDef {
            kind: PolarKind::Control,
            range: Range {
                from: 0.0,
                to: 0.0,
                step: 1.0,
            },
            ..d.clone()
        },
    );
    println!(
        "T1/T6 equivalent relaxed settings: CL={:.9}/{:.9}, wake_lines={}/{}",
        fixed.cl,
        t6.cl,
        fixed.wake.len(),
        t6.wake.len()
    );
    let invalid = PlanePolarDef {
        kind: PolarKind::List,
        free_wake: 0,
        list: vec![(5.0, 0.0, 0.0)],
        ..d.clone()
    };
    println!(
        "T8 zero speed: validate={:?}, CL={}",
        invalid.validate(),
        run(&p, &db, &invalid).cl
    );

    let pa = polar(&f.name, "N9", 1e6, 9.0, 0.01);
    let pb = polar(&f.name, "N5", 1e6, 5.0, 0.03);
    let mut ab = db.clone();
    ab.add_polar(&pa);
    ab.add_polar(&pb);
    let mut ba = db.clone();
    ba.add_polar(&pb);
    ba.add_polar(&pa);
    let n9 = foiley_plane::FamilySelector { ncrit: Some(9.0), ..Default::default() };
    let cd = |db: &FoilDb| db.at_re_cl(db.resolve(&f.name, &n9).unwrap().unwrap(), 1e6, 0.5).unwrap().0.cd;
    println!("Polar insertion order (Ncrit 9 family): CD={}/{}", cd(&ab), cd(&ba));

    let mut ldb = db.clone();
    for re in [1e4, 1e8] {
        ldb.add_polar(&polar(&f.name, &format!("Re{re}"), re, 9.0, 0.3));
    }
    let llt = run(
        &p,
        &ldb,
        &PlanePolarDef {
            method: Method::Llt,
            kind: PolarKind::Glide,
            free_wake: 0,
            viscous: true,
            ..d.clone()
        },
    );
    println!(
        "LLT T3 constant-CD fixture: resultant/weight={:.9}",
        llt.lift.hypot(llt.drag) / (inertia.mass * G)
    );

    let f9 = Foil::naca("0009").unwrap();
    db.add_foil(&f9);
    let glider = Plane::default_plane(&f.name, &f9.name);
    let sd = PlanePolarDef {
        viscous: false,
        derivatives: true,
        range: Range {
            from: 4.0,
            to: 4.0,
            step: 1.0,
        },
        ..Default::default()
    };
    let clean = run(&glider, &db, &sd);
    let draggy = run(
        &glider,
        &db,
        &PlanePolarDef {
            extra_drag: vec![ExtraDrag {
                name: "audit drag".into(),
                area: 1.0,
                cd: 0.1,
            }],
            ..sd
        },
    );
    println!(
        "Extra drag: CD={:.9}/{:.9}, Xu={:.9}/{:.9}",
        clean.cd,
        draggy.cd,
        clean.stability.unwrap().dim.xu,
        draggy.stability.unwrap().dim.xu
    );

    let mut pods = p.clone();
    pods.fuses.push(Fuse::default_pod());
    let opts = foiley_plane::mesh::MeshOptions {
        thick: true,
        include_fuse: true,
        ..Default::default()
    };
    let count_fuse = |plane: &Plane| {
        foiley_plane::mesh::build(plane, &db, &opts)
            .unwrap()
            .panels
            .iter()
            .filter(|p| p.kind == foiley_plane::mesh::PanelKind::Fuse)
            .count()
    };
    let one = count_fuse(&pods);
    let mut second = Fuse::default_pod();
    second.position.y = 3.0;
    pods.fuses.push(second);
    println!(
        "One/two fuselages: body_panel_count={}/{}",
        one,
        count_fuse(&pods)
    );

    let mut asym = p.clone();
    asym.wings[0].symmetric = false;
    for sec in &mut asym.wings[0].sections {
        sec.foil_left = f9.name.clone();
    }
    let avl = foiley_plane::export::avl(&asym, &db, 0.0);
    println!(
        "Asymmetric AVL export: has_right_0012={}, has_left_0009={}, duplicates_right={}",
        avl.contains("NACA\n0012"),
        avl.contains("NACA\n0009"),
        avl.contains("YDUPLICATE")
    );
}
