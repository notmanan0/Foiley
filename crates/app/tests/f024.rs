//! F024: CAD import (STL / STEP reading and tessellation) runs as a job in
//! the job manager, with owned input bytes, progress, cancel, and errors and
//! notes reported through the job instead of being swallowed on the UI thread.

use foiley_app::jobs::{CadFormat, CadImport, JobKind, JobManager, JobOutput, JobStatus};
use foiley_core::{LogLevel, Vec3};
use foiley_plane::{Fuse, Plane};

const BOX_STEP: &str = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(0.,1.,0.));#3=CARTESIAN_POINT('',(1.,0.,0.));#4=CARTESIAN_POINT('',(1.,1.,0.5));\n\
#5=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#1,#2),(#3,#4)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.);\n\
ENDSEC;\nEND-ISO-10303-21;\n";

fn request(format: CadFormat, bytes: &[u8], body: usize) -> CadImport {
    CadImport { name: "part".into(), format, bytes: bytes.to_vec(), scale: 1.0, target_triangles: 200, body }
}

#[test]
fn f024_cad_import_is_job() {
    let mut m = JobManager::new(2);

    // a bad file fails visibly: terminal Failed status with the reader's message, logged as an error
    let bad = m.submit("import garbage.stl", JobKind::CadImport(request(CadFormat::Stl, b"garbage", 0)));
    let st = m.wait(bad.id, |_| {});
    let JobStatus::Failed(msg) = st else { panic!("bad file must fail, got {st:?}") };
    assert!(msg.contains("STL") || msg.contains("triangle"), "error names the problem: {msg}");
    let log = bad.progress.drain_log();
    assert!(log.iter().any(|l| l.level == LogLevel::Error && l.text.contains(&msg)), "error reaches the log: {log:?}");
    assert!(bad.take_output().is_none());

    // an unreadable STEP file fails too (no surfaces)
    let bad2 = m.submit("import empty.step", JobKind::CadImport(request(CadFormat::Step, b"ISO-10303-21;\nDATA;\nENDSEC;\n", 0)));
    assert!(matches!(m.wait(bad2.id, |_| {}), JobStatus::Failed(_)));

    // a good file finishes with triangles and notes, aimed at body 1
    let good = m.submit("import part.step", JobKind::CadImport(request(CadFormat::Step, BOX_STEP.as_bytes(), 1)));
    assert_eq!(m.wait(good.id, |_| {}), JobStatus::Done);
    let log = good.progress.drain_log();
    assert!(log.iter().any(|l| l.level == LogLevel::Info && l.text.contains("triangles")), "notes are logged: {log:?}");
    let Some(JobOutput::CadImport(r)) = good.take_output() else { panic!("expected a CAD import output") };
    assert_eq!(r.body, 1);
    assert!(r.triangles.len() >= 2);
    assert!(!r.notes.is_empty());

    // applying it to a plane fills the requested body and keeps body 0
    let mut p = Plane::default_plane("n0012", "n0012");
    p.fuses = vec![Fuse::default_pod()];
    let pod = p.fuses[0].clone();
    let summary = r.apply(&mut p).unwrap();
    assert!(summary.contains("triangles"), "{summary}");
    assert_eq!(p.fuses.len(), 2);
    assert_eq!(p.fuses[0], pod);
    assert!(!p.fuses[1].triangles.is_empty());
    assert_eq!(p.fuses[1].name, "part");

    // re-importing into body 0 keeps its structural mass
    let mut req = request(CadFormat::Step, BOX_STEP.as_bytes(), 0);
    req.scale = 2.0;
    let j = m.submit("import part.step", JobKind::CadImport(req));
    assert_eq!(m.wait(j.id, |_| {}), JobStatus::Done);
    let Some(JobOutput::CadImport(r)) = j.take_output() else { panic!() };
    p.fuses[0].mass = 3.5;
    r.apply(&mut p).unwrap();
    assert_eq!(p.fuses.len(), 2);
    assert_eq!(p.fuses[0].mass, 3.5);
    assert!(!p.fuses[0].triangles.is_empty());
    let len = p.fuses[0].length();
    assert!((len - 2.0).abs() < 1e-9, "scale applied: length {len}");

    // a body index beyond the end appends; one far beyond is refused
    let mut q = Plane::default_plane("n0012", "n0012");
    q.fuses.clear();
    let j = m.submit("import part.step", JobKind::CadImport(request(CadFormat::Step, BOX_STEP.as_bytes(), 5)));
    assert_eq!(m.wait(j.id, |_| {}), JobStatus::Done);
    let Some(JobOutput::CadImport(r)) = j.take_output() else { panic!() };
    assert!(r.apply(&mut q).is_err());
    assert!(q.fuses.is_empty());

    // a cancelled import ends Cancelled, not Done
    let big: Vec<u8> = {
        let t: Vec<[Vec3; 3]> = (0..200_000).map(|i| {
            let x = i as f64 * 1e-3;
            [Vec3::new(x, 0.0, 0.0), Vec3::new(x + 1e-3, 0.0, 0.0), Vec3::new(x, 1.0, 0.0)]
        }).collect();
        foiley_plane::cad::write_stl(&t, "big")
    };
    let j = m.submit("import big.stl", JobKind::CadImport(request(CadFormat::Stl, &big, 0)));
    j.cancel();
    assert_eq!(m.wait(j.id, |_| {}), JobStatus::Cancelled);

    // the format comes from the file extension
    assert_eq!(CadFormat::from_path(std::path::Path::new("a/b.STP")), Some(CadFormat::Step));
    assert_eq!(CadFormat::from_path(std::path::Path::new("b.stl")), Some(CadFormat::Stl));
    assert_eq!(CadFormat::from_path(std::path::Path::new("b.obj")), None);
}
