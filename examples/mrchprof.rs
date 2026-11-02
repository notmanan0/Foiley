use foiley::xfoil::{XFoil, DTOR};
use std::time::Instant;
fn main() {
    let mut xf = XFoil::new();
    xf.naca(2412);
    xf.lvisc = true; xf.reinf1 = 3e5; let mut r = 0.0; xf.mrcl(1.0, &mut r);
    xf.alfa = 4.0 * DTOR; xf.lalfa = true;
    xf.specal();
    xf.viscal(20);
    let base = xf.clone();
    let reps = 2000;
    let mut x = base.clone();
    let t = Instant::now();
    for _ in 0..reps { x.mrchdu(); }
    println!("mrchdu {:.1} us", t.elapsed().as_secs_f64() / reps as f64 * 1e6);
    // blsys micro
    let mut b = x.bl.clone();
    let t = Instant::now();
    for _ in 0..200000 { b.blvar(2); }
    println!("blvar(2) {:.1} ns", t.elapsed().as_secs_f64() / 200000.0 * 1e9);
    let t = Instant::now();
    for _ in 0..200000 { b.blvar(1); }
    println!("blvar(1) {:.1} ns", t.elapsed().as_secs_f64() / 200000.0 * 1e9);
    let t = Instant::now();
    for _ in 0..200000 { b.blkin(); }
    println!("blkin {:.1} ns", t.elapsed().as_secs_f64() / 200000.0 * 1e9);
    let t = Instant::now();
    for _ in 0..200000 { b.bldif(2); }
    println!("bldif(2) {:.1} ns", t.elapsed().as_secs_f64() / 200000.0 * 1e9);
    let t = Instant::now();
    for _ in 0..200000 { b.bldif(1); }
    println!("bldif(1) {:.1} ns", t.elapsed().as_secs_f64() / 200000.0 * 1e9);
    let t = Instant::now();
    let mut s = 0.0; for i in 0..1000000 { s += (1.0 + i as f64 * 1e-7f64).powf(0.43); }
    println!("powf {:.1} ns ({})", t.elapsed().as_secs_f64() / 1e6 * 1e9, s);
    let t = Instant::now();
    let mut s = 0.0; for i in 0..1000000 { s += (1.0 + i as f64 * 1e-7f64).ln(); }
    println!("ln {:.1} ns ({})", t.elapsed().as_secs_f64() / 1e6 * 1e9, s);
    let t = Instant::now();
    let mut s = 0.0; for i in 0..1000000 { s += (0.1 + i as f64 * 1e-7f64).tanh(); }
    println!("tanh {:.1} ns ({})", t.elapsed().as_secs_f64() / 1e6 * 1e9, s);
    let t = Instant::now();
    let mut s = 0.0; for i in 0..1000000 { s += (0.1 + i as f64 * 1e-7f64).exp(); }
    println!("exp {:.1} ns ({})", t.elapsed().as_secs_f64() / 1e6 * 1e9, s);
    let t = Instant::now();
    let mut s = 0.0; for i in 0..1000000 { s += (0.1 + i as f64 * 1e-7f64).atan2(0.3); }
    println!("atan2 {:.1} ns ({})", t.elapsed().as_secs_f64() / 1e6 * 1e9, s);
}
