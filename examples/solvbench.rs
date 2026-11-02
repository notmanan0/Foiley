use foiley::xfoil::{XFoil, DTOR};
use std::time::Instant;
fn main() {
    let mut xf = XFoil::new();
    xf.naca(2412);
    xf.lvisc = true; xf.reinf1 = 3e5; let mut r = 0.0; xf.mrcl(1.0, &mut r);
    xf.alfa = 4.0 * DTOR; xf.lalfa = true;
    xf.specal();
    xf.viscal(3);
    xf.setbl();
    let base = xf.clone();
    {
        let mut a = base.clone(); a.blsolv_ref();
        let mut b = base.clone(); b.blsolv();
        let mut nd = 0; let mut mx: f64 = 0.0;
        for iv in 1..=base.nsys { for k in 0..3 { for j in 0..2 {
            let (x, y) = (a.vdel[iv][k][j], b.vdel[iv][k][j]);
            if x.to_bits() != y.to_bits() { nd += 1; mx = mx.max((x-y).abs()/(x.abs()+1e-300)); }
        }}}
        println!("bitwise mismatches: {} / {}  max rel diff {:e}", nd, base.nsys*6, mx);
    }
    let reps = 300;
    for (name, which) in [("ref", 0), ("blocked", 1), ("blocked-1thread", 2)] {
        let mut best = 1e9;
        let mut out = None;
        for _ in 0..reps {
            let mut x = base.clone();
            let t = Instant::now();
            match which {
                0 => x.blsolv_ref(),
                1 => x.blsolv(),
                _ => { let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap(); pool.install(|| x.blsolv()); }
            }
            let e = t.elapsed().as_secs_f64();
            if e < best { best = e; }
            out = Some(x.vdel.clone());
        }
        let v = out.unwrap();
        println!("{:16} {:8.1} us   vdel[1]={:?}", name, best * 1e6, v[1][0]);
    }
    // timing of setbl
    let mut best = 1e9;
    for _ in 0..50 { let mut x = base.clone(); let t = Instant::now(); x.setbl(); best = f64::min(best, t.elapsed().as_secs_f64()); }
    println!("setbl {:8.1} us nsys {}", best * 1e6, base.nsys);
    let mut best = 1e9;
    for _ in 0..50 { let mut x = base.clone(); x.blsolv(); let t = Instant::now(); x.update(); best = f64::min(best, t.elapsed().as_secs_f64()); }
    println!("update {:8.1} us", best * 1e6);
}
