//! Phase-level profiler: replicates OPER/ASEQ + VISCAL with timers.
use foiley::xfoil::{XFoil, DTOR};
use std::time::Instant;

fn main() {
    foiley::set_quiet(true);
    let args: Vec<String> = std::env::args().collect();
    let naca: i32 = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(2412);
    let re: f64 = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(3e5);
    let mut t = [0f64; 12];
    let names = ["ggcalc+specal", "xywake", "qwcalc", "qdcalc", "setbl", "blsolv", "update", "rest", "mrchdu(in setbl)", "ueset", "", ""];
    let t0 = Instant::now();
    let mut xf = XFoil::new();
    xf.naca(naca);
    xf.lvisc = true;
    xf.reinf1 = re;
    let mut r = 0.0;
    xf.mrcl(1.0, &mut r);
    let mut iters = 0;
    let mut npts = 0;
    let mut a = -5.0;
    while a <= 14.0001 {
        xf.alfa = a * DTOR;
        xf.lalfa = true;
        if (xf.alfa - xf.awake).abs() > 1e-5 { xf.lwake = false; }
        if (xf.alfa - xf.avisc).abs() > 1e-5 { xf.lvconv = false; }
        let s = Instant::now();
        xf.specal();
        t[0] += s.elapsed().as_secs_f64();
        let s = Instant::now();
        if !xf.lwake { xf.xywake(); }
        t[1] += s.elapsed().as_secs_f64();
        let s = Instant::now();
        xf.qwcalc(); xf.qiset();
        if !xf.lipan { if xf.lblini { xf.gamqv(); } xf.stfind(); xf.iblpan(); xf.xicalc(); xf.iblsys(); }
        xf.uicalc();
        if !xf.lblini { for is in 1..=2 { for ibl in 1..=xf.nbl[is] { xf.uedg[is][ibl] = xf.uinv[is][ibl]; } } }
        t[2] += s.elapsed().as_secs_f64();
        let s = Instant::now();
        if !xf.lwdij || !xf.ladij { xf.qdcalc(); }
        t[3] += s.elapsed().as_secs_f64();
        for _ in 0..(xf.itmax + 5) {
            iters += 1;
            let s = Instant::now(); xf.setbl(); t[4] += s.elapsed().as_secs_f64();
            let s = Instant::now(); xf.blsolv(); t[5] += s.elapsed().as_secs_f64();
            let s = Instant::now(); xf.update(); t[6] += s.elapsed().as_secs_f64();
            let s = Instant::now();
            let cl = xf.cl; let mut rr = 0.0; xf.mrcl(cl, &mut rr); xf.comset();
            xf.qvfue(); xf.gamqv(); xf.stmove(); xf.clcalc(); xf.cdcalc();
            t[7] += s.elapsed().as_secs_f64();
            if xf.rmsbl < 1e-4 { xf.lvconv = true; xf.avisc = xf.alfa; xf.mvisc = xf.minf; break; }
        }
        npts += 1;
        a += 0.25;
    }
    // micro: mrchdu and ueset alone
    let s = Instant::now(); for _ in 0..50 { xf.mrchdu(); } t[8] = s.elapsed().as_secs_f64() / 50.0 * iters as f64;
    let s = Instant::now(); for _ in 0..50 { xf.ueset(); } t[9] = s.elapsed().as_secs_f64() / 50.0 * iters as f64;
    let tot = t0.elapsed().as_secs_f64();
    println!("points {} newton iters {} total {:.1} ms  nsys {}", npts, iters, tot * 1e3, xf.nsys);
    for i in 0..10 { println!("  {:18} {:8.1} ms  {:5.1}%", names[i], t[i] * 1e3, 100.0 * t[i] / tot); }
}
#[allow(dead_code)]
fn _x() {}
