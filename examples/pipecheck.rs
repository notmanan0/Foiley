//! Check that pipelined sequential sweeps are bit-identical to the plain loop.
use foiley::sweep::{run_sequential, sequence_specs};
use foiley::xfoil::XFoil;
use std::time::Instant;

fn setup(naca: i32, re: f64) -> XFoil {
    let mut xf = XFoil::new();
    xf.naca(naca);
    xf.lvisc = true;
    xf.reinf1 = re;
    let mut r = 0.0;
    xf.mrcl(1.0, &mut r);
    xf.itmax = 100;
    xf
}

fn main() {
    foiley::set_quiet(true);
    let pipe = std::env::var("FOILEY_PIPELINE").map(|v| v != "0").unwrap_or(true);
    let mut tot = 0.0;
    let mut sig = 0u64;
    for (naca, re, a, b, d) in [(12, 1e6, 0.0, 16.0, 0.5), (2412, 3e5, -5.0, 14.0, 0.25), (4412, 1e6, -10.0, 18.0, 0.25)] {
        let mut xf = setup(naca, re);
        let specs = sequence_specs(true, a, b, d);
        let t0 = Instant::now();
        let p = run_sequential(&mut xf, true, &specs, 105, |_, _| {});
        tot += t0.elapsed().as_secs_f64();
        for q in &p {
            for v in [q.point.cl, q.point.cd, q.point.cm, q.point.cdp, q.point.xtr_top, q.point.xtr_bot] {
                sig = sig.wrapping_mul(1099511628211).wrapping_add(v.to_bits());
            }
        }
    }
    println!("pipeline={} time {:.1} ms  signature {:016x}", pipe, tot * 1e3, sig);
}
