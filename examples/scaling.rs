#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use foiley::sweep::{run_sequential, sequence_specs};
use foiley::xfoil::XFoil;
use std::time::Instant;
fn one() {
    let mut xf = XFoil::new();
    xf.par_inner = false;
    xf.naca(2412);
    xf.lvisc = true; xf.reinf1 = 3e5; let mut r = 0.0; xf.mrcl(1.0, &mut r);
    let specs = sequence_specs(true, -5.0, 14.0, 0.25);
    run_sequential(&mut xf, true, &specs, 25, |_, _| {});
}
fn main() {
    foiley::set_quiet(true);
    for nt in [1usize, 2, 4, 8, 12, 16] {
        let t0 = Instant::now();
        std::thread::scope(|s| { for _ in 0..nt { s.spawn(one); } });
        let t = t0.elapsed().as_secs_f64();
        println!("{:2} threads: {:7.1} ms  (per-polar {:6.1} ms, throughput x{:.2})", nt, t * 1e3, t * 1e3, nt as f64 / t);
    }
}
