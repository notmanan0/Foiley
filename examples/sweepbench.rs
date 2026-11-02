//! Compare sequential (XFoil-exact) and parallel polar sweeps.
use foiley::sweep::{run_parallel, run_sequential, sequence_specs};
use foiley::xfoil::XFoil;
use std::time::Instant;

fn setup(naca: i32, re: f64, mach: f64) -> XFoil {
    let mut xf = XFoil::new();
    xf.naca(naca);
    xf.lvisc = true;
    xf.reinf1 = re;
    xf.minf1 = mach;
    let mut r = 0.0;
    xf.mrcl(1.0, &mut r);
    xf.comset();
    xf.itmax = 100;
    xf
}

fn main() {
    foiley::set_quiet(true);
    let threads: usize = std::env::args().nth(1).map(|s| s.parse().unwrap()).unwrap_or(8);
    let cases: Vec<(i32, f64, f64, bool, f64, f64, f64)> = vec![
        (12, 1e6, 0.0, true, 0.0, 16.0, 0.5),
        (2412, 3e5, 0.0, true, -5.0, 14.0, 0.25),
        (4415, 6e6, 0.3, true, -4.0, 12.0, 0.5),
        (23012, 2e6, 0.0, false, 0.0, 1.2, 0.1),
        (4412, 1e6, 0.0, true, -10.0, 18.0, 0.25),
    ];
    let (mut ts, mut tp) = (0.0, 0.0);
    for (naca, re, m, lalfa, a, b, d) in cases {
        let specs = sequence_specs(lalfa, a, b, d);
        let mut x1 = setup(naca, re, m);
        let t0 = Instant::now();
        let s = run_sequential(&mut x1, lalfa, &specs, 105, |_, _| {});
        let t_seq = t0.elapsed().as_secs_f64();
        let mut x2 = setup(naca, re, m);
        let t0 = Instant::now();
        let p = run_parallel(&mut x2, lalfa, &specs, 105, threads, |_, _| {});
        let t_par = t0.elapsed().as_secs_f64();
        ts += t_seq;
        tp += t_par;
        let mut mx = [0.0f64; 3];
        let mut nboth = 0;
        let mut nmis = 0;
        for (u, v) in s.iter().zip(p.iter()) {
            if u.converged && v.converged {
                nboth += 1;
                mx[0] = mx[0].max((u.point.cl - v.point.cl).abs());
                mx[1] = mx[1].max((u.point.cd - v.point.cd).abs());
                mx[2] = mx[2].max((u.point.cm - v.point.cm).abs());
            } else if u.converged != v.converged {
                nmis += 1;
            }
        }
        println!(
            "NACA {:5} Re {:8.0} M {:.2} pts {:3}/{:3}  seq {:7.1} ms  par {:7.1} ms  x{:4.1}  conv-both {:3} conv-mismatch {:2}  max|dCL| {:.5} |dCD| {:.6} |dCM| {:.5}",
            naca, re, m, s.len(), p.len(), t_seq * 1e3, t_par * 1e3, t_seq / t_par, nboth, nmis, mx[0], mx[1], mx[2]
        );
    }
    println!("TOTAL seq {:.1} ms par {:.1} ms  speedup {:.2}x", ts * 1e3, tp * 1e3, ts / tp);
}
