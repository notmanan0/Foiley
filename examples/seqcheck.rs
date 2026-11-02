use foiley::sweep::{run_sequential, sequence_specs};
use foiley::xfoil::XFoil;
use std::time::Instant;
fn main() {
    foiley::set_quiet(true);
    for par in [true, false] {
        let mut xf = XFoil::new();
        xf.par_inner = par;
        xf.verbose = std::env::var_os("VV").is_some();
        xf.naca(12);
        xf.lvisc = true; xf.reinf1 = 1e6; let mut r = 0.0; xf.mrcl(1.0, &mut r);
        xf.itmax = 100;
        let specs = sequence_specs(true, 0.0, 16.0, 0.5);
        let t0 = Instant::now();
        let p = run_sequential(&mut xf, true, &specs, 105, |_, _| {});
        println!("par_inner {} : {:.1} ms  pts {}", par, t0.elapsed().as_secs_f64() * 1e3, p.len());
    }
}
