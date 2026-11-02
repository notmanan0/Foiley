//! Operating-point sequences (ASEQ / CSEQ), sequential and parallel.
//!
//! The sequential path reproduces XFoil exactly.  The parallel path splits
//! the sequence into contiguous chunks that are solved concurrently:
//!
//! * chunk 0 continues from the caller's state, exactly like XFoil;
//! * every other chunk starts from a cold boundary-layer initialisation at
//!   its first point (what XFoil does for the first point of a session);
//! * if a chunk's first point fails to converge, the chunk is re-solved
//!   sequentially starting from the previous chunk's final state, i.e. along
//!   XFoil's own continuation path.
//!
//! Converged points satisfy the same discrete equations to the same Newton
//! tolerance, so they agree with the sequential results to within that
//! tolerance (typically identical to all printed digits).

use crate::polar::PolarPoint;
use crate::xfoil::{WakePack, XFoil, DTOR};

/// Result of one point of a sequence.
#[derive(Clone, Copy, Debug)]
pub struct SeqPoint {
    pub point: PolarPoint,
    pub converged: bool,
}

/// Snapshot the current operating point as a polar point (PLRADD semantics).
pub fn polar_point(xf: &XFoil) -> PolarPoint {
    let (cd, re) = if xf.lvisc { (xf.cd, xf.reinf) } else { (0.0, 0.0) };
    PolarPoint {
        alpha: xf.adeg,
        cl: xf.cl,
        cd,
        cdp: xf.cdp,
        cm: xf.cm,
        cdv: cd,
        mach: xf.minf,
        re,
        ncrit: xf.acrit,
        xtr_top: if xf.lvisc { xf.xoctr[1] } else { 0.0 },
        xtr_bot: if xf.lvisc { xf.xoctr[2] } else { 0.0 },
        xtrip_top: xf.xstrip[1],
        xtrip_bot: xf.xstrip[2],
        chinge: xf.hmom,
        cpmin: xf.cpmn,
        converged: xf.lvconv || !xf.lvisc,
    }
}

/// Solve one point of an ASEQ/CSEQ sequence (body of the OPER sequence loop).
pub fn solve_point(xf: &mut XFoil, lalfa: bool, spec: f64, itmaxs: usize) -> SeqPoint {
    solve_point_with(xf, lalfa, spec, itmaxs, None)
}

/// `solve_point` with optionally precomputed wake data for this alpha.
pub fn solve_point_with(xf: &mut XFoil, lalfa: bool, spec: f64, itmaxs: usize, pack: Option<WakePack>) -> SeqPoint {
    if let Some(p) = pack {
        xf.install_wake(p);
    }
    xf.lalfa = lalfa;
    xf.qinf = 1.0;
    if lalfa {
        xf.alfa = spec;
    } else {
        xf.clspec = spec;
        xf.speccl();
    }
    if (xf.alfa - xf.awake).abs() > 1.0e-5 {
        xf.lwake = false;
    }
    if (xf.alfa - xf.avisc).abs() > 1.0e-5 {
        xf.lvconv = false;
    }
    if (xf.minf - xf.mvisc).abs() > 1.0e-5 {
        xf.lvconv = false;
    }
    xf.specal();
    if xf.lvisc {
        xf.viscal(itmaxs);
    }
    xf.adeg = xf.alfa / DTOR;
    xf.fcpmin();
    let converged = !xf.lvisc || xf.lvconv;
    SeqPoint { point: polar_point(xf), converged }
}

/// Specification values (alpha in radians, or CL) for a sequence.
pub fn sequence_specs(lalfa: bool, first: f64, last: f64, step: f64) -> Vec<f64> {
    let step = if last < first { -step.abs() } else { step.abs() };
    if lalfa {
        let (a1, a2, da) = (first * DTOR, last * DTOR, step * DTOR);
        let np = if da != 0.0 { ((a2 - a1) / da + 0.5) as i64 + 1 } else { 1 };
        (0..np.max(0)).map(|i| a1 + da * i as f64).collect()
    } else {
        let np = if step != 0.0 { ((last - first) / step + 0.5) as i64 + 1 } else { 1 };
        (0..np.max(0)).map(|i| first + step * i as f64).collect()
    }
}

/// Run a sequence sequentially, calling `each` after every point (for polar
/// accumulation / printing).  Stops after `nseqex` consecutive failures,
/// like XFoil.  Returns the points computed.
pub fn run_sequential(
    xf: &mut XFoil,
    lalfa: bool,
    specs: &[f64],
    itmaxs: usize,
    mut each: impl FnMut(&mut XFoil, &SeqPoint),
) -> Vec<SeqPoint> {
    // Alpha sweeps: precompute each next point's wake + source-influence
    // matrix on a helper thread (it depends only on alpha and geometry), so
    // the main thread only runs the boundary-layer Newton iterations.
    // Results are bit-identical to the plain sequential loop.
    let pipeline = lalfa && xf.lvisc && specs.len() > 1 && pipeline_enabled();
    if !pipeline {
        return run_seq_inner(xf, lalfa, specs, itmaxs, &mut each, None);
    }
    if !xf.lgamu || !xf.lqaij {
        let v = xf.verbose;
        xf.verbose = false;
        xf.ggcalc();
        xf.verbose = v;
    }
    let mut twin = xf.clone();
    twin.verbose = false;
    twin.par_inner = false;
    let (tx, rx) = std::sync::mpsc::sync_channel::<WakePack>(2);
    let rest: Vec<f64> = specs[1..].to_vec();
    std::thread::scope(|sc| {
        sc.spawn(move || {
            crate::set_quiet(true);
            for a in rest {
                twin.alfa = a;
                twin.lalfa = true;
                twin.qinf = 1.0;
                if (twin.alfa - twin.awake).abs() > 1.0e-5 {
                    twin.lwake = false;
                }
                twin.specal();
                if !twin.lwake {
                    twin.xywake();
                }
                twin.qwcalc();
                if !twin.lwdij || !twin.ladij {
                    twin.qdcalc();
                }
                if tx.send(twin.wake_pack()).is_err() {
                    break;
                }
            }
        });
        let r = run_seq_inner(xf, lalfa, specs, itmaxs, &mut each, Some(&rx));
        drop(rx);
        r
    })
}

/// Whether the bit-exact helper-thread pipelining is enabled
/// (set FOILEY_PIPELINE=0 to disable).
pub fn pipeline_enabled() -> bool {
    static F: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *F.get_or_init(|| {
        std::env::var("FOILEY_PIPELINE").map(|v| v != "0").unwrap_or(true)
            && std::thread::available_parallelism().map(|n| n.get() > 1).unwrap_or(false)
    })
}

fn run_seq_inner(
    xf: &mut XFoil,
    lalfa: bool,
    specs: &[f64],
    itmaxs: usize,
    each: &mut impl FnMut(&mut XFoil, &SeqPoint),
    packs: Option<&std::sync::mpsc::Receiver<WakePack>>,
) -> Vec<SeqPoint> {
    let mut out = Vec::with_capacity(specs.len());
    let mut iseqex = 0;
    let mut alast = xf.adeg;
    let mut clast = xf.cl;
    for (k, &sp) in specs.iter().enumerate() {
        let pack = match packs {
            Some(rx) if k > 0 => rx.recv().ok(),
            _ => None,
        };
        let p = solve_point_with(xf, lalfa, sp, itmaxs, pack);
        each(xf, &p);
        out.push(p);
        if xf.lvisc && !p.converged {
            iseqex += 1;
            if iseqex >= xf.nseqex {
                if xf.verbose {
                    crate::outln!(
                        "\n Sequence halted since previous{:3} points did not converge\n Last-converged  alpha ={:8.3}    CL ={:10.5}",
                        iseqex,
                        alast,
                        clast
                    );
                }
                break;
            }
        } else {
            iseqex = 0;
            alast = xf.adeg;
            clast = xf.cl;
        }
    }
    out
}

/// Number of worker chunks for a parallel sequence of `npts` points.
fn chunk_count(npts: usize, threads: usize) -> usize {
    // each extra chunk pays for a cold BL start (~2-3 point solves), so keep
    // chunks reasonably long
    let max_by_len = (npts / 3).max(1);
    threads.min(max_by_len).max(1)
}

/// Run a viscous sequence in parallel chunks (see module docs).  On return
/// `xf` holds the state after the last computed point.  The `each` callback
/// is invoked afterwards for every point in sequence order (so polar files
/// are written in order).
pub fn run_parallel(
    xf: &mut XFoil,
    lalfa: bool,
    specs: &[f64],
    itmaxs: usize,
    threads: usize,
    mut each: impl FnMut(&mut XFoil, &SeqPoint),
) -> Vec<SeqPoint> {
    let npts = specs.len();
    let nch = chunk_count(npts, threads);
    if !xf.lvisc || nch <= 1 {
        return run_sequential(xf, lalfa, specs, itmaxs, each);
    }
    // contiguous chunk boundaries
    let bounds: Vec<(usize, usize)> = (0..nch).map(|c| (c * npts / nch, (c + 1) * npts / nch)).collect();

    let verbose = xf.verbose;
    // Build the alpha-independent influence matrices once (bit-identical to
    // what the first point would compute) so chunks don't each redo them.
    {
        let v = xf.verbose;
        xf.verbose = false;
        if !xf.lgamu || !xf.lqaij {
            xf.ggcalc();
        }
        if lalfa && (!xf.ladij || !xf.lwake || (xf.alfa - specs[0]).abs() > 1.0e-5) {
            xf.alfa = specs[0];
            xf.lalfa = true;
            xf.qinf = 1.0;
            if (xf.alfa - xf.awake).abs() > 1.0e-5 {
                xf.lwake = false;
            }
            xf.specal();
            if !xf.lwake {
                xf.xywake();
            }
            if !xf.lwdij || !xf.ladij {
                xf.qdcalc();
            }
        }
        xf.verbose = v;
    }
    let base = {
        let mut b = xf.clone();
        b.verbose = false;
        b
    };
    let solve_chunk = |c: usize, start: &XFoil, cold: bool| -> (XFoil, Vec<SeqPoint>) {
        let mut w = start.clone();
        w.verbose = false;
        w.par_inner = false;
        if cold {
            w.lblini = false;
            w.lipan = false;
            w.lvconv = false;
        }
        let (a, b) = bounds[c];
        let pts = if chunk_pipeline() {
            // pipelined (wake precompute on a helper thread); halts like XFoil
            run_sequential(&mut w, lalfa, &specs[a..b], itmaxs, |_, _| {})
        } else {
            let mut pts = Vec::with_capacity(b - a);
            let mut nfail = 0;
            for &sp in &specs[a..b] {
                let p = solve_point(&mut w, lalfa, sp, itmaxs);
                nfail = if p.converged { 0 } else { nfail + 1 };
                pts.push(p);
                if nfail >= w.nseqex && !(cold && pts.len() == nfail) {
                    break;
                }
            }
            pts
        };
        (w, pts)
    };

    // dedicated worker threads pulling chunks from a shared counter (no
    // work-stealing interplay with nested parallel regions)
    let next = std::sync::atomic::AtomicUsize::new(0);
    let slots: Vec<std::sync::Mutex<Option<(XFoil, Vec<SeqPoint>)>>> = (0..nch).map(|_| std::sync::Mutex::new(None)).collect();
    std::thread::scope(|sc| {
        for _ in 0..threads.min(nch) {
            sc.spawn(|| loop {
                crate::set_quiet(true);
                let c = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if c >= nch {
                    break;
                }
                let r = solve_chunk(c, &base, c > 0);
                *slots[c].lock().unwrap() = Some(r);
            });
        }
    });
    let mut results: Vec<(XFoil, Vec<SeqPoint>)> = slots.into_iter().map(|m| m.into_inner().unwrap().unwrap()).collect();

    // repair chunks whose cold start failed: continue from the previous chunk
    for c in 1..nch {
        if !results[c].1[0].converged {
            let prev = results[c - 1].0.clone();
            let (mut st, pts) = solve_chunk(c, &prev, false);
            st.par_inner = true;
            results[c] = (st, pts);
        }
    }

    // emit in order, applying XFoil's "halt after NSEQEX failures" rule
    let mut out = Vec::with_capacity(npts);
    let mut iseqex = 0;
    let mut last_state: Option<usize> = None;
    'outer: for (c, (st, pts)) in results.iter_mut().enumerate() {
        for p in pts {
            each(st, p);
            out.push(*p);
            if !p.converged {
                iseqex += 1;
                if iseqex >= xf.nseqex {
                    if verbose {
                        crate::outln!(
                            "\n Sequence halted since previous{:3} points did not converge",
                            iseqex
                        );
                    }
                    last_state = Some(c);
                    break 'outer;
                }
            } else {
                iseqex = 0;
            }
        }
        last_state = Some(c);
    }
    if let Some(c) = last_state {
        let mut st = results.swap_remove(c).0;
        st.verbose = verbose;
        st.par_inner = true;
        *xf = st;
    }
    out
}

/// Number of parallel sweep jobs: `FOILEY_JOBS` if set, otherwise
/// min(8, half the logical CPUs) (hyperthreads and the L3 working set make
/// more threads counter-productive for this solver).
pub fn default_jobs() -> usize {
    if let Ok(v) = std::env::var("FOILEY_JOBS") {
        if let Ok(n) = v.trim().parse::<usize>() {
            return n.max(1);
        }
    }
    let n = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    (n / 2).clamp(1, 8)
}

fn chunk_pipeline() -> bool {
    static F: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *F.get_or_init(|| std::env::var("FOILEY_CHUNK_PIPE").map(|v| v != "0").unwrap_or(true))
}
