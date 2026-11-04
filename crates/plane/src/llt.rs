//! Nonlinear lifting-line theory (main wing only).
//!
//! Stations are cosine-spaced (y = b/2 cos theta).  The induced angle is
//! iterated with under-relaxation using the 2D polars (lift from the section
//! polar at the effective angle of attack), so stall and nonlinear lift are
//! captured.  Planform only: sweep and dihedral are ignored, as in xflr5.

use crate::analysis::{Ctx, PlaneOpp, PolarKind, StripResult, G};
use crate::geom::Wing;
use crate::PlaneError;
use foiley_core::{Progress, Vec3};
use std::f64::consts::PI;

struct Station {
    y: f64,
    chord: f64,
    x_c4: f64,
    z: f64,
    twist: f64,
    foil_a: String,
    foil_b: String,
    tau: f64,
    width: f64,
}

fn stations(w: &Wing, n: usize) -> Vec<Station> {
    let b = w.planform_span();
    let half = 0.5 * b;
    let mut out = Vec::new();
    for k in 1..n {
        let th = k as f64 * PI / n as f64;
        let y = half * th.cos();
        let ya = y.abs();
        let left = y < 0.0;
        // find section bracket
        let s = &w.sections;
        let mut i = 0;
        while i + 2 < s.len() && s[i + 1].y < ya {
            i += 1;
        }
        let (a, bb) = (&s[i], &s[i + 1]);
        let t = if bb.y > a.y { ((ya - a.y) / (bb.y - a.y)).clamp(0.0, 1.0) } else { 0.0 };
        let chord = a.chord + t * (bb.chord - a.chord);
        let off = a.offset + t * (bb.offset - a.offset);
        let twist = a.twist + t * (bb.twist - a.twist) + w.ry;
        // dihedral height (approximate, developed)
        let mut z = 0.0;
        for j in 0..i {
            z += (s[j + 1].y - s[j].y) * (s[j].dihedral.to_radians()).sin();
        }
        z += (ya - a.y).max(0.0) * a.dihedral.to_radians().sin();
        let w_th = half * (PI / n as f64) * th.sin();
        out.push(Station {
            y,
            chord,
            x_c4: w.position.x + off + 0.25 * chord,
            z: w.position.z + z,
            twist,
            foil_a: w.foil(i, left).to_string(),
            foil_b: w.foil(i + 1, left).to_string(),
            tau: t,
            width: w_th,
        });
    }
    out
}

/// Run an LLT polar (types 1-4).
pub fn run(ctx: &Ctx, prog: &Progress, on_point: &(dyn Fn(&PlaneOpp) + Sync), failed: &mut Vec<String>) -> Result<Vec<PlaneOpp>, PlaneError> {
    let def = ctx.def;
    let wing = ctx.plane.main_wing().ok_or_else(|| PlaneError::Geometry("no main wing".into()))?;
    let n = def.llt.stations.clamp(8, 400);
    let st = stations(wing, n);
    let b = wing.planform_span();
    let s_w = wing.planform_area();
    if b <= 0.0 || s_w <= 0.0 {
        return Err(PlaneError::Geometry("main wing has zero span or area".into()));
    }
    for s in &st {
        for f in [&s.foil_a, &s.foil_b] {
            if ctx.db.re_range(f).is_none() {
                return Err(PlaneError::Input(format!("LLT needs type-1 viscous 2D polars for foil '{f}' (run a batch polar mesh first)")));
            }
        }
    }
    let m = st.len();
    let theta: Vec<f64> = (1..n).map(|k| k as f64 * PI / n as f64).collect();
    // influence (degrees per unit Cl*c/b)
    let mut beta = vec![vec![0.0; m]; m];
    for (k, row) in beta.iter_mut().enumerate() {
        let kk = k + 1;
        for (j, v) in row.iter_mut().enumerate() {
            let mm = j + 1;
            *v = if mm == kk {
                180.0 * n as f64 / (8.0 * PI * theta[k].sin())
            } else if (mm + kk) % 2 == 0 {
                0.0
            } else {
                let c = 180.0 / (4.0 * PI * n as f64 * theta[k].sin());
                c * (1.0 / (1.0 - (((kk + mm) as f64) * PI / n as f64).cos()) - 1.0 / (1.0 - (((kk as f64) - (mm as f64)) * PI / n as f64).cos()))
            };
        }
    }
    let rho = def.fluid.rho;
    let nu = def.fluid.nu();
    let mass = ctx.inertia.mass;
    let cog = if mass > 0.0 { ctx.inertia.cog } else { Vec3::ZERO };
    let ops: Vec<(f64, f64)> = match def.kind {
        PolarKind::FixedAlpha => def.range.values().into_iter().map(|v| (def.alpha, v)).collect(),
        _ => def.range.values().into_iter().map(|a| (a, def.velocity)).collect(),
    };
    let interp = |s: &Station, re: f64, a: f64| -> Option<(f64, f64, f64, f64, f64, bool)> {
        let (pa, oa) = ctx.db.at_re_alpha(&s.foil_a, re, a)?;
        let (pb, ob) = ctx.db.at_re_alpha(&s.foil_b, re, a)?;
        let l = |x: f64, y: f64| x + s.tau * (y - x);
        Some((l(pa.cl, pb.cl), l(pa.cd, pb.cd), l(pa.cm, pb.cm), l(pa.xtr_top, pb.xtr_top), l(pa.xtr_bot, pb.xtr_bot), oa || ob))
    };
    let mut ai = vec![0.0; m];
    let mut out = Vec::new();
    // Per-station relaxation from the local feedback gain d(ai)/d(ai) =
    // beta_kk c_k / b * dCl/dalpha: the fixed xflr5 factor (1/20) diverges
    // near the tips when many stations are used; this one does not.
    let omega: Vec<f64> = (0..m)
        .map(|k| {
            let re = def.velocity.max(1.0) * st[k].chord / nu;
            let slope = ctx.db.linear_lift(&st[k].foil_a, re).map(|v| v.1).unwrap_or(0.1).clamp(0.02, 0.2);
            let lam = beta[k][k] * st[k].chord / b * slope;
            1.0 / (def.llt.relax.max(1.0) * (1.0 + lam))
        })
        .collect();
    let total = ops.len();
    let ar = b * b / s_w;
    for (ip, &(alpha, v_spec)) in ops.iter().enumerate() {
        prog.check()?;
        let mut v = v_spec;
        if matches!(def.kind, PolarKind::FixedLift | PolarKind::Glide) {
            if mass <= 0.0 {
                return Err(PlaneError::Input("T2/T3 need the plane mass".into()));
            }
            v = (2.0 * mass * G / (rho * s_w * 0.5)).sqrt();
        }
        let mut converged = false;
        let mut cl = vec![0.0; m];
        let mut extra = vec![(0.0, 0.0, 1.0, 1.0, false); m];
        let mut ok = true;
        for _it in 0..def.llt.max_iter {
            let mut maxd: f64 = 0.0;
            let mut ai_new = vec![0.0; m];
            for k in 0..m {
                let mut s = 0.0;
                for j in 0..m {
                    s += beta[k][j] * cl[j] * st[j].chord / b;
                }
                ai_new[k] = -s;
            }
            for k in 0..m {
                maxd = maxd.max((ai_new[k] - ai[k]).abs());
                ai[k] += (ai_new[k] - ai[k]) * omega[k];
            }
            for k in 0..m {
                let re = v * st[k].chord / nu;
                match interp(&st[k], re, alpha + st[k].twist + ai[k]) {
                    Some((c, cd, cm, xt, xb, o)) => {
                        cl[k] = c;
                        extra[k] = (cd, cm, xt, xb, o);
                    }
                    None => ok = false,
                }
            }
            if !ok {
                break;
            }
            if matches!(def.kind, PolarKind::FixedLift | PolarKind::Glide) {
                let clw: f64 = (0..m).map(|k| cl[k] * st[k].chord * st[k].width).sum::<f64>() / s_w;
                if clw <= 1e-6 {
                    ok = false;
                    break;
                }
                v = (2.0 * mass * G / (rho * s_w * clw)).sqrt();
            }
            if maxd < def.llt.tolerance && _it > 2 {
                converged = true;
                break;
            }
        }
        prog.set((ip + 1) as f64 / total as f64);
        if !ok || !converged {
            prog.count_failed();
            failed.push(format!("alpha={alpha:.2} V={v:.2}: LLT {}", if ok { "did not converge (stall?)" } else { "left the 2D polar range" }));
            ai.iter_mut().for_each(|a| *a = 0.0);
            continue;
        }
        let q = 0.5 * rho * v * v;
        let mut l_tot = 0.0;
        let mut di = 0.0;
        let mut dv = 0.0;
        let mut my = 0.0;
        let mut strips = Vec::with_capacity(m);
        let mut n_out = 0;
        for k in 0..m {
            let s = &st[k];
            let (cd, cm, xt, xb, o) = extra[k];
            if o {
                n_out += 1;
            }
            let dy = s.width;
            let lk = q * cl[k] * s.chord * dy;
            let cdi = -cl[k] * ai[k].to_radians();
            l_tot += lk;
            di += q * cdi * s.chord * dy;
            dv += q * cd * s.chord * dy;
            my += q * s.chord * s.chord * dy * cm - (s.x_c4 - cog.x) * lk + (s.z - cog.z) * q * (cd + cdi) * s.chord * dy;
            strips.push(StripResult {
                wing: 0,
                span_pos: s.y,
                chord: s.chord,
                cl: cl[k],
                cdi,
                cdv: cd,
                ai: ai[k],
                re: v * s.chord / nu,
                xtr_top: xt,
                xtr_bot: xb,
                gamma: 0.5 * v * s.chord * cl[k],
                lift_per_span: q * cl[k] * s.chord,
                twist: s.twist,
                alpha_eff: alpha + s.twist + ai[k],
                cm_c4: cm,
                bending: 0.0,
                out_of_polar: o,
            });
        }
        strips.reverse();
        // bending: integrate outboard of each station on its side
        for k in 0..m {
            let yk = strips[k].span_pos;
            let mut bm = 0.0;
            for j in 0..m {
                let yj = strips[j].span_pos;
                if (yk >= 0.0 && yj > yk) || (yk < 0.0 && yj < yk) {
                    bm += strips[j].lift_per_span * st[m - 1 - j].width * (yj - yk).abs();
                }
            }
            strips[k].bending = bm;
        }
        let qs = q * ctx.s_ref;
        let mut d_extra: f64 = def.extra_drag.iter().map(|e| q * e.area * e.cd).sum();
        if !def.viscous {
            dv = 0.0;
        }
        d_extra += dv;
        let cl_t = l_tot / qs;
        let cdi_t = di / qs;
        let cdv_t = d_extra / qs;
        let mut warnings = vec!["LLT: main wing only; sweep and dihedral ignored".to_string()];
        if n_out > 0 {
            warnings.push(format!("{n_out} station(s) outside the 2D polar envelope"));
        }
        let opp = PlaneOpp {
            alpha,
            beta: 0.0,
            velocity: v,
            ctrl: 0.0,
            cl: cl_t,
            cd: cdi_t + cdv_t,
            cdi: cdi_t,
            cdv: cdv_t,
            cy: 0.0,
            cm: my / (qs * ctx.c_ref),
            cr: 0.0,
            cn: 0.0,
            lift: l_tot,
            drag: (cdi_t + cdv_t) * qs,
            x_cp: if l_tot.abs() > 1e-9 { cog.x - my / l_tot } else { cog.x },
            efficiency: if cdi_t.abs() > 1e-12 { cl_t * cl_t / (PI * ar * cdi_t) } else { 0.0 },
            mass,
            cog,
            root_bending: strips.iter().map(|s| s.bending.abs()).fold(0.0, f64::max),
            strips,
            cp: vec![],
            warnings,
            stability: None,
            wake: vec![],
            balance_residual: None,
        };
        prog.count_ok();
        on_point(&opp);
        out.push(opp);
    }
    Ok(out)
}
