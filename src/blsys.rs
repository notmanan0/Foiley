//! Boundary-layer equation system: correlations, secondary variables and
//! the local 2-station Newton system (port of XFoil xblsys.f).
//!
//! The Fortran COMMON blocks V_VAR1/V_VAR2 ("1" and "2" station variables
//! with all their sensitivities) are represented by [`Station`]; the rest of
//! XBL.INC lives in [`BlSys`].  The VS1/VS2/... matrices use 0-based
//! indices: `vs1[k][l]` == `VS1(K+1,L+1)`.

/// BL calibration constants (BLPAR.INC).
#[derive(Clone, Copy, Debug)]
pub struct BlPar {
    pub sccon: f64,
    pub gacon: f64,
    pub gbcon: f64,
    pub gccon: f64,
    pub dlcon: f64,
    pub ctrcon: f64,
    pub ctrcex: f64,
    pub duxcon: f64,
    pub ctcon: f64,
    pub cffac: f64,
}

impl Default for BlPar {
    fn default() -> Self {
        let gacon = 6.70;
        let gbcon = 0.75;
        BlPar {
            sccon: 5.6,
            gacon,
            gbcon,
            gccon: 18.0,
            dlcon: 0.9,
            ctrcon: 1.8,
            ctrcex: 3.3,
            duxcon: 1.0,
            ctcon: 0.5 / (gacon * gacon * gbcon),
            cffac: 1.0,
        }
    }
}

/// Primary and secondary BL variables at one station, with sensitivities.
#[derive(Clone, Copy, Debug, Default)]
pub struct Station {
    pub x: f64,
    pub u: f64,
    pub t: f64,
    pub d: f64,
    pub s: f64,
    pub ampl: f64,
    pub u_uei: f64,
    pub u_ms: f64,
    pub dw: f64,
    pub h: f64,
    pub h_t: f64,
    pub h_d: f64,
    pub m: f64,
    pub m_u: f64,
    pub m_ms: f64,
    pub r: f64,
    pub r_u: f64,
    pub r_ms: f64,
    pub v: f64,
    pub v_u: f64,
    pub v_ms: f64,
    pub v_re: f64,
    pub hk: f64,
    pub hk_u: f64,
    pub hk_t: f64,
    pub hk_d: f64,
    pub hk_ms: f64,
    pub hs: f64,
    pub hs_u: f64,
    pub hs_t: f64,
    pub hs_d: f64,
    pub hs_ms: f64,
    pub hs_re: f64,
    pub hc: f64,
    pub hc_u: f64,
    pub hc_t: f64,
    pub hc_d: f64,
    pub hc_ms: f64,
    pub rt: f64,
    pub rt_u: f64,
    pub rt_t: f64,
    pub rt_ms: f64,
    pub rt_re: f64,
    pub cf: f64,
    pub cf_u: f64,
    pub cf_t: f64,
    pub cf_d: f64,
    pub cf_ms: f64,
    pub cf_re: f64,
    pub di: f64,
    pub di_u: f64,
    pub di_t: f64,
    pub di_d: f64,
    pub di_s: f64,
    pub di_ms: f64,
    pub di_re: f64,
    pub us: f64,
    pub us_u: f64,
    pub us_t: f64,
    pub us_d: f64,
    pub us_ms: f64,
    pub us_re: f64,
    pub cq: f64,
    pub cq_u: f64,
    pub cq_t: f64,
    pub cq_d: f64,
    pub cq_ms: f64,
    pub cq_re: f64,
    pub de: f64,
    pub de_u: f64,
    pub de_t: f64,
    pub de_d: f64,
    pub de_ms: f64,
}

/// Global BL state (XBL.INC) and the local 2-station system.
#[derive(Clone, Debug, Default)]
pub struct BlSys {
    pub c1: Station,
    pub c2: Station,
    c1sav: Station,
    c2sav: Station,

    // V_VARA
    pub cfm: f64,
    pub cfm_ms: f64,
    pub cfm_re: f64,
    pub cfm_u1: f64,
    pub cfm_t1: f64,
    pub cfm_d1: f64,
    pub cfm_u2: f64,
    pub cfm_t2: f64,
    pub cfm_d2: f64,
    pub xt: f64,
    pub xt_a1: f64,
    pub xt_ms: f64,
    pub xt_re: f64,
    pub xt_xf: f64,
    pub xt_x1: f64,
    pub xt_t1: f64,
    pub xt_d1: f64,
    pub xt_u1: f64,
    pub xt_x2: f64,
    pub xt_t2: f64,
    pub xt_d2: f64,
    pub xt_u2: f64,

    // V_VAR
    pub dwte: f64,
    pub qinfbl: f64,
    pub tkbl: f64,
    pub tkbl_ms: f64,
    pub rstbl: f64,
    pub rstbl_ms: f64,
    pub hstinv: f64,
    pub hstinv_ms: f64,
    pub reybl: f64,
    pub reybl_ms: f64,
    pub reybl_re: f64,
    pub gambl: f64,
    pub gm1bl: f64,
    pub hvrat: f64,
    pub bule: f64,
    pub xiforc: f64,
    pub amcrit: f64,

    // V_INT
    pub simi: bool,
    pub tran: bool,
    pub turb: bool,
    pub wake: bool,
    pub trforc: bool,
    pub trfree: bool,
    pub idampv: i32,

    // V_SYS
    pub vs1: [[f64; 5]; 4],
    pub vs2: [[f64; 5]; 4],
    pub vsrez: [f64; 4],
    pub vsr: [f64; 4],
    pub vsm: [f64; 4],
    pub vsx: [f64; 4],

    pub par: BlPar,
}

// ---------------------------------------------------------------------------
// Correlations (pure functions)
// ---------------------------------------------------------------------------

/// Kinematic shape parameter (Whitfield). Returns (hk, hk_h, hk_msq).
#[inline]
pub fn hkin(h: f64, msq: f64) -> (f64, f64, f64) {
    let den = 1.0 + 0.113 * msq;
    let hk = (h - 0.29 * msq) / den;
    let hk_h = 1.0 / den;
    let hk_msq = (-0.29 - 0.113 * hk) / den;
    (hk, hk_h, hk_msq)
}

/// Laminar dissipation function. Returns (di, di_hk, di_rt).
#[inline]
pub fn dil(hk: f64, rt: f64) -> (f64, f64, f64) {
    let (di, di_hk);
    if hk < 4.0 {
        di = (0.00205 * (4.0 - hk).powf(5.5) + 0.207) / rt;
        di_hk = (-0.00205 * 5.5 * (4.0 - hk).powf(4.5)) / rt;
    } else {
        let hkb = hk - 4.0;
        let den = 1.0 + 0.02 * hkb * hkb;
        di = (-0.0016 * hkb * hkb / den + 0.207) / rt;
        di_hk = (-0.0016 * 2.0 * hkb * (1.0 / den - 0.02 * hkb * hkb / (den * den))) / rt;
    }
    (di, di_hk, -di / rt)
}

/// Laminar wake dissipation function.
#[inline]
pub fn dilw(hk: f64, rt: f64) -> (f64, f64, f64) {
    let (hs, hs_hk, hs_rt, _) = hsl(hk, rt, 0.0);
    let rcd = 1.10 * (1.0 - 1.0 / hk).powi(2) / hk;
    let rcd_hk = -1.10 * (1.0 - 1.0 / hk) * 2.0 / (hk * hk * hk) - rcd / hk;
    let di = 2.0 * rcd / (hs * rt);
    let di_hk = 2.0 * rcd_hk / (hs * rt) - (di / hs) * hs_hk;
    let di_rt = -di / rt - (di / hs) * hs_rt;
    (di, di_hk, di_rt)
}

/// Laminar H* correlation. Returns (hs, hs_hk, hs_rt, hs_msq).
#[inline]
pub fn hsl(hk: f64, _rt: f64, _msq: f64) -> (f64, f64, f64, f64) {
    let (hs, hs_hk);
    if hk < 4.35 {
        let tmp = hk - 4.35;
        let hp1 = hk + 1.0;
        hs = 0.0111 * tmp * tmp / hp1 - 0.0278 * tmp * tmp * tmp / hp1 + 1.528
            - 0.0002 * (tmp * hk) * (tmp * hk);
        hs_hk = 0.0111 * (2.0 * tmp - tmp * tmp / hp1) / hp1
            - 0.0278 * (3.0 * tmp * tmp - tmp * tmp * tmp / hp1) / hp1
            - 0.0002 * 2.0 * tmp * hk * (tmp + hk);
    } else {
        hs = 0.015 * (hk - 4.35).powi(2) / hk + 1.528;
        hs_hk = 0.015 * 2.0 * (hk - 4.35) / hk - 0.015 * (hk - 4.35).powi(2) / (hk * hk);
    }
    (hs, hs_hk, 0.0, 0.0)
}

/// Laminar skin friction (Falkner-Skan). Returns (cf, cf_hk, cf_rt, cf_msq).
#[inline]
pub fn cfl(hk: f64, rt: f64, _msq: f64) -> (f64, f64, f64, f64) {
    let (cf, cf_hk);
    if hk < 5.5 {
        let tmp = (5.5 - hk).powi(3) / (hk + 1.0);
        cf = (0.0727 * tmp - 0.07) / rt;
        cf_hk = (-0.0727 * tmp * 3.0 / (5.5 - hk) - 0.0727 * tmp / (hk + 1.0)) / rt;
    } else {
        let tmp = 1.0 - 1.0 / (hk - 4.5);
        cf = (0.015 * tmp * tmp - 0.07) / rt;
        cf_hk = (0.015 * tmp * 2.0 / (hk - 4.5).powi(2)) / rt;
    }
    (cf, cf_hk, -cf / rt, 0.0)
}

/// Turbulent H* correlation. Returns (hs, hs_hk, hs_rt, hs_msq).
#[inline]
pub fn hst(hk: f64, rt: f64, msq: f64) -> (f64, f64, f64, f64) {
    const HSMIN: f64 = 1.500;
    const DHSINF: f64 = 0.015;
    let (ho, ho_rt) = if rt > 400.0 {
        (3.0 + 400.0 / rt, -400.0 / (rt * rt))
    } else {
        (4.0, 0.0)
    };
    let (rtz, rtz_rt) = if rt > 200.0 { (rt, 1.0) } else { (200.0, 0.0) };
    let (mut hs, mut hs_hk, mut hs_rt);
    if hk < ho {
        let hr = (ho - hk) / (ho - 1.0);
        let hr_hk = -1.0 / (ho - 1.0);
        let hr_rt = (1.0 - hr) / (ho - 1.0) * ho_rt;
        let a = 2.0 - HSMIN - 4.0 / rtz;
        hs = a * hr * hr * 1.5 / (hk + 0.5) + HSMIN + 4.0 / rtz;
        hs_hk = -a * hr * hr * 1.5 / ((hk + 0.5) * (hk + 0.5)) + a * hr * 2.0 * 1.5 / (hk + 0.5) * hr_hk;
        hs_rt = a * hr * 2.0 * 1.5 / (hk + 0.5) * hr_rt
            + (hr * hr * 1.5 / (hk + 0.5) - 1.0) * 4.0 / (rtz * rtz) * rtz_rt;
    } else {
        let grt = rtz.ln();
        let hdif = hk - ho;
        let rtmp = hk - ho + 4.0 / grt;
        let htmp = 0.007 * grt / (rtmp * rtmp) + DHSINF / hk;
        let htmp_hk = -0.014 * grt / (rtmp * rtmp * rtmp) - DHSINF / (hk * hk);
        let htmp_rt = -0.014 * grt / (rtmp * rtmp * rtmp) * (-ho_rt - 4.0 / (grt * grt) / rtz * rtz_rt)
            + 0.007 / (rtmp * rtmp) / rtz * rtz_rt;
        hs = hdif * hdif * htmp + HSMIN + 4.0 / rtz;
        hs_hk = hdif * 2.0 * htmp + hdif * hdif * htmp_hk;
        hs_rt = hdif * hdif * htmp_rt - 4.0 / (rtz * rtz) * rtz_rt + hdif * 2.0 * htmp * (-ho_rt);
    }
    let fm = 1.0 + 0.014 * msq;
    hs = (hs + 0.028 * msq) / fm;
    hs_hk /= fm;
    hs_rt /= fm;
    let hs_msq = 0.028 / fm - 0.014 * hs / fm;
    (hs, hs_hk, hs_rt, hs_msq)
}

/// Turbulent skin friction (Coles). Returns (cf, cf_hk, cf_rt, cf_msq).
#[inline]
pub fn cft(hk: f64, rt: f64, msq: f64, cffac: f64) -> (f64, f64, f64, f64) {
    const GAM: f64 = 1.4;
    let gm1 = GAM - 1.0;
    let fc = (1.0 + 0.5 * gm1 * msq).sqrt();
    let grt = (rt / fc).ln().max(3.0);
    let gex = -1.74 - 0.31 * hk;
    let arg = (-1.33 * hk).max(-20.0);
    let thk = (4.0 - hk / 0.875).tanh();
    let cfo = cffac * 0.3 * arg.exp() * (grt / 2.3026).powf(gex);
    let cf = (cfo + 1.1e-4 * (thk - 1.0)) / fc;
    let cf_hk = (-1.33 * cfo - 0.31 * (grt / 2.3026).ln() * cfo - 1.1e-4 * (1.0 - thk * thk) / 0.875) / fc;
    let cf_rt = gex * cfo / (fc * grt) / rt;
    let cf_msq = gex * cfo / (fc * grt) * (-0.25 * gm1 / (fc * fc)) - 0.25 * gm1 * cf / (fc * fc);
    (cf, cf_hk, cf_rt, cf_msq)
}

/// Density shape parameter (Whitfield). Returns (hc, hc_hk, hc_msq).
#[inline]
pub fn hct(hk: f64, msq: f64) -> (f64, f64, f64) {
    let hc = msq * (0.064 / (hk - 0.8) + 0.251);
    let hc_hk = msq * (-0.064 / ((hk - 0.8) * (hk - 0.8)));
    let hc_msq = 0.064 / (hk - 0.8) + 0.251;
    (hc, hc_hk, hc_msq)
}

/// Envelope e^n amplification rate. Returns (ax, ax_hk, ax_th, ax_rt).
#[inline]
pub fn dampl(hk: f64, th: f64, rt: f64) -> (f64, f64, f64, f64) {
    const DGR: f64 = 0.08;
    let hmi = 1.0 / (hk - 1.0);
    let hmi_hk = -hmi * hmi;
    let aa = 2.492 * hmi.powf(0.43);
    let aa_hk = (aa / hmi) * 0.43 * hmi_hk;
    let bb = (14.0 * hmi - 9.24).tanh();
    let bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;
    let grcrit = aa + 0.7 * (bb + 1.0);
    let grc_hk = aa_hk + 0.7 * bb_hk;
    let gr = rt.log10();
    let gr_rt = 1.0 / (2.3025851 * rt);
    if gr < grcrit - DGR {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let rnorm = (gr - (grcrit - DGR)) / (2.0 * DGR);
    let rn_hk = -grc_hk / (2.0 * DGR);
    let rn_rt = gr_rt / (2.0 * DGR);
    let (rfac, rfac_hk, rfac_rt) = if rnorm >= 1.0 {
        (1.0, 0.0, 0.0)
    } else {
        let rfac = 3.0 * rnorm * rnorm - 2.0 * rnorm * rnorm * rnorm;
        let rfac_rn = 6.0 * rnorm - 6.0 * rnorm * rnorm;
        (rfac, rfac_rn * rn_hk, rfac_rn * rn_rt)
    };
    let arg = 3.87 * hmi - 2.52;
    let arg_hk = 3.87 * hmi_hk;
    let ex = (-arg * arg).exp();
    let ex_hk = ex * (-2.0 * arg * arg_hk);
    let dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    let dadr_hk = 0.028 - 0.0345 * ex_hk;
    let af = -0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi * hmi * hmi;
    let af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi * hmi;
    let af_hk = af_hmi * hmi_hk;
    let ax = (af * dadr / th) * rfac;
    let ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk;
    let ax_th = -ax / th;
    let ax_rt = (af * dadr / th) * rfac_rt;
    (ax, ax_hk, ax_th, ax_rt)
}

/// Modified envelope e^n amplification rate (with separated-profile correction).
#[inline]
pub fn dampl2(hk: f64, th: f64, rt: f64) -> (f64, f64, f64, f64) {
    const DGR: f64 = 0.08;
    const HK1: f64 = 3.5;
    const HK2: f64 = 4.0;
    let hmi = 1.0 / (hk - 1.0);
    let hmi_hk = -hmi * hmi;
    let aa = 2.492 * hmi.powf(0.43);
    let aa_hk = (aa / hmi) * 0.43 * hmi_hk;
    let bb = (14.0 * hmi - 9.24).tanh();
    let bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;
    let grc = aa + 0.7 * (bb + 1.0);
    let grc_hk = aa_hk + 0.7 * bb_hk;
    let gr = rt.log10();
    let gr_rt = 1.0 / (2.3025851 * rt);
    let (mut ax, mut ax_hk, mut ax_th, mut ax_rt);
    if gr < grc - DGR {
        ax = 0.0;
        ax_hk = 0.0;
        ax_th = 0.0;
        ax_rt = 0.0;
    } else {
        let rnorm = (gr - (grc - DGR)) / (2.0 * DGR);
        let rn_hk = -grc_hk / (2.0 * DGR);
        let rn_rt = gr_rt / (2.0 * DGR);
        let (rfac, rfac_hk, rfac_rt) = if rnorm >= 1.0 {
            (1.0, 0.0, 0.0)
        } else {
            let rfac = 3.0 * rnorm * rnorm - 2.0 * rnorm * rnorm * rnorm;
            let rfac_rn = 6.0 * rnorm - 6.0 * rnorm * rnorm;
            (rfac, rfac_rn * rn_hk, rfac_rn * rn_rt)
        };
        let arg = 3.87 * hmi - 2.52;
        let arg_hk = 3.87 * hmi_hk;
        let ex = (-arg * arg).exp();
        let ex_hk = ex * (-2.0 * arg * arg_hk);
        let dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
        let dadr_hk = 0.028 - 0.0345 * ex_hk;
        let brg = -20.0 * hmi;
        let af = -0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi * hmi * hmi + 0.1 * brg.exp();
        let af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi * hmi - 2.0 * brg.exp();
        let af_hk = af_hmi * hmi_hk;
        ax = (af * dadr / th) * rfac;
        ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk;
        ax_th = -ax / th;
        ax_rt = (af * dadr / th) * rfac_rt;
    }
    if hk < HK1 {
        return (ax, ax_hk, ax_th, ax_rt);
    }
    let hnorm = (hk - HK1) / (HK2 - HK1);
    let hn_hk = 1.0 / (HK2 - HK1);
    let (hfac, hf_hk) = if hnorm >= 1.0 {
        (1.0, 0.0)
    } else {
        (3.0 * hnorm * hnorm - 2.0 * hnorm * hnorm * hnorm, (6.0 * hnorm - 6.0 * hnorm * hnorm) * hn_hk)
    };
    let (ax1, ax1_hk, ax1_th, ax1_rt) = (ax, ax_hk, ax_th, ax_rt);
    let gr0 = 0.30 + 0.35 * (-0.15 * (hk - 5.0)).exp();
    let gr0_hk = -0.35 * (-0.15 * (hk - 5.0)).exp() * 0.15;
    let tnr = (1.2 * (gr - gr0)).tanh();
    let tnr_rt = (1.0 - tnr * tnr) * 1.2 * gr_rt;
    let tnr_hk = -(1.0 - tnr * tnr) * 1.2 * gr0_hk;
    let mut ax2 = (0.086 * tnr - 0.25 / (hk - 1.0).powf(1.5)) / th;
    let mut ax2_hk = (0.086 * tnr_hk + 1.5 * 0.25 / (hk - 1.0).powf(2.5)) / th;
    let mut ax2_rt = (0.086 * tnr_rt) / th;
    let mut ax2_th = -ax2 / th;
    if ax2 < 0.0 {
        ax2 = 0.0;
        ax2_hk = 0.0;
        ax2_rt = 0.0;
        ax2_th = 0.0;
    }
    ax = hfac * ax2 + (1.0 - hfac) * ax1;
    ax_hk = hfac * ax2_hk + (1.0 - hfac) * ax1_hk + hf_hk * (ax2 - ax1);
    ax_rt = hfac * ax2_rt + (1.0 - hfac) * ax1_rt;
    ax_th = hfac * ax2_th + (1.0 - hfac) * ax1_th;
    (ax, ax_hk, ax_th, ax_rt)
}

/// Result of AXSET.
#[derive(Clone, Copy, Default)]
struct Ax {
    ax: f64,
    ax_hk1: f64,
    ax_t1: f64,
    ax_rt1: f64,
    ax_a1: f64,
    ax_hk2: f64,
    ax_t2: f64,
    ax_rt2: f64,
    ax_a2: f64,
}

/// Average amplification rate over interval 1..2.
#[allow(clippy::too_many_arguments)]
#[inline]
fn axset(hk1: f64, t1: f64, rt1: f64, a1: f64, hk2: f64, t2: f64, rt2: f64, a2: f64, acrit: f64, idampv: i32) -> Ax {
    let (ax1, ax1_hk1, ax1_t1, ax1_rt1, ax2, ax2_hk2, ax2_t2, ax2_rt2);
    if idampv == 0 {
        (ax1, ax1_hk1, ax1_t1, ax1_rt1) = dampl(hk1, t1, rt1);
        (ax2, ax2_hk2, ax2_t2, ax2_rt2) = dampl(hk2, t2, rt2);
    } else {
        (ax1, ax1_hk1, ax1_t1, ax1_rt1) = dampl2(hk1, t1, rt1);
        (ax2, ax2_hk2, ax2_t2, ax2_rt2) = dampl2(hk2, t2, rt2);
    }
    let axsq = 0.5 * (ax1 * ax1 + ax2 * ax2);
    let (axa, axa_ax1, axa_ax2) = if axsq <= 0.0 {
        (0.0, 0.0, 0.0)
    } else {
        let axa = axsq.sqrt();
        (axa, 0.5 * ax1 / axa, 0.5 * ax2 / axa)
    };
    let arg = (20.0 * (acrit - 0.5 * (a1 + a2))).min(20.0);
    let (exn, exn_a1, exn_a2) = if arg <= 0.0 {
        (1.0, 0.0, 0.0)
    } else {
        let exn = (-arg).exp();
        (exn, 20.0 * 0.5 * exn, 20.0 * 0.5 * exn)
    };
    let dax = exn * 0.002 / (t1 + t2);
    let dax_a1 = exn_a1 * 0.002 / (t1 + t2);
    let dax_a2 = exn_a2 * 0.002 / (t1 + t2);
    let dax_t1 = -dax / (t1 + t2);
    let dax_t2 = -dax / (t1 + t2);
    Ax {
        ax: axa + dax,
        ax_hk1: axa_ax1 * ax1_hk1,
        ax_t1: axa_ax1 * ax1_t1 + dax_t1,
        ax_rt1: axa_ax1 * ax1_rt1,
        ax_a1: dax_a1,
        ax_hk2: axa_ax2 * ax2_hk2,
        ax_t2: axa_ax2 * ax2_t2 + dax_t2,
        ax_rt2: axa_ax2 * ax2_rt2,
        ax_a2: dax_a2,
    }
}

// ---------------------------------------------------------------------------
// Station-level routines
// ---------------------------------------------------------------------------

impl BlSys {
    pub fn new() -> Self {
        BlSys { par: BlPar::default(), ..Default::default() }
    }

    #[inline]
    fn zero_sys(&mut self) {
        self.vs1 = [[0.0; 5]; 4];
        self.vs2 = [[0.0; 5]; 4];
        self.vsrez = [0.0; 4];
        self.vsm = [0.0; 4];
        self.vsr = [0.0; 4];
        self.vsx = [0.0; 4];
    }

    /// Set primary "2" variables (BLPRV).
    #[allow(clippy::too_many_arguments)]
    pub fn blprv(&mut self, xsi: f64, ami: f64, cti: f64, thi: f64, dsi: f64, dswaki: f64, uei: f64) {
        let tk = self.tkbl;
        let q = self.qinfbl;
        let c2 = &mut self.c2;
        c2.x = xsi;
        c2.ampl = ami;
        c2.s = cti;
        c2.t = thi;
        c2.d = dsi - dswaki;
        c2.dw = dswaki;
        let uq2 = (uei / q) * (uei / q);
        c2.u = uei * (1.0 - tk) / (1.0 - tk * uq2);
        c2.u_uei = (1.0 + tk * (2.0 * c2.u * uei / (q * q) - 1.0)) / (1.0 - tk * uq2);
        c2.u_ms = (c2.u * uq2 - uei) * self.tkbl_ms / (1.0 - tk * uq2);
    }

    /// Turbulence-independent secondary "2" variables (BLKIN).
    pub fn blkin(&mut self) {
        let hstinv = self.hstinv;
        let hstinv_ms = self.hstinv_ms;
        let gm1 = self.gm1bl;
        let c = &mut self.c2;
        let u2 = c.u;
        c.m = u2 * u2 * hstinv / (gm1 * (1.0 - 0.5 * u2 * u2 * hstinv));
        let tr2 = 1.0 + 0.5 * gm1 * c.m;
        c.m_u = 2.0 * c.m * tr2 / u2;
        c.m_ms = u2 * u2 * tr2 / (gm1 * (1.0 - 0.5 * u2 * u2 * hstinv)) * hstinv_ms;

        let trp = tr2.powf(-1.0 / gm1);
        c.r = self.rstbl * trp;
        c.r_u = -c.r / tr2 * 0.5 * c.m_u;
        c.r_ms = -c.r / tr2 * 0.5 * c.m_ms + self.rstbl_ms * trp;

        c.h = c.d / c.t;
        c.h_d = 1.0 / c.t;
        c.h_t = -c.h / c.t;

        let herat = 1.0 - 0.5 * u2 * u2 * hstinv;
        let he_u2 = -u2 * hstinv;
        let he_ms = -0.5 * u2 * u2 * hstinv_ms;

        let hvrat = self.hvrat;
        c.v = (herat * herat * herat).sqrt() * (1.0 + hvrat) / (herat + hvrat) / self.reybl;
        let v2_he = c.v * (1.5 / herat - 1.0 / (herat + hvrat));
        c.v_u = v2_he * he_u2;
        c.v_ms = -c.v / self.reybl * self.reybl_ms + v2_he * he_ms;
        c.v_re = -c.v / self.reybl * self.reybl_re;

        let (hk, hk_h2, hk_m2) = hkin(c.h, c.m);
        c.hk = hk;
        c.hk_u = hk_m2 * c.m_u;
        c.hk_t = hk_h2 * c.h_t;
        c.hk_d = hk_h2 * c.h_d;
        c.hk_ms = hk_m2 * c.m_ms;

        c.rt = c.r * u2 * c.t / c.v;
        c.rt_u = c.rt * (1.0 / u2 + c.r_u / c.r - c.v_u / c.v);
        c.rt_t = c.rt / c.t;
        c.rt_ms = c.rt * (c.r_ms / c.r - c.v_ms / c.v);
        c.rt_re = c.rt * (-c.v_re / c.v);
    }

    /// All secondary "2" variables (BLVAR). ityp: 1 lam, 2 turb, 3 wake.
    pub fn blvar(&mut self, ityp: i32) {
        let p = self.par;
        let c = &mut self.c2;
        if ityp == 3 {
            c.hk = c.hk.max(1.00005);
        } else {
            c.hk = c.hk.max(1.05000);
        }
        let hk2 = c.hk;
        let h2 = c.h;
        let rt2 = c.rt;
        let m2 = c.m;

        let (hc, hc_hk2, hc_m2) = hct(hk2, m2);
        c.hc = hc;
        c.hc_u = hc_hk2 * c.hk_u + hc_m2 * c.m_u;
        c.hc_t = hc_hk2 * c.hk_t;
        c.hc_d = hc_hk2 * c.hk_d;
        c.hc_ms = hc_hk2 * c.hk_ms + hc_m2 * c.m_ms;

        let (hs, hs_hk2, hs_rt2, hs_m2) = if ityp == 1 { hsl(hk2, rt2, m2) } else { hst(hk2, rt2, m2) };
        c.hs = hs;
        c.hs_u = hs_hk2 * c.hk_u + hs_rt2 * c.rt_u + hs_m2 * c.m_u;
        c.hs_t = hs_hk2 * c.hk_t + hs_rt2 * c.rt_t;
        c.hs_d = hs_hk2 * c.hk_d;
        c.hs_ms = hs_hk2 * c.hk_ms + hs_rt2 * c.rt_ms + hs_m2 * c.m_ms;
        c.hs_re = hs_rt2 * c.rt_re;

        // normalized slip velocity Us
        c.us = 0.5 * hs * (1.0 - (hk2 - 1.0) / (p.gbcon * h2));
        let us2_hs2 = 0.5 * (1.0 - (hk2 - 1.0) / (p.gbcon * h2));
        let us2_hk2 = 0.5 * hs * (-1.0 / (p.gbcon * h2));
        let us2_h2 = 0.5 * hs * (hk2 - 1.0) / (p.gbcon * h2 * h2);
        c.us_u = us2_hs2 * c.hs_u + us2_hk2 * c.hk_u;
        c.us_t = us2_hs2 * c.hs_t + us2_hk2 * c.hk_t + us2_h2 * c.h_t;
        c.us_d = us2_hs2 * c.hs_d + us2_hk2 * c.hk_d + us2_h2 * c.h_d;
        c.us_ms = us2_hs2 * c.hs_ms + us2_hk2 * c.hk_ms;
        c.us_re = us2_hs2 * c.hs_re;

        if ityp <= 2 && c.us > 0.95 {
            c.us = 0.98;
            c.us_u = 0.0;
            c.us_t = 0.0;
            c.us_d = 0.0;
            c.us_ms = 0.0;
            c.us_re = 0.0;
        }
        if ityp == 3 && c.us > 0.99995 {
            c.us = 0.99995;
            c.us_u = 0.0;
            c.us_t = 0.0;
            c.us_d = 0.0;
            c.us_ms = 0.0;
            c.us_re = 0.0;
        }

        // equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2
        let mut hkc = hk2 - 1.0;
        let mut hkc_hk2 = 1.0;
        let mut hkc_rt2 = 0.0;
        if ityp == 2 {
            let gcc = p.gccon;
            hkc = hk2 - 1.0 - gcc / rt2;
            hkc_hk2 = 1.0;
            hkc_rt2 = gcc / (rt2 * rt2);
            if hkc < 0.01 {
                hkc = 0.01;
                hkc_hk2 = 0.0;
                hkc_rt2 = 0.0;
            }
        }
        let hkb = hk2 - 1.0;
        let usb = 1.0 - c.us;
        let ctc = p.ctcon;
        let den = usb * h2 * hk2 * hk2;
        c.cq = (ctc * hs * hkb * hkc * hkc / den).sqrt();
        let cq = c.cq;
        let cq2_hs2 = ctc * hkb * hkc * hkc / den * 0.5 / cq;
        let cq2_us2 = ctc * hs * hkb * hkc * hkc / den / usb * 0.5 / cq;
        let cq2_hk2 = ctc * hs * hkc * hkc / den * 0.5 / cq
            - ctc * hs * hkb * hkc * hkc / (usb * h2 * hk2 * hk2 * hk2) * 2.0 * 0.5 / cq
            + ctc * hs * hkb * hkc / den * 2.0 * 0.5 / cq * hkc_hk2;
        let cq2_rt2 = ctc * hs * hkb * hkc / den * 2.0 * 0.5 / cq * hkc_rt2;
        let cq2_h2 = -ctc * hs * hkb * hkc * hkc / den / h2 * 0.5 / cq;

        c.cq_u = cq2_hs2 * c.hs_u + cq2_us2 * c.us_u + cq2_hk2 * c.hk_u;
        c.cq_t = cq2_hs2 * c.hs_t + cq2_us2 * c.us_t + cq2_hk2 * c.hk_t;
        c.cq_d = cq2_hs2 * c.hs_d + cq2_us2 * c.us_d + cq2_hk2 * c.hk_d;
        c.cq_ms = cq2_hs2 * c.hs_ms + cq2_us2 * c.us_ms + cq2_hk2 * c.hk_ms;
        c.cq_re = cq2_hs2 * c.hs_re + cq2_us2 * c.us_re;

        c.cq_u += cq2_rt2 * c.rt_u;
        c.cq_t += cq2_h2 * c.h_t + cq2_rt2 * c.rt_t;
        c.cq_d += cq2_h2 * c.h_d;
        c.cq_ms += cq2_rt2 * c.rt_ms;
        c.cq_re += cq2_rt2 * c.rt_re;

        // skin friction (turbulent CFT is also reused for the wall dissipation
        // below: same arguments, so the result is identical to recomputing it)
        let cft_t = if ityp == 2 { cft(hk2, rt2, m2, p.cffac) } else { (0.0, 0.0, 0.0, 0.0) };
        let (cf2, cf2_hk2, cf2_rt2, cf2_m2);
        if ityp == 3 {
            cf2 = 0.0;
            cf2_hk2 = 0.0;
            cf2_rt2 = 0.0;
            cf2_m2 = 0.0;
        } else if ityp == 1 {
            (cf2, cf2_hk2, cf2_rt2, cf2_m2) = cfl(hk2, rt2, m2);
        } else {
            let (a, b, cc, d) = cft_t;
            let (la, lb, lc, ld) = cfl(hk2, rt2, m2);
            if la > a {
                (cf2, cf2_hk2, cf2_rt2, cf2_m2) = (la, lb, lc, ld);
            } else {
                (cf2, cf2_hk2, cf2_rt2, cf2_m2) = (a, b, cc, d);
            }
        }
        c.cf = cf2;
        c.cf_u = cf2_hk2 * c.hk_u + cf2_rt2 * c.rt_u + cf2_m2 * c.m_u;
        c.cf_t = cf2_hk2 * c.hk_t + cf2_rt2 * c.rt_t;
        c.cf_d = cf2_hk2 * c.hk_d;
        c.cf_ms = cf2_hk2 * c.hk_ms + cf2_rt2 * c.rt_ms + cf2_m2 * c.m_ms;
        c.cf_re = cf2_rt2 * c.rt_re;

        // dissipation function 2 CD / H*
        if ityp == 1 {
            let (di, di_hk2, di_rt2) = dil(hk2, rt2);
            c.di = di;
            c.di_u = di_hk2 * c.hk_u + di_rt2 * c.rt_u;
            c.di_t = di_hk2 * c.hk_t + di_rt2 * c.rt_t;
            c.di_d = di_hk2 * c.hk_d;
            c.di_s = 0.0;
            c.di_ms = di_hk2 * c.hk_ms + di_rt2 * c.rt_ms;
            c.di_re = di_rt2 * c.rt_re;
        } else if ityp == 2 {
            let (cf2t, cf2t_hk2, cf2t_rt2, cf2t_m2) = cft_t;
            let cf2t_u2 = cf2t_hk2 * c.hk_u + cf2t_rt2 * c.rt_u + cf2t_m2 * c.m_u;
            let cf2t_t2 = cf2t_hk2 * c.hk_t + cf2t_rt2 * c.rt_t;
            let cf2t_d2 = cf2t_hk2 * c.hk_d;
            let cf2t_ms = cf2t_hk2 * c.hk_ms + cf2t_rt2 * c.rt_ms + cf2t_m2 * c.m_ms;
            let cf2t_re = cf2t_rt2 * c.rt_re;

            c.di = (0.5 * cf2t * c.us) * 2.0 / hs;
            let di2_hs2 = -(0.5 * cf2t * c.us) * 2.0 / (hs * hs);
            let di2_us2 = (0.5 * cf2t) * 2.0 / hs;
            let di2_cf2t = (0.5 * c.us) * 2.0 / hs;

            c.di_s = 0.0;
            c.di_u = di2_hs2 * c.hs_u + di2_us2 * c.us_u + di2_cf2t * cf2t_u2;
            c.di_t = di2_hs2 * c.hs_t + di2_us2 * c.us_t + di2_cf2t * cf2t_t2;
            c.di_d = di2_hs2 * c.hs_d + di2_us2 * c.us_d + di2_cf2t * cf2t_d2;
            c.di_ms = di2_hs2 * c.hs_ms + di2_us2 * c.us_ms + di2_cf2t * cf2t_ms;
            c.di_re = di2_hs2 * c.hs_re + di2_us2 * c.us_re + di2_cf2t * cf2t_re;

            let grt = rt2.ln();
            let hmin = 1.0 + 2.1 / grt;
            let hm_rt2 = -(2.1 / (grt * grt)) / rt2;
            let fl = (hk2 - 1.0) / (hmin - 1.0);
            let fl_hk2 = 1.0 / (hmin - 1.0);
            let fl_rt2 = (-fl / (hmin - 1.0)) * hm_rt2;
            let tfl = fl.tanh();
            let dfac = 0.5 + 0.5 * tfl;
            let df_fl = 0.5 * (1.0 - tfl * tfl);
            let df_hk2 = df_fl * fl_hk2;
            let df_rt2 = df_fl * fl_rt2;

            c.di_s *= dfac;
            c.di_u = c.di_u * dfac + c.di * (df_hk2 * c.hk_u + df_rt2 * c.rt_u);
            c.di_t = c.di_t * dfac + c.di * (df_hk2 * c.hk_t + df_rt2 * c.rt_t);
            c.di_d = c.di_d * dfac + c.di * (df_hk2 * c.hk_d);
            c.di_ms = c.di_ms * dfac + c.di * (df_hk2 * c.hk_ms + df_rt2 * c.rt_ms);
            c.di_re = c.di_re * dfac + c.di * (df_rt2 * c.rt_re);
            c.di *= dfac;
        } else {
            c.di = 0.0;
            c.di_s = 0.0;
            c.di_u = 0.0;
            c.di_t = 0.0;
            c.di_d = 0.0;
            c.di_ms = 0.0;
            c.di_re = 0.0;
        }

        // turbulent outer layer contribution
        if ityp != 1 {
            let s2 = c.s;
            let us = c.us;
            let dd = s2 * s2 * (0.995 - us) * 2.0 / hs;
            let dd_hs2 = -s2 * s2 * (0.995 - us) * 2.0 / (hs * hs);
            let dd_us2 = -s2 * s2 * 2.0 / hs;
            let dd_s2 = s2 * 2.0 * (0.995 - us) * 2.0 / hs;
            c.di += dd;
            c.di_s = dd_s2;
            c.di_u += dd_hs2 * c.hs_u + dd_us2 * c.us_u;
            c.di_t += dd_hs2 * c.hs_t + dd_us2 * c.us_t;
            c.di_d += dd_hs2 * c.hs_d + dd_us2 * c.us_d;
            c.di_ms += dd_hs2 * c.hs_ms + dd_us2 * c.us_ms;
            c.di_re += dd_hs2 * c.hs_re + dd_us2 * c.us_re;

            let dd = 0.15 * (0.995 - us) * (0.995 - us) / rt2 * 2.0 / hs;
            let dd_us2 = -0.15 * (0.995 - us) * 2.0 / rt2 * 2.0 / hs;
            let dd_hs2 = -dd / hs;
            let dd_rt2 = -dd / rt2;
            c.di += dd;
            c.di_u += dd_hs2 * c.hs_u + dd_us2 * c.us_u + dd_rt2 * c.rt_u;
            c.di_t += dd_hs2 * c.hs_t + dd_us2 * c.us_t + dd_rt2 * c.rt_t;
            c.di_d += dd_hs2 * c.hs_d + dd_us2 * c.us_d;
            c.di_ms += dd_hs2 * c.hs_ms + dd_us2 * c.us_ms + dd_rt2 * c.rt_ms;
            c.di_re += dd_hs2 * c.hs_re + dd_us2 * c.us_re + dd_rt2 * c.rt_re;
        }

        if ityp == 2 {
            let (di2l, di2l_hk2, di2l_rt2) = dil(hk2, rt2);
            if di2l > c.di {
                c.di = di2l;
                c.di_s = 0.0;
                c.di_u = di2l_hk2 * c.hk_u + di2l_rt2 * c.rt_u;
                c.di_t = di2l_hk2 * c.hk_t + di2l_rt2 * c.rt_t;
                c.di_d = di2l_hk2 * c.hk_d;
                c.di_ms = di2l_hk2 * c.hk_ms + di2l_rt2 * c.rt_ms;
                c.di_re = di2l_rt2 * c.rt_re;
            }
        }

        if ityp == 3 {
            let (di2l, di2l_hk2, di2l_rt2) = dilw(hk2, rt2);
            if di2l > c.di {
                c.di = di2l;
                c.di_s = 0.0;
                c.di_u = di2l_hk2 * c.hk_u + di2l_rt2 * c.rt_u;
                c.di_t = di2l_hk2 * c.hk_t + di2l_rt2 * c.rt_t;
                c.di_d = di2l_hk2 * c.hk_d;
                c.di_ms = di2l_hk2 * c.hk_ms + di2l_rt2 * c.rt_ms;
                c.di_re = di2l_rt2 * c.rt_re;
            }
            c.di *= 2.0;
            c.di_s *= 2.0;
            c.di_u *= 2.0;
            c.di_t *= 2.0;
            c.di_d *= 2.0;
            c.di_ms *= 2.0;
            c.di_re *= 2.0;
        }

        // BL thickness (Delta) from simplified Green's correlation
        c.de = (3.15 + 1.72 / (hk2 - 1.0)) * c.t + c.d;
        let de2_hk2 = (-1.72 / ((hk2 - 1.0) * (hk2 - 1.0))) * c.t;
        c.de_u = de2_hk2 * c.hk_u;
        c.de_t = de2_hk2 * c.hk_t + (3.15 + 1.72 / (hk2 - 1.0));
        c.de_d = de2_hk2 * c.hk_d + 1.0;
        c.de_ms = de2_hk2 * c.hk_ms;
        let hdmax = 12.0;
        if c.de > hdmax * c.t {
            c.de = hdmax * c.t;
            c.de_u = 0.0;
            c.de_t = hdmax;
            c.de_d = 0.0;
            c.de_ms = 0.0;
        }
    }

    /// Midpoint skin friction CFM (BLMID).
    pub fn blmid(&mut self, ityp: i32) {
        if self.simi {
            let c2 = self.c2;
            let c1 = &mut self.c1;
            c1.hk = c2.hk;
            c1.hk_t = c2.hk_t;
            c1.hk_d = c2.hk_d;
            c1.hk_u = c2.hk_u;
            c1.hk_ms = c2.hk_ms;
            c1.rt = c2.rt;
            c1.rt_t = c2.rt_t;
            c1.rt_u = c2.rt_u;
            c1.rt_ms = c2.rt_ms;
            c1.rt_re = c2.rt_re;
            c1.m = c2.m;
            c1.m_u = c2.m_u;
            c1.m_ms = c2.m_ms;
        }
        let c1 = &self.c1;
        let c2 = &self.c2;
        let hka = 0.5 * (c1.hk + c2.hk);
        let rta = 0.5 * (c1.rt + c2.rt);
        let ma = 0.5 * (c1.m + c2.m);
        let (cfm, cfm_hka, cfm_rta, cfm_ma);
        if ityp == 3 {
            cfm = 0.0;
            cfm_hka = 0.0;
            cfm_rta = 0.0;
            cfm_ma = 0.0;
        } else if ityp == 1 {
            (cfm, cfm_hka, cfm_rta, cfm_ma) = cfl(hka, rta, ma);
        } else {
            let t = cft(hka, rta, ma, self.par.cffac);
            let l = cfl(hka, rta, ma);
            if l.0 > t.0 {
                (cfm, cfm_hka, cfm_rta, cfm_ma) = l;
            } else {
                (cfm, cfm_hka, cfm_rta, cfm_ma) = t;
            }
        }
        self.cfm = cfm;
        self.cfm_u1 = 0.5 * (cfm_hka * c1.hk_u + cfm_ma * c1.m_u + cfm_rta * c1.rt_u);
        self.cfm_t1 = 0.5 * (cfm_hka * c1.hk_t + cfm_rta * c1.rt_t);
        self.cfm_d1 = 0.5 * (cfm_hka * c1.hk_d);
        self.cfm_u2 = 0.5 * (cfm_hka * c2.hk_u + cfm_ma * c2.m_u + cfm_rta * c2.rt_u);
        self.cfm_t2 = 0.5 * (cfm_hka * c2.hk_t + cfm_rta * c2.rt_t);
        self.cfm_d2 = 0.5 * (cfm_hka * c2.hk_d);
        self.cfm_ms = 0.5
            * (cfm_hka * c1.hk_ms + cfm_ma * c1.m_ms + cfm_rta * c1.rt_ms + cfm_hka * c2.hk_ms + cfm_ma * c2.m_ms + cfm_rta * c2.rt_ms);
        self.cfm_re = 0.5 * (cfm_rta * c1.rt_re + cfm_rta * c2.rt_re);
    }

    /// Transition check, 2nd-order implicit amplification (TRCHEK2).
    #[allow(unused_assignments)]
    pub fn trchek(&mut self) {
        const DAEPS: f64 = 5.0e-5;
        self.c2sav = self.c2;
        let amcrit = self.amcrit;
        let idampv = self.idampv;
        let xiforc = self.xiforc;

        let (x1, t1, d1, u1, ampl1, hk1, rt1) = (self.c1.x, self.c1.t, self.c1.d, self.c1.u, self.c1.ampl, self.c1.hk, self.c1.rt);
        let (x2, t2, d2, u2) = (self.c2.x, self.c2.t, self.c2.d, self.c2.u);

        let a = axset(hk1, t1, rt1, ampl1, self.c2.hk, t2, self.c2.rt, self.c2.ampl, amcrit, idampv);
        let mut ampl2 = ampl1 + a.ax * (x2 - x1);

        let mut ax = 0.0;
        let mut axr = Ax::default();
        let (mut wf1, mut wf2) = (0.0, 0.0);
        let (mut wf1_a1, mut wf2_a1) = (0.0, 0.0);
        let (mut wf1_x1, mut wf2_x1, mut wf1_x2, mut wf2_x2, mut wf1_xf, mut wf2_xf) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let (mut xt, mut tt, mut dt, mut ut) = (0.0, 0.0, 0.0, 0.0);
        let (mut xt_a2, mut tt_a2, mut dt_a2, mut ut_a2) = (0.0, 0.0, 0.0, 0.0);
        let (mut hkt_tt, mut hkt_dt, mut hkt_ut, mut hkt_ms) = (0.0, 0.0, 0.0, 0.0);
        let (mut rtt_tt, mut rtt_ut, mut rtt_ms, mut rtt_re) = (0.0, 0.0, 0.0, 0.0);
        let mut amplt_a2 = 0.0;

        let mut converged = false;
        for _itam in 0..30 {
            let (amplt, sfa, sfa_a1, sfa_a2);
            if ampl2 <= amcrit {
                amplt = ampl2;
                amplt_a2 = 1.0;
                sfa = 1.0;
                sfa_a1 = 0.0;
                sfa_a2 = 0.0;
            } else {
                amplt = amcrit;
                amplt_a2 = 0.0;
                sfa = (amplt - ampl1) / (ampl2 - ampl1);
                sfa_a1 = (sfa - 1.0) / (ampl2 - ampl1);
                sfa_a2 = (-sfa) / (ampl2 - ampl1);
            }
            let (sfx, sfx_x1, sfx_x2, sfx_xf) = if xiforc < x2 {
                let sfx = (xiforc - x1) / (x2 - x1);
                (sfx, (sfx - 1.0) / (x2 - x1), (-sfx) / (x2 - x1), 1.0 / (x2 - x1))
            } else {
                (1.0, 0.0, 0.0, 0.0)
            };
            let wf2_a2;
            if sfa < sfx {
                wf2 = sfa;
                wf2_a1 = sfa_a1;
                wf2_a2 = sfa_a2;
                wf2_x1 = 0.0;
                wf2_x2 = 0.0;
                wf2_xf = 0.0;
            } else {
                wf2 = sfx;
                wf2_a1 = 0.0;
                wf2_a2 = 0.0;
                wf2_x1 = sfx_x1;
                wf2_x2 = sfx_x2;
                wf2_xf = sfx_xf;
            }
            wf1 = 1.0 - wf2;
            wf1_a1 = -wf2_a1;
            let wf1_a2 = -wf2_a2;
            wf1_x1 = -wf2_x1;
            wf1_x2 = -wf2_x2;
            wf1_xf = -wf2_xf;

            xt = x1 * wf1 + x2 * wf2;
            tt = t1 * wf1 + t2 * wf2;
            dt = d1 * wf1 + d2 * wf2;
            ut = u1 * wf1 + u2 * wf2;
            xt_a2 = x1 * wf1_a2 + x2 * wf2_a2;
            tt_a2 = t1 * wf1_a2 + t2 * wf2_a2;
            dt_a2 = d1 * wf1_a2 + d2 * wf2_a2;
            ut_a2 = u1 * wf1_a2 + u2 * wf2_a2;

            // temporarily set "2" variables from "T" for BLKIN
            self.c2.x = xt;
            self.c2.t = tt;
            self.c2.d = dt;
            self.c2.u = ut;
            self.blkin();
            let hkt = self.c2.hk;
            hkt_tt = self.c2.hk_t;
            hkt_dt = self.c2.hk_d;
            hkt_ut = self.c2.hk_u;
            hkt_ms = self.c2.hk_ms;
            let rtt = self.c2.rt;
            rtt_tt = self.c2.rt_t;
            rtt_ut = self.c2.rt_u;
            rtt_ms = self.c2.rt_ms;
            rtt_re = self.c2.rt_re;

            // restore clobbered "2" variables, except for AMPL2
            self.c2 = self.c2sav;
            self.c2.ampl = ampl2;

            axr = axset(hk1, t1, rt1, ampl1, hkt, tt, rtt, amplt, amcrit, idampv);
            ax = axr.ax;
            if ax <= 0.0 {
                converged = true;
                break;
            }
            let ax_a2 = (axr.ax_hk2 * hkt_tt + axr.ax_t2 + axr.ax_rt2 * rtt_tt) * tt_a2
                + (axr.ax_hk2 * hkt_dt) * dt_a2
                + (axr.ax_hk2 * hkt_ut + axr.ax_rt2 * rtt_ut) * ut_a2
                + axr.ax_a2 * amplt_a2;
            let res = ampl2 - ampl1 - ax * (x2 - x1);
            let res_a2 = 1.0 - ax_a2 * (x2 - x1);
            let da2 = -res / res_a2;
            let mut rlx = 1.0;
            let dxt = xt_a2 * da2;
            if rlx * (dxt / (x2 - x1)).abs() > 0.05 {
                rlx = 0.05 * ((x2 - x1) / dxt).abs();
            }
            if rlx * da2.abs() > 1.0 {
                rlx = 1.0 * (1.0 / da2).abs();
            }
            if da2.abs() < DAEPS {
                converged = true;
                break;
            }
            if (ampl2 > amcrit && ampl2 + rlx * da2 < amcrit) || (ampl2 < amcrit && ampl2 + rlx * da2 > amcrit) {
                ampl2 = amcrit;
            } else {
                ampl2 += rlx * da2;
            }
            self.c2.ampl = ampl2;
        }
        if !converged {
            if crate::verbose() {
                crate::outln!(" TRCHEK2: N2 convergence failed.");
            }
        }
        self.c2.ampl = ampl2;
        self.xt = xt;

        self.trfree = ampl2 >= amcrit;
        self.trforc = xiforc > x1 && xiforc <= x2;
        self.tran = self.trforc || self.trfree;
        if !self.tran {
            return;
        }
        if self.trfree && self.trforc {
            self.trforc = xiforc < xt;
            self.trfree = xiforc >= xt;
        }
        if self.trforc {
            self.xt = xiforc;
            self.xt_a1 = 0.0;
            self.xt_x1 = 0.0;
            self.xt_t1 = 0.0;
            self.xt_d1 = 0.0;
            self.xt_u1 = 0.0;
            self.xt_x2 = 0.0;
            self.xt_t2 = 0.0;
            self.xt_d2 = 0.0;
            self.xt_u2 = 0.0;
            self.xt_ms = 0.0;
            self.xt_re = 0.0;
            self.xt_xf = 1.0;
            return;
        }

        // free transition: set sensitivities of XT
        let mut xt_x1 = wf1;
        let tt_t1 = wf1;
        let dt_d1 = wf1;
        let ut_u1 = wf1;
        let mut xt_x2 = wf2;
        let tt_t2 = wf2;
        let dt_d2 = wf2;
        let ut_u2 = wf2;

        let mut xt_a1 = x1 * wf1_a1 + x2 * wf2_a1;
        let tt_a1 = t1 * wf1_a1 + t2 * wf2_a1;
        let dt_a1 = d1 * wf1_a1 + d2 * wf2_a1;
        let ut_a1 = u1 * wf1_a1 + u2 * wf2_a1;

        xt_x1 += x1 * wf1_x1 + x2 * wf2_x1;
        let tt_x1 = t1 * wf1_x1 + t2 * wf2_x1;
        let dt_x1 = d1 * wf1_x1 + d2 * wf2_x1;
        let ut_x1 = u1 * wf1_x1 + u2 * wf2_x1;

        xt_x2 += x1 * wf1_x2 + x2 * wf2_x2;
        let tt_x2 = t1 * wf1_x2 + t2 * wf2_x2;
        let dt_x2 = d1 * wf1_x2 + d2 * wf2_x2;
        let ut_x2 = u1 * wf1_x2 + u2 * wf2_x2;

        let _xt_xf = x1 * wf1_xf + x2 * wf2_xf;
        let tt_xf = t1 * wf1_xf + t2 * wf2_xf;
        let dt_xf = d1 * wf1_xf + d2 * wf2_xf;
        let ut_xf = u1 * wf1_xf + u2 * wf2_xf;

        let c1 = &self.c1;
        let gt = axr.ax_hk2 * hkt_tt + axr.ax_t2 + axr.ax_rt2 * rtt_tt;
        let gd = axr.ax_hk2 * hkt_dt;
        let gu = axr.ax_hk2 * hkt_ut + axr.ax_rt2 * rtt_ut;

        let ax_t1 = axr.ax_hk1 * c1.hk_t + axr.ax_t1 + axr.ax_rt1 * c1.rt_t + gt * tt_t1;
        let ax_d1 = axr.ax_hk1 * c1.hk_d + gd * dt_d1;
        let ax_u1 = axr.ax_hk1 * c1.hk_u + axr.ax_rt1 * c1.rt_u + gu * ut_u1;
        let ax_a1 = axr.ax_a1 + gt * tt_a1 + gd * dt_a1 + gu * ut_a1;
        let ax_x1 = gt * tt_x1 + gd * dt_x1 + gu * ut_x1;

        let ax_t2 = gt * tt_t2;
        let ax_d2 = gd * dt_d2;
        let ax_u2 = gu * ut_u2;
        let ax_a2 = axr.ax_a2 * amplt_a2 + gt * tt_a2 + gd * dt_a2 + gu * ut_a2;
        let ax_x2 = gt * tt_x2 + gd * dt_x2 + gu * ut_x2;
        let _ax_xf = gt * tt_xf + gd * dt_xf + gu * ut_xf;

        let ax_ms = axr.ax_hk2 * hkt_ms + axr.ax_rt2 * rtt_ms + axr.ax_hk1 * c1.hk_ms + axr.ax_rt1 * c1.rt_ms;
        let ax_re = axr.ax_rt2 * rtt_re + axr.ax_rt1 * c1.rt_re;

        let z_ax = -(x2 - x1);
        let z_a1 = z_ax * ax_a1 - 1.0;
        let z_t1 = z_ax * ax_t1;
        let z_d1 = z_ax * ax_d1;
        let z_u1 = z_ax * ax_u1;
        let z_x1 = z_ax * ax_x1 + ax;
        let z_a2 = z_ax * ax_a2 + 1.0;
        let z_t2 = z_ax * ax_t2;
        let z_d2 = z_ax * ax_d2;
        let z_u2 = z_ax * ax_u2;
        let z_x2 = z_ax * ax_x2 - ax;
        let z_ms = z_ax * ax_ms;
        let z_re = z_ax * ax_re;

        let f = xt_a2 / z_a2;
        xt_a1 -= f * z_a1;
        self.xt_a1 = xt_a1;
        self.xt_t1 = -f * z_t1;
        self.xt_d1 = -f * z_d1;
        self.xt_u1 = -f * z_u1;
        self.xt_x1 = xt_x1 - f * z_x1;
        self.xt_t2 = -f * z_t2;
        self.xt_d2 = -f * z_d2;
        self.xt_u2 = -f * z_u2;
        self.xt_x2 = xt_x2 - f * z_x2;
        self.xt_ms = -f * z_ms;
        self.xt_re = -f * z_re;
        self.xt_xf = 0.0;
    }

    /// Set up BL Newton system for the current interval (BLSYS).
    pub fn blsys(&mut self) {
        if self.wake {
            self.blvar(3);
            self.blmid(3);
        } else if self.turb || self.tran {
            self.blvar(2);
            self.blmid(2);
        } else {
            self.blvar(1);
            self.blmid(1);
        }
        if self.simi {
            self.c1 = self.c2;
        }
        if self.tran {
            self.trdif();
        } else if self.simi {
            self.bldif(0);
        } else if !self.turb {
            self.bldif(1);
        } else if self.wake {
            self.bldif(3);
        } else {
            self.bldif(2);
        }
        if self.simi {
            for k in 0..4 {
                for l in 0..5 {
                    self.vs2[k][l] += self.vs1[k][l];
                    self.vs1[k][l] = 0.0;
                }
            }
        }
        // change system over into incompressible Uei and Mach
        let (u1_uei, u1_ms, u2_uei, u2_ms) = (self.c1.u_uei, self.c1.u_ms, self.c2.u_uei, self.c2.u_ms);
        for k in 0..4 {
            let res_u1 = self.vs1[k][3];
            let res_u2 = self.vs2[k][3];
            let res_ms = self.vsm[k];
            self.vs1[k][3] = res_u1 * u1_uei;
            self.vs2[k][3] = res_u2 * u2_uei;
            self.vsm[k] = res_u1 * u1_ms + res_u2 * u2_ms + res_ms;
        }
    }

    /// Dummy system between TE and first wake point (TESYS).
    pub fn tesys(&mut self, cte: f64, tte: f64, dte: f64) {
        self.zero_sys();
        self.blvar(3);
        self.vs1[0][0] = -1.0;
        self.vs2[0][0] = 1.0;
        self.vsrez[0] = cte - self.c2.s;
        self.vs1[1][1] = -1.0;
        self.vs2[1][1] = 1.0;
        self.vsrez[1] = tte - self.c2.t;
        self.vs1[2][2] = -1.0;
        self.vs2[2][2] = 1.0;
        self.vsrez[2] = dte - self.c2.d - self.c2.dw;
    }

    /// Transition interval system (TRDIF).
    fn trdif(&mut self) {
        self.c1sav = self.c1;
        self.c2sav = self.c2;
        let (x1, t1, d1, u1) = (self.c1.x, self.c1.t, self.c1.d, self.c1.u);
        let (x2, t2, d2, u2) = (self.c2.x, self.c2.t, self.c2.d, self.c2.u);
        let xt = self.xt;

        let wf2 = (xt - x1) / (x2 - x1);
        let wf2_xt = 1.0 / (x2 - x1);
        let wf2_a1 = wf2_xt * self.xt_a1;
        let wf2_x1 = wf2_xt * self.xt_x1 + (wf2 - 1.0) / (x2 - x1);
        let wf2_x2 = wf2_xt * self.xt_x2 - wf2 / (x2 - x1);
        let wf2_t1 = wf2_xt * self.xt_t1;
        let wf2_t2 = wf2_xt * self.xt_t2;
        let wf2_d1 = wf2_xt * self.xt_d1;
        let wf2_d2 = wf2_xt * self.xt_d2;
        let wf2_u1 = wf2_xt * self.xt_u1;
        let wf2_u2 = wf2_xt * self.xt_u2;
        let wf2_ms = wf2_xt * self.xt_ms;
        let wf2_re = wf2_xt * self.xt_re;
        let wf2_xf = wf2_xt * self.xt_xf;

        let wf1 = 1.0 - wf2;
        let (wf1_a1, wf1_x1, wf1_x2, wf1_t1, wf1_t2, wf1_d1, wf1_d2, wf1_u1, wf1_u2, wf1_ms, wf1_re, wf1_xf) =
            (-wf2_a1, -wf2_x1, -wf2_x2, -wf2_t1, -wf2_t2, -wf2_d1, -wf2_d2, -wf2_u1, -wf2_u2, -wf2_ms, -wf2_re, -wf2_xf);

        // interpolated "T" quantities and their sensitivities
        macro_rules! interp {
            ($v1:expr, $v2:expr) => {
                [
                    $v1 * wf1_a1 + $v2 * wf2_a1, // a1
                    $v1 * wf1_x1 + $v2 * wf2_x1, // x1
                    $v1 * wf1_x2 + $v2 * wf2_x2, // x2
                    $v1 * wf1_t1 + $v2 * wf2_t1, // t1
                    $v1 * wf1_t2 + $v2 * wf2_t2, // t2
                    $v1 * wf1_d1 + $v2 * wf2_d1, // d1
                    $v1 * wf1_d2 + $v2 * wf2_d2, // d2
                    $v1 * wf1_u1 + $v2 * wf2_u1, // u1
                    $v1 * wf1_u2 + $v2 * wf2_u2, // u2
                    $v1 * wf1_ms + $v2 * wf2_ms, // ms
                    $v1 * wf1_re + $v2 * wf2_re, // re
                    $v1 * wf1_xf + $v2 * wf2_xf, // xf
                ]
            };
        }
        const A1: usize = 0;
        const X1: usize = 1;
        const X2: usize = 2;
        const T1: usize = 3;
        const T2: usize = 4;
        const D1: usize = 5;
        const D2: usize = 6;
        const U1: usize = 7;
        const U2: usize = 8;
        const MS: usize = 9;
        const RE: usize = 10;
        const XF: usize = 11;

        let tt = t1 * wf1 + t2 * wf2;
        let mut tt_ = interp!(t1, t2);
        tt_[T1] += wf1;
        tt_[T2] += wf2;
        let dt = d1 * wf1 + d2 * wf2;
        let mut dt_ = interp!(d1, d2);
        dt_[D1] += wf1;
        dt_[D2] += wf2;
        let ut = u1 * wf1 + u2 * wf2;
        let mut ut_ = interp!(u1, u2);
        ut_[U1] += wf1;
        ut_[U2] += wf2;
        // XT sensitivities in same ordering
        let xt_ = [
            self.xt_a1, self.xt_x1, self.xt_x2, self.xt_t1, self.xt_t2, self.xt_d1, self.xt_d2, self.xt_u1, self.xt_u2, self.xt_ms,
            self.xt_re, self.xt_xf,
        ];

        // **** FIRST, laminar part between X1 and XT
        self.c2.x = xt;
        self.c2.t = tt;
        self.c2.d = dt;
        self.c2.u = ut;
        self.c2.ampl = self.amcrit;
        self.c2.s = 0.0;
        self.blkin();
        self.blvar(1);
        self.blmid(1);
        self.bldif(1);

        let mut bl1 = [[0.0; 5]; 4];
        let mut bl2 = [[0.0; 5]; 4];
        let mut blrez = [0.0; 4];
        let mut blm = [0.0; 4];
        let mut blr = [0.0; 4];
        let mut blx = [0.0; 4];
        for k in 1..3 {
            let v2 = self.vs2[k];
            let comb = |i: usize| v2[1] * tt_[i] + v2[2] * dt_[i] + v2[3] * ut_[i] + v2[4] * xt_[i];
            blrez[k] = self.vsrez[k];
            blm[k] = self.vsm[k] + comb(MS);
            blr[k] = self.vsr[k] + comb(RE);
            blx[k] = self.vsx[k] + comb(XF);
            bl1[k][0] = self.vs1[k][0] + comb(A1);
            bl1[k][1] = self.vs1[k][1] + comb(T1);
            bl1[k][2] = self.vs1[k][2] + comb(D1);
            bl1[k][3] = self.vs1[k][3] + comb(U1);
            bl1[k][4] = self.vs1[k][4] + comb(X1);
            bl2[k][0] = 0.0;
            bl2[k][1] = comb(T2);
            bl2[k][2] = comb(D2);
            bl2[k][3] = comb(U2);
            bl2[k][4] = comb(X2);
        }

        // **** SECOND, turbulent part between XT and X2
        self.blvar(2);
        let p = self.par;
        let hk2 = self.c2.hk;
        let ctr = p.ctrcon * (-p.ctrcex / (hk2 - 1.0)).exp();
        let ctr_hk2 = ctr * p.ctrcex / ((hk2 - 1.0) * (hk2 - 1.0));
        let c = &self.c2;
        let st = ctr * c.cq;
        let st_tt = ctr * c.cq_t + c.cq * ctr_hk2 * c.hk_t;
        let st_dt = ctr * c.cq_d + c.cq * ctr_hk2 * c.hk_d;
        let st_ut = ctr * c.cq_u + c.cq * ctr_hk2 * c.hk_u;
        let st_ms0 = ctr * c.cq_ms + c.cq * ctr_hk2 * c.hk_ms;
        let st_re0 = ctr * c.cq_re;
        let mut st_ = [0.0; 12];
        for i in 0..12 {
            st_[i] = st_tt * tt_[i] + st_dt * dt_[i] + st_ut * ut_[i];
        }
        st_[MS] += st_ms0;
        st_[RE] += st_re0;

        self.c2.ampl = 0.0;
        self.c2.s = st;
        self.blvar(2);

        self.c1 = self.c2;
        self.c2 = self.c2sav;
        self.blmid(2);
        self.bldif(2);

        let mut bt1 = [[0.0; 5]; 4];
        let mut bt2 = [[0.0; 5]; 4];
        let mut btrez = [0.0; 4];
        let mut btm = [0.0; 4];
        let mut btr = [0.0; 4];
        let mut btx = [0.0; 4];
        for k in 0..3 {
            let v1 = self.vs1[k];
            let comb = |i: usize| v1[0] * st_[i] + v1[1] * tt_[i] + v1[2] * dt_[i] + v1[3] * ut_[i] + v1[4] * xt_[i];
            btrez[k] = self.vsrez[k];
            btm[k] = self.vsm[k] + comb(MS);
            btr[k] = self.vsr[k] + comb(RE);
            btx[k] = self.vsx[k] + comb(XF);
            bt1[k][0] = comb(A1);
            bt1[k][1] = comb(T1);
            bt1[k][2] = comb(D1);
            bt1[k][3] = comb(U1);
            bt1[k][4] = comb(X1);
            bt2[k][0] = self.vs2[k][0];
            bt2[k][1] = self.vs2[k][1] + comb(T2);
            bt2[k][2] = self.vs2[k][2] + comb(D2);
            bt2[k][3] = self.vs2[k][3] + comb(U2);
            bt2[k][4] = self.vs2[k][4] + comb(X2);
        }

        self.vsrez[0] = btrez[0];
        self.vsrez[1] = blrez[1] + btrez[1];
        self.vsrez[2] = blrez[2] + btrez[2];
        self.vsm[0] = btm[0];
        self.vsm[1] = blm[1] + btm[1];
        self.vsm[2] = blm[2] + btm[2];
        self.vsr[0] = btr[0];
        self.vsr[1] = blr[1] + btr[1];
        self.vsr[2] = blr[2] + btr[2];
        self.vsx[0] = btx[0];
        self.vsx[1] = blx[1] + btx[1];
        self.vsx[2] = blx[2] + btx[2];
        for l in 0..5 {
            self.vs1[0][l] = bt1[0][l];
            self.vs2[0][l] = bt2[0][l];
            self.vs1[1][l] = bl1[1][l] + bt1[1][l];
            self.vs2[1][l] = bl2[1][l] + bt2[1][l];
            self.vs1[2][l] = bl1[2][l] + bt1[2][l];
            self.vs2[2][l] = bl2[2][l] + bt2[2][l];
        }
        self.c1 = self.c1sav;
    }

    /// Newton system coefficients for an interval (BLDIF).
    /// ityp: 0 similarity, 1 laminar, 2 turbulent, 3 wake.
    pub fn bldif(&mut self, ityp: i32) {
        let p = self.par;
        let c1 = self.c1;
        let c2 = self.c2;
        let (xlog, ulog, tlog, hlog, ddlog);
        if ityp == 0 {
            xlog = 1.0;
            ulog = self.bule;
            tlog = 0.5 * (1.0 - self.bule);
            hlog = 0.0;
            ddlog = 0.0;
        } else {
            xlog = (c2.x / c1.x).ln();
            ulog = (c2.u / c1.u).ln();
            tlog = (c2.t / c1.t).ln();
            hlog = (c2.hs / c1.hs).ln();
            ddlog = 1.0;
        }
        self.zero_sys();

        let hupwt = 1.0;
        let mut hdcon = 5.0 * hupwt / (c2.hk * c2.hk);
        let hd_hk1 = 0.0;
        let mut hd_hk2 = -hdcon * 2.0 / c2.hk;
        if ityp == 3 {
            hdcon = hupwt / (c2.hk * c2.hk);
            hd_hk2 = -hdcon * 2.0 / c2.hk;
        }
        let arg = ((c2.hk - 1.0) / (c1.hk - 1.0)).abs();
        let hl = arg.ln();
        let hl_hk1 = -1.0 / (c1.hk - 1.0);
        let hl_hk2 = 1.0 / (c2.hk - 1.0);
        let hlsq = (hl * hl).min(15.0);
        let ehh = (-hlsq * hdcon).exp();
        let upw = 1.0 - 0.5 * ehh;
        let upw_hl = ehh * hl * hdcon;
        let upw_hd = 0.5 * ehh * hlsq;
        let upw_hk1 = upw_hl * hl_hk1 + upw_hd * hd_hk1;
        let upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;
        let upw_u1 = upw_hk1 * c1.hk_u;
        let upw_t1 = upw_hk1 * c1.hk_t;
        let upw_d1 = upw_hk1 * c1.hk_d;
        let upw_u2 = upw_hk2 * c2.hk_u;
        let upw_t2 = upw_hk2 * c2.hk_t;
        let upw_d2 = upw_hk2 * c2.hk_d;
        let upw_ms = upw_hk1 * c1.hk_ms + upw_hk2 * c2.hk_ms;

        if ityp == 0 {
            self.vs2[0][0] = 1.0;
            self.vsr[0] = 0.0;
            self.vsrez[0] = -c2.ampl;
        } else if ityp == 1 {
            let a = axset(c1.hk, c1.t, c1.rt, c1.ampl, c2.hk, c2.t, c2.rt, c2.ampl, self.amcrit, self.idampv);
            let rezc = c2.ampl - c1.ampl - a.ax * (c2.x - c1.x);
            let z_ax = -(c2.x - c1.x);
            self.vs1[0][0] = z_ax * a.ax_a1 - 1.0;
            self.vs1[0][1] = z_ax * (a.ax_hk1 * c1.hk_t + a.ax_t1 + a.ax_rt1 * c1.rt_t);
            self.vs1[0][2] = z_ax * (a.ax_hk1 * c1.hk_d);
            self.vs1[0][3] = z_ax * (a.ax_hk1 * c1.hk_u + a.ax_rt1 * c1.rt_u);
            self.vs1[0][4] = a.ax;
            self.vs2[0][0] = z_ax * a.ax_a2 + 1.0;
            self.vs2[0][1] = z_ax * (a.ax_hk2 * c2.hk_t + a.ax_t2 + a.ax_rt2 * c2.rt_t);
            self.vs2[0][2] = z_ax * (a.ax_hk2 * c2.hk_d);
            self.vs2[0][3] = z_ax * (a.ax_hk2 * c2.hk_u + a.ax_rt2 * c2.rt_u);
            self.vs2[0][4] = -a.ax;
            self.vsm[0] = z_ax * (a.ax_hk1 * c1.hk_ms + a.ax_rt1 * c1.rt_ms + a.ax_hk2 * c2.hk_ms + a.ax_rt2 * c2.rt_ms);
            self.vsr[0] = z_ax * (a.ax_rt1 * c1.rt_re + a.ax_rt2 * c2.rt_re);
            self.vsx[0] = 0.0;
            self.vsrez[0] = -rezc;
        } else {
            // turbulent part --> shear lag equation
            let sa = (1.0 - upw) * c1.s + upw * c2.s;
            let cqa = (1.0 - upw) * c1.cq + upw * c2.cq;
            let cfa = (1.0 - upw) * c1.cf + upw * c2.cf;
            let hka = (1.0 - upw) * c1.hk + upw * c2.hk;
            let usa = 0.5 * (c1.us + c2.us);
            let rta = 0.5 * (c1.rt + c2.rt);
            let dea = 0.5 * (c1.de + c2.de);
            let da = 0.5 * (c1.d + c2.d);
            let ald = if ityp == 3 { p.dlcon } else { 1.0 };

            let (hkc, hkc_hka, hkc_rta);
            if ityp == 2 {
                let gcc = p.gccon;
                let h = hka - 1.0 - gcc / rta;
                if h < 0.01 {
                    hkc = 0.01;
                    hkc_hka = 0.0;
                    hkc_rta = 0.0;
                } else {
                    hkc = h;
                    hkc_hka = 1.0;
                    hkc_rta = gcc / (rta * rta);
                }
            } else {
                hkc = hka - 1.0;
                hkc_hka = 1.0;
                hkc_rta = 0.0;
            }
            let hr = hkc / (p.gacon * ald * hka);
            let hr_hka = hkc_hka / (p.gacon * ald * hka) - hr / hka;
            let hr_rta = hkc_rta / (p.gacon * ald * hka);

            let uq = (0.5 * cfa - hr * hr) / (p.gbcon * da);
            let uq_hka = -2.0 * hr * hr_hka / (p.gbcon * da);
            let uq_rta = -2.0 * hr * hr_rta / (p.gbcon * da);
            let uq_cfa = 0.5 / (p.gbcon * da);
            let uq_da = -uq / da;
            let uq_upw = uq_cfa * (c2.cf - c1.cf) + uq_hka * (c2.hk - c1.hk);

            let mut uq_t1 = (1.0 - upw) * (uq_cfa * c1.cf_t + uq_hka * c1.hk_t) + uq_upw * upw_t1;
            let mut uq_d1 = (1.0 - upw) * (uq_cfa * c1.cf_d + uq_hka * c1.hk_d) + uq_upw * upw_d1;
            let mut uq_u1 = (1.0 - upw) * (uq_cfa * c1.cf_u + uq_hka * c1.hk_u) + uq_upw * upw_u1;
            let mut uq_t2 = upw * (uq_cfa * c2.cf_t + uq_hka * c2.hk_t) + uq_upw * upw_t2;
            let mut uq_d2 = upw * (uq_cfa * c2.cf_d + uq_hka * c2.hk_d) + uq_upw * upw_d2;
            let mut uq_u2 = upw * (uq_cfa * c2.cf_u + uq_hka * c2.hk_u) + uq_upw * upw_u2;
            let mut uq_ms = (1.0 - upw) * (uq_cfa * c1.cf_ms + uq_hka * c1.hk_ms)
                + uq_upw * upw_ms
                + upw * (uq_cfa * c2.cf_ms + uq_hka * c2.hk_ms);
            let mut uq_re = (1.0 - upw) * uq_cfa * c1.cf_re + upw * uq_cfa * c2.cf_re;

            uq_t1 += 0.5 * uq_rta * c1.rt_t;
            uq_d1 += 0.5 * uq_da;
            uq_u1 += 0.5 * uq_rta * c1.rt_u;
            uq_t2 += 0.5 * uq_rta * c2.rt_t;
            uq_d2 += 0.5 * uq_da;
            uq_u2 += 0.5 * uq_rta * c2.rt_u;
            uq_ms += 0.5 * uq_rta * c1.rt_ms + 0.5 * uq_rta * c2.rt_ms;
            uq_re += 0.5 * uq_rta * c1.rt_re + 0.5 * uq_rta * c2.rt_re;
            let _ = (uq_t1, uq_d1, uq_u1, uq_t2, uq_d2, uq_u2, uq_ms, uq_re);

            let scc = p.sccon * 1.333 / (1.0 + usa);
            let scc_usa = -scc / (1.0 + usa);

            let slog = (c2.s / c1.s).ln();
            let dxi = c2.x - c1.x;

            let rezc = scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog + dea * 2.0 * (uq * dxi - ulog) * p.duxcon;

            let z_cfa = dea * 2.0 * uq_cfa * dxi * p.duxcon;
            let z_hka = dea * 2.0 * uq_hka * dxi * p.duxcon;
            let z_da = dea * 2.0 * uq_da * dxi * p.duxcon;
            let z_sl = -dea * 2.0;
            let z_ul = -dea * 2.0 * p.duxcon;
            let z_dxi = scc * (cqa - sa * ald) + dea * 2.0 * uq * p.duxcon;
            let z_usa = scc_usa * (cqa - sa * ald) * dxi;
            let z_cqa = scc * dxi;
            let z_sa = -scc * dxi * ald;
            let z_dea = 2.0 * ((uq * dxi - ulog) * p.duxcon - slog);

            let z_upw = z_cqa * (c2.cq - c1.cq) + z_sa * (c2.s - c1.s) + z_cfa * (c2.cf - c1.cf) + z_hka * (c2.hk - c1.hk);
            let z_de1 = 0.5 * z_dea;
            let z_de2 = 0.5 * z_dea;
            let z_us1 = 0.5 * z_usa;
            let z_us2 = 0.5 * z_usa;
            let z_d1 = 0.5 * z_da;
            let z_d2 = 0.5 * z_da;
            let z_u1 = -z_ul / c1.u;
            let z_u2 = z_ul / c2.u;
            let z_x1 = -z_dxi;
            let z_x2 = z_dxi;
            let z_s1 = (1.0 - upw) * z_sa - z_sl / c1.s;
            let z_s2 = upw * z_sa + z_sl / c2.s;
            let z_cq1 = (1.0 - upw) * z_cqa;
            let z_cq2 = upw * z_cqa;
            let z_cf1 = (1.0 - upw) * z_cfa;
            let z_cf2 = upw * z_cfa;
            let z_hk1 = (1.0 - upw) * z_hka;
            let z_hk2 = upw * z_hka;

            self.vs1[0][0] = z_s1;
            self.vs1[0][1] = z_upw * upw_t1 + z_de1 * c1.de_t + z_us1 * c1.us_t;
            self.vs1[0][2] = z_d1 + z_upw * upw_d1 + z_de1 * c1.de_d + z_us1 * c1.us_d;
            self.vs1[0][3] = z_u1 + z_upw * upw_u1 + z_de1 * c1.de_u + z_us1 * c1.us_u;
            self.vs1[0][4] = z_x1;
            self.vs2[0][0] = z_s2;
            self.vs2[0][1] = z_upw * upw_t2 + z_de2 * c2.de_t + z_us2 * c2.us_t;
            self.vs2[0][2] = z_d2 + z_upw * upw_d2 + z_de2 * c2.de_d + z_us2 * c2.us_d;
            self.vs2[0][3] = z_u2 + z_upw * upw_u2 + z_de2 * c2.de_u + z_us2 * c2.us_u;
            self.vs2[0][4] = z_x2;
            self.vsm[0] = z_upw * upw_ms + z_de1 * c1.de_ms + z_us1 * c1.us_ms + z_de2 * c2.de_ms + z_us2 * c2.us_ms;

            self.vs1[0][1] += z_cq1 * c1.cq_t + z_cf1 * c1.cf_t + z_hk1 * c1.hk_t;
            self.vs1[0][2] += z_cq1 * c1.cq_d + z_cf1 * c1.cf_d + z_hk1 * c1.hk_d;
            self.vs1[0][3] += z_cq1 * c1.cq_u + z_cf1 * c1.cf_u + z_hk1 * c1.hk_u;
            self.vs2[0][1] += z_cq2 * c2.cq_t + z_cf2 * c2.cf_t + z_hk2 * c2.hk_t;
            self.vs2[0][2] += z_cq2 * c2.cq_d + z_cf2 * c2.cf_d + z_hk2 * c2.hk_d;
            self.vs2[0][3] += z_cq2 * c2.cq_u + z_cf2 * c2.cf_u + z_hk2 * c2.hk_u;

            self.vsm[0] += z_cq1 * c1.cq_ms + z_cf1 * c1.cf_ms + z_hk1 * c1.hk_ms + z_cq2 * c2.cq_ms + z_cf2 * c2.cf_ms + z_hk2 * c2.hk_ms;
            self.vsr[0] = z_cq1 * c1.cq_re + z_cf1 * c1.cf_re + z_cq2 * c2.cq_re + z_cf2 * c2.cf_re;
            self.vsx[0] = 0.0;
            self.vsrez[0] = -rezc;
        }

        // **** momentum equation
        let ha = 0.5 * (c1.h + c2.h);
        let ma = 0.5 * (c1.m + c2.m);
        let xa = 0.5 * (c1.x + c2.x);
        let ta = 0.5 * (c1.t + c2.t);
        let hwa = 0.5 * (c1.dw / c1.t + c2.dw / c2.t);

        let cfm = self.cfm;
        let cfx = 0.50 * cfm * xa / ta + 0.25 * (c1.cf * c1.x / c1.t + c2.cf * c2.x / c2.t);
        let cfx_xa = 0.50 * cfm / ta;
        let cfx_ta = -0.50 * cfm * xa / (ta * ta);
        let cfx_x1 = 0.25 * c1.cf / c1.t + cfx_xa * 0.5;
        let cfx_x2 = 0.25 * c2.cf / c2.t + cfx_xa * 0.5;
        let cfx_t1 = -0.25 * c1.cf * c1.x / (c1.t * c1.t) + cfx_ta * 0.5;
        let cfx_t2 = -0.25 * c2.cf * c2.x / (c2.t * c2.t) + cfx_ta * 0.5;
        let cfx_cf1 = 0.25 * c1.x / c1.t;
        let cfx_cf2 = 0.25 * c2.x / c2.t;
        let cfx_cfm = 0.50 * xa / ta;

        let btmp = ha + 2.0 - ma + hwa;
        let rezt = tlog + btmp * ulog - xlog * 0.5 * cfx;
        let z_cfx = -xlog * 0.5;
        let z_ha = ulog;
        let z_hwa = ulog;
        let z_ma = -ulog;
        let z_xl = -ddlog * 0.5 * cfx;
        let z_ul = ddlog * btmp;
        let z_tl = ddlog;

        let z_cfm = z_cfx * cfx_cfm;
        let z_cf1 = z_cfx * cfx_cf1;
        let z_cf2 = z_cfx * cfx_cf2;

        let z_t1 = -z_tl / c1.t + z_cfx * cfx_t1 + z_hwa * 0.5 * (-c1.dw / (c1.t * c1.t));
        let z_t2 = z_tl / c2.t + z_cfx * cfx_t2 + z_hwa * 0.5 * (-c2.dw / (c2.t * c2.t));
        let z_x1 = -z_xl / c1.x + z_cfx * cfx_x1;
        let z_x2 = z_xl / c2.x + z_cfx * cfx_x2;
        let z_u1 = -z_ul / c1.u;
        let z_u2 = z_ul / c2.u;

        self.vs1[1][1] = 0.5 * z_ha * c1.h_t + z_cfm * self.cfm_t1 + z_cf1 * c1.cf_t + z_t1;
        self.vs1[1][2] = 0.5 * z_ha * c1.h_d + z_cfm * self.cfm_d1 + z_cf1 * c1.cf_d;
        self.vs1[1][3] = 0.5 * z_ma * c1.m_u + z_cfm * self.cfm_u1 + z_cf1 * c1.cf_u + z_u1;
        self.vs1[1][4] = z_x1;
        self.vs2[1][1] = 0.5 * z_ha * c2.h_t + z_cfm * self.cfm_t2 + z_cf2 * c2.cf_t + z_t2;
        self.vs2[1][2] = 0.5 * z_ha * c2.h_d + z_cfm * self.cfm_d2 + z_cf2 * c2.cf_d;
        self.vs2[1][3] = 0.5 * z_ma * c2.m_u + z_cfm * self.cfm_u2 + z_cf2 * c2.cf_u + z_u2;
        self.vs2[1][4] = z_x2;
        self.vsm[1] = 0.5 * z_ma * c1.m_ms + z_cfm * self.cfm_ms + z_cf1 * c1.cf_ms + 0.5 * z_ma * c2.m_ms + z_cf2 * c2.cf_ms;
        self.vsr[1] = z_cfm * self.cfm_re + z_cf1 * c1.cf_re + z_cf2 * c2.cf_re;
        self.vsx[1] = 0.0;
        self.vsrez[1] = -rezt;

        // **** shape parameter equation
        let xot1 = c1.x / c1.t;
        let xot2 = c2.x / c2.t;
        let ha = 0.5 * (c1.h + c2.h);
        let hsa = 0.5 * (c1.hs + c2.hs);
        let hca = 0.5 * (c1.hc + c2.hc);
        let hwa = 0.5 * (c1.dw / c1.t + c2.dw / c2.t);

        let dix = (1.0 - upw) * c1.di * xot1 + upw * c2.di * xot2;
        let cfx = (1.0 - upw) * c1.cf * xot1 + upw * c2.cf * xot2;
        let dix_upw = c2.di * xot2 - c1.di * xot1;
        let cfx_upw = c2.cf * xot2 - c1.cf * xot1;

        let btmp = 2.0 * hca / hsa + 1.0 - ha - hwa;
        let rezh = hlog + btmp * ulog + xlog * (0.5 * cfx - dix);
        let z_cfx = xlog * 0.5;
        let z_dix = -xlog;
        let z_hca = 2.0 * ulog / hsa;
        let z_ha = -ulog;
        let z_hwa = -ulog;
        let z_xl = ddlog * (0.5 * cfx - dix);
        let z_ul = ddlog * btmp;
        let z_hl = ddlog;

        let z_upw = z_cfx * cfx_upw + z_dix * dix_upw;
        let z_hs1 = -hca * ulog / (hsa * hsa) - z_hl / c1.hs;
        let z_hs2 = -hca * ulog / (hsa * hsa) + z_hl / c2.hs;
        let z_cf1 = (1.0 - upw) * z_cfx * xot1;
        let z_cf2 = upw * z_cfx * xot2;
        let z_di1 = (1.0 - upw) * z_dix * xot1;
        let z_di2 = upw * z_dix * xot2;

        let mut z_t1 = (1.0 - upw) * (z_cfx * c1.cf + z_dix * c1.di) * (-xot1 / c1.t);
        let mut z_t2 = upw * (z_cfx * c2.cf + z_dix * c2.di) * (-xot2 / c2.t);
        let z_x1 = (1.0 - upw) * (z_cfx * c1.cf + z_dix * c1.di) / c1.t - z_xl / c1.x;
        let z_x2 = upw * (z_cfx * c2.cf + z_dix * c2.di) / c2.t + z_xl / c2.x;
        let z_u1 = -z_ul / c1.u;
        let z_u2 = z_ul / c2.u;
        z_t1 += z_hwa * 0.5 * (-c1.dw / (c1.t * c1.t));
        z_t2 += z_hwa * 0.5 * (-c2.dw / (c2.t * c2.t));

        self.vs1[2][0] = z_di1 * c1.di_s;
        self.vs1[2][1] = z_hs1 * c1.hs_t + z_cf1 * c1.cf_t + z_di1 * c1.di_t + z_t1;
        self.vs1[2][2] = z_hs1 * c1.hs_d + z_cf1 * c1.cf_d + z_di1 * c1.di_d;
        self.vs1[2][3] = z_hs1 * c1.hs_u + z_cf1 * c1.cf_u + z_di1 * c1.di_u + z_u1;
        self.vs1[2][4] = z_x1;
        self.vs2[2][0] = z_di2 * c2.di_s;
        self.vs2[2][1] = z_hs2 * c2.hs_t + z_cf2 * c2.cf_t + z_di2 * c2.di_t + z_t2;
        self.vs2[2][2] = z_hs2 * c2.hs_d + z_cf2 * c2.cf_d + z_di2 * c2.di_d;
        self.vs2[2][3] = z_hs2 * c2.hs_u + z_cf2 * c2.cf_u + z_di2 * c2.di_u + z_u2;
        self.vs2[2][4] = z_x2;
        self.vsm[2] = z_hs1 * c1.hs_ms + z_cf1 * c1.cf_ms + z_di1 * c1.di_ms + z_hs2 * c2.hs_ms + z_cf2 * c2.cf_ms + z_di2 * c2.di_ms;
        self.vsr[2] = z_hs1 * c1.hs_re + z_cf1 * c1.cf_re + z_di1 * c1.di_re + z_hs2 * c2.hs_re + z_cf2 * c2.cf_re + z_di2 * c2.di_re;

        self.vs1[2][1] += 0.5 * (z_hca * c1.hc_t + z_ha * c1.h_t) + z_upw * upw_t1;
        self.vs1[2][2] += 0.5 * (z_hca * c1.hc_d + z_ha * c1.h_d) + z_upw * upw_d1;
        self.vs1[2][3] += 0.5 * (z_hca * c1.hc_u) + z_upw * upw_u1;
        self.vs2[2][1] += 0.5 * (z_hca * c2.hc_t + z_ha * c2.h_t) + z_upw * upw_t2;
        self.vs2[2][2] += 0.5 * (z_hca * c2.hc_d + z_ha * c2.h_d) + z_upw * upw_d2;
        self.vs2[2][3] += 0.5 * (z_hca * c2.hc_u) + z_upw * upw_u2;

        self.vsm[2] += 0.5 * (z_hca * c1.hc_ms) + z_upw * upw_ms + 0.5 * (z_hca * c2.hc_ms);
        self.vsx[2] = 0.0;
        self.vsrez[2] = -rezh;
    }
}

/// Solve small dense system in place (GAUSS with NRHS=1). z is n x n, r length n.
pub fn gauss4(z: &mut [[f64; 5]; 4], r: &mut [f64; 4]) {
    let nn = 4;
    for np in 0..nn - 1 {
        let np1 = np + 1;
        let mut nx = np;
        for n in np1..nn {
            if z[n][np].abs() - z[nx][np].abs() > 0.0 {
                nx = n;
            }
        }
        let pivot = 1.0 / z[nx][np];
        z[nx][np] = z[np][np];
        for l in np1..nn {
            let temp = z[nx][l] * pivot;
            z[nx][l] = z[np][l];
            z[np][l] = temp;
        }
        let temp = r[nx] * pivot;
        r[nx] = r[np];
        r[np] = temp;
        for k in np1..nn {
            let ztmp = z[k][np];
            for l in np1..nn {
                z[k][l] -= ztmp * z[np][l];
            }
            r[k] -= ztmp * r[np];
        }
    }
    r[nn - 1] /= z[nn - 1][nn - 1];
    for np in (0..nn - 1).rev() {
        for k in np + 1..nn {
            r[np] -= z[np][k] * r[k];
        }
    }
}
