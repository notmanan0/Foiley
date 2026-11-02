//! QDES mixed-inverse surface-speed design menu (port of xqdes.f).
//!
//! `QdesState` also holds the XFOIL.INC inverse-design variables that are
//! shared between QDES and MDES (QSPEC, SSPEC, QGAMM, ALQSP, IQ1/IQ2, ...),
//! exactly as the Fortran shares them through COMMON.
//!
//! Graphics are not available: plot-only commands are accepted and ignored,
//! cursor-driven commands (MODI, WIND, MARK without arguments) print a note.
//! Extension: `MARK i1 i2 [k]` sets the target segment from node indices
//! (the Fortran gets these by mouse clicks), which makes `eXec` scriptable.

use crate::cli::{args, Session};
use crate::fmt::{e_fmt, f_fmt, fixed};
use crate::geom::{atanc, lefind};
use crate::panel::ncalc;
use crate::spline::{scalc, seval, splind, trisol};
use crate::xfoil::{clcalc, cpcalc, PsiOut, XFoil, DTOR, IBX};

/// Max number of Qspec distributions (IPX).
pub const IPX: usize = 5;

/// Persistent inverse-design state (XFOIL.INC design variables).
#[derive(Clone, Debug)]
pub struct QdesState {
    pub nsp: usize,
    pub nqsp: usize,
    pub iq1: usize,
    pub iq2: usize,
    pub kqtarg: usize,
    pub iacqsp: i32,
    pub liqset: bool,
    pub lqsym: bool,
    pub lqslop: bool,
    pub lcpxx: bool,
    pub lqvdes: bool,
    pub lqrefl: bool,
    pub sspec: Vec<f64>,
    pub xspoc: Vec<f64>,
    pub yspoc: Vec<f64>,
    /// `qspec[k][i]`, k = 1..IPX.  `qspec[0]` is QGAMM: in the Fortran COMMON
    /// block QGAMM immediately precedes QSPEC, so QSPEC(i,0) == QGAMM(i).
    pub qspec: Vec<Vec<f64>>,
    pub qspecp: Vec<Vec<f64>>,
    pub ssple: f64,
    pub algam: f64,
    pub clgam: f64,
    pub cmgam: f64,
    pub alqsp: [f64; IPX + 1],
    pub clqsp: [f64; IPX + 1],
    pub cmqsp: [f64; IPX + 1],
    pub qdof0: f64,
    pub qdof1: f64,
    pub qdof2: f64,
    pub qdof3: f64,
    /// last iteration count given to eXec
    pub niterq: i64,
    /// plot-object size (plotting only, kept for SIZE)
    pub size: f64,
}

impl Default for QdesState {
    fn default() -> Self {
        QdesState {
            nsp: 0,
            nqsp: 0,
            iq1: 0,
            iq2: 0,
            kqtarg: 0,
            iacqsp: 0,
            liqset: false,
            lqsym: false,
            lqslop: true,
            lcpxx: false,
            lqvdes: false,
            lqrefl: false,
            sspec: vec![0.0; IBX + 1],
            xspoc: vec![0.0; IBX + 1],
            yspoc: vec![0.0; IBX + 1],
            qspec: vec![vec![0.0; IBX + 1]; IPX + 1],
            qspecp: vec![vec![0.0; IBX + 1]; IPX + 1],
            ssple: 0.0,
            algam: 0.0,
            clgam: 0.0,
            cmgam: 0.0,
            alqsp: [0.0; IPX + 1],
            clqsp: [0.0; IPX + 1],
            cmqsp: [0.0; IPX + 1],
            qdof0: 0.0,
            qdof1: 0.0,
            qdof2: 0.0,
            qdof3: 0.0,
            niterq: 10,
            size: 10.0,
        }
    }
}

impl QdesState {
    /// QGAMM array (aliases QSPEC(.,0)).
    #[inline]
    pub fn qgamm(&self) -> &Vec<f64> {
        &self.qspec[0]
    }

    /// SPLQSP: spline Qspec(s), end intervals splined separately.
    pub fn splqsp(&mut self, kqsp: usize) {
        let nsp = self.nsp;
        let sspec = &self.sspec;
        let qs = &self.qspec[kqsp];
        let qp = &mut self.qspecp[kqsp];
        splind(&qs[1..], &mut qp[1..], &sspec[1..], nsp - 2, -999.0, -999.0);
        let i = 1;
        let xs2 = qp[i + 1];
        splind(&qs[i - 1..], &mut qp[i - 1..], &sspec[i - 1..], 2, -999.0, xs2);
        let i = nsp - 1;
        let xs1 = qp[i];
        splind(&qs[i - 1..], &mut qp[i - 1..], &sspec[i - 1..], 2, xs1, -999.0);
    }

    /// GAMQSP: set Qspec(s,k) from current speed Q(s).
    pub fn gamqsp(&mut self, kqsp: usize) {
        self.alqsp[kqsp] = self.algam;
        self.clqsp[kqsp] = self.clgam;
        self.cmqsp[kqsp] = self.cmgam;
        for i in 1..=self.nsp {
            self.qspec[kqsp][i] = self.qspec[0][i];
        }
        self.qdof0 = 0.0;
        self.qdof1 = 0.0;
        self.qdof2 = 0.0;
        self.qdof3 = 0.0;
        self.splqsp(kqsp);
        if !self.liqset {
            self.iq1 = 1;
            self.iq2 = self.nsp;
        }
    }

    /// SYMQSP: force symmetry of Qspec(kqsp) (unused by the menus, kept for completeness).
    pub fn symqsp(&mut self, kqsp: usize) {
        let nsp = self.nsp;
        self.alqsp[kqsp] = 0.0;
        self.clqsp[kqsp] = 0.0;
        self.cmqsp[kqsp] = 0.0;
        let sspmid = 0.5 * (self.sspec[nsp] - self.sspec[1]);
        for i in 1..=(nsp + 1) / 2 {
            self.sspec[i] = sspmid + 0.5 * (self.sspec[i] - self.sspec[nsp - i + 1]);
            self.qspec[kqsp][i] = 0.5 * (self.qspec[kqsp][i] - self.qspec[kqsp][nsp - i + 1]);
        }
        for i in (nsp + 1) / 2 + 1..=nsp {
            self.sspec[i] = -self.sspec[nsp - i + 1] + 2.0 * sspmid;
            self.qspec[kqsp][i] = -self.qspec[kqsp][nsp - i + 1];
        }
        self.qdof0 = 0.0;
        self.qdof1 = 0.0;
        self.qdof2 = 0.0;
        self.qdof3 = 0.0;
        self.splqsp(kqsp);
        crate::outln!("\n Qspec{:2}  made symmetric", kqsp);
    }

    /// SMOOQ: smooth Qspec(s) inside target segment kq1..kq2.
    pub fn smooq(&mut self, kq1: usize, kq2: usize, kqsp: usize) {
        let nsp = self.nsp;
        let mut w1 = vec![0.0; nsp + 2];
        let mut w2 = vec![0.0; nsp + 2];
        let mut w3 = vec![0.0; nsp + 2];
        let w8: Vec<f64> = self.sspec[..=nsp].to_vec();

        if kq1 < 1 || kq2 < kq1 + 2 {
            crate::outln!(" Segment is too short.  No smoothing possible.");
            return;
        }

        let smool = 0.002 * (w8[nsp] - w8[1]);
        let smoosq = smool * smool;
        for i in kq1 + 1..kq2 {
            let dsm = w8[i] - w8[i - 1];
            let dsp = w8[i + 1] - w8[i];
            let dso = 0.5 * (w8[i + 1] - w8[i - 1]);
            w1[i] = smoosq * (-1.0 / dsm) / dso;
            w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            w3[i] = smoosq * (-1.0 / dsp) / dso;
        }

        w2[kq1] = 1.0;
        w3[kq1] = 0.0;
        w1[kq2] = 0.0;
        w2[kq2] = 1.0;

        let qs = &mut self.qspec[kqsp];
        if self.lqslop {
            let i = kq1 + 1;
            let dsm = w8[i] - w8[i - 1];
            let dsp = w8[i + 1] - w8[i];
            let ds = w8[i + 1] - w8[i - 1];
            w1[i] = -1.0 / dsm - (dsm / ds) / dsm;
            w2[i] = 1.0 / dsm + (dsm / ds) / dsm + (dsm / ds) / dsp;
            w3[i] = -(dsm / ds) / dsp;
            let qspp1 = w1[i] * qs[i - 1] + w2[i] * qs[i] + w3[i] * qs[i + 1];

            let i = kq2 - 1;
            let dsm = w8[i] - w8[i - 1];
            let dsp = w8[i + 1] - w8[i];
            let ds = w8[i + 1] - w8[i - 1];
            w1[i] = (dsp / ds) / dsm;
            w2[i] = -1.0 / dsp - (dsp / ds) / dsp - (dsp / ds) / dsm;
            w3[i] = 1.0 / dsp + (dsp / ds) / dsp;
            let qspp2 = w1[i] * qs[i - 1] + w2[i] * qs[i] + w3[i] * qs[i + 1];

            qs[kq1 + 1] = qspp1;
            qs[kq2 - 1] = qspp2;
        }

        let o = kq1 - 1;
        trisol(&mut w2[o..], &w1[o..], &mut w3[o..], &mut qs[o..], kq2 - kq1 + 1);
    }
}

/// QINCOM: incompressible speed from Karman-Tsien compressible speed.
pub fn qincom(qc: f64, qinf: f64, tklam: f64) -> f64 {
    if tklam < 1.0e-4 || qc.abs() < 1.0e-4 {
        qc / (1.0 - tklam)
    } else {
        let tmp = 0.5 * (1.0 - tklam) * qinf / (qc * tklam);
        qinf * tmp * ((1.0 + 1.0 / (tklam * tmp * tmp)).sqrt() - 1.0)
    }
}

/// GAUSS: dense solve with partial pivoting (Z row-major, 1-based, leading dim `ld`).
pub fn gauss(nn: usize, z: &mut [f64], ld: usize, r: &mut [f64]) {
    for np in 1..nn {
        let np1 = np + 1;
        let mut nx = np;
        for n in np1..=nn {
            if z[n * ld + np].abs() - z[nx * ld + np].abs() > 0.0 {
                nx = n;
            }
        }
        let pivot = 1.0 / z[nx * ld + np];
        z[nx * ld + np] = z[np * ld + np];
        for l in np1..=nn {
            let temp = z[nx * ld + l] * pivot;
            z[nx * ld + l] = z[np * ld + l];
            z[np * ld + l] = temp;
        }
        let temp = r[nx] * pivot;
        r[nx] = r[np];
        r[np] = temp;
        for k in np1..=nn {
            let ztmp = z[k * ld + np];
            for l in np1..=nn {
                z[k * ld + l] -= ztmp * z[np * ld + l];
            }
            r[k] -= ztmp * r[np];
        }
    }
    r[nn] /= z[nn * ld + nn];
    for np in (1..nn).rev() {
        for k in np + 1..=nn {
            r[np] -= z[np * ld + k] * r[k];
        }
    }
}

/// Target-segment marking.  The Fortran IQSGET uses cursor clicks; here
/// node indices may be given as arguments (`MARK i1 i2 [kqsp]`).
pub fn iqsget(q: &mut QdesState, iv: &[i64]) {
    if iv.len() < 2 {
        crate::outln!(" MARK needs cursor input (not available).  Use  MARK i1 i2 [kqsp]  with node indices.");
        return;
    }
    let nsp = q.nsp as i64;
    let clamp = |v: i64| v.max(1).min(nsp.max(1)) as usize;
    let i1 = clamp(iv[0]);
    let i2 = clamp(iv[1]);
    let k = if iv.len() >= 3 { iv[2].max(1).min(q.nqsp.max(1) as i64) as usize } else { 1 };
    q.kqtarg = k;
    if i1 == i2 {
        crate::outln!(" ***  Endpoints must be distinct  ***");
        crate::outln!(" ***  NEW SEGMENT NOT MARKED OFF  ***");
        return;
    }
    q.iq1 = i1.min(i2);
    q.iq2 = i1.max(i2);
    q.liqset = true;
}

/// GAMLIN: Jacobian entry for row i due to node speed GAM(j).
#[inline]
#[allow(clippy::too_many_arguments)]
fn gamlin(qm: &mut [f64], ld: usize, n: usize, xf: &XFoil, iq1: usize, iq2: usize, i: usize, j: usize, coef: f64) {
    if j >= iq1 && j <= iq2 {
        qm[i * ld + n + 2] += coef * xf.qf0[j];
        qm[i * ld + n + 3] += coef * xf.qf1[j];
        qm[i * ld + n + 4] += coef * xf.qf2[j];
        qm[i * ld + n + 5] += coef * xf.qf3[j];
    } else {
        qm[i * ld + j] += coef;
    }
}

/// MIXED: mixed-inverse calculation using QSPEC(.,kqsp).
pub fn mixed(xf: &mut XFoil, q: &mut QdesState, kqsp: usize, niterq: i64) {
    let bwt = 0.1;
    let n = xf.n;
    let iq1 = q.iq1;
    let iq2 = q.iq2;

    xf.cosa = xf.alfa.cos();
    xf.sina = xf.alfa.sin();
    scalc(&xf.x, &xf.y, &mut xf.s, n);

    for i in 1..=n {
        xf.qf0[i] = 0.0;
        xf.qf1[i] = 0.0;
        xf.qf2[i] = 0.0;
        xf.qf3[i] = 0.0;
    }
    for i in iq1..=iq2 {
        let fs = (xf.s[i] - xf.s[iq1]) / (xf.s[iq2] - xf.s[iq1]);
        xf.qf0[i] = 1.0 - fs;
        xf.qf1[i] = fs;
        if q.lcpxx {
            xf.qf2[i] = (-5.0 * fs).exp();
            xf.qf3[i] = (-5.0 * (1.0 - fs)).exp();
        } else {
            xf.qf2[i] = 0.0;
            xf.qf3[i] = 0.0;
        }
        xf.gam[i] = q.qspec[kqsp][i] + q.qdof0 * xf.qf0[i] + q.qdof1 * xf.qf1[i] + q.qdof2 * xf.qf2[i] + q.qdof3 * xf.qf3[i];
    }

    let ld = n + 6;
    let mut qm = vec![0.0; (n + 6) * ld];
    let mut dq = vec![0.0; n + 6];
    let mut o = PsiOut::new();

    for _iter in 1..=niterq {
        for v in qm.iter_mut() {
            *v = 0.0;
        }

        ncalc(&xf.x, &xf.y, &xf.s, n, &mut xf.nx, &mut xf.ny);

        for i in 1..=n {
            xf.psilin(i, xf.x[i], xf.y[i], xf.nx[i], xf.ny[i], true, false, &mut o);
            let psi = o.psi;
            let psi_n = o.psi_ni;
            o.dzdn[i] += psi_n;

            let row = i * ld;
            for j in 1..iq1 {
                qm[row + j] += o.dzdg[j];
            }
            for j in iq1..=iq2 {
                qm[row + j] += o.dzdn[j];
            }
            for j in iq2 + 1..=n {
                qm[row + j] += o.dzdg[j];
            }

            dq[i] = xf.psio - psi;

            qm[row + n + 1] -= 1.0;
            qm[row + n + 2] += o.z_qdof0;
            qm[row + n + 3] += o.z_qdof1;
            qm[row + n + 4] += o.z_qdof2;
            qm[row + n + 5] += o.z_qdof3;
        }

        // Kutta condition
        dq[n + 1] = -(xf.gam[1] + xf.gam[n]);
        gamlin(&mut qm, ld, n, xf, iq1, iq2, n + 1, 1, 1.0);
        gamlin(&mut qm, ld, n, xf, iq1, iq2, n + 1, n, 1.0);

        if xf.sharp {
            // zero internal velocity in TE corner
            let ag1 = (-xf.yp[1]).atan2(-xf.xp[1]);
            let ag2 = atanc(xf.yp[n], xf.xp[n], ag1);
            let abis = 0.5 * (ag1 + ag2);
            let cbis = abis.cos();
            let sbis = abis.sin();

            let ds1 = ((xf.x[1] - xf.x[2]).powi(2) + (xf.y[1] - xf.y[2]).powi(2)).sqrt();
            let ds2 = ((xf.x[n] - xf.x[n - 1]).powi(2) + (xf.y[n] - xf.y[n - 1]).powi(2)).sqrt();
            let dsmin = ds1.min(ds2);

            let xbis = xf.xte - bwt * dsmin * cbis;
            let ybis = xf.yte - bwt * dsmin * sbis;

            xf.psilin(0, xbis, ybis, -sbis, cbis, false, true, &mut o);
            let qbis = o.psi_ni;
            let res = qbis;

            for j in 1..=n + 5 {
                qm[n * ld + j] = 0.0;
            }
            for j in 1..=n {
                gamlin(&mut qm, ld, n, xf, iq1, iq2, n, j, o.dqdg[j]);
                qm[n * ld + j] = o.dqdg[j];
            }
            qm[n * ld + n + 1] = 0.0;
            dq[n] = -res;
        }

        // pinned IQ1, IQ2 point conditions
        qm[(n + 2) * ld + iq1] = 1.0;
        dq[n + 2] = 0.0;
        qm[(n + 3) * ld + iq2] = 1.0;
        dq[n + 3] = 0.0;

        let qs = &q.qspec[kqsp];
        if iq1 > 1 && q.lcpxx {
            let res = xf.gam[iq1 - 1] - 2.0 * xf.gam[iq1] + xf.gam[iq1 + 1] - (qs[iq1 - 1] - 2.0 * qs[iq1] + qs[iq1 + 1]);
            gamlin(&mut qm, ld, n, xf, iq1, iq2, n + 4, iq1 - 1, 1.0);
            gamlin(&mut qm, ld, n, xf, iq1, iq2, n + 4, iq1, -2.0);
            gamlin(&mut qm, ld, n, xf, iq1, iq2, n + 4, iq1 + 1, 1.0);
            dq[n + 4] = -res;
        } else {
            qm[(n + 4) * ld + n + 4] = 1.0;
            dq[n + 4] = -q.qdof2;
        }

        if iq2 < n && q.lcpxx {
            let res = xf.gam[iq2 - 1] - 2.0 * xf.gam[iq2] + xf.gam[iq2 + 1] - (qs[iq2 - 1] - 2.0 * qs[iq2] + qs[iq2 + 1]);
            gamlin(&mut qm, ld, n, xf, iq1, iq2, n + 5, iq2 - 1, 1.0);
            gamlin(&mut qm, ld, n, xf, iq1, iq2, n + 5, iq2, -2.0);
            gamlin(&mut qm, ld, n, xf, iq1, iq2, n + 5, iq2 + 1, 1.0);
            dq[n + 5] = -res;
        } else {
            qm[(n + 5) * ld + n + 5] = 1.0;
            dq[n + 5] = -q.qdof3;
        }

        gauss(n + 5, &mut qm, ld, &mut dq);

        let mut inmax = 0usize;
        let mut igmax = 0usize;
        let mut dnmax = 0.0f64;
        let mut dgmax = 0.0f64;

        for i in 1..iq1 {
            xf.gam[i] += dq[i];
            if dq[i].abs() > dgmax.abs() {
                dgmax = dq[i];
                igmax = i;
            }
        }
        for i in iq1..=iq2 {
            xf.x[i] += xf.nx[i] * dq[i];
            xf.y[i] += xf.ny[i] * dq[i];
            if dq[i].abs() > dnmax.abs() {
                dnmax = dq[i];
                inmax = i;
            }
        }
        for i in iq2 + 1..=n {
            xf.gam[i] += dq[i];
            if dq[i].abs() > dgmax.abs() {
                dgmax = dq[i];
                igmax = i;
            }
        }

        xf.psio += dq[n + 1];
        q.qdof0 += dq[n + 2];
        q.qdof1 += dq[n + 3];
        q.qdof2 += dq[n + 4];
        q.qdof3 += dq[n + 5];

        xf.cosa = xf.alfa.cos();
        xf.sina = xf.alfa.sin();
        scalc(&xf.x, &xf.y, &mut xf.s, n);

        for i in iq1..=iq2 {
            xf.gam[i] = q.qspec[kqsp][i] + q.qdof0 * xf.qf0[i] + q.qdof1 * xf.qf1[i] + q.qdof2 * xf.qf2[i] + q.qdof3 * xf.qf3[i];
        }

        xf.tecalc();
        xf.clcalc();
        crate::outln!(
            "\n dNmax ={}{:4}   dQmax ={}{:4}    CL ={}\n dQf1  ={}       dQf2  ={}\n dQf3  ={}       dQf4  ={}",
            e_fmt(dnmax, 10, 3),
            inmax,
            e_fmt(dgmax, 10, 3),
            igmax,
            f_fmt(xf.cl, 7, 4),
            e_fmt(dq[n + 2], 10, 3),
            e_fmt(dq[n + 3], 10, 3),
            e_fmt(dq[n + 4], 10, 3),
            e_fmt(dq[n + 5], 10, 3)
        );

        if dnmax.abs() < 5.0e-5 && dgmax.abs() < 5.0e-4 {
            crate::outln!();
            crate::outln!(" New current airfoil generated");
            crate::outln!(" Old buffer  airfoil unchanged");
            return;
        }
    }
    crate::outln!(" Not quite converged.  Can EXEC again if necessary.");
}

/// Re-spline the current airfoil after an in-place node change (EXEC/REST).
fn respline_current(xf: &mut XFoil) {
    let n = xf.n;
    scalc(&xf.x, &xf.y, &mut xf.s, n);
    splind(&xf.x, &mut xf.xp, &xf.s, n, -999.0, -999.0);
    splind(&xf.y, &mut xf.yp, &xf.s, n, -999.0, -999.0);
    ncalc(&xf.x, &xf.y, &xf.s, n, &mut xf.nx, &mut xf.ny);
    xf.sle = lefind(&xf.x, &xf.xp, &xf.y, &xf.yp, &xf.s, n);
    xf.xle = seval(xf.sle, &xf.x, &xf.xp, &xf.s, n);
    xf.yle = seval(xf.sle, &xf.y, &xf.yp, &xf.s, n);
    xf.chord = ((0.5 * (xf.x[1] + xf.x[n]) - xf.xle).powi(2) + (0.5 * (xf.y[1] + xf.y[n]) - xf.yle).powi(2)).sqrt();
    xf.tecalc();
    xf.apcalc();
}

/// gfortran list-directed REAL(8) item (field incl. separator blank, 26 chars).
pub fn ld_real(v: f64) -> String {
    let a = v.abs();
    if v == 0.0 || (a >= 0.1 && a < 1.0e17) {
        let k = if v == 0.0 { 1 } else { a.log10().floor() as i32 + 1 };
        let mut dec = (17 - k).max(0) as usize;
        let mut s = format!("{:.*}", dec, v);
        // rounding may add a digit
        if v != 0.0 && s.trim_start_matches('-').split('.').next().unwrap_or("").len() > k.max(1) as usize && k >= 1 {
            dec = dec.saturating_sub(1);
            s = format!("{:.*}", dec, v);
        }
        format!(" {:>20}     ", s)
    } else {
        let mut e = a.log10().floor() as i32;
        let mut m = v / 10f64.powi(e);
        let mut ms = format!("{:.16}", m);
        if ms.trim_start_matches('-').starts_with("10") {
            e += 1;
            m = v / 10f64.powi(e);
            ms = format!("{:.16}", m);
        }
        let es = format!("E{}{:03}", if e < 0 { '-' } else { '+' }, e.abs());
        format!(" {:>25}", format!("{}{}", ms, es))
    }
}

/// gfortran list-directed default INTEGER item (separator blank + I11).
pub fn ld_int(i: i64) -> String {
    format!(" {:>11}", i)
}

/// Statement function QCOMP (Karman-Tsien compressible speed).
#[inline]
pub fn qcomp(g: f64, tklam: f64, qinf: f64) -> f64 {
    g * (1.0 - tklam) / (1.0 - tklam * (g / qinf) * (g / qinf))
}

fn help() {
    crate::outln!(
        "\n   <cr>   Return to Top Level\n\n   QSET   Reset Qspec <== Q\n\n   Modi   Modify Qspec\n   MARK   Mark off target segment\n   SMOO   Smooth Qspec inside target segment\n   SLOP   Toggle modified-Qspec slope matching flag\n\n   eXec i Execute mixed-inverse calculation\n   REST   Restore geometry from buffer airfoil\n   CPXX   CPxx endpoint constraint toggle\n\n   Visc   Qvis overlay toggle\n   REFL   Reflected Qspec overlay toggle\n\n   Plot   Plot Qspec (line) and Q (symbols)\n   Blow   Blowup plot region\n   Rese   Reset plot scale and origin\n   Wind     Plot window adjust via cursor and keys\n\n   SIZE r Change absolute plot-object size\n  .ANNO   Annotate plot\n   HARD   Hardcopy current plot\n\n   (Foiley: no cursor input -- use  MARK i1 i2  to set the target segment\n    by node indices; Modi and Wind are not available)"
    );
}

/// SIZE command shared by QDES and MDES.
pub fn size_cmd(s: &mut Session, r: &[f64], msg: &str, prompt: &str) {
    if !r.is_empty() {
        s.qdes_state.size = r[0];
    } else {
        crate::outln!(" {}{}", msg, ld_real(s.qdes_state.size));
        s.qdes_state.size = s.input.askr(prompt, s.qdes_state.size);
    }
}

/// QDES menu entry point.
pub fn qdes(s: &mut Session) {
    if s.xf.n == 0 {
        crate::outln!();
        crate::outln!(" ***  No airfoil available  ***");
        return;
    }

    {
        let xf = &mut s.xf;
        let q = &mut s.qdes_state;

        // make sure a current solution exists
        xf.specal();

        // see if current Qspec, if any, didn't come from Full-Inverse
        if q.nsp != xf.n {
            xf.lqspec = false;
            q.liqset = false;
        }

        q.algam = xf.alfa;
        q.clgam = xf.cl;
        q.cmgam = xf.cm;

        let n = xf.n;
        let chx = xf.xte - xf.xle;
        let chy = xf.yte - xf.yle;
        let chsq = chx * chx + chy * chy;
        q.nsp = n;
        for i in 1..=n {
            q.qspec[0][i] = xf.gam[i];
            q.sspec[i] = xf.s[i] / xf.s[n];
            q.xspoc[i] = ((xf.x[i] - xf.xle) * chx + (xf.y[i] - xf.yle) * chy) / chsq;
            q.yspoc[i] = ((xf.y[i] - xf.yle) * chx - (xf.x[i] - xf.xle) * chy) / chsq;
        }
        q.ssple = xf.sle / xf.s[n];

        crate::outln!(
            "\n Current Q operating condition:\n alpha = {} deg.      CL = {}\n",
            f_fmt(q.algam / DTOR, 8, 3),
            f_fmt(q.clgam, 8, 4)
        );

        if !xf.lqspec {
            q.nqsp = 1;
            q.kqtarg = 1;
            q.gamqsp(1);
            crate::outln!("\n Qspec initialized to current Q.\n");
            xf.lqspec = true;
        }
    }

    let mut comand = String::from("****");
    let mut comarg = String::new();

    loop {
        // label 500
        let comold = comand.clone();
        let argold = comarg.clone();

        // label 501
        loop {
            let (c, a) = s.input.askc(".QDES");
            if s.input.eof {
                s.qdes_state.lqsym = false;
                return;
            }
            if c.starts_with('!') {
                if comold == "****" {
                    crate::outln!(" Previous .QDES command not valid");
                    continue;
                }
                comand = comold.clone();
                comarg = argold.clone();
            } else {
                comand = c;
                comarg = a;
            }
            break;
        }

        if comand.is_empty() {
            s.qdes_state.lqsym = false;
            return;
        }

        let (r, iv) = args(&comarg);
        let ninput = r.len();

        match comand.as_str() {
            "?" => help(),
            "Z" | "U" => {}
            "QSET" => {
                s.qdes_state.gamqsp(1);
            }
            "VISC" | "V" => {
                let q = &mut s.qdes_state;
                q.lqvdes = !q.lqvdes;
                if q.lqvdes {
                    crate::outln!(" Qspec & Qvis will be plotted");
                } else {
                    crate::outln!(" Only Qspec will be plotted");
                }
            }
            "REFL" => {
                let q = &mut s.qdes_state;
                q.lqrefl = !q.lqrefl;
                if q.lqrefl {
                    crate::outln!(" Reflected Qspec will be plotted");
                } else {
                    crate::outln!(" Reflected Qspec will not be plotted");
                }
            }
            "MARK" => iqsget(&mut s.qdes_state, &iv),
            "MODI" | "M" => {
                crate::outln!(" MODI needs cursor input (not available in Foiley)");
            }
            "SMOO" => {
                let xf = &mut s.xf;
                let q = &mut s.qdes_state;
                let kqsp = 1;
                let (iq1, iq2) = (q.iq1, q.iq2);
                q.smooq(iq1, iq2, kqsp);
                q.splqsp(kqsp);
                qspec_clcm_print(xf, q, kqsp);
            }
            "SLOP" => {
                let q = &mut s.qdes_state;
                q.lqslop = !q.lqslop;
                if q.lqslop {
                    crate::outln!(" Modified Qspec piece will be made tangent at endpoints");
                } else {
                    crate::outln!(" Modified Qspec piece will not be made tangent at endpoints");
                }
            }
            "HARD" | "PLOT" | "P" | "BLOW" | "B" | "RESE" | "R" => {}
            "WIND" | "W" => {
                crate::outln!(" ");
                crate::outln!(" Type I,O,P to In,Out,Pan with cursor...");
                crate::outln!(" (cursor input not available)");
            }
            "ANNO" => {
                crate::outln!(" No active plot to annotate");
            }
            "SIZE" => size_cmd(s, &r, "Current plot size =", "Enter new plot size^"),
            "CPXX" => {
                let q = &mut s.qdes_state;
                q.lcpxx = !q.lcpxx;
                if q.lcpxx {
                    crate::outln!(" CPxx will be constrained");
                } else {
                    crate::outln!(" CPxx will not be constrained");
                }
            }
            "EXEC" | "X" => {
                if !s.qdes_state.liqset {
                    crate::outln!(" ***  Must mark off target segment first  ***");
                    continue;
                }
                {
                    let q = &s.qdes_state;
                    let mut ist = 0;
                    for i in q.iq1..q.iq2 {
                        if q.qspec[0][i] >= 0.0 && q.qspec[0][i + 1] < 0.0 {
                            ist = i;
                        }
                    }
                    if ist != 0 {
                        crate::outln!();
                        crate::outln!(" Target segment cannot include stagnation point in mixed-inverse.");
                        continue;
                    }
                }
                let kqsp = 1;
                s.xf.clspec = s.qdes_state.clqsp[kqsp];

                // save current coordinates for restoration if requested
                {
                    let xf = &mut s.xf;
                    for i in 1..=xf.n {
                        xf.xb[i] = xf.x[i];
                        xf.yb[i] = xf.y[i];
                        xf.sb[i] = xf.s[i];
                        xf.xbp[i] = xf.xp[i];
                        xf.ybp[i] = xf.yp[i];
                    }
                    xf.nb = xf.n;
                    xf.lgsame = true;
                }
                crate::outln!();
                crate::outln!(" Current airfoil saved in buffer airfoil");

                let niterq = if ninput >= 1 {
                    iv[0]
                } else {
                    let d = s.qdes_state.niterq;
                    s.input.aski("Enter max number of iterations^", d)
                };
                s.qdes_state.niterq = niterq;

                let xf = &mut s.xf;
                let q = &mut s.qdes_state;
                mixed(xf, q, kqsp, niterq);
                xf.adeg = xf.alfa / DTOR;

                respline_current(xf);

                q.algam = xf.alfa;
                let n = xf.n;
                q.nsp = n;
                for i in 1..=n {
                    q.qspec[0][i] = xf.gam[i];
                    q.sspec[i] = xf.s[i] / xf.s[n];
                }
                q.ssple = xf.sle / xf.s[n];

                for i in 1..=n {
                    xf.qinv[i] = xf.gam[i];
                }
                cpcalc(n, &xf.qinv, xf.qinf, xf.minf, &mut xf.cpi);

                xf.lgamu = false;
                xf.lqinu = false;
                xf.lwake = false;
                xf.lqaij = false;
                xf.ladij = false;
                xf.lwdij = false;
                xf.lipan = false;
                xf.lvconv = false;
                xf.lscini = false;
                xf.lgsame = false;
            }
            "REST" => {
                let xf = &mut s.xf;
                for i in 1..=xf.n {
                    xf.x[i] = xf.xb[i];
                    xf.y[i] = xf.yb[i];
                }
                respline_current(xf);
                xf.lgamu = false;
                xf.lqinu = false;
                xf.lgsame = true;
            }
            _ => {
                crate::outln!(" Command {} not recognized.  Type a \" ? \" for list.", fixed(&comand, 4));
                comand = "****".into();
            }
        }
    }
}

/// CLCALC on Qspec(.,kqsp) and print Q / Qspec CL, CM (format 1200 of QDES).
fn qspec_clcm_print(xf: &XFoil, q: &mut QdesState, kqsp: usize) {
    let w1 = vec![0.0; xf.n + 1];
    let r = clcalc(xf.n, &xf.x, &xf.y, &q.qspec[kqsp], &w1, xf.alfa, xf.minf, xf.qinf, xf.xcmref, xf.ycmref);
    q.clqsp[kqsp] = r.0;
    q.cmqsp[kqsp] = r.1;
    crate::outln!(
        "\n Q    :   CL ={}    CM ={}\n Qspec:   CL ={}    CM ={}",
        f_fmt(xf.cl, 11, 6),
        f_fmt(xf.cm, 11, 6),
        f_fmt(q.clqsp[kqsp], 11, 6),
        f_fmt(q.cmqsp[kqsp], 11, 6)
    );
}
