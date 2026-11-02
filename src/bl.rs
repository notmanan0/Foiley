//! Viscous/inviscid coupling: BL marching, global Newton system setup,
//! block solver and update (port of XFoil xbl.f and BLSOLV from xsolve.f).

use crate::blsys::{gauss4, hkin};
use crate::spline::{seval, sinvrt, splind};
use crate::xfoil::{XFoil, DTOR};

/// Limit DSTR so that Hk >= HKLIM (DSLIM).
#[inline]
pub fn dslim(dstr: &mut f64, thet: f64, msq: f64, hklim: f64) {
    let h = *dstr / thet;
    let (hk, hk_h, _) = hkin(h, msq);
    let dh = (hklim - hk).max(0.0) / hk_h;
    *dstr += dh * thet;
}

/// How the "1"-station D sensitivities D1_M are formed for a station.
#[derive(Clone, Copy, Debug, Default)]
pub enum D1Kind {
    #[default]
    Zero,
    Prev { row: usize, d_u: f64, d_m: f64 },
    Te { ute1: f64, ute2: f64, mte1: f64, mte2: f64 },
}

/// Per-station coefficients of the VM rows (filled by SETBL's sequential pass).
#[derive(Clone, Copy, Debug, Default)]
pub struct StCoef {
    pub co: [[f64; 5]; 3],
    pub xi1: f64,
    pub xi2: f64,
    pub d2_u2: f64,
    pub d2_m2: f64,
    pub d1: D1Kind,
    pub u1row: usize,
}

/// out[k][iv] = sum_{jv=1..nsys} M[iv][jv] * v[k][jv], accumulated left to
/// right from 0.0 exactly like the Fortran loops, with 4 rows interleaved so
/// the independent accumulation chains overlap.
fn matvec_rows(m: &[f64], ld: usize, nsys: usize, v: &[&[f64]], out: &mut [&mut [f64]]) {
    let nv = v.len();
    let mut iv = 1;
    while iv + 3 <= nsys {
        let r: [&[f64]; 4] = [&m[iv * ld..], &m[(iv + 1) * ld..], &m[(iv + 2) * ld..], &m[(iv + 3) * ld..]];
        for k in 0..nv {
            let vk = v[k];
            let mut a = [0.0f64; 4];
            for jv in 1..=nsys {
                let x = vk[jv];
                a[0] += r[0][jv] * x;
                a[1] += r[1][jv] * x;
                a[2] += r[2][jv] * x;
                a[3] += r[3][jv] * x;
            }
            for j in 0..4 {
                out[k][iv + j] = a[j];
            }
        }
        iv += 4;
    }
    while iv <= nsys {
        let r = &m[iv * ld..];
        for k in 0..nv {
            let mut a = 0.0;
            for jv in 1..=nsys {
                a += r[jv] * v[k][jv];
            }
            out[k][iv] = a;
        }
        iv += 1;
    }
}

impl XFoil {
    /// BL station -> Newton system line pointers (IBLSYS).
    pub fn iblsys(&mut self) {
        let mut iv = 0;
        for is in 1..=2 {
            for ibl in 2..=self.nbl[is] {
                iv += 1;
                self.isys[is][ibl] = iv;
            }
        }
        self.nsys = iv;
        let nsys = self.nsys;
        if self.va.len() < nsys + 2 {
            self.va = vec![[[0.0; 2]; 3]; nsys + 2];
            self.vb = vec![[[0.0; 2]; 3]; nsys + 2];
            self.vdel = vec![[[0.0; 2]; 3]; nsys + 2];
        }
        // two extra columns (nsys+1, nsys+2) carry VDEL in the blocked solver
        let ldm = (nsys + 3 + 7) & !7;
        if self.ldm != ldm {
            self.ldm = ldm;
            self.vm = vec![0.0; (nsys + 2) * 3 * ldm];
        }
    }

    /// Forced transition arc length location on side `is` (XIFSET).
    pub fn xifset(&mut self, is: usize) {
        if self.xstrip[is] >= 1.0 {
            self.bl.xiforc = self.xssi[is][self.iblte[is]];
            return;
        }
        let n = self.n;
        let chx = self.xte - self.xle;
        let chy = self.yte - self.yle;
        let chsq = chx * chx + chy * chy;
        let mut w1 = vec![0.0; n + 1];
        let mut w2 = vec![0.0; n + 1];
        let mut w3 = vec![0.0; n + 1];
        let mut w4 = vec![0.0; n + 1];
        for i in 1..=n {
            w1[i] = ((self.x[i] - self.xle) * chx + (self.y[i] - self.yle) * chy) / chsq;
            w2[i] = ((self.y[i] - self.yle) * chx - (self.x[i] - self.xle) * chy) / chsq;
        }
        splind(&w1, &mut w3, &self.s, n, -999.0, -999.0);
        splind(&w2, &mut w4, &self.s, n, -999.0, -999.0);
        if is == 1 {
            let mut str_ = self.sle + (self.s[1] - self.sle) * self.xstrip[is];
            sinvrt(&mut str_, self.xstrip[is], &w1, &w3, &self.s, n);
            self.bl.xiforc = (self.sst - str_).min(self.xssi[is][self.iblte[is]]);
        } else {
            let mut str_ = self.sle + (self.s[n] - self.sle) * self.xstrip[is];
            sinvrt(&mut str_, self.xstrip[is], &w1, &w3, &self.s, n);
            self.bl.xiforc = (str_ - self.sst).min(self.xssi[is][self.iblte[is]]);
        }
        if self.bl.xiforc < 0.0 {
            if self.verbose {
                crate::outln!("\n ***  Stagnation point is past trip on side{:2}  ***", is);
            }
            self.bl.xiforc = self.xssi[is][self.iblte[is]];
        }
    }

    /// Store station results into BL arrays (shared tail of MRCHUE/MRCHDU).
    fn store_station(&mut self, is: usize, ibl: usize, ami: f64, cti: f64, thi: f64, dsi: f64, uei: f64) {
        if ibl < self.itran[is] {
            self.ctau[is][ibl] = ami;
        }
        if ibl >= self.itran[is] {
            self.ctau[is][ibl] = cti;
        }
        self.thet[is][ibl] = thi;
        self.dstr[is][ibl] = dsi;
        self.uedg[is][ibl] = uei;
        self.mass[is][ibl] = dsi * uei;
        let c2 = &self.bl.c2;
        self.tau[is][ibl] = 0.5 * c2.r * c2.u * c2.u * c2.cf;
        self.dis[is][ibl] = c2.r * c2.u * c2.u * c2.u * c2.di * c2.hs * 0.5;
        self.ctq[is][ibl] = c2.cq;
        self.delt[is][ibl] = c2.de;
        self.tstr[is][ibl] = c2.hs * c2.t;
    }

    fn te_quantities(&self) -> (f64, f64, f64) {
        let (i1, i2) = (self.iblte[1], self.iblte[2]);
        let tte = self.thet[1][i1] + self.thet[2][i2];
        let dte = self.dstr[1][i1] + self.dstr[2][i2] + self.ante;
        let cte = (self.ctau[1][i1] * self.thet[1][i1] + self.ctau[2][i2] * self.thet[2][i2]) / tte;
        (cte, tte, dte)
    }

    /// Fallback extrapolation & secondary vars when station Newton fails badly.
    #[allow(clippy::too_many_arguments)]
    fn station_recover(
        &mut self,
        is: usize,
        ibl: usize,
        dmax: f64,
        xsi: f64,
        ami: &mut f64,
        cti: &mut f64,
        thi: &mut f64,
        dsi: &mut f64,
        uei: &mut f64,
        dswaki: f64,
        cte: f64,
        tte: f64,
        dte: f64,
        mrchdu: bool,
    ) {
        let ibm = ibl - 1;
        if dmax > 0.1 && ibl > 3 {
            if ibl <= self.iblte[is] {
                *thi = self.thet[is][ibm] * (self.xssi[is][ibl] / self.xssi[is][ibm]).sqrt();
                *dsi = self.dstr[is][ibm] * (self.xssi[is][ibl] / self.xssi[is][ibm]).sqrt();
                if mrchdu {
                    *uei = self.uedg[is][ibm];
                }
            } else if ibl == self.iblte[is] + 1 {
                *cti = cte;
                *thi = tte;
                *dsi = dte;
                if mrchdu {
                    *uei = self.uedg[is][ibm];
                }
            } else {
                *thi = self.thet[is][ibm];
                let ratlen = (self.xssi[is][ibl] - self.xssi[is][ibm]) / (10.0 * self.dstr[is][ibm]);
                *dsi = (self.dstr[is][ibm] + *thi * ratlen) / (1.0 + ratlen);
                if mrchdu {
                    *uei = self.uedg[is][ibm];
                }
            }
            if ibl == self.itran[is] {
                *cti = 0.05;
            }
            if ibl > self.itran[is] {
                *cti = self.ctau[is][ibm];
            }
            if !mrchdu {
                *uei = self.uedg[is][ibl];
                if ibl > 2 && ibl < self.nbl[is] {
                    *uei = 0.5 * (self.uedg[is][ibl - 1] + self.uedg[is][ibl + 1]);
                }
            }
        }
        // label 109
        self.bl.blprv(xsi, *ami, *cti, *thi, *dsi, dswaki, *uei);
        self.bl.blkin();
        if !self.bl.simi && !self.bl.turb {
            self.bl.trchek();
            *ami = self.bl.c2.ampl;
            if self.bl.tran {
                self.itran[is] = ibl;
            } else {
                self.itran[is] = ibl + 2;
            }
        }
        if ibl < self.itran[is] {
            self.bl.blvar(1);
        }
        if ibl >= self.itran[is] {
            self.bl.blvar(2);
        }
        if self.bl.wake {
            self.bl.blvar(3);
        }
        if ibl < self.itran[is] {
            self.bl.blmid(1);
        }
        if ibl >= self.itran[is] {
            self.bl.blmid(2);
        }
        if self.bl.wake {
            self.bl.blmid(3);
        }
    }

    /// Direct-mode BL march with Ue, inverse at separation (MRCHUE).
    pub fn mrchue(&mut self) {
        let hlmax = 3.8;
        let htmax = 2.5;
        for is in 1..=2 {
            if self.verbose {
                crate::outln!("    side  {} ...", is);
            }
            self.xifset(is);
            let ibl = 2;
            let mut xsi = self.xssi[is][ibl];
            let mut uei = self.uedg[is][ibl];
            let bule = 1.0;
            let ucon = uei / xsi.powf(bule);
            let tsq = 0.45 / (ucon * (5.0 * bule + 1.0) * self.bl.reybl) * xsi.powf(1.0 - bule);
            let mut thi = tsq.sqrt();
            let mut dsi = 2.2 * thi;
            let mut ami = 0.0;
            let mut cti = 0.03;
            self.bl.tran = false;
            self.bl.turb = false;
            self.itran[is] = self.iblte[is];

            let (mut cte, mut tte, mut dte) = (0.0, 0.0, 0.0);
            for ibl in 2..=self.nbl[is] {
                self.bl.simi = ibl == 2;
                self.bl.wake = ibl > self.iblte[is];
                xsi = self.xssi[is][ibl];
                uei = self.uedg[is][ibl];
                let dswaki = if self.bl.wake { self.wgap[ibl - self.iblte[is]] } else { 0.0 };
                let mut direct = true;
                let mut dmax = 0.0;
                let mut htarg = 0.0;
                let mut hmax = 0.0;
                let mut converged = false;
                for _itbl in 1..=25 {
                    self.bl.blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
                    self.bl.blkin();
                    if !self.bl.simi && !self.bl.turb {
                        self.bl.trchek();
                        ami = self.bl.c2.ampl;
                        if self.bl.tran {
                            self.itran[is] = ibl;
                            if cti <= 0.0 {
                                cti = 0.03;
                                self.bl.c2.s = cti;
                            }
                        } else {
                            self.itran[is] = ibl + 2;
                        }
                    }
                    if ibl == self.iblte[is] + 1 {
                        (cte, tte, dte) = self.te_quantities();
                        self.bl.tesys(cte, tte, dte);
                    } else {
                        self.bl.blsys();
                    }

                    if direct {
                        self.bl.vs2[3] = [0.0, 0.0, 0.0, 1.0, self.bl.vs2[3][4]];
                        self.bl.vsrez[3] = 0.0;
                        let mut z = self.bl.vs2;
                        let mut r = self.bl.vsrez;
                        gauss4(&mut z, &mut r);
                        self.bl.vsrez = r;
                        let vr = r;
                        dmax = (vr[1] / thi).abs().max((vr[2] / dsi).abs());
                        if ibl < self.itran[is] {
                            dmax = dmax.max((vr[0] / 10.0).abs());
                        }
                        if ibl >= self.itran[is] {
                            dmax = dmax.max((vr[0] / cti).abs());
                        }
                        let mut rlx = 1.0;
                        if dmax > 0.3 {
                            rlx = 0.3 / dmax;
                        }
                        if ibl != self.iblte[is] + 1 {
                            let msq = uei * uei * self.bl.hstinv / (self.bl.gm1bl * (1.0 - 0.5 * uei * uei * self.bl.hstinv));
                            let htest = (dsi + rlx * vr[2]) / (thi + rlx * vr[1]);
                            let (hktest, _, _) = hkin(htest, msq);
                            if ibl < self.itran[is] {
                                hmax = hlmax;
                            }
                            if ibl >= self.itran[is] {
                                hmax = htmax;
                            }
                            direct = hktest < hmax;
                        }
                        if direct {
                            if ibl >= self.itran[is] {
                                cti += rlx * vr[0];
                            }
                            thi += rlx * vr[1];
                            dsi += rlx * vr[2];
                        } else {
                            let c1 = &self.bl.c1;
                            let (x1, x2, t1, hk1) = (c1.x, self.bl.c2.x, c1.t, c1.hk);
                            if ibl < self.itran[is] {
                                htarg = hk1 + 0.03 * (x2 - x1) / t1;
                            } else if ibl == self.itran[is] {
                                htarg = hk1 + (0.03 * (self.bl.xt - x1) - 0.15 * (x2 - self.bl.xt)) / t1;
                            } else if self.bl.wake {
                                let cnst = 0.03 * (x2 - x1) / t1;
                                let mut hk2 = hk1;
                                for _ in 0..3 {
                                    hk2 -= (hk2 + cnst * (hk2 - 1.0).powi(3) - hk1) / (1.0 + 3.0 * cnst * (hk2 - 1.0).powi(2));
                                }
                                htarg = hk2;
                            } else {
                                htarg = hk1 - 0.15 * (x2 - x1) / t1;
                            }
                            if self.bl.wake {
                                htarg = htarg.max(1.01);
                            } else {
                                htarg = htarg.max(hmax);
                            }
                            if self.verbose {
                                crate::outln!(" MRCHUE: Inverse mode at{:4}     Hk ={:8.3}", ibl, htarg);
                            }
                            continue;
                        }
                    } else {
                        let c2 = self.bl.c2;
                        self.bl.vs2[3] = [0.0, c2.hk_t, c2.hk_d, c2.hk_u, self.bl.vs2[3][4]];
                        self.bl.vsrez[3] = htarg - c2.hk;
                        let mut z = self.bl.vs2;
                        let mut r = self.bl.vsrez;
                        gauss4(&mut z, &mut r);
                        self.bl.vsrez = r;
                        dmax = (r[1] / thi).abs().max((r[2] / dsi).abs()).max((r[3] / uei).abs());
                        if ibl >= self.itran[is] {
                            dmax = dmax.max((r[0] / cti).abs());
                        }
                        let mut rlx = 1.0;
                        if dmax > 0.3 {
                            rlx = 0.3 / dmax;
                        }
                        if ibl >= self.itran[is] {
                            cti += rlx * r[0];
                        }
                        thi += rlx * r[1];
                        dsi += rlx * r[2];
                        uei += rlx * r[3];
                    }
                    if ibl >= self.itran[is] {
                        cti = cti.min(0.30).max(0.0000001);
                    }
                    let hklim = if ibl <= self.iblte[is] { 1.02 } else { 1.00005 };
                    let msq = uei * uei * self.bl.hstinv / (self.bl.gm1bl * (1.0 - 0.5 * uei * uei * self.bl.hstinv));
                    let mut dsw = dsi - dswaki;
                    dslim(&mut dsw, thi, msq, hklim);
                    dsi = dsw + dswaki;
                    if dmax <= 1.0e-5 {
                        converged = true;
                        break;
                    }
                }
                if !converged {
                    if self.verbose {
                        crate::outln!(" MRCHUE: Convergence failed at{:4}  side{:2}    Res ={}", ibl, is, crate::fmt::e_fmt(dmax, 12, 4));
                    }
                    self.station_recover(
                        is, ibl, dmax, xsi, &mut ami, &mut cti, &mut thi, &mut dsi, &mut uei, dswaki, cte, tte, dte, false,
                    );
                }
                // label 110
                self.store_station(is, ibl, ami, cti, thi, dsi, uei);
                self.bl.blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
                self.bl.blkin();
                self.bl.c1 = self.bl.c2;
                if self.bl.tran || ibl == self.iblte[is] {
                    self.bl.turb = true;
                    self.tforce[is] = self.bl.trforc;
                    self.xssitr[is] = self.bl.xt;
                }
                self.bl.tran = false;
                if ibl == self.iblte[is] {
                    thi = self.thet[1][self.iblte[1]] + self.thet[2][self.iblte[2]];
                    dsi = self.dstr[1][self.iblte[1]] + self.dstr[2][self.iblte[2]] + self.ante;
                }
            }
        }
    }

    /// Mixed-mode BL march with Ue-Hk relation (MRCHDU).
    pub fn mrchdu(&mut self) {
        const DEPS: f64 = 5.0e-6;
        let senswt = 1000.0;
        let mut sens = 0.0;
        let mut sennew = 0.0;
        for is in 1..=2 {
            self.xifset(is);
            self.bl.bule = 1.0;
            let itrold = self.itran[is];
            self.bl.tran = false;
            self.bl.turb = false;
            self.itran[is] = self.iblte[is];
            let (mut cte, mut tte, mut dte) = (0.0, 0.0, 0.0);
            let mut ami = 0.0;
            for ibl in 2..=self.nbl[is] {
                self.bl.simi = ibl == 2;
                self.bl.wake = ibl > self.iblte[is];
                let xsi = self.xssi[is][ibl];
                let mut uei = self.uedg[is][ibl];
                let mut thi = self.thet[is][ibl];
                let mut dsi = self.dstr[is][ibl];
                let mut cti;
                if ibl < itrold {
                    ami = self.ctau[is][ibl];
                    cti = 0.03;
                } else {
                    cti = self.ctau[is][ibl];
                    if cti <= 0.0 {
                        cti = 0.03;
                    }
                }
                let dswaki = if self.bl.wake { self.wgap[ibl - self.iblte[is]] } else { 0.0 };
                if ibl <= self.iblte[is] {
                    dsi = (dsi - dswaki).max(1.02000 * thi) + dswaki;
                }
                if ibl > self.iblte[is] {
                    dsi = (dsi - dswaki).max(1.00005 * thi) + dswaki;
                }
                let mut ueref = 0.0;
                let mut hkref = 0.0;
                let mut dmax = 0.0;
                let mut converged = false;
                for itbl in 1..=25 {
                    self.bl.blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
                    self.bl.blkin();
                    if !self.bl.simi && !self.bl.turb {
                        self.bl.trchek();
                        ami = self.bl.c2.ampl;
                        if self.bl.tran {
                            self.itran[is] = ibl;
                        } else {
                            self.itran[is] = ibl + 2;
                        }
                    }
                    if ibl == self.iblte[is] + 1 {
                        (cte, tte, dte) = self.te_quantities();
                        self.bl.tesys(cte, tte, dte);
                    } else {
                        self.bl.blsys();
                    }
                    if itbl == 1 {
                        ueref = self.bl.c2.u;
                        hkref = self.bl.c2.hk;
                        if ibl < self.itran[is] && ibl >= itrold {
                            let uem = self.uedg[is][ibl - 1];
                            let dsm = self.dstr[is][ibl - 1];
                            let thm = self.thet[is][ibl - 1];
                            let msq = uem * uem * self.bl.hstinv / (self.bl.gm1bl * (1.0 - 0.5 * uem * uem * self.bl.hstinv));
                            hkref = hkin(dsm / thm, msq).0;
                        }
                        if ibl < itrold {
                            if self.bl.tran {
                                self.ctau[is][ibl] = 0.03;
                            }
                            if self.bl.turb {
                                self.ctau[is][ibl] = self.ctau[is][ibl - 1];
                            }
                            if self.bl.tran || self.bl.turb {
                                cti = self.ctau[is][ibl];
                                self.bl.c2.s = cti;
                            }
                        }
                    }
                    let c2 = self.bl.c2;
                    if self.bl.simi || ibl == self.iblte[is] + 1 {
                        self.bl.vs2[3] = [0.0, 0.0, 0.0, c2.u_uei, self.bl.vs2[3][4]];
                        self.bl.vsrez[3] = ueref - c2.u;
                    } else {
                        let mut vtmp = self.bl.vs2;
                        let mut vztmp = self.bl.vsrez;
                        vtmp[3] = [0.0, c2.hk_t, c2.hk_d, c2.hk_u * c2.u_uei, vtmp[3][4]];
                        vztmp[3] = 1.0;
                        gauss4(&mut vtmp, &mut vztmp);
                        sennew = senswt * vztmp[3] * hkref / ueref;
                        if itbl <= 5 {
                            sens = sennew;
                        } else if itbl <= 15 {
                            sens = 0.5 * (sens + sennew);
                        }
                        self.bl.vs2[3] = [
                            0.0,
                            c2.hk_t * hkref,
                            c2.hk_d * hkref,
                            (c2.hk_u * hkref + sens / ueref) * c2.u_uei,
                            self.bl.vs2[3][4],
                        ];
                        self.bl.vsrez[3] = -(hkref * hkref) * (c2.hk / hkref - 1.0) - sens * (c2.u / ueref - 1.0);
                    }
                    let mut z = self.bl.vs2;
                    let mut r = self.bl.vsrez;
                    gauss4(&mut z, &mut r);
                    self.bl.vsrez = r;
                    dmax = (r[1] / thi).abs().max((r[2] / dsi).abs()).max((r[3] / uei).abs());
                    if ibl >= self.itran[is] {
                        dmax = dmax.max((r[0] / (10.0 * cti)).abs());
                    }
                    let mut rlx = 1.0;
                    if dmax > 0.3 {
                        rlx = 0.3 / dmax;
                    }
                    if ibl < self.itran[is] {
                        ami += rlx * r[0];
                    }
                    if ibl >= self.itran[is] {
                        cti += rlx * r[0];
                    }
                    thi += rlx * r[1];
                    dsi += rlx * r[2];
                    uei += rlx * r[3];
                    if ibl >= self.itran[is] {
                        cti = cti.min(0.30).max(0.0000001);
                    }
                    let hklim = if ibl <= self.iblte[is] { 1.02 } else { 1.00005 };
                    let msq = uei * uei * self.bl.hstinv / (self.bl.gm1bl * (1.0 - 0.5 * uei * uei * self.bl.hstinv));
                    let mut dsw = dsi - dswaki;
                    dslim(&mut dsw, thi, msq, hklim);
                    dsi = dsw + dswaki;
                    if dmax <= DEPS {
                        converged = true;
                        break;
                    }
                }
                if !converged {
                    if self.verbose {
                        crate::outln!(" MRCHDU: Convergence failed at{:4}  side{:2}    Res ={}", ibl, is, crate::fmt::e_fmt(dmax, 12, 4));
                    }
                    self.station_recover(
                        is, ibl, dmax, xsi, &mut ami, &mut cti, &mut thi, &mut dsi, &mut uei, dswaki, cte, tte, dte, true,
                    );
                }
                sens = sennew;
                self.store_station(is, ibl, ami, cti, thi, dsi, uei);
                self.bl.blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
                self.bl.blkin();
                self.bl.c1 = self.bl.c2;
                if self.bl.tran || ibl == self.iblte[is] {
                    self.bl.turb = true;
                    self.tforce[is] = self.bl.trforc;
                    self.xssitr[is] = self.bl.xt;
                }
                self.bl.tran = false;
            }
        }
    }

    /// Build the mass-influence matrix M[iv][jv] = dUe(iv)/dm(jv)
    /// (the UE_M rows of SETBL/UESET/UPDATE), rows in parallel.
    pub fn build_uem(&mut self) {
        let nsys = self.nsys;
        let ld = self.ldm;
        if self.uem.len() < (nsys + 1) * ld {
            self.uem = vec![0.0; (nsys + 1) * ld];
        }
        let mut jpan = vec![0usize; nsys + 1];
        let mut jvti = vec![0.0f64; nsys + 1];
        for is in 1..=2 {
            for ibl in 2..=self.nbl[is] {
                let iv = self.isys[is][ibl];
                jpan[iv] = self.ipan[is][ibl];
                jvti[iv] = self.vti[is][ibl];
            }
        }
        let dij = &self.dij;
        let fill = |iv: usize, row: &mut [f64]| {
            let drow = dij.row(jpan[iv]);
            let vt = jvti[iv];
            for jv in 1..=nsys {
                row[jv] = -vt * jvti[jv] * drow[jpan[jv]];
            }
        };
        let m = &mut self.uem[..(nsys + 1) * ld];
        if crate::par_fine() && nsys >= 64 {
            use rayon::prelude::*;
            m.par_chunks_mut(ld).enumerate().skip(1).for_each(|(iv, row)| fill(iv, row));
        } else {
            for (iv, row) in m.chunks_mut(ld).enumerate().skip(1) {
                fill(iv, row);
            }
        }
    }

    /// UESET using the precomputed M matrix (same summation order).
    fn ueset_m(&mut self) {
        let nsys = self.nsys;
        let mut mv = vec![0.0; nsys + 1];
        for js in 1..=2 {
            for jbl in 2..=self.nbl[js] {
                mv[self.isys[js][jbl]] = self.mass[js][jbl];
            }
        }
        let mut dui = vec![0.0; nsys + 1];
        matvec_rows(&self.uem, self.ldm, nsys, &[&mv], &mut [&mut dui]);
        for is in 1..=2 {
            for ibl in 2..=self.nbl[is] {
                let iv = self.isys[is][ibl];
                self.uedg[is][ibl] = self.uinv[is][ibl] + dui[iv];
            }
        }
    }

    /// Set up global Newton system (SETBL).
    pub fn setbl(&mut self) {
        let clmr = if self.lalfa { self.cl } else { self.clspec };
        let mut re_clmr = 0.0;
        let ma_clmr = self.mrcl(clmr, &mut re_clmr);
        let msq_clmr = 2.0 * self.minf * ma_clmr;
        self.comset();

        {
            let b = &mut self.bl;
            b.gambl = self.gamma;
            b.gm1bl = self.gamm1;
            b.qinfbl = self.qinf;
            b.tkbl = self.tklam;
            b.tkbl_ms = self.tkl_msq;
            let gm1 = b.gm1bl;
            let minf = self.minf;
            let f = 1.0 + 0.5 * gm1 * minf * minf;
            b.rstbl = f.powf(1.0 / gm1);
            b.rstbl_ms = 0.5 * b.rstbl / f;
            b.hstinv = gm1 * (minf / b.qinfbl).powi(2) / f;
            b.hstinv_ms = gm1 * (1.0 / b.qinfbl).powi(2) / f - 0.5 * gm1 * b.hstinv / f;
            let herat = 1.0 - 0.5 * b.qinfbl * b.qinfbl * b.hstinv;
            let herat_ms = -0.5 * b.qinfbl * b.qinfbl * b.hstinv_ms;
            let hv = b.hvrat;
            b.reybl = self.reinf * (herat * herat * herat).sqrt() * (1.0 + hv) / (herat + hv);
            b.reybl_re = (herat * herat * herat).sqrt() * (1.0 + hv) / (herat + hv);
            b.reybl_ms = b.reybl * (1.5 / herat - 1.0 / (herat + hv)) * herat_ms;
            b.amcrit = self.acrit;
            b.idampv = self.idamp;
            b.dwte = self.wgap[1];
        }

        if !self.lblini {
            if self.verbose {
                crate::outln!("\n Initializing BL ...");
            }
            self.mrchue();
            self.lblini = true;
        }
        self.mrchdu();

        self.build_uem();

        let nsys = self.nsys;
        let mut usav = [vec![0.0; self.nbl[1] + 1], vec![0.0; self.nbl[1] + 1], vec![0.0; self.nbl[2] + 1]];
        for is in 1..=2 {
            for ibl in 2..=self.nbl[is] {
                usav[is][ibl] = self.uedg[is][ibl];
            }
        }
        self.ueset_m();
        for is in 1..=2 {
            for ibl in 2..=self.nbl[is] {
                std::mem::swap(&mut usav[is][ibl], &mut self.uedg[is][ibl]);
            }
        }

        let jvle1 = self.isys[1][2];
        let jvle2 = self.isys[2][2];
        let jvte1 = self.isys[1][self.iblte[1]];
        let jvte2 = self.isys[2][self.iblte[2]];
        let dule1 = self.uedg[1][2] - usav[1][2];
        let dule2 = self.uedg[2][2] - usav[2][2];
        let ule1_a = self.uinv_a[1][2];
        let ule2_a = self.uinv_a[2][2];

        // per-station coefficients for the parallel VM fill
        let mut sc = std::mem::take(&mut self.stcoef);
        sc.clear();
        sc.resize(nsys + 1, StCoef::default());

        let mut ami = 0.0;
        let mut cti = 0.0;
        for is in 1..=2 {
            let mut u1_a = 0.0;
            let mut d1_a = 0.0;
            let mut due1 = 0.0;
            let mut dds1 = 0.0;
            let mut d1: D1Kind = D1Kind::Zero;
            let mut u1row: usize = 0;
            self.bl.bule = 1.0;
            self.xifset(is);
            self.bl.tran = false;
            self.bl.turb = false;

            let (mut tte_tte1, mut tte_tte2, mut cte_cte1, mut cte_cte2, mut cte_tte1, mut cte_tte2) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            for ibl in 2..=self.nbl[is] {
                let iv = self.isys[is][ibl];
                self.bl.simi = ibl == 2;
                self.bl.wake = ibl > self.iblte[is];
                self.bl.tran = ibl == self.itran[is];
                self.bl.turb = ibl > self.itran[is];

                let xsi = self.xssi[is][ibl];
                if ibl < self.itran[is] {
                    ami = self.ctau[is][ibl];
                }
                if ibl >= self.itran[is] {
                    cti = self.ctau[is][ibl];
                }
                let uei = self.uedg[is][ibl];
                let thi = self.thet[is][ibl];
                let mdi = self.mass[is][ibl];
                let dsi = mdi / uei;
                let dswaki = if self.bl.wake { self.wgap[ibl - self.iblte[is]] } else { 0.0 };

                let d2_m2 = 1.0 / uei;
                let d2_u2 = -dsi / uei;
                let u2_a = self.uinv_a[is][ibl];
                let d2_a = d2_u2 * u2_a;
                let due2 = self.uedg[is][ibl] - usav[is][ibl];
                let dds2 = d2_u2 * due2;

                self.bl.blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
                self.bl.blkin();
                if self.bl.tran {
                    self.bl.trchek();
                    ami = self.bl.c2.ampl;
                }
                if ibl == self.itran[is] && !self.bl.tran && self.verbose {
                    crate::outln!(" SETBL: Xtr???  n1 n2: {} {}", self.bl.c1.ampl, self.bl.c2.ampl);
                }

                if ibl == self.iblte[is] + 1 {
                    let (i1, i2) = (self.iblte[1], self.iblte[2]);
                    let (cte, tte, dte) = self.te_quantities();
                    self.bl.tesys(cte, tte, dte);
                    tte_tte1 = 1.0;
                    tte_tte2 = 1.0;
                    let dte_mte1 = 1.0 / self.uedg[1][i1];
                    let dte_ute1 = -self.dstr[1][i1] / self.uedg[1][i1];
                    let dte_mte2 = 1.0 / self.uedg[2][i2];
                    let dte_ute2 = -self.dstr[2][i2] / self.uedg[2][i2];
                    cte_cte1 = self.thet[1][i1] / tte;
                    cte_cte2 = self.thet[2][i2] / tte;
                    cte_tte1 = (self.ctau[1][i1] - cte) / tte;
                    cte_tte2 = (self.ctau[2][i2] - cte) / tte;
                    d1 = D1Kind::Te { ute1: dte_ute1, ute2: dte_ute2, mte1: dte_mte1, mte2: dte_mte2 };
                    due1 = 0.0;
                    dds1 = dte_ute1 * (self.uedg[1][i1] - usav[1][i1]) + dte_ute2 * (self.uedg[2][i2] - usav[2][i2]);
                } else {
                    self.bl.blsys();
                }

                {
                    let c2 = &self.bl.c2;
                    self.tau[is][ibl] = 0.5 * c2.r * c2.u * c2.u * c2.cf;
                    self.dis[is][ibl] = c2.r * c2.u * c2.u * c2.u * c2.di * c2.hs * 0.5;
                    self.ctq[is][ibl] = c2.cq;
                    self.delt[is][ibl] = c2.de;
                    self.uslp[is][ibl] = 1.60 / (1.0 + c2.us);
                }

                let (xi_ule1, xi_ule2) = if is == 1 { (self.sst_go, -self.sst_gp) } else { (-self.sst_go, self.sst_gp) };

                let vs1 = self.bl.vs1;
                let vs2 = self.bl.vs2;
                let vsx = self.bl.vsx;
                let mut co = [[0.0; 5]; 3];
                for k in 0..3 {
                    let e = vs1[k][4] + vs2[k][4] + vsx[k];
                    co[k] = [vs1[k][2], vs1[k][3], vs2[k][2], vs2[k][3], e];
                    self.vb[iv][k][0] = vs1[k][0];
                    self.vb[iv][k][1] = vs1[k][1];
                    self.va[iv][k][0] = vs2[k][0];
                    self.va[iv][k][1] = vs2[k][1];
                    if self.lalfa {
                        self.vdel[iv][k][1] = self.bl.vsr[k] * re_clmr + self.bl.vsm[k] * msq_clmr;
                    } else {
                        self.vdel[iv][k][1] = (vs1[k][3] * u1_a + vs1[k][2] * d1_a)
                            + (vs2[k][3] * u2_a + vs2[k][2] * d2_a)
                            + e * (xi_ule1 * ule1_a + xi_ule2 * ule2_a);
                    }
                    self.vdel[iv][k][0] = self.bl.vsrez[k]
                        + (vs1[k][3] * due1 + vs1[k][2] * dds1)
                        + (vs2[k][3] * due2 + vs2[k][2] * dds2)
                        + e * (xi_ule1 * dule1 + xi_ule2 * dule2);
                }
                sc[iv] = StCoef { co, xi1: xi_ule1, xi2: xi_ule2, d2_u2, d2_m2, d1, u1row };

                if ibl == self.iblte[is] + 1 {
                    for k in 0..3 {
                        self.vz[k][0] = vs1[k][0] * cte_cte1;
                        self.vz[k][1] = vs1[k][0] * cte_tte1 + vs1[k][1] * tte_tte1;
                        self.vb[iv][k][0] = vs1[k][0] * cte_cte2;
                        self.vb[iv][k][1] = vs1[k][0] * cte_tte2 + vs1[k][1] * tte_tte2;
                    }
                }

                if self.bl.tran {
                    self.bl.turb = true;
                    self.itran[is] = ibl;
                    self.tforce[is] = self.bl.trforc;
                    self.xssitr[is] = self.bl.xt;
                    let str_ = if is == 1 { self.sst - self.bl.xt } else { self.sst + self.bl.xt };
                    let chx = self.xte - self.xle;
                    let chy = self.yte - self.yle;
                    let chsq = chx * chx + chy * chy;
                    let n = self.n;
                    let xtr = seval(str_, &self.x, &self.xp, &self.s, n);
                    let ytr = seval(str_, &self.y, &self.yp, &self.s, n);
                    self.xoctr[is] = ((xtr - self.xle) * chx + (ytr - self.yle) * chy) / chsq;
                    self.yoctr[is] = ((ytr - self.yle) * chx - (xtr - self.xle) * chy) / chsq;
                }
                self.bl.tran = false;

                if ibl == self.iblte[is] {
                    self.bl.turb = true;
                    self.bl.wake = true;
                    self.bl.blvar(3);
                    self.bl.blmid(3);
                }

                // "2" -> "1" for the next station
                d1 = D1Kind::Prev { row: iv, d_u: d2_u2, d_m: d2_m2 };
                u1row = iv;
                u1_a = u2_a;
                d1_a = d2_a;
                due1 = due2;
                dds1 = dds2;
                self.bl.c1 = self.bl.c2;
            }
            if self.verbose {
                if self.tforce[is] {
                    crate::outln!(" Side{:2} forced transition at x/c = {:7.4}{:5}", is, self.xoctr[is], self.itran[is]);
                } else {
                    crate::outln!(" Side{:2}  free  transition at x/c = {:7.4}{:5}", is, self.xoctr[is], self.itran[is]);
                }
            }
        }

        // ---- parallel fill of VM rows (bit-identical to the sequential Fortran formulas)
        let ldm = self.ldm;
        let mut vm = std::mem::take(&mut self.vm);
        {
            let uem = &self.uem;
            let row = |iv: usize| &uem[iv * ldm..iv * ldm + nsys + 1];
            let zero = vec![0.0f64; nsys + 1];
            let sc = &sc;
            let fill = |iv: usize, blk: &mut [f64], d1v: &mut Vec<f64>, d2v: &mut Vec<f64>| {
                let c = &sc[iv];
                let u2 = row(iv);
                let u1: &[f64] = if c.u1row > 0 { row(c.u1row) } else { &zero };
                let l1 = row(jvle1);
                let l2 = row(jvle2);
                d1v.clear();
                d1v.resize(nsys + 1, 0.0);
                match c.d1 {
                    D1Kind::Zero => {}
                    D1Kind::Prev { row: pr, d_u, d_m } => {
                        let p = row(pr);
                        for (o, &pv) in d1v.iter_mut().zip(p.iter()) {
                            *o = d_u * pv;
                        }
                        d1v[pr] += d_m;
                    }
                    D1Kind::Te { ute1, ute2, mte1, mte2 } => {
                        let t1 = row(jvte1);
                        let t2 = row(jvte2);
                        for ((o, &a), &b) in d1v.iter_mut().zip(t1.iter()).zip(t2.iter()) {
                            *o = ute1 * a + ute2 * b;
                        }
                        d1v[jvte1] += mte1;
                        d1v[jvte2] += mte2;
                    }
                }
                d2v.clear();
                d2v.resize(nsys + 1, 0.0);
                for (o, &u) in d2v.iter_mut().zip(u2.iter()) {
                    *o = c.d2_u2 * u;
                }
                d2v[iv] += c.d2_m2;
                for k in 0..3 {
                    let [a, b, cc, d, e] = c.co[k];
                    let (x1, x2) = (c.xi1, c.xi2);
                    let r = &mut blk[k * ldm + 1..k * ldm + nsys + 1];
                    let it = r
                        .iter_mut()
                        .zip(&d1v[1..])
                        .zip(&u1[1..])
                        .zip(&d2v[1..])
                        .zip(&u2[1..])
                        .zip(&l1[1..])
                        .zip(&l2[1..]);
                    for ((((((o, &dd1), &uu1), &dd2), &uu2), &ll1), &ll2) in it {
                        *o = a * dd1 + b * uu1 + cc * dd2 + d * uu2 + e * (x1 * ll1 + x2 * ll2);
                    }
                }
            };
            let blocks = &mut vm[..(nsys + 1) * 3 * ldm];
            if crate::par_fine() && nsys >= 64 {
                use rayon::prelude::*;
                blocks
                    .par_chunks_mut(3 * ldm)
                    .enumerate()
                    .skip(1)
                    .for_each_init(|| (Vec::new(), Vec::new()), |(a, b), (iv, blk)| fill(iv, blk, a, b));
            } else {
                let (mut a, mut b) = (Vec::new(), Vec::new());
                for (iv, blk) in blocks.chunks_mut(3 * ldm).enumerate().skip(1) {
                    fill(iv, blk, &mut a, &mut b);
                }
            }
        }
        self.vm = vm;
        self.stcoef = sc;
    }

    /// Custom block solver for the coupled Newton system (BLSOLV).
    /// Uses the cache-blocked, bit-exact reformulation in `blsolv.rs`.
    pub fn blsolv(&mut self) {
        let nsys = self.nsys;
        let ivte1 = self.isys[1][self.iblte[1]];
        let ivz = self.isys[2][self.iblte[2] + 1];
        let n = self.n;
        let vacc = [self.vaccel, self.vaccel * 2.0 / (self.s[n] - self.s[1]), self.vaccel * 2.0 / (self.s[n] - self.s[1])];
        crate::blsolv::blsolv_blocked(&mut self.vm, self.ldm, &mut self.va, &self.vb, &self.vz, &mut self.vdel, nsys, ivte1, ivz, vacc);
    }

    /// Reference (unblocked) BLSOLV, a direct transcription of the Fortran.
    pub fn blsolv_ref(&mut self) {
        let nsys = self.nsys;
        let ivte1 = self.isys[1][self.iblte[1]];
        let n = self.n;
        let vacc1 = self.vaccel;
        let vacc2 = self.vaccel * 2.0 / (self.s[n] - self.s[1]);
        let vacc3 = self.vaccel * 2.0 / (self.s[n] - self.s[1]);
        let ldm = self.ldm;
        let vm = &mut self.vm;
        let va = &mut self.va;
        let vb = &self.vb;
        let vdel = &mut self.vdel;
        let idx = |k: usize, l: usize, iv: usize| (iv * 3 + k) * ldm + l;
        let ivz = self.isys[2][self.iblte[2] + 1];

        for iv in 1..=nsys {
            let ivp = iv + 1;
            // normalize first row
            let pivot = 1.0 / va[iv][0][0];
            va[iv][0][1] *= pivot;
            for l in iv..=nsys {
                vm[idx(0, l, iv)] *= pivot;
            }
            vdel[iv][0][0] *= pivot;
            vdel[iv][0][1] *= pivot;
            // eliminate lower first column in VA block
            for k in 1..3 {
                let vtmp = va[iv][k][0];
                va[iv][k][1] -= vtmp * va[iv][0][1];
                let (r0, rk) = (idx(0, 0, iv), idx(k, 0, iv));
                for l in iv..=nsys {
                    vm[rk + l] -= vtmp * vm[r0 + l];
                }
                vdel[iv][k][0] -= vtmp * vdel[iv][0][0];
                vdel[iv][k][1] -= vtmp * vdel[iv][0][1];
            }
            // normalize second row
            let pivot = 1.0 / va[iv][1][1];
            for l in iv..=nsys {
                vm[idx(1, l, iv)] *= pivot;
            }
            vdel[iv][1][0] *= pivot;
            vdel[iv][1][1] *= pivot;
            // eliminate lower second column
            let vtmp = va[iv][2][1];
            {
                let (r1, r2) = (idx(1, 0, iv), idx(2, 0, iv));
                for l in iv..=nsys {
                    vm[r2 + l] -= vtmp * vm[r1 + l];
                }
            }
            vdel[iv][2][0] -= vtmp * vdel[iv][1][0];
            vdel[iv][2][1] -= vtmp * vdel[iv][1][1];
            // normalize third row
            let pivot = 1.0 / vm[idx(2, iv, iv)];
            for l in ivp..=nsys {
                vm[idx(2, l, iv)] *= pivot;
            }
            vdel[iv][2][0] *= pivot;
            vdel[iv][2][1] *= pivot;
            // eliminate upper third column
            let vtmp1 = vm[idx(0, iv, iv)];
            let vtmp2 = vm[idx(1, iv, iv)];
            {
                let (r0, r1, r2) = (idx(0, 0, iv), idx(1, 0, iv), idx(2, 0, iv));
                for l in ivp..=nsys {
                    let v3 = vm[r2 + l];
                    vm[r0 + l] -= vtmp1 * v3;
                    vm[r1 + l] -= vtmp2 * v3;
                }
            }
            vdel[iv][0][0] -= vtmp1 * vdel[iv][2][0];
            vdel[iv][1][0] -= vtmp2 * vdel[iv][2][0];
            vdel[iv][0][1] -= vtmp1 * vdel[iv][2][1];
            vdel[iv][1][1] -= vtmp2 * vdel[iv][2][1];
            // eliminate upper second column
            let vtmp = va[iv][0][1];
            {
                let (r0, r1) = (idx(0, 0, iv), idx(1, 0, iv));
                for l in ivp..=nsys {
                    vm[r0 + l] -= vtmp * vm[r1 + l];
                }
            }
            vdel[iv][0][0] -= vtmp * vdel[iv][1][0];
            vdel[iv][0][1] -= vtmp * vdel[iv][1][1];

            if iv == nsys {
                continue;
            }

            // eliminate VB(IV+1) block
            {
                let (head, tail) = vm.split_at_mut(idx(0, 0, ivp));
                let p0 = &head[idx(0, 0, iv)..idx(0, 0, iv) + ldm];
                let p1 = &head[idx(1, 0, iv)..idx(1, 0, iv) + ldm];
                let p2 = &head[idx(2, 0, iv)..idx(2, 0, iv) + ldm];
                for k in 0..3 {
                    let vtmp1 = vb[ivp][k][0];
                    let vtmp2 = vb[ivp][k][1];
                    let vtmp3 = tail[k * ldm + iv];
                    let row = &mut tail[k * ldm..(k + 1) * ldm];
                    for l in ivp..=nsys {
                        row[l] -= vtmp1 * p0[l] + vtmp2 * p1[l] + vtmp3 * p2[l];
                    }
                    let d0 = vdel[iv];
                    vdel[ivp][k][0] -= vtmp1 * d0[0][0] + vtmp2 * d0[1][0] + vtmp3 * d0[2][0];
                    vdel[ivp][k][1] -= vtmp1 * d0[0][1] + vtmp2 * d0[1][1] + vtmp3 * d0[2][1];
                }
            }

            if iv == ivte1 {
                // eliminate VZ block
                let (head, tail) = vm.split_at_mut(idx(0, 0, ivz));
                let p0 = &head[idx(0, 0, iv)..idx(0, 0, iv) + ldm];
                let p1 = &head[idx(1, 0, iv)..idx(1, 0, iv) + ldm];
                for k in 0..3 {
                    let vtmp1 = self.vz[k][0];
                    let vtmp2 = self.vz[k][1];
                    let row = &mut tail[k * ldm..(k + 1) * ldm];
                    for l in ivp..=nsys {
                        row[l] -= vtmp1 * p0[l] + vtmp2 * p1[l];
                    }
                    let d0 = vdel[iv];
                    vdel[ivz][k][0] -= vtmp1 * d0[0][0] + vtmp2 * d0[1][0];
                    vdel[ivz][k][1] -= vtmp1 * d0[0][1] + vtmp2 * d0[1][1];
                }
            }

            if ivp == nsys {
                continue;
            }

            // eliminate lower VM column
            let (head, tail) = vm.split_at_mut(idx(0, 0, iv + 2));
            let p2 = &head[idx(2, 0, iv)..idx(2, 0, iv) + ldm];
            let d2 = vdel[iv][2];
            let blk = 3 * ldm;
            for (kk, rows) in tail[..(nsys - iv - 1) * blk].chunks_exact_mut(blk).enumerate() {
                let kv = iv + 2 + kk;
                let vacc = [vacc1, vacc2, vacc3];
                for k in 0..3 {
                    let row = &mut rows[k * ldm..(k + 1) * ldm];
                    let vtmp = row[iv];
                    if vtmp.abs() > vacc[k] {
                        for l in ivp..=nsys {
                            row[l] -= vtmp * p2[l];
                        }
                        vdel[kv][k][0] -= vtmp * d2[0];
                        vdel[kv][k][1] -= vtmp * d2[1];
                    }
                }
            }
        }

        // back substitution
        for iv in (2..=nsys).rev() {
            let vtmp = vdel[iv][2][0];
            for kv in (1..iv).rev() {
                vdel[kv][0][0] -= vm[idx(0, iv, kv)] * vtmp;
                vdel[kv][1][0] -= vm[idx(1, iv, kv)] * vtmp;
                vdel[kv][2][0] -= vm[idx(2, iv, kv)] * vtmp;
            }
            let vtmp = vdel[iv][2][1];
            for kv in (1..iv).rev() {
                vdel[kv][0][1] -= vm[idx(0, iv, kv)] * vtmp;
                vdel[kv][1][1] -= vm[idx(1, iv, kv)] * vtmp;
                vdel[kv][2][1] -= vm[idx(2, iv, kv)] * vtmp;
            }
        }
    }

    /// Apply Newton deltas with under-relaxation (UPDATE).
    pub fn update(&mut self) {
        let dalmax = 0.5 * DTOR;
        let dalmin = -0.5 * DTOR;
        let dclmax = 0.5;
        let mut dclmin = -0.5;
        if self.matyp != 1 {
            dclmin = (-0.5f64).max(-0.9 * self.cl);
        }
        let hstinv = self.gamm1 * (self.minf / self.qinf).powi(2) / (1.0 + 0.5 * self.gamm1 * self.minf * self.minf);

        let nsys = self.nsys;
        // flatten mass+vdel per jv to speed the inner loop
        let mut mpd = vec![0.0; nsys + 1];
        let mut vd2 = vec![0.0; nsys + 1];
        for js in 1..=2 {
            for jbl in 2..=self.nbl[js] {
                let jv = self.isys[js][jbl];
                mpd[jv] = self.mass[js][jbl] + self.vdel[jv][2][0];
                vd2[jv] = -self.vdel[jv][2][1];
            }
        }
        let mut unew = [vec![0.0; self.nbl[1] + 1], vec![0.0; self.nbl[1] + 1], vec![0.0; self.nbl[2] + 1]];
        let mut u_ac = [vec![0.0; self.nbl[1] + 1], vec![0.0; self.nbl[1] + 1], vec![0.0; self.nbl[2] + 1]];
        let mut duiv = vec![0.0; nsys + 1];
        let mut duacv = vec![0.0; nsys + 1];
        matvec_rows(&self.uem, self.ldm, nsys, &[&mpd, &vd2], &mut [&mut duiv, &mut duacv]);
        for is in 1..=2 {
            for ibl in 2..=self.nbl[is] {
                let iv = self.isys[is][ibl];
                let dui = duiv[iv];
                let dui_ac = duacv[iv];
                let uinv_ac = if self.lalfa { 0.0 } else { self.uinv_a[is][ibl] };
                unew[is][ibl] = self.uinv[is][ibl] + dui;
                u_ac[is][ibl] = uinv_ac + dui_ac;
            }
        }
        let n = self.n;
        let mut qnew = vec![0.0; n + 1];
        let mut q_ac = vec![0.0; n + 1];
        for is in 1..=2 {
            for ibl in 2..=self.iblte[is] {
                let i = self.ipan[is][ibl];
                qnew[i] = self.vti[is][ibl] * unew[is][ibl];
                q_ac[i] = self.vti[is][ibl] * u_ac[is][ibl];
            }
        }
        let sa = self.alfa.sin();
        let ca = self.alfa.cos();
        let minf = self.minf;
        let qinf = self.qinf;
        let beta = (1.0 - minf * minf).sqrt();
        let beta_msq = -0.5 / beta;
        let bfac = 0.5 * minf * minf / (1.0 + beta);
        let bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;
        let mut clnew = 0.0;
        let mut cl_a = 0.0;
        let mut cl_ms = 0.0;
        let mut cl_ac = 0.0;
        let i = 1;
        let cginc = 1.0 - (qnew[i] / qinf).powi(2);
        let mut cpg1 = cginc / (beta + bfac * cginc);
        let mut cpg1_ms = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);
        let cpi_q = -2.0 * qnew[i] / (qinf * qinf);
        let cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
        let mut cpg1_ac = cpc_cpi * cpi_q * q_ac[i];
        for i in 1..=n {
            let ip = if i == n { 1 } else { i + 1 };
            let cginc = 1.0 - (qnew[ip] / qinf).powi(2);
            let cpg2 = cginc / (beta + bfac * cginc);
            let cpg2_ms = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);
            let cpi_q = -2.0 * qnew[ip] / (qinf * qinf);
            let cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
            let cpg2_ac = cpc_cpi * cpi_q * q_ac[ip];
            let dx = (self.x[ip] - self.x[i]) * ca + (self.y[ip] - self.y[i]) * sa;
            let dx_a = -(self.x[ip] - self.x[i]) * sa + (self.y[ip] - self.y[i]) * ca;
            let ag = 0.5 * (cpg2 + cpg1);
            let ag_ms = 0.5 * (cpg2_ms + cpg1_ms);
            let ag_ac = 0.5 * (cpg2_ac + cpg1_ac);
            clnew += dx * ag;
            cl_a += dx_a * ag;
            cl_ms += dx * ag_ms;
            cl_ac += dx * ag_ac;
            cpg1 = cpg2;
            cpg1_ms = cpg2_ms;
            cpg1_ac = cpg2_ac;
        }

        let mut rlx = 1.0;
        let dac;
        if self.lalfa {
            dac = (clnew - self.cl) / (1.0 - cl_ac - cl_ms * 2.0 * minf * self.minf_cl);
            if rlx * dac > dclmax {
                rlx = dclmax / dac;
            }
            if rlx * dac < dclmin {
                rlx = dclmin / dac;
            }
        } else {
            dac = (clnew - self.clspec) / (0.0 - cl_ac - cl_a);
            if rlx * dac > dalmax {
                rlx = dalmax / dac;
            }
            if rlx * dac < dalmin {
                rlx = dalmin / dac;
            }
        }

        let mut rmsbl = 0.0;
        let mut rmxbl: f64 = 0.0;
        let dhi = 1.5;
        let dlo = -0.5;
        for is in 1..=2 {
            for ibl in 2..=self.nbl[is] {
                let iv = self.isys[is][ibl];
                let vd = self.vdel[iv];
                let dctau = vd[0][0] - dac * vd[0][1];
                let dthet = vd[1][0] - dac * vd[1][1];
                let dmass = vd[2][0] - dac * vd[2][1];
                let duedg = unew[is][ibl] + dac * u_ac[is][ibl] - self.uedg[is][ibl];
                let ddstr = (dmass - self.dstr[is][ibl] * duedg) / self.uedg[is][ibl];
                let dn1 = if ibl < self.itran[is] { dctau / 10.0 } else { dctau / self.ctau[is][ibl] };
                let dn2 = dthet / self.thet[is][ibl];
                let dn3 = ddstr / self.dstr[is][ibl];
                let dn4 = duedg.abs() / 0.25;
                rmsbl += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;

                let rdn1 = rlx * dn1;
                if dn1.abs() > rmxbl.abs() {
                    rmxbl = dn1;
                    self.vmxbl = if ibl < self.itran[is] { 'n' } else { 'C' };
                    self.imxbl = ibl;
                    self.ismxbl = is;
                }
                if rdn1 > dhi {
                    rlx = dhi / dn1;
                }
                if rdn1 < dlo {
                    rlx = dlo / dn1;
                }
                let rdn2 = rlx * dn2;
                if dn2.abs() > rmxbl.abs() {
                    rmxbl = dn2;
                    self.vmxbl = 'T';
                    self.imxbl = ibl;
                    self.ismxbl = is;
                }
                if rdn2 > dhi {
                    rlx = dhi / dn2;
                }
                if rdn2 < dlo {
                    rlx = dlo / dn2;
                }
                let rdn3 = rlx * dn3;
                if dn3.abs() > rmxbl.abs() {
                    rmxbl = dn3;
                    self.vmxbl = 'D';
                    self.imxbl = ibl;
                    self.ismxbl = is;
                }
                if rdn3 > dhi {
                    rlx = dhi / dn3;
                }
                if rdn3 < dlo {
                    rlx = dlo / dn3;
                }
                let rdn4 = rlx * dn4;
                if dn4.abs() > rmxbl.abs() {
                    rmxbl = duedg;
                    self.vmxbl = 'U';
                    self.imxbl = ibl;
                    self.ismxbl = is;
                }
                if rdn4 > dhi {
                    rlx = dhi / dn4;
                }
                if rdn4 < dlo {
                    rlx = dlo / dn4;
                }
            }
        }
        self.rmsbl = (rmsbl / (4.0 * (self.nbl[1] + self.nbl[2]) as f64)).sqrt();
        self.rmxbl = rmxbl;
        self.rlx = rlx;

        if self.lalfa {
            self.cl += rlx * dac;
        } else {
            self.alfa += rlx * dac;
            self.adeg = self.alfa / DTOR;
        }

        for is in 1..=2 {
            for ibl in 2..=self.nbl[is] {
                let iv = self.isys[is][ibl];
                let vd = self.vdel[iv];
                let dctau = vd[0][0] - dac * vd[0][1];
                let dthet = vd[1][0] - dac * vd[1][1];
                let dmass = vd[2][0] - dac * vd[2][1];
                let duedg = unew[is][ibl] + dac * u_ac[is][ibl] - self.uedg[is][ibl];
                let ddstr = (dmass - self.dstr[is][ibl] * duedg) / self.uedg[is][ibl];
                self.ctau[is][ibl] += rlx * dctau;
                self.thet[is][ibl] += rlx * dthet;
                self.dstr[is][ibl] += rlx * ddstr;
                self.uedg[is][ibl] += rlx * duedg;
                let dswaki = if ibl > self.iblte[is] { self.wgap[ibl - self.iblte[is]] } else { 0.0 };
                if ibl >= self.itran[is] {
                    self.ctau[is][ibl] = self.ctau[is][ibl].min(0.25);
                }
                let hklim = if ibl <= self.iblte[is] { 1.02 } else { 1.00005 };
                let ue = self.uedg[is][ibl];
                let msq = ue * ue * hstinv / (self.gamm1 * (1.0 - 0.5 * ue * ue * hstinv));
                let mut dsw = self.dstr[is][ibl] - dswaki;
                dslim(&mut dsw, self.thet[is][ibl], msq, hklim);
                self.dstr[is][ibl] = dsw + dswaki;
                self.mass[is][ibl] = self.dstr[is][ibl] * self.uedg[is][ibl];
            }
            for ibl in 3..=self.iblte[is] {
                if self.uedg[is][ibl - 1] > 0.0 && self.uedg[is][ibl] <= 0.0 {
                    self.uedg[is][ibl] = self.uedg[is][ibl - 1];
                    self.mass[is][ibl] = self.dstr[is][ibl] * self.uedg[is][ibl];
                }
            }
        }
        for kbl in 1..=self.nbl[2] - self.iblte[2] {
            let (a, b) = (self.iblte[1] + kbl, self.iblte[2] + kbl);
            self.ctau[1][a] = self.ctau[2][b];
            self.thet[1][a] = self.thet[2][b];
            self.dstr[1][a] = self.dstr[2][b];
            self.uedg[1][a] = self.uedg[2][b];
            self.tau[1][a] = self.tau[2][b];
            self.dis[1][a] = self.dis[2][b];
            self.ctq[1][a] = self.ctq[2][b];
            self.delt[1][a] = self.delt[2][b];
            self.tstr[1][a] = self.tstr[2][b];
        }
    }
}
