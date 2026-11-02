//! MDES complex-mapping full-inverse design menu (port of xmdes.f and the
//! circle-plane machinery of CIRCLE.INC).
//!
//! The Qspec arrays and other design variables shared with QDES live in
//! `Session::qdes_state` (see `crate::qdes::QdesState`); the circle-plane
//! (CIRCLE.INC) arrays live in `MdesState`.
//!
//! Graphics are not available: plot-only commands are accepted and ignored,
//! cursor-driven commands (MODI, WIND, MARK without arguments) print a note.

use crate::cli::{args, Session};
use crate::fmt::{f_fmt, fixed};
use crate::geom::geopar;
use crate::qdes::{iqsget, ld_int, ld_real, qincom, size_cmd, QdesState, IPX};
use crate::spline::{curv, d2val, deval, scalc, seval, sinvrt, splind};
use crate::xfoil::{XFoil, DTOR, IQX};
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

/// Number of circle-plane points (2^n + 1).
pub const ICX: usize = 257;
/// Number of complex mapping coefficients Cn.
pub const IMX: usize = (ICX - 1) / 4;
const IMX4: usize = IMX / 4;
const PI: f64 = std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Minimal complex arithmetic with gfortran (-fcx-fortran-rules) semantics.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct C64 {
    pub re: f64,
    pub im: f64,
}

#[inline]
pub const fn c(re: f64, im: f64) -> C64 {
    C64 { re, im }
}

impl C64 {
    #[inline]
    pub fn conj(self) -> C64 {
        c(self.re, -self.im)
    }
    #[inline]
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    #[inline]
    pub fn exp(self) -> C64 {
        let e = self.re.exp();
        c(e * self.im.cos(), e * self.im.sin())
    }
    #[inline]
    pub fn ln(self) -> C64 {
        c(self.abs().ln(), self.im.atan2(self.re))
    }
    /// complex * real
    #[inline]
    pub fn sc(self, r: f64) -> C64 {
        c(self.re * r, self.im * r)
    }
    /// complex / real
    #[inline]
    pub fn dr(self, r: f64) -> C64 {
        c(self.re / r, self.im / r)
    }
}

impl Add for C64 {
    type Output = C64;
    #[inline]
    fn add(self, b: C64) -> C64 {
        c(self.re + b.re, self.im + b.im)
    }
}
impl Sub for C64 {
    type Output = C64;
    #[inline]
    fn sub(self, b: C64) -> C64 {
        c(self.re - b.re, self.im - b.im)
    }
}
impl Neg for C64 {
    type Output = C64;
    #[inline]
    fn neg(self) -> C64 {
        c(-self.re, -self.im)
    }
}
impl Mul for C64 {
    type Output = C64;
    #[inline]
    fn mul(self, b: C64) -> C64 {
        c(self.re * b.re - self.im * b.im, self.re * b.im + self.im * b.re)
    }
}
impl Div for C64 {
    type Output = C64;
    /// Smith's algorithm, as GCC expands complex division for Fortran.
    #[inline]
    fn div(self, b: C64) -> C64 {
        let (ar, ai, br, bi) = (self.re, self.im, b.re, b.im);
        if br.abs() < bi.abs() {
            let ratio = br / bi;
            let div = br * ratio + bi;
            c((ar * ratio + ai) / div, (ai * ratio - ar) / div)
        } else {
            let ratio = bi / br;
            let div = bi * ratio + br;
            c((ai * ratio + ar) / div, (ai - ar * ratio) / div)
        }
    }
}
impl AddAssign for C64 {
    #[inline]
    fn add_assign(&mut self, b: C64) {
        *self = *self + b;
    }
}
impl SubAssign for C64 {
    #[inline]
    fn sub_assign(&mut self, b: C64) {
        *self = *self - b;
    }
}

// ---------------------------------------------------------------------------
// Circle-plane state (CIRCLE.INC)
// ---------------------------------------------------------------------------

/// Persistent MDES / circle-plane state.
#[derive(Clone, Debug)]
pub struct MdesState {
    /// unit-circle exp(inw) array initialized
    pub leiw: bool,
    /// requested number of circle-plane points
    pub nc1: usize,
    pub nc: usize,
    pub mc: usize,
    pub mct: usize,
    pub agte: f64,
    pub ag0: f64,
    pub qim0: f64,
    pub qimold: f64,
    pub dwc: f64,
    pub wc: Vec<f64>,
    pub sc: Vec<f64>,
    pub scold: Vec<f64>,
    pub xcold: Vec<f64>,
    pub ycold: Vec<f64>,
    pub dzte: C64,
    pub chordz: C64,
    pub zleold: C64,
    pub zcoldw: Vec<C64>,
    pub zc: Vec<C64>,
    /// zc_cn[ic][m], m = 1..IMX/4
    pub zc_cn: Vec<[C64; IMX4 + 1]>,
    pub piq: Vec<C64>,
    /// cn[m], m = 0..IMX
    pub cn: Vec<C64>,
    /// eiw[ic*(IMX+1) + m]
    pub eiw: Vec<C64>,
}

impl Default for MdesState {
    fn default() -> Self {
        // circle plane array size (257, or largest 2^n + 1 that fits)
        let ann = ((2 * IQX - 1) as f64).ln() / 2.0f64.ln();
        let nn = (ann + 0.00001) as u32;
        let nc1 = (2usize.pow(nn) + 1).min(257);
        MdesState {
            leiw: false,
            nc1,
            nc: 0,
            mc: 0,
            mct: 0,
            agte: 0.0,
            ag0: 0.0,
            qim0: 0.0,
            qimold: 0.0,
            dwc: 0.0,
            wc: vec![0.0; ICX + 1],
            sc: vec![0.0; ICX + 1],
            scold: vec![0.0; ICX + 1],
            xcold: vec![0.0; ICX + 1],
            ycold: vec![0.0; ICX + 1],
            dzte: C64::default(),
            chordz: C64::default(),
            zleold: C64::default(),
            zcoldw: vec![C64::default(); ICX + 1],
            zc: vec![C64::default(); ICX + 1],
            zc_cn: vec![[C64::default(); IMX4 + 1]; ICX + 1],
            piq: vec![C64::default(); ICX + 1],
            cn: vec![C64::default(); IMX + 1],
            eiw: vec![C64::default(); (ICX + 1) * (IMX + 1)],
        }
    }
}

/// ZLEFIND: leading edge of the complex-plane airfoil.
fn zlefind(zc: &[C64], wc: &[f64], nc: usize, piq: &[C64], agte: f64) -> C64 {
    const NTX: usize = 33;
    let zte = (zc[1] + zc[nc]).sc(0.5);

    let mut dmax = 0.0;
    let mut icle = 1usize;
    for ic in 1..=nc {
        let dist = (zc[ic] - zte).abs();
        if dist > dmax {
            dmax = dist;
            icle = ic;
        }
    }

    let ic1 = if icle > (NTX - 1) / 2 { icle - (NTX - 1) / 2 } else { 1 }.max(1);
    let ic2 = (icle + (NTX - 1) / 2).min(nc);

    let dzdw_at = |ic: usize| -> C64 {
        let sinw = 2.0 * (0.5 * wc[ic]).sin();
        let sinwe = sinw.powf(1.0 - agte);
        let hwc = 0.5 * (wc[ic] - PI) * (1.0 + agte) - 0.5 * PI;
        (piq[ic] + c(0.0, hwc)).exp().sc(sinwe)
    };
    let dzdw1 = dzdw_at(ic1);
    let dzdw2 = dzdw_at(ic2);

    let mut xc = [0.0f64; NTX + 1];
    let mut yc = [0.0f64; NTX + 1];
    let mut xcw = [0.0f64; NTX + 1];
    let mut ycw = [0.0f64; NTX + 1];
    for ic in ic1..=ic2 {
        let i = ic - ic1 + 1;
        xc[i] = zc[ic].re;
        yc[i] = zc[ic].im;
    }

    let nic = ic2 - ic1 + 1;
    let w = &wc[ic1 - 1..];
    splind(&xc, &mut xcw, w, nic, dzdw1.re, dzdw2.re);
    splind(&yc, &mut ycw, w, nic, dzdw1.im, dzdw2.im);

    let xcte = 0.5 * (zc[1] + zc[nc]).re;
    let ycte = 0.5 * (zc[1] + zc[nc]).im;

    let mut wcle = wc[icle];
    let mut ok = false;
    for _ in 0..10 {
        let xcle = seval(wcle, &xc, &xcw, w, nic);
        let ycle = seval(wcle, &yc, &ycw, w, nic);
        let dxdw = deval(wcle, &xc, &xcw, w, nic);
        let dydw = deval(wcle, &yc, &ycw, w, nic);
        let dxdd = d2val(wcle, &xc, &xcw, w, nic);
        let dydd = d2val(wcle, &yc, &ycw, w, nic);

        let xchord = xcle - xcte;
        let ychord = ycle - ycte;

        let res = xchord * dxdw + ychord * dydw;
        let resw = dxdw * dxdw + dydw * dydw + xchord * dxdd + ychord * dydd;

        let dwcle = -res / resw;
        wcle += dwcle;

        if dwcle.abs() < 1.0e-5 {
            ok = true;
            break;
        }
    }
    if !ok {
        crate::outln!(" ZLEFIND: LE location failed.");
        wcle = wc[icle];
    }

    let xcle = seval(wcle, &xc, &xcw, w, nic);
    let ycle = seval(wcle, &yc, &ycw, w, nic);
    c(xcle, ycle)
}

impl MdesState {
    #[inline]
    fn eiw(&self, ic: usize, m: usize) -> C64 {
        self.eiw[ic * (IMX + 1) + m]
    }

    /// EIWSET: circle-plane coordinates and exp(inw) array.
    pub fn eiwset(&mut self, nc1: usize) {
        self.nc = nc1;
        self.mc = nc1 / 4;
        self.mct = nc1 / 16;
        if self.nc > ICX {
            panic!("EIWSET: Array overflow. Increase ICX.");
        }
        let nc = self.nc;
        self.dwc = 2.0 * PI / ((nc - 1) as f64);
        for ic in 1..=nc {
            self.wc[ic] = self.dwc * ((ic - 1) as f64);
        }
        for ic in 1..=nc {
            self.eiw[ic * (IMX + 1)] = c(1.0, 0.0);
        }
        for ic in 1..=nc {
            self.eiw[ic * (IMX + 1) + 1] = c(0.0, self.wc[ic]).exp();
        }
        for m in 2..=self.mc {
            for ic in 1..=nc {
                let ic1 = (m * (ic - 1)) % (nc - 1) + 1;
                self.eiw[ic * (IMX + 1) + m] = self.eiw[ic1 * (IMX + 1) + 1];
            }
        }
    }

    /// FTP: slow Fourier transform of P(w) (trapezoidal integration).
    pub fn ftp(&mut self) {
        let nc = self.nc;
        for m in 0..=self.mc {
            let mut zsum = c(0.0, 0.0);
            for ic in 2..nc {
                zsum += self.piq[ic] * self.eiw(ic, m);
            }
            self.cn[m] = ((self.piq[1] * self.eiw(1, m) + self.piq[nc] * self.eiw(nc, m)).sc(0.5) + zsum).sc(self.dwc).dr(PI);
        }
        self.cn[0] = self.cn[0].sc(0.5);
    }

    /// PIQSUM: inverse transform to get P(w) + iQ(w).
    pub fn piqsum(&mut self) {
        for ic in 1..=self.nc {
            let mut zsum = c(0.0, 0.0);
            for m in 0..=self.mc {
                zsum += self.cn[m] * self.eiw(ic, m).conj();
            }
            self.piq[ic] = zsum;
        }
    }

    /// CNSYMM: eliminate imaginary (camber) parts of mapping coefficients.
    pub fn cnsymm(&mut self) {
        for m in 1..=self.mc {
            self.cn[m] = c(self.cn[m].re, 0.0);
        }
        self.piqsum();
    }

    /// CNFILT: modified Hanning filter on Cn.
    pub fn cnfilt(&mut self, ffilt: f64) {
        if ffilt == 0.0 {
            return;
        }
        for m in 0..=self.mc {
            let freq = (m as f64) / (self.mc as f64);
            let cwt = 0.5 * (1.0 + (PI * freq).cos());
            let mut cwtx = cwt;
            if ffilt > 0.0 {
                cwtx = cwt.powf(ffilt);
            }
            self.cn[m] = self.cn[m].sc(cwtx);
        }
    }

    /// ZCCALC: airfoil geometry z(w) from P(w) + iQ(w), with dz/dCn.
    pub fn zccalc(&mut self, mtest: usize) {
        let nc = self.nc;
        let agte = self.agte;
        let dwc = self.dwc;
        let dzdw = |wc: f64, piq: C64| -> C64 {
            let sinw = 2.0 * (0.5 * wc).sin();
            let mut sinwe = 0.0;
            if sinw > 0.0 {
                sinwe = sinw.powf(1.0 - agte);
            }
            let hwc = 0.5 * (wc - PI) * (1.0 + agte) - 0.5 * PI;
            (piq + c(0.0, hwc)).exp().sc(sinwe)
        };

        let ic = 1;
        self.zc[ic] = c(4.0, 0.0);
        for m in 1..=mtest {
            self.zc_cn[ic][m] = c(0.0, 0.0);
        }
        let mut dzdw1 = dzdw(self.wc[ic], self.piq[ic]);
        for ic in 2..=nc {
            let dzdw2 = dzdw(self.wc[ic], self.piq[ic]);

            self.zc[ic] = (dzdw1 + dzdw2).sc(0.5).sc(dwc) + self.zc[ic - 1];
            let dz_piq1 = dzdw1.sc(0.5).sc(dwc);
            let dz_piq2 = dzdw2.sc(0.5).sc(dwc);

            for m in 1..=mtest {
                self.zc_cn[ic][m] =
                    dz_piq1 * self.eiw(ic - 1, m).conj() + dz_piq2 * self.eiw(ic, m).conj() + self.zc_cn[ic - 1][m];
            }
            dzdw1 = dzdw2;
        }

        self.sc[1] = 0.0;
        for ic in 2..=nc {
            self.sc[ic] = self.sc[ic - 1] + (self.zc[ic] - self.zc[ic - 1]).abs();
        }
        let scn = self.sc[nc];
        for ic in 1..=nc {
            self.sc[ic] /= scn;
        }
    }

    /// ZCNORM: normalize z(w) to the old chord and angle, reset dz/dCn.
    pub fn zcnorm(&mut self, mtest: usize) {
        let nc = self.nc;
        let zle = zlefind(&self.zc, &self.wc, nc, &self.piq, self.agte);

        for ic in 1..=nc {
            self.zc[ic] -= zle;
        }

        let zte = (self.zc[1] + self.zc[nc]).sc(0.5);
        let mut zte_cn = [c(0.0, 0.0); IMX4 + 1];
        for m in 1..=mtest {
            zte_cn[m] = (self.zc_cn[1][m] + self.zc_cn[nc][m]).sc(0.5);
        }

        let chordz = self.chordz;
        for ic in 1..=nc {
            let zcnew = chordz * self.zc[ic] / zte;
            let zc_zte = -(zcnew / zte);
            self.zc[ic] = zcnew;
            for m in 1..=mtest {
                self.zc_cn[ic][m] = chordz * self.zc_cn[ic][m] / zte + zc_zte * zte_cn[m];
            }
        }

        let qimoff = -(chordz / zte).ln().im;
        self.cn[0] = self.cn[0] - c(0.0, qimoff);

        for ic in 1..=nc {
            self.zc[ic] += self.zleold;
        }
    }

    /// SCINIT: circle-plane coordinate s(w) for the current airfoil; also
    /// produces the Cn coefficients for its geometry.
    #[allow(clippy::too_many_arguments)]
    pub fn scinit(&mut self, n: usize, x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64], sle: f64) {
        const CEPS: f64 = 1.0e-7;
        const SEPS: f64 = 5.0e-7;
        let nc = self.nc;

        self.agte = (xp[n].atan2(-yp[n]) - xp[1].atan2(-yp[1])) / PI - 1.0;
        self.ag0 = xp[1].atan2(-yp[1]);
        self.qim0 = self.ag0 + 0.5 * PI * (1.0 + self.agte);

        let xle = seval(sle, x, xp, s, n);
        let yle = seval(sle, y, yp, s, n);

        let dxte = x[1] - x[n];
        let dyte = y[1] - y[n];
        self.dzte = c(dxte, dyte);

        let chordx = 0.5 * (x[1] + x[n]) - xle;
        let chordy = 0.5 * (y[1] + y[n]) - yle;
        self.chordz = c(chordx, chordy);
        self.zleold = c(xle, yle);

        crate::outln!(
            "\n Current TE gap  dx dy ={}{}    TE angle ={} deg.\n",
            f_fmt(self.dzte.re, 7, 4),
            f_fmt(self.dzte.im, 7, 4),
            f_fmt(self.agte * 180.0, 7, 3)
        );
        crate::outln!(" Initializing mapping coordinate ...");

        let cvle = curv(sle, x, xp, y, yp, s, n) * s[n];
        let cvabs = cvle.abs();
        let dsdwle = (0.5 / cvabs).max(1.0e-3);

        let tops = sle / s[n];
        let bots = (s[n] - sle) / s[n];

        let wwt = 1.0 - 2.0 * dsdwle / tops;
        for ic in 1..=(nc - 1) / 2 + 1 {
            self.sc[ic] = tops * (1.0 - (wwt * self.wc[ic]).cos()) / (1.0 - (wwt * PI).cos());
        }
        let wwt = 1.0 - 2.0 * dsdwle / bots;
        for ic in (nc - 1) / 2 + 2..=nc {
            self.sc[ic] = 1.0 - bots * (1.0 - (wwt * (self.wc[nc] - self.wc[ic])).cos()) / (1.0 - (wwt * PI).cos());
        }

        for ipass in 1..=30 {
            for ic in 1..=nc {
                let sic = s[1] + (s[n] - s[1]) * self.sc[ic];
                let dxds = deval(sic, x, xp, s, n);
                let dyds = deval(sic, y, yp, s, n);
                let qim = dxds.atan2(-dyds) - 0.5 * (self.wc[ic] - PI) * (1.0 + self.agte) - self.qim0;
                self.piq[ic] = c(0.0, qim);
            }

            self.ftp();
            self.cn[0] = c(0.0, self.cn[0].im + self.qim0);
            self.piqsum();

            for ic in 1..=nc {
                self.scold[ic] = self.sc[ic];
            }

            for _itgap in 1..=5 {
                self.zccalc(1);
                let zle = zlefind(&self.zc, &self.wc, nc, &self.piq, self.agte);
                let zte = (self.zc[1] + self.zc[nc]).sc(0.5);
                let dzwt = (zte - zle).abs() / self.chordz.abs();
                let dcn = -((self.zc[1] - self.zc[nc] - self.dzte.sc(dzwt)) / (self.zc_cn[1][1] - self.zc_cn[nc][1]));
                self.cn[1] += dcn;
                self.piqsum();
                if dcn.abs() < CEPS {
                    break;
                }
            }

            let mut dscmax = 0.0f64;
            for ic in 1..=nc {
                dscmax = dscmax.max((self.sc[ic] - self.scold[ic]).abs());
            }
            crate::outln!("{}      max(dw) ={}", ld_int(ipass), ld_real(dscmax));
            if dscmax < SEPS {
                break;
            }
        }

        self.zcnorm(1);

        for ic in 1..=nc {
            self.scold[ic] = self.sc[ic];
            self.xcold[ic] = self.zc[ic].re;
            self.ycold[ic] = self.zc[ic].im;
        }

        for ic in 1..=nc {
            let sinw = 2.0 * (0.5 * self.wc[ic]).sin();
            let mut sinwe = 0.0;
            if sinw > 0.0 {
                sinwe = sinw.powf(1.0 - self.agte);
            }
            let hwc = 0.5 * (self.wc[ic] - PI) * (1.0 + self.agte) - 0.5 * PI;
            self.zcoldw[ic] = (self.piq[ic] + c(0.0, hwc)).exp().sc(sinwe);
        }

        self.qimold = self.cn[0].im;
    }

    /// CNCALC: Cn coefficients from surface speed q(w).
    pub fn cncalc(&mut self, qc: &[f64], lsymm: bool) {
        let nc = self.nc;
        if nc > ICX {
            panic!("CNCALC: Array overflow.");
        }
        let mut qcw = vec![0.0; ICX + 1];
        splind(qc, &mut qcw, &self.wc, nc, -999.0, -999.0);

        let mut ic = 2;
        while ic <= nc {
            if qc[ic] < 0.0 {
                break;
            }
            ic += 1;
        }
        let mut wcle = self.wc[ic.min(nc)];
        sinvrt(&mut wcle, 0.0, qc, &qcw, &self.wc, nc);

        let alfcir = 0.5 * (wcle - PI);

        for ic in 2..nc {
            let cosw = 2.0 * (0.5 * self.wc[ic] - alfcir).cos();
            let sinw = 2.0 * (0.5 * self.wc[ic]).sin();
            let sinwe = sinw.powf(self.agte);
            let pfun = if cosw.abs() < 1.0e-4 { (sinwe / qcw[ic]).abs() } else { (cosw * sinwe / qc[ic]).abs() };
            self.piq[ic] = c(pfun.ln(), 0.0);
        }

        self.piq[1] = self.piq[2].sc(3.0) - self.piq[3].sc(3.0) + self.piq[4];
        self.piq[nc] = self.piq[nc - 1].sc(3.0) - self.piq[nc - 2].sc(3.0) + self.piq[nc - 3];

        let cnsav: Vec<C64> = self.cn[..=self.mc].to_vec();

        self.ftp();
        self.cn[0] = c(0.0, self.qimold);

        if lsymm {
            for m in 1..=self.mc {
                let cnr = 2.0 * self.cn[m].re - cnsav[m].re;
                self.cn[m] = c(cnr, 0.0);
            }
        }

        self.piqsum();
    }

    /// QCCALC: surface speed from the complex speed function for a
    /// prescribed alpha (ispec=1) or CL (otherwise).  Returns NCIR.
    #[allow(clippy::too_many_arguments)]
    pub fn qccalc(
        &self,
        ispec: i32,
        alfa: &mut f64,
        cl: &mut f64,
        cm: &mut f64,
        minf: f64,
        qinf: f64,
        xcir: &mut [f64],
        ycir: &mut [f64],
        scir: &mut [f64],
        qcir: &mut [f64],
    ) -> usize {
        const AEPS: f64 = 5.0e-7;
        let nc = self.nc;
        let beta = (1.0 - minf * minf).sqrt();
        let bfac = 0.5 * minf * minf / (1.0 + beta);
        let mut qc_a = vec![0.0; ICX + 1];
        let mut dalfa = 0.0;

        for _ipass in 1..=10 {
            let alfcir = *alfa - self.cn[0].im;

            let mut cmt = c(0.0, 0.0);
            let mut cft = c(0.0, 0.0);
            let mut cft_a = c(0.0, 0.0);

            for ic in 1..=nc {
                let ppp = self.piq[ic].re;
                let eppp = (-ppp).exp();
                let sinw = 2.0 * (0.5 * self.wc[ic]).sin();
                let sinwe = if self.agte == 0.0 {
                    1.0
                } else if sinw > 0.0 {
                    sinw.powf(self.agte)
                } else {
                    0.0
                };
                qcir[ic] = 2.0 * (0.5 * self.wc[ic] - alfcir).cos() * sinwe * eppp;
                qc_a[ic] = 2.0 * (0.5 * self.wc[ic] - alfcir).sin() * sinwe * eppp;
                xcir[ic] = self.zc[ic].re;
                ycir[ic] = self.zc[ic].im;
                scir[ic] = self.sc[ic];
            }

            let ic = 1;
            let cpinc1 = 1.0 - (qcir[ic] / qinf) * (qcir[ic] / qinf);
            let cpi_q1 = -2.0 * qcir[ic] / (qinf * qinf);
            let mut cpcom1 = cpinc1 / (beta + bfac * cpinc1);
            let cpc_q1 = (1.0 - bfac * cpcom1) / (beta + bfac * cpinc1) * cpi_q1;
            let mut cpc_a1 = cpc_q1 * qc_a[ic];
            for ic in 1..=nc {
                let icp = if ic == nc { 1 } else { ic + 1 };
                let cpinc2 = 1.0 - (qcir[icp] / qinf) * (qcir[icp] / qinf);
                let cpi_q2 = -2.0 * qcir[icp] / (qinf * qinf);
                let cpcom2 = cpinc2 / (beta + bfac * cpinc2);
                let cpc_q2 = (1.0 - bfac * cpcom2) / (beta + bfac * cpinc2) * cpi_q2;
                let cpc_a2 = cpc_q2 * qc_a[icp];

                let za = (self.zc[icp] + self.zc[ic]).sc(0.5) - c(0.25, 0.0);
                let dz = self.zc[icp] - self.zc[ic];

                cmt = cmt - dz.sc(0.5 * (cpcom1 + cpcom2)) * za.conj() + (dz.sc(cpcom1 - cpcom2) * dz.conj()).dr(12.0);
                cft += dz.sc(0.5 * (cpcom1 + cpcom2));
                cft_a += dz.sc(0.5 * (cpc_a1 + cpc_a2));

                cpcom1 = cpcom2;
                cpc_a1 = cpc_a2;
            }

            let eia = c(0.0, -*alfa).exp();
            cft = cft * eia;
            cft_a = cft_a * eia + cft * c(0.0, -1.0);

            let clt = cft.re;
            let clt_a = cft_a.re;
            *cm = cmt.re;

            if ispec == 1 {
                *cl = clt;
                return nc;
            } else {
                dalfa = (*cl - clt) / clt_a;
                *alfa += dalfa;
                if dalfa.abs() < AEPS {
                    return nc;
                }
            }
        }
        crate::outln!(" QCCALC: CL convergence failed.  dAlpha ={}", ld_real(dalfa));
        nc
    }

    /// QSPINT: integrate circle-plane surface pressures for CL, CM.
    pub fn qspint(&self, alqsp: f64, qspec: &[f64], qinf: f64, minf: f64) -> (f64, f64) {
        let nc = self.nc;
        let sa = alqsp.sin();
        let ca = alqsp.cos();
        let beta = (1.0 - minf * minf).sqrt();
        let bfac = 0.5 * minf * minf / (1.0 + beta);

        let mut clqsp = 0.0;
        let mut cmqsp = 0.0;

        let i = 1;
        let cqinc = 1.0 - (qspec[i] / qinf) * (qspec[i] / qinf);
        let mut cpq1 = cqinc / (beta + bfac * cqinc);
        let (xc, yc) = (&self.xcold, &self.ycold);
        for i in 1..=nc {
            let ip = if i == nc { 1 } else { i + 1 };
            let cqinc = 1.0 - (qspec[ip] / qinf) * (qspec[ip] / qinf);
            let cpq2 = cqinc / (beta + bfac * cqinc);

            let dx = (xc[ip] - xc[i]) * ca + (yc[ip] - yc[i]) * sa;
            let dy = (yc[ip] - yc[i]) * ca - (xc[ip] - xc[i]) * sa;
            let du = cpq2 - cpq1;

            let ax = 0.5 * (xc[ip] + xc[i]) * ca + 0.5 * (yc[ip] + yc[i]) * sa;
            let ay = 0.5 * (yc[ip] + yc[i]) * ca - 0.5 * (xc[ip] + xc[i]) * sa;
            let aq = 0.5 * (cpq2 + cpq1);

            clqsp += dx * aq;
            cmqsp = cmqsp - dx * (aq * (ax - 0.25) + du * dx / 12.0) - dy * (aq * ay + du * dy / 12.0);

            cpq1 = cpq2;
        }
        (clqsp, cmqsp)
    }

    /// Newton iteration on CN(1) to fix the TE gap (shared by MAPGEN and PERT).
    fn tegap_newton(&mut self, reset_qq: bool, failmsg: &str) {
        let nc = self.nc;
        let mct = self.mct;
        let mut converged = false;
        for itercn in 1..=10 {
            // single free coefficient (NCN = 1): CGAUSS reduces to one division
            let _ = reset_qq;
            let mut dcn = self.zc[1] - self.zc[nc] - self.dzte;
            let qq = self.zc_cn[1][1] - self.zc_cn[nc][1];
            dcn = dcn / qq;

            self.cn[1] -= dcn;
            let dcnmax = dcn.abs().max(0.0);

            self.piqsum();
            self.zccalc(mct);
            self.zcnorm(mct);

            crate::outln!("{}{}", ld_int(itercn), ld_real(dcnmax));
            if dcnmax <= 5.0e-5 {
                converged = true;
                break;
            }
        }
        if !converged {
            crate::outln!("{}", failmsg);
        }
    }

    /// MAPGEN: geometry from the Cn coefficients, adjusting them to meet
    /// the TE-gap constraint.  Writes the new airfoil into x,y and returns n.
    pub fn mapgen(&mut self, x: &mut [f64], y: &mut [f64]) -> usize {
        let dx = self.xcold[2] - self.xcold[1];
        let dy = self.ycold[2] - self.ycold[1];
        self.qim0 = dx.atan2(-dy) + 0.5 * PI * (1.0 + self.agte);
        let qimoff = self.qim0 - self.cn[0].im;
        self.cn[0] = self.cn[0] + c(0.0, qimoff);

        self.piqsum();
        let mct = self.mct;
        self.zccalc(mct);
        self.zcnorm(mct);

        self.cn[0] = c(0.0, self.cn[0].im);

        self.tegap_newton(true, "\n MAPGEN: Geometric constraints not fully converged");

        let nc = self.nc;
        for i in 1..=nc {
            x[i] = self.zc[i].re;
            y[i] = self.zc[i].im;
        }
        nc
    }

    /// CNDUMP: write Fourier coefficients.
    pub fn cndump(&self) -> String {
        let mut out = String::new();
        for m in 0..=self.mc {
            out += &format!(
                " {:3}{}{}{}{}\n",
                m,
                f_fmt(self.cn[m].re, 11, 6),
                f_fmt(self.cn[m].im, 11, 6),
                f_fmt(self.piq[m + 1].re, 11, 6),
                f_fmt(self.piq[m + 1].im, 11, 6)
            );
        }
        for m in self.mc + 1..self.nc {
            out += &format!(
                " {:3}{}{}{}{}\n",
                m,
                f_fmt(0.0, 11, 6),
                f_fmt(0.0, 11, 6),
                f_fmt(self.piq[m + 1].re, 11, 6),
                f_fmt(self.piq[m + 1].im, 11, 6)
            );
        }
        out
    }
}

// ---------------------------------------------------------------------------
// routines touching both the design state and the circle-plane state
// ---------------------------------------------------------------------------

/// MAPGAM: mapped Q for the current airfoil (IAC=1: alpha, 2: CL).
fn mapgam(m: &MdesState, q: &mut QdesState, xf: &XFoil, iac: i32) {
    let mut w1 = vec![0.0; ICX + 1];
    let mut w2 = vec![0.0; ICX + 1];
    let mut w5 = vec![0.0; ICX + 1];
    let mut w6 = vec![0.0; ICX + 1];
    let (mut alg, mut clg, mut cmg) = (q.algam, q.clgam, q.cmgam);
    q.nsp = m.qccalc(iac, &mut alg, &mut clg, &mut cmg, xf.minf, xf.qinf, &mut w1, &mut w2, &mut w5, &mut w6);
    q.algam = alg;
    q.clgam = clg;
    q.cmgam = cmg;

    let n = xf.n;
    let chx = xf.xte - xf.xle;
    let chy = xf.yte - xf.yle;
    let chsq = chx * chx + chy * chy;
    for i in 1..=q.nsp {
        q.qspec[0][i] = w6[i];
        q.sspec[i] = w5[i];
        let xic = seval(xf.s[n] * q.sspec[i], &xf.x, &xf.xp, &xf.s, n);
        let yic = seval(xf.s[n] * q.sspec[i], &xf.y, &xf.yp, &xf.s, n);
        q.xspoc[i] = ((xic - xf.xle) * chx + (yic - xf.yle) * chy) / chsq;
        q.yspoc[i] = ((yic - xf.yle) * chx - (xic - xf.xle) * chy) / chsq;
    }
    q.ssple = xf.sle / xf.s[n];
}

/// QSPCIR: Qspec arrays for all design alphas or CLs.
fn qspcir(m: &MdesState, q: &mut QdesState, xf: &mut XFoil) {
    let mut w1 = vec![0.0; ICX + 1];
    let mut w2 = vec![0.0; ICX + 1];
    let mut w5 = vec![0.0; ICX + 1];
    for k in 1..=q.nqsp {
        let iacqsp = q.iacqsp;
        let QdesState { alqsp, clqsp, cmqsp, qspec, .. } = &mut *q;
        let nsp = m.qccalc(
            iacqsp,
            &mut alqsp[k],
            &mut clqsp[k],
            &mut cmqsp[k],
            xf.minf,
            xf.qinf,
            &mut w1,
            &mut w2,
            &mut w5,
            &mut qspec[k],
        );
        q.nsp = nsp;
        q.splqsp(k);
    }
    xf.lqspec = true;
}

/// GETVOV: read Qspec(1) from a V/Vinf vs s file.
fn getvov(s: &mut Session, kqsp: usize) {
    let fname = s.input.ask("Enter V/Vinf vs s data filename^");
    let fname = fname.trim().to_string();
    let txt = match std::fs::read_to_string(&fname) {
        Ok(t) => t,
        Err(_) => {
            crate::outln!(" GETVOV: File OPEN error.");
            return;
        }
    };
    let mut w1 = vec![0.0; 1];
    let mut w2 = vec![0.0; 1];
    let mut toks: Vec<f64> = Vec::new();
    for line in txt.lines() {
        if line.trim().is_empty() {
            continue;
        }
        match crate::geom::getflt(line, 20) {
            Some(v) if !v.is_empty() => {
                toks.extend_from_slice(&v);
                while toks.len() >= 2 {
                    w1.push(toks[0]);
                    w2.push(toks[1]);
                    toks.clear();
                }
            }
            _ => {
                crate::outln!(" GETVOV: File READ error.");
                return;
            }
        }
        if w1.len() > IQX {
            break;
        }
    }
    let mut kk = w1.len() - 1;
    if kk < 2 {
        crate::outln!(" GETVOV: File READ error.");
        return;
    }

    let sspan = w1[kk] - w1[1];
    let sstart = w1[1];
    for k in 1..=kk {
        w1[k] = 1.0 - (w1[k] - sstart) / sspan;
    }

    sort(&mut kk, &mut w1, &mut w2);

    let mut w3 = vec![0.0; kk + 1];
    splind(&w2, &mut w3, &w1, kk, -999.0, -999.0);

    let xf = &s.xf;
    let q = &mut s.qdes_state;
    for i in 1..=q.nsp {
        let ss = q.sspec[i];
        let qsnew = seval(ss, &w2, &w3, &w1, kk);
        q.qspec[kqsp][i] = qincom(qsnew, xf.qinf, xf.tklam);
    }
    q.splqsp(kqsp);
}

/// SORT: bubble sort of s(.), w(.) and removal of duplicate pairs.
fn sort(kk: &mut usize, s: &mut [f64], w: &mut [f64]) {
    let mut done = false;
    for _ipass in 1..=1234 {
        done = true;
        for n in 1..*kk {
            let np = n + 1;
            if s[np] >= s[n] {
                continue;
            }
            s.swap(np, n);
            w.swap(np, n);
            done = false;
        }
        if done {
            break;
        }
    }
    if !done {
        crate::outln!(" Sort failed");
    }
    let kks = *kk;
    for k in 1..=kks {
        if k >= *kk {
            return;
        }
        if s[k] != s[k + 1] {
            continue;
        }
        *kk -= 2;
        for kt in k..=*kk {
            s[kt] = s[kt + 2];
            w[kt] = w[kt + 2];
        }
    }
}

/// PERT: perturb one Cn coefficient (user input) and regenerate geometry.
fn pert(s: &mut Session) {
    {
        let m = &mut s.mdes_state;
        let q = &s.qdes_state;
        let qs = q.qspec[1].clone();
        m.cncalc(&qs, false);

        let dx = m.xcold[2] - m.xcold[1];
        let dy = m.ycold[2] - m.ycold[1];
        m.qim0 = dx.atan2(-dy) + 0.5 * PI * (1.0 + m.agte);
        let qimoff = m.qim0 - m.cn[0].im;
        m.cn[0] = m.cn[0] + c(0.0, qimoff);

        crate::outln!();
        crate::outln!(" Current mapping coefficients...");
        crate::outln!("       n    Re(Cn)      Im(Cn)");
        for k in 1..=m.nc.min(32) {
            crate::outln!("    {:4}{}{}", k, f_fmt(m.cn[k].re, 12, 6), f_fmt(m.cn[k].im, 12, 6));
        }
    }

    // read  n, delta(Cnr), delta(Cni)  (list-directed: may span lines)
    let (mm, dcnr, dcni) = loop {
        crate::outln!();
        let mut toks: Vec<String> = Vec::new();
        while toks.len() < 3 {
            let l = s.input.ask("   Enter  n, delta(Cnr), delta(Cni):");
            if s.input.eof {
                return;
            }
            toks.extend(l.split(|ch: char| ch == ',' || ch.is_whitespace()).filter(|t| !t.is_empty()).map(|t| t.to_string()));
        }
        let mm = match toks[0].parse::<i64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let pf = |t: &str| t.replace(['d', 'D'], "E").parse::<f64>();
        let (dr, di) = match (pf(&toks[1]), pf(&toks[2])) {
            (Ok(a), Ok(b)) => (a, b),
            _ => continue,
        };
        if mm <= 0 {
            continue;
        }
        if mm as usize > s.mdes_state.nc {
            crate::outln!(" Max number of modes is{}", ld_int(s.mdes_state.nc as i64));
            continue;
        }
        break (mm as usize, dr, di);
    };

    let m = &mut s.mdes_state;
    if mm <= IMX {
        m.cn[mm] += c(dcnr, dcni);
    } else {
        crate::outln!(" (Foiley: Cn index beyond {} has no effect)", IMX);
    }

    m.piqsum();
    let mct = m.mct;
    m.zccalc(mct);
    m.zcnorm(mct);
    m.cn[0] = c(0.0, m.cn[0].im);

    m.tegap_newton(false, " TE gap,chord did not converge");
}

/// MDES EXEC: generate new buffer airfoil from the mapping.
fn exec(s: &mut Session) {
    let xf = &mut s.xf;
    let m = &mut s.mdes_state;
    let nb = m.mapgen(&mut xf.xb, &mut xf.yb);
    xf.nb = nb;

    scalc(&xf.xb, &xf.yb, &mut xf.sb, nb);
    splind(&xf.xb, &mut xf.xbp, &xf.sb, nb, -999.0, -999.0);
    splind(&xf.yb, &mut xf.ybp, &xf.sb, nb, -999.0, -999.0);

    let g = geopar(&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, nb, true);
    xf.sble = g.sle;
    xf.chordb = g.chord;
    xf.areab = g.area;
    xf.radble = g.radle;
    xf.angbte = g.angte;
    xf.ei11ba = g.ei11a;
    xf.ei22ba = g.ei22a;
    xf.apx1ba = g.apx1a;
    xf.apx2ba = g.apx2a;
    xf.ei11bt = g.ei11t;
    xf.ei22bt = g.ei22t;
    xf.apx1bt = g.apx1t;
    xf.apx2bt = g.apx2t;
    xf.thickb = g.thick;
    xf.cambrb = g.cambr;

    // (Cp / Qspec plotting omitted)
    xf.lgsame = false;

    crate::outln!("\n\n New buffer airfoil generated\n Execute PANE at Top Level to set new current airfoil\n");
}

/// Print current / Qspec alpha, CL, CM after a Qspec change (formats 1200/1210).
fn print_qspec_forces(m: &MdesState, q: &mut QdesState, xf: &XFoil) {
    crate::outln!(
        "\n Current :  alpha ={}    CL ={}    CM ={}",
        f_fmt(q.algam / DTOR, 9, 4),
        f_fmt(q.clgam, 11, 6),
        f_fmt(q.cmgam, 11, 6)
    );
    for k in 1..=q.nqsp {
        let (clq, cmq) = m.qspint(q.alqsp[k], &q.qspec[k], xf.qinf, xf.minf);
        q.cmqsp[k] = cmq;
        if q.iacqsp == 1 {
            q.clqsp[k] = clq;
        }
        crate::outln!(
            " Qspec{:2} :  alpha ={}    CL ={}    CM ={}",
            k,
            f_fmt(q.alqsp[k] / DTOR, 9, 4),
            f_fmt(q.clqsp[k], 11, 6),
            f_fmt(q.cmqsp[k], 11, 6)
        );
    }
}

fn help() {
    crate::outln!(
        "\n   <cr>   Return to Top Level\n   !      Redo previous command\n\n   INIT   Re-initialize mapping\n   QSET   Reset Qspec <== Q\n   AQ r.. Show/select alpha(s) for Qspec\n   CQ r.. Show/select  CL(s)   for Qspec\n\n   Symm   Toggle symmetry flag\n   TGAP r Set new TE gap\n   TANG r Set new TE angle\n\n   Modi   Modify Qspec\n   MARK   Mark off target segment for smoothing\n   SMOO   Smooth Qspec inside target segment\n   Filt   Apply Hanning filter to entire Qspec\n   SLOP   Toggle modified-Qspec slope matching flag\n\n   eXec   Execute  full-inverse calculation\n\n   Visc   Qvis overlay toggle\n   REFL   Reflected Qspec overlay toggle\n   SPEC   Plot mapping coefficient spectrum\n\n   Plot   Replot Qspec (line) and Q (symbols)\n   Blow   Blowup plot region\n   Rese   Reset plot scale and origin\n   Wind     Plot window adjust via cursor and keys\n\n   SIZE r Change absolute plot-object size\n  .ANNO   Annotate plot\n   HARD   Hardcopy current plot\n\n   PERT   Perturb one Cn and generate geometry\n\n   (Foiley: no cursor input -- use  MARK i1 i2  to set the smoothing\n    segment by circle-plane point indices; Modi and Wind are not available)"
    );
}

/// Read a list of reals for AQ/CQ when none were given on the command line.
/// Returns None if just <return> was typed.
fn read_list(s: &mut Session, prompt: &str) -> Option<Vec<f64>> {
    loop {
        let line = s.input.ask(prompt);
        if s.input.eof {
            return None;
        }
        match crate::geom::getflt(&line, IPX) {
            Some(v) => {
                if v.is_empty() {
                    return None;
                }
                let n = v.len().min(IPX);
                return Some(v[..n].to_vec());
            }
            None => continue,
        }
    }
}

/// MDES menu entry point.
#[allow(unused_assignments)]
pub fn mdes(s: &mut Session) {
    if s.xf.n == 0 {
        crate::outln!();
        crate::outln!(" ***  No airfoil available  ***");
        return;
    }

    let mut lcnpl;

    'init: loop {
        // label 1
        lcnpl = false;
        {
            let xf = &mut s.xf;
            let q = &mut s.qdes_state;
            let m = &mut s.mdes_state;

            if q.nsp != m.nc1 {
                xf.lqspec = false;
                q.iq1 = 1;
                q.iq2 = m.nc1;
            }

            if !m.leiw {
                let nc1 = m.nc1;
                m.eiwset(nc1);
            }
            m.leiw = true;

            if q.nqsp == 0 {
                q.iacqsp = 1;
                q.alqsp[1] = xf.alfa;
                q.nqsp = 1;
            }

            if !xf.lscini {
                let n = xf.n;
                m.scinit(n, &xf.x, &xf.xp, &xf.y, &xf.yp, &xf.s, xf.sle);
                xf.lscini = true;
                xf.lqspec = false;
            }

            q.algam = xf.alfa;
            mapgam(m, q, xf, 1);
            crate::outln!(
                "\n Current Q operating condition:   alpha = {}     CL = {}",
                f_fmt(q.algam / DTOR, 7, 3),
                f_fmt(q.clgam, 8, 4)
            );

            if !xf.lqspec {
                let qg = q.qspec[0].clone();
                m.cncalc(&qg, false);
                qspcir(m, q, xf);
                crate::outln!("\n Qspec initialized to current Q");
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
                let prompt = if s.qdes_state.lqsym { ".MDESs" } else { ".MDES" };
                let (cmd, a) = s.input.askc(prompt);
                if s.input.eof {
                    s.qdes_state.lqsym = false;
                    return;
                }
                if cmd.starts_with('!') {
                    if comold == "****" {
                        crate::outln!(" Previous .MDES command not valid");
                        continue;
                    }
                    comand = comold.clone();
                    comarg = argold.clone();
                } else {
                    comand = cmd;
                    comarg = a;
                }
                break;
            }

            // label 505 (re-entered by PERT)
            let mut redo = true;
            while redo {
                redo = false;

                if comand.is_empty() {
                    s.qdes_state.lqsym = false;
                    return;
                }

                let (r, iv) = args(&comarg);
                let ninput = r.len();

                match comand.as_str() {
                    "?" => help(),
                    "Z" | "U" => {}
                    "INIT" => {
                        s.xf.lqspec = false;
                        s.xf.lscini = false;
                        continue 'init;
                    }
                    "QSET" => {
                        let xf = &mut s.xf;
                        let q = &mut s.qdes_state;
                        let m = &mut s.mdes_state;
                        let qg = q.qspec[0].clone();
                        m.cncalc(&qg, false);
                        if q.lqsym {
                            m.cnsymm();
                        }
                        qspcir(m, q, xf);
                        lcnpl = false;
                    }
                    "AQ" | "CQ" => {
                        let is_aq = comand == "AQ";
                        let vals = if ninput >= 1 {
                            Some(r[..ninput.min(IPX)].to_vec())
                        } else {
                            let q = &s.qdes_state;
                            crate::outln!(
                                "\n Current Q operating condition:   alpha = {}     CL = {}",
                                f_fmt(q.algam / DTOR, 7, 3),
                                f_fmt(q.clgam, 8, 4)
                            );
                            let mut l = String::new();
                            for k in 1..=q.nqsp {
                                if is_aq {
                                    l += &f_fmt(q.alqsp[k] / DTOR, 9, 3);
                                } else {
                                    l += &f_fmt(q.clqsp[k], 8, 4);
                                }
                            }
                            if is_aq {
                                crate::outln!("\n Current Qspec alphas  ={}", l);
                                read_list(s, "New alphas or <return>:")
                            } else {
                                crate::outln!("\n Current Qspec CLs  ={}", l);
                                read_list(s, "New CLs or <return>:")
                            }
                        };
                        let vals = match vals {
                            Some(v) => v,
                            None => continue, // just <return>: GO TO 500
                        };
                        let xf = &mut s.xf;
                        let q = &mut s.qdes_state;
                        let m = &s.mdes_state;
                        q.nqsp = vals.len();
                        for (k, v) in vals.iter().enumerate() {
                            if is_aq {
                                q.alqsp[k + 1] = v * DTOR;
                            } else {
                                q.clqsp[k + 1] = *v;
                            }
                        }
                        q.iacqsp = if is_aq { 1 } else { 2 };
                        qspcir(m, q, xf);
                        lcnpl = false;
                    }
                    "SYMM" | "S" => {
                        let xf = &mut s.xf;
                        let q = &mut s.qdes_state;
                        let m = &mut s.mdes_state;
                        q.lqsym = !q.lqsym;
                        if q.lqsym {
                            crate::outln!(" Qspec symmetry forcing enabled.");
                            m.cnsymm();
                            qspcir(m, q, xf);
                            lcnpl = false;
                        } else {
                            crate::outln!(" Qspec symmetry forcing disabled.");
                        }
                    }
                    "TGAP" => {
                        let (dxnew, dynew) = if ninput >= 2 {
                            (r[0], r[1])
                        } else {
                            let dz = s.mdes_state.dzte;
                            crate::outln!("\n Current TE gap  dx/c dy/c ={}{}", f_fmt(dz.re, 7, 4), f_fmt(dz.im, 7, 4));
                            let a = s.input.askr("Enter new TE gap dx/c^", dz.re);
                            let b = s.input.askr("Enter new TE gap dy/c^", dz.im);
                            (a, b)
                        };
                        s.mdes_state.dzte = c(dxnew, dynew);
                    }
                    "TANG" => {
                        // (as in the Fortran, the value is taken from the command
                        //  line only if at least two numbers are given)
                        let agted = if ninput >= 2 {
                            r[0]
                        } else {
                            let a = s.mdes_state.agte * 180.0;
                            crate::outln!("\n Current TE angle ={} deg.", f_fmt(a, 7, 3));
                            s.input.askr("Enter new TE angle (deg)^", a)
                        };
                        s.mdes_state.agte = agted / 180.0;
                    }
                    "VISC" | "V" => {
                        let q = &mut s.qdes_state;
                        q.lqvdes = !q.lqvdes;
                        if q.lqvdes {
                            crate::outln!(" Qspec & Qvis will be plotted");
                        } else {
                            crate::outln!(" Only Qspec will be plotted");
                        }
                        lcnpl = false;
                    }
                    "REFL" => {
                        let q = &mut s.qdes_state;
                        q.lqrefl = !q.lqrefl;
                        if q.lqrefl {
                            crate::outln!(" Reflected Qspec will be plotted");
                        } else {
                            crate::outln!(" Reflected Qspec will not be plotted");
                        }
                        lcnpl = false;
                    }
                    "MODI" | "M" => {
                        crate::outln!(" MODI needs cursor input (not available in Foiley)");
                    }
                    "MARK" => {
                        iqsget(&mut s.qdes_state, &iv);
                        lcnpl = false;
                    }
                    "READ" => {
                        let kqsp = 1;
                        getvov(s, kqsp);
                        let xf = &mut s.xf;
                        let q = &mut s.qdes_state;
                        let m = &mut s.mdes_state;
                        let qs = q.qspec[kqsp].clone();
                        m.cncalc(&qs, false);
                        if q.lqsym {
                            m.cnsymm();
                        }
                        lcnpl = false;
                        let (clq, cmq) = m.qspint(q.alqsp[kqsp], &q.qspec[kqsp], xf.qinf, xf.minf);
                        q.clqsp[kqsp] = clq;
                        q.cmqsp[kqsp] = cmq;
                        crate::outln!(
                            "\n Current :  alpha ={}    CL ={}    CM ={}",
                            f_fmt(q.algam / DTOR, 9, 4),
                            f_fmt(q.clgam, 11, 6),
                            f_fmt(q.cmgam, 11, 6)
                        );
                        crate::outln!(
                            " Qspec{:2} :  alpha ={}    CL ={}    CM ={}",
                            kqsp,
                            f_fmt(q.alqsp[kqsp] / DTOR, 9, 4),
                            f_fmt(q.clqsp[kqsp], 11, 6),
                            f_fmt(q.cmqsp[kqsp], 11, 6)
                        );
                    }
                    "SMOO" => {
                        let xf = &mut s.xf;
                        let q = &mut s.qdes_state;
                        let m = &mut s.mdes_state;
                        let kqsp = q.kqtarg;
                        let (iq1, iq2) = (q.iq1, q.iq2);
                        q.smooq(iq1, iq2, kqsp);
                        let qs = q.qspec[kqsp].clone();
                        m.cncalc(&qs, q.lqsym);
                        qspcir(m, q, xf);
                        print_qspec_forces(m, q, xf);
                    }
                    "FILT" | "F" => {
                        let xf = &mut s.xf;
                        let q = &mut s.qdes_state;
                        let m = &mut s.mdes_state;
                        let cfilt = 0.2;
                        m.cnfilt(cfilt);
                        m.piqsum();
                        qspcir(m, q, xf);
                        print_qspec_forces(m, q, xf);
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
                    "HARD" | "BLOW" | "B" | "RESE" | "R" => {
                        lcnpl = false;
                    }
                    "PLOT" | "P" => {
                        lcnpl = false;
                    }
                    "SPEC" => {
                        lcnpl = true;
                    }
                    "WIND" | "W" => {
                        crate::outln!(" ");
                        crate::outln!(" Type I,O,P to In,Out,Pan with cursor...");
                        crate::outln!(" (cursor input not available)");
                    }
                    "SIZE" => {
                        size_cmd(s, &r, "Current plot-object size =", "Enter new plot-object size^");
                        lcnpl = false;
                    }
                    "ANNO" => {
                        crate::outln!(" No active plot to annotate");
                    }
                    "DUMP" => {
                        let mut fname = comarg.trim().to_string();
                        if fname.is_empty() {
                            fname = s.input.ask("Enter Cn output filename^").trim().to_string();
                        }
                        let txt = s.mdes_state.cndump();
                        if std::fs::write(&fname, txt).is_err() {
                            crate::outln!(" Cannot write file {}", fname);
                        }
                    }
                    "EXEC" | "X" => {
                        exec(s);
                        lcnpl = false;
                    }
                    "PERT" => {
                        pert(s);
                        if s.input.eof {
                            return;
                        }
                        {
                            let xf = &mut s.xf;
                            let q = &mut s.qdes_state;
                            let m = &s.mdes_state;
                            qspcir(m, q, xf);
                        }
                        comand = "EXEC".into();
                        comarg = String::new();
                        redo = true;
                    }
                    _ => {
                        crate::outln!(" Command {} not recognized.  Type a \" ? \" for list.", fixed(&comand, 4));
                        comand = "****".into();
                    }
                }
            }
            let _ = lcnpl;
        }
    }
}
