//! Linear-vorticity panel method (port of XFoil xpanel.f) plus the dense
//! LU routines from xsolve.f.

use crate::spline::{segspl, setexp};
use crate::xfoil::{Mat, PsiOut, XFoil, HOPI, PI, QOPI};
use rayon::prelude::*;

/// Fortran SIGN(1.0, x)
#[inline(always)]
fn sign1(x: f64) -> f64 {
    if x >= 0.0 {
        1.0
    } else {
        -1.0
    }
}

/// LU factorization with implicit (scaled) partial pivoting (LUDCMP).
/// `a` is row-major with leading dimension `a.ld`, 1-based, size n x n.
pub fn ludcmp(n: usize, a: &mut Mat, indx: &mut [usize]) {
    let ld = a.ld;
    let d = &mut a.d;
    let mut vv = vec![0.0; n + 1];
    for i in 1..=n {
        let mut aamax: f64 = 0.0;
        for j in 1..=n {
            aamax = aamax.max(d[i * ld + j].abs());
        }
        vv[i] = 1.0 / aamax;
    }
    // Right-looking (outer-product) elimination: identical pivot choices to the
    // Crout version in exact arithmetic, but contiguous, vectorizable updates.
    for j in 1..=n {
        // pivot search on column j (rows j..n), using already-updated values
        let mut aamax = 0.0;
        let mut imax = j;
        for i in j..=n {
            let dum = vv[i] * d[i * ld + j].abs();
            if dum >= aamax {
                imax = i;
                aamax = dum;
            }
        }
        if j != imax {
            for k in 1..=n {
                d.swap(imax * ld + k, j * ld + k);
            }
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if j != n {
            let piv = d[j * ld + j];
            let dum = 1.0 / piv;
            let (top, bot) = d.split_at_mut((j + 1) * ld);
            let rowj = &top[j * ld + j + 1..j * ld + n + 1];
            for i in j + 1..=n {
                let r = &mut bot[(i - j - 1) * ld..(i - j) * ld];
                let lij = r[j] * dum;
                r[j] = lij;
                let rr = &mut r[j + 1..n + 1];
                for (x, y) in rr.iter_mut().zip(rowj.iter()) {
                    *x -= lij * *y;
                }
            }
        }
    }
}

/// Back-substitution for LUDCMP factors (BAKSUB). b is 1-based.
pub fn baksub(n: usize, a: &Mat, indx: &[usize], b: &mut [f64]) {
    let ld = a.ld;
    let d = &a.d;
    let mut ii = 0usize;
    for i in 1..=n {
        let ll = indx[i];
        let mut sum = b[ll];
        b[ll] = b[i];
        if ii != 0 {
            let row = &d[i * ld..];
            for j in ii..i {
                sum -= row[j] * b[j];
            }
        } else if sum != 0.0 {
            ii = i;
        }
        b[i] = sum;
    }
    for i in (1..=n).rev() {
        let mut sum = b[i];
        let row = &d[i * ld..];
        for j in i + 1..=n {
            sum -= row[j] * b[j];
        }
        b[i] = sum / row[i];
    }
}

impl XFoil {
    /// Panel angles (APCALC).
    pub fn apcalc(&mut self) {
        let n = self.n;
        for i in 1..n {
            let sx = self.x[i + 1] - self.x[i];
            let sy = self.y[i + 1] - self.y[i];
            if sx == 0.0 && sy == 0.0 {
                self.apanel[i] = (-self.ny[i]).atan2(-self.nx[i]);
            } else {
                self.apanel[i] = sx.atan2(-sy);
            }
        }
        let i = n;
        let ip = 1;
        if self.sharp {
            self.apanel[i] = PI;
        } else {
            let sx = self.x[ip] - self.x[i];
            let sy = self.y[ip] - self.y[i];
            self.apanel[i] = (-sx).atan2(sy) + PI;
        }
    }

    /// Streamfunction & sensitivities at point (xi,yi) due to airfoil
    /// vorticity/sources and freestream (PSILIN).
    #[allow(clippy::too_many_arguments)]
    pub fn psilin(&self, io: usize, xi: f64, yi: f64, nxi: f64, nyi: f64, geolin: bool, siglin: bool, o: &mut PsiOut) {
        let n = self.n;
        let x = &self.x;
        let y = &self.y;
        let seps = (self.s[n] - self.s[1]) * 1.0e-5;
        let cosa = self.alfa.cos();
        let sina = self.alfa.sin();

        for jo in 1..=n {
            o.dzdg[jo] = 0.0;
            o.dzdn[jo] = 0.0;
            o.dqdg[jo] = 0.0;
            o.dzdm[jo] = 0.0;
            o.dqdm[jo] = 0.0;
        }
        o.z_qinf = 0.0;
        o.z_alfa = 0.0;
        o.z_qdof0 = 0.0;
        o.z_qdof1 = 0.0;
        o.z_qdof2 = 0.0;
        o.z_qdof3 = 0.0;
        let mut psi = 0.0;
        let mut psi_ni = 0.0;
        let mut qtan1 = 0.0;
        let mut qtan2 = 0.0;

        let (scs, sds) = if self.sharp { (1.0, 0.0) } else { (self.ante / self.dste, self.aste / self.dste) };
        let gamu1 = &self.gamu[1];
        let gamu2 = &self.gamu[2];
        let gam = &self.gam;
        let sig = &self.sig;
        let on_surface = io >= 1 && io <= n;

        for jo in 1..=n {
            let mut jp = jo + 1;
            let mut jm = jo.wrapping_sub(1);
            let mut jq = jp + 1;
            if jo == 1 {
                jm = jo;
            } else if jo == n - 1 {
                jq = jp;
            } else if jo == n {
                jp = 1;
                if (x[jo] - x[jp]).powi(2) + (y[jo] - y[jp]).powi(2) < seps * seps {
                    break; // GO TO 12: no TE panel
                }
            }
            let dso = ((x[jo] - x[jp]).powi(2) + (y[jo] - y[jp]).powi(2)).sqrt();
            if dso == 0.0 {
                continue;
            }
            let dsio = 1.0 / dso;
            let apan = self.apanel[jo];
            let rx1 = xi - x[jo];
            let ry1 = yi - y[jo];
            let rx2 = xi - x[jp];
            let ry2 = yi - y[jp];
            let sx = (x[jp] - x[jo]) * dsio;
            let sy = (y[jp] - y[jo]) * dsio;
            let x1 = sx * rx1 + sy * ry1;
            let x2 = sx * rx2 + sy * ry2;
            let yy = sx * ry1 - sy * rx1;
            let rs1 = rx1 * rx1 + ry1 * ry1;
            let rs2 = rx2 * rx2 + ry2 * ry2;
            let sgn = if on_surface { 1.0 } else { sign1(yy) };
            let (g1, t1) = if io != jo && rs1 > 0.0 {
                (rs1.ln(), (sgn * x1).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI)
            } else {
                (0.0, 0.0)
            };
            let (g2, t2) = if io != jp && rs2 > 0.0 {
                (rs2.ln(), (sgn * x2).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI)
            } else {
                (0.0, 0.0)
            };
            let x1i = sx * nxi + sy * nyi;
            let x2i = sx * nxi + sy * nyi;
            let yyi = sx * nyi - sy * nxi;

            let (mut x1o, mut x1p, mut x2o, mut x2p, mut yyo, mut yyp) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            if geolin {
                let nxo = self.nx[jo];
                let nyo = self.ny[jo];
                let nxp = self.nx[jp];
                let nyp = self.ny[jp];
                x1o = -((rx1 - x1 * sx) * nxo + (ry1 - x1 * sy) * nyo) * dsio - (sx * nxo + sy * nyo);
                x1p = ((rx1 - x1 * sx) * nxp + (ry1 - x1 * sy) * nyp) * dsio;
                x2o = -((rx2 - x2 * sx) * nxo + (ry2 - x2 * sy) * nyo) * dsio;
                x2p = ((rx2 - x2 * sx) * nxp + (ry2 - x2 * sy) * nyp) * dsio - (sx * nxp + sy * nyp);
                yyo = ((rx1 + x1 * sy) * nyo - (ry1 - x1 * sx) * nxo) * dsio - (sx * nyo - sy * nxo);
                yyp = -((rx1 - x1 * sy) * nyp - (ry1 + x1 * sx) * nxp) * dsio;
            }

            if jo == n {
                // TE panel (label 11)
                let psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan);
                let pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
                let psigx1 = -(t1 - apan);
                let psigx2 = t2 - apan;
                let psigyy = 0.5 * (g1 - g2);
                let pgamx1 = 0.5 * g1;
                let pgamx2 = -0.5 * g2;
                let pgamyy = t1 - t2;
                let psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
                let pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;

                let sigte1 = 0.5 * scs * (gamu1[jp] - gamu1[jo]);
                let sigte2 = 0.5 * scs * (gamu2[jp] - gamu2[jo]);
                let gamte1 = -0.5 * sds * (gamu1[jp] - gamu1[jo]);
                let gamte2 = -0.5 * sds * (gamu2[jp] - gamu2[jo]);
                let sigte = 0.5 * scs * (gam[jp] - gam[jo]);
                let gamte = -0.5 * sds * (gam[jp] - gam[jo]);

                psi += HOPI * (psig * sigte + pgam * gamte);
                o.dzdg[jo] += -HOPI * psig * scs * 0.5;
                o.dzdg[jp] += HOPI * psig * scs * 0.5;
                o.dzdg[jo] += HOPI * pgam * sds * 0.5;
                o.dzdg[jp] += -HOPI * pgam * sds * 0.5;
                psi_ni += HOPI * (psigni * sigte + pgamni * gamte);
                qtan1 += HOPI * (psigni * sigte1 + pgamni * gamte1);
                qtan2 += HOPI * (psigni * sigte2 + pgamni * gamte2);
                o.dqdg[jo] -= HOPI * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
                o.dqdg[jp] += HOPI * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
                if geolin {
                    o.dzdn[jo] += HOPI * (psigx1 * x1o + psigx2 * x2o + psigyy * yyo) * sigte
                        + HOPI * (pgamx1 * x1o + pgamx2 * x2o + pgamyy * yyo) * gamte;
                    o.dzdn[jp] += HOPI * (psigx1 * x1p + psigx2 * x2p + psigyy * yyp) * sigte
                        + HOPI * (pgamx1 * x1p + pgamx2 * x2p + pgamyy * yyp) * gamte;
                    o.z_qdof0 += HOPI * psig * 0.5 * (self.qf0[jp] - self.qf0[jo]) * scs
                        - HOPI * pgam * 0.5 * (self.qf0[jp] - self.qf0[jo]) * sds;
                    o.z_qdof1 += HOPI * psig * 0.5 * (self.qf1[jp] - self.qf1[jo]) * scs
                        - HOPI * pgam * 0.5 * (self.qf1[jp] - self.qf1[jo]) * sds;
                    o.z_qdof2 += HOPI * psig * 0.5 * (self.qf2[jp] - self.qf2[jo]) * scs
                        - HOPI * pgam * 0.5 * (self.qf2[jp] - self.qf2[jo]) * sds;
                    o.z_qdof3 += HOPI * psig * 0.5 * (self.qf3[jp] - self.qf3[jo]) * scs
                        - HOPI * pgam * 0.5 * (self.qf3[jp] - self.qf3[jo]) * sds;
                }
                break;
            }

            if siglin {
                let x0 = 0.5 * (x1 + x2);
                let rs0 = x0 * x0 + yy * yy;
                let g0 = rs0.ln();
                let t0 = (sgn * x0).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI;

                // 1-0 half panel
                let dxinv = 1.0 / (x1 - x0);
                let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
                let pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv;
                let psx1 = -(t1 - apan);
                let psx0 = t0 - apan;
                let psyy = 0.5 * (g1 - g0);
                let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
                let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
                let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;
                let dsm = ((x[jp] - x[jm]).powi(2) + (y[jp] - y[jm]).powi(2)).sqrt();
                let dsim = 1.0 / dsm;
                let ssum = (sig[jp] - sig[jo]) * dsio + (sig[jp] - sig[jm]) * dsim;
                let sdif = (sig[jp] - sig[jo]) * dsio - (sig[jp] - sig[jm]) * dsim;
                psi += QOPI * (psum * ssum + pdif * sdif);
                o.dzdm[jm] += QOPI * (-psum * dsim + pdif * dsim);
                o.dzdm[jo] += QOPI * (-psum * dsio - pdif * dsio);
                o.dzdm[jp] += QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim));
                let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
                let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
                psi_ni += QOPI * (psni * ssum + pdni * sdif);
                o.dqdm[jm] += QOPI * (-psni * dsim + pdni * dsim);
                o.dqdm[jo] += QOPI * (-psni * dsio - pdni * dsio);
                o.dqdm[jp] += QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim));

                // 0-2 half panel
                let dxinv = 1.0 / (x0 - x2);
                let psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
                let pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv;
                let psx0 = -(t0 - apan);
                let psx2 = t2 - apan;
                let psyy = 0.5 * (g0 - g2);
                let pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
                let pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
                let pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;
                let dsp = ((x[jq] - x[jo]).powi(2) + (y[jq] - y[jo]).powi(2)).sqrt();
                let dsip = 1.0 / dsp;
                let ssum = (sig[jq] - sig[jo]) * dsip + (sig[jp] - sig[jo]) * dsio;
                let sdif = (sig[jq] - sig[jo]) * dsip - (sig[jp] - sig[jo]) * dsio;
                psi += QOPI * (psum * ssum + pdif * sdif);
                o.dzdm[jo] += QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
                o.dzdm[jp] += QOPI * (psum * dsio - pdif * dsio);
                o.dzdm[jq] += QOPI * (psum * dsip + pdif * dsip);
                let psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
                let pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
                psi_ni += QOPI * (psni * ssum + pdni * sdif);
                o.dqdm[jo] += QOPI * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
                o.dqdm[jp] += QOPI * (psni * dsio - pdni * dsio);
                o.dqdm[jq] += QOPI * (psni * dsip + pdni * dsip);
            }

            // vortex panel contribution
            let dxinv = 1.0 / (x1 - x2);
            let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
            let psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;
            let psx1 = 0.5 * g1;
            let psx2 = -0.5 * g2;
            let psyy = t1 - t2;
            let pdx1 = ((x1 + x2) * psx1 + psis - x1 * g1 - psid) * dxinv;
            let pdx2 = ((x1 + x2) * psx2 + psis + x2 * g2 + psid) * dxinv;
            let pdyy = ((x1 + x2) * psyy - yy * (g1 - g2)) * dxinv;

            let gsum1 = gamu1[jp] + gamu1[jo];
            let gsum2 = gamu2[jp] + gamu2[jo];
            let gdif1 = gamu1[jp] - gamu1[jo];
            let gdif2 = gamu2[jp] - gamu2[jo];
            let gsum = gam[jp] + gam[jo];
            let gdif = gam[jp] - gam[jo];

            psi += QOPI * (psis * gsum + psid * gdif);
            o.dzdg[jo] += QOPI * (psis - psid);
            o.dzdg[jp] += QOPI * (psis + psid);
            let psni = psx1 * x1i + psx2 * x2i + psyy * yyi;
            let pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi;
            psi_ni += QOPI * (gsum * psni + gdif * pdni);
            qtan1 += QOPI * (gsum1 * psni + gdif1 * pdni);
            qtan2 += QOPI * (gsum2 * psni + gdif2 * pdni);
            o.dqdg[jo] += QOPI * (psni - pdni);
            o.dqdg[jp] += QOPI * (psni + pdni);

            if geolin {
                o.dzdn[jo] += QOPI * gsum * (psx1 * x1o + psx2 * x2o + psyy * yyo) + QOPI * gdif * (pdx1 * x1o + pdx2 * x2o + pdyy * yyo);
                o.dzdn[jp] += QOPI * gsum * (psx1 * x1p + psx2 * x2p + psyy * yyp) + QOPI * gdif * (pdx1 * x1p + pdx2 * x2p + pdyy * yyp);
                o.z_qdof0 += QOPI * ((psis - psid) * self.qf0[jo] + (psis + psid) * self.qf0[jp]);
                o.z_qdof1 += QOPI * ((psis - psid) * self.qf1[jo] + (psis + psid) * self.qf1[jp]);
                o.z_qdof2 += QOPI * ((psis - psid) * self.qf2[jo] + (psis + psid) * self.qf2[jp]);
                o.z_qdof3 += QOPI * ((psis - psid) * self.qf3[jo] + (psis + psid) * self.qf3[jp]);
            }
        }

        // freestream terms
        let qinf = self.qinf;
        psi += qinf * (cosa * yi - sina * xi);
        psi_ni += qinf * (cosa * nyi - sina * nxi);
        qtan1 += qinf * nyi;
        qtan2 -= qinf * nxi;
        o.z_qinf += cosa * yi - sina * xi;
        o.z_alfa -= qinf * (sina * yi + cosa * xi);

        o.psi = psi;
        o.psi_ni = psi_ni;
        o.qtan1 = qtan1;
        o.qtan2 = qtan2;
    }

    /// Streamfunction due to wake sources (PSWLIN).
    pub fn pswlin(&self, io: usize, xi: f64, yi: f64, nxi: f64, nyi: f64, o: &mut PsiOut) {
        let n = self.n;
        let nw = self.nw;
        let x = &self.x;
        let y = &self.y;
        let sig = &self.sig;
        for jo in n + 1..=n + nw {
            o.dzdm[jo] = 0.0;
            o.dqdm[jo] = 0.0;
        }
        let mut psi = 0.0;
        let mut psi_ni = 0.0;
        let in_wake = (n + 1..=n + nw).contains(&io);
        for jo in n + 1..n + nw {
            let jp = jo + 1;
            let mut jm = jo - 1;
            let mut jq = jp + 1;
            if jo == n + 1 {
                jm = jo;
            } else if jo == n + nw - 1 {
                jq = jp;
            }
            let dso = ((x[jo] - x[jp]).powi(2) + (y[jo] - y[jp]).powi(2)).sqrt();
            let dsio = 1.0 / dso;
            let apan = self.apanel[jo];
            let rx1 = xi - x[jo];
            let ry1 = yi - y[jo];
            let rx2 = xi - x[jp];
            let ry2 = yi - y[jp];
            let sx = (x[jp] - x[jo]) * dsio;
            let sy = (y[jp] - y[jo]) * dsio;
            let x1 = sx * rx1 + sy * ry1;
            let x2 = sx * rx2 + sy * ry2;
            let yy = sx * ry1 - sy * rx1;
            let rs1 = rx1 * rx1 + ry1 * ry1;
            let rs2 = rx2 * rx2 + ry2 * ry2;
            let sgn = if in_wake { 1.0 } else { sign1(yy) };
            let (g1, t1) = if io != jo && rs1 > 0.0 {
                (rs1.ln(), (sgn * x1).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * PI)
            } else {
                (0.0, 0.0)
            };
            let (g2, t2) = if io != jp && rs2 > 0.0 {
                (rs2.ln(), (sgn * x2).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * PI)
            } else {
                (0.0, 0.0)
            };
            let x1i = sx * nxi + sy * nyi;
            let x2i = sx * nxi + sy * nyi;
            let yyi = sx * nyi - sy * nxi;

            let x0 = 0.5 * (x1 + x2);
            let rs0 = x0 * x0 + yy * yy;
            let g0 = rs0.ln();
            let t0 = (sgn * x0).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * PI;

            let dxinv = 1.0 / (x1 - x0);
            let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
            let pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv;
            let psx1 = -(t1 - apan);
            let psx0 = t0 - apan;
            let psyy = 0.5 * (g1 - g0);
            let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
            let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
            let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;
            let dsm = ((x[jp] - x[jm]).powi(2) + (y[jp] - y[jm]).powi(2)).sqrt();
            let dsim = 1.0 / dsm;
            let ssum = (sig[jp] - sig[jo]) * dsio + (sig[jp] - sig[jm]) * dsim;
            let sdif = (sig[jp] - sig[jo]) * dsio - (sig[jp] - sig[jm]) * dsim;
            psi += QOPI * (psum * ssum + pdif * sdif);
            o.dzdm[jm] += QOPI * (-psum * dsim + pdif * dsim);
            o.dzdm[jo] += QOPI * (-psum * dsio - pdif * dsio);
            o.dzdm[jp] += QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim));
            let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
            let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
            psi_ni += QOPI * (psni * ssum + pdni * sdif);
            o.dqdm[jm] += QOPI * (-psni * dsim + pdni * dsim);
            o.dqdm[jo] += QOPI * (-psni * dsio - pdni * dsio);
            o.dqdm[jp] += QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim));

            let dxinv = 1.0 / (x0 - x2);
            let psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
            let pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv;
            let psx0 = -(t0 - apan);
            let psx2 = t2 - apan;
            let psyy = 0.5 * (g0 - g2);
            let pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
            let pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
            let pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;
            let dsp = ((x[jq] - x[jo]).powi(2) + (y[jq] - y[jo]).powi(2)).sqrt();
            let dsip = 1.0 / dsp;
            let ssum = (sig[jq] - sig[jo]) * dsip + (sig[jp] - sig[jo]) * dsio;
            let sdif = (sig[jq] - sig[jo]) * dsip - (sig[jp] - sig[jo]) * dsio;
            psi += QOPI * (psum * ssum + pdif * sdif);
            o.dzdm[jo] += QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
            o.dzdm[jp] += QOPI * (psum * dsio - pdif * dsio);
            o.dzdm[jq] += QOPI * (psum * dsip + pdif * dsip);
            let psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
            let pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
            psi_ni += QOPI * (psni * ssum + pdni * sdif);
            o.dqdm[jo] += QOPI * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
            o.dqdm[jp] += QOPI * (psni * dsio - pdni * dsio);
            o.dqdm[jq] += QOPI * (psni * dsip + pdni * dsip);
        }
        o.psi = psi;
        o.psi_ni = psi_ni;
    }

    fn ensure_bij(&mut self) {
        let ldb = self.n + 2;
        let nc = self.n + self.nw.max(crate::xfoil::IWX) + 1;
        if self.ldb != ldb || self.bijc.len() < ldb * (nc + 1) {
            self.ldb = ldb;
            self.bijc = vec![0.0; ldb * (nc + 1)];
        }
    }

    /// Unit vorticity distributions for alpha = 0, 90 (GGCALC).
    pub fn ggcalc(&mut self) {
        const BWT: f64 = 0.1;
        let n = self.n;
        if self.verbose {
            crate::outln!(" Calculating unit vorticity distributions ...");
        }
        for i in 1..=n {
            self.gam[i] = 0.0;
            self.gamu[1][i] = 0.0;
            self.gamu[2][i] = 0.0;
        }
        self.psio = 0.0;
        self.qw_valid = false;
        self.aij.ensure(n + 1, n + 1);
        self.ensure_bij();
        let ldb = self.ldb;
        let qinf = self.qinf;

        // rows of AIJ and BIJ in parallel
        let rows: Vec<(Vec<f64>, Vec<f64>)> = {
            let me = &*self;
            let f = |o: &mut PsiOut, i: usize| {
                me.psilin(i, me.x[i], me.y[i], me.nx[i], me.ny[i], false, true, o);
                (o.dzdg[1..=n].to_vec(), o.dzdm[1..=n].to_vec())
            };
            if me.par_inner {
                (1..=n).into_par_iter().map_init(PsiOut::new, f).collect()
            } else {
                let mut o = PsiOut::new();
                (1..=n).map(|i| f(&mut o, i)).collect()
            }
        };
        for (k, (dzdg, dzdm)) in rows.into_iter().enumerate() {
            let i = k + 1;
            let row = self.aij.row_mut(i);
            row[1..=n].copy_from_slice(&dzdg);
            row[n + 1] = -1.0;
            for j in 1..=n {
                self.bijc[j * ldb + i] = -dzdm[j - 1];
            }
            self.gamu[1][i] = -(qinf * self.y[i]);
            self.gamu[2][i] = -(-qinf * self.x[i]);
        }
        // Kutta condition
        let res = 0.0;
        for j in 1..=n + 1 {
            self.aij.set(n + 1, j, 0.0);
        }
        self.aij.set(n + 1, 1, 1.0);
        self.aij.set(n + 1, n, 1.0);
        self.gamu[1][n + 1] = -res;
        self.gamu[2][n + 1] = -res;
        for j in 1..=n {
            self.bijc[j * ldb + n + 1] = 0.0;
        }

        if self.sharp {
            let ag1 = (-self.yp[1]).atan2(-self.xp[1]);
            let ag2 = crate::geom::atanc(self.yp[n], self.xp[n], ag1);
            let abis = 0.5 * (ag1 + ag2);
            let cbis = abis.cos();
            let sbis = abis.sin();
            let ds1 = ((self.x[1] - self.x[2]).powi(2) + (self.y[1] - self.y[2]).powi(2)).sqrt();
            let ds2 = ((self.x[n] - self.x[n - 1]).powi(2) + (self.y[n] - self.y[n - 1]).powi(2)).sqrt();
            let dsmin = ds1.min(ds2);
            let xbis = self.xte - BWT * dsmin * cbis;
            let ybis = self.yte - BWT * dsmin * sbis;
            let mut o = std::mem::take(&mut self.pso);
            self.psilin(0, xbis, ybis, -sbis, cbis, false, true, &mut o);
            for j in 1..=n {
                self.aij.set(n, j, o.dqdg[j]);
                self.bijc[j * ldb + n] = -o.dqdm[j];
            }
            self.pso = o;
            self.aij.set(n, n + 1, 0.0);
            self.gamu[1][n] = -cbis;
            self.gamu[2][n] = -sbis;
        }

        let mut piv = std::mem::take(&mut self.aijpiv);
        ludcmp(n + 1, &mut self.aij, &mut piv);
        self.aijpiv = piv;
        self.lqaij = true;
        let mut g1 = std::mem::take(&mut self.gamu[1]);
        let mut g2 = std::mem::take(&mut self.gamu[2]);
        baksub(n + 1, &self.aij, &self.aijpiv, &mut g1);
        baksub(n + 1, &self.aij, &self.aijpiv, &mut g2);
        self.gamu[1] = g1;
        self.gamu[2] = g2;
        for i in 1..=n {
            self.qinvu[1][i] = self.gamu[1][i];
            self.qinvu[2][i] = self.gamu[2][i];
        }
        self.lgamu = true;
    }

    /// Inviscid wake tangential velocities for alpha = 0, 90 (QWCALC).
    pub fn qwcalc(&mut self) {
        let n = self.n;
        self.qinvu[1][n + 1] = self.qinvu[1][n];
        self.qinvu[2][n + 1] = self.qinvu[2][n];
        if self.qw_valid {
            return; // already computed by XYWAKE's fused pass
        }
        let mut o = std::mem::take(&mut self.pso);
        for i in n + 2..=n + self.nw {
            self.psilin(i, self.x[i], self.y[i], self.nx[i], self.ny[i], false, false, &mut o);
            self.qinvu[1][i] = o.qtan1;
            self.qinvu[2][i] = o.qtan2;
        }
        self.pso = o;
        self.qw_valid = true;
    }

    /// Source influence matrix DIJ (QDCALC).
    pub fn qdcalc(&mut self) {
        let n = self.n;
        let nw = self.nw;
        let nt = n + nw;
        if self.verbose {
            crate::outln!(" Calculating source influence matrix ...");
        }
        self.ensure_bij();
        let ldb = self.ldb;
        if self.dij.nr != nt {
            self.dij = Mat::new(nt, nt);
            self.ladij = false;
        }
        self.cij.ensure(nw, n);

        if !self.ladij {
            // dGam/dSig for airfoil sources: back-substitute each column (parallel)
            let aij = &self.aij;
            let piv = &self.aijpiv;
            if self.par_inner {
                self.bijc[ldb..(n + 1) * ldb].par_chunks_mut(ldb).for_each(|col| {
                    baksub(n + 1, aij, piv, col);
                });
            } else {
                self.bijc[ldb..(n + 1) * ldb].chunks_mut(ldb).for_each(|col| {
                    baksub(n + 1, aij, piv, col);
                });
            }
            for j in 1..=n {
                for i in 1..=n {
                    let v = self.bijc[j * ldb + i];
                    self.dij.set(i, j, v);
                }
            }
            self.ladij = true;
        }

        // dPsi/dm for wake sources on airfoil surface
        {
            let me = &*self;
            let f = |o: &mut PsiOut, i: usize| {
                me.pswlin(i, me.x[i], me.y[i], me.nx[i], me.ny[i], o);
                o.dzdm[n + 1..=n + nw].to_vec()
            };
            let rows: Vec<Vec<f64>> = if me.par_inner {
                (1..=n).into_par_iter().map_init(PsiOut::new, f).collect()
            } else {
                let mut o = PsiOut::new();
                (1..=n).map(|i| f(&mut o, i)).collect()
            };
            for (k, r) in rows.into_iter().enumerate() {
                let i = k + 1;
                for (jj, val) in r.into_iter().enumerate() {
                    let j = n + 1 + jj;
                    self.bijc[j * ldb + i] = -val;
                }
            }
        }
        for j in n + 1..=nt {
            self.bijc[j * ldb + n + 1] = 0.0;
        }
        if self.sharp {
            for j in n + 1..=nt {
                self.bijc[j * ldb + n] = 0.0;
            }
        }
        {
            let aij = &self.aij;
            let piv = &self.aijpiv;
            if self.par_inner {
                self.bijc[(n + 1) * ldb..(nt + 1) * ldb].par_chunks_mut(ldb).for_each(|col| {
                    baksub(n + 1, aij, piv, col);
                });
            } else {
                self.bijc[(n + 1) * ldb..(nt + 1) * ldb].chunks_mut(ldb).for_each(|col| {
                    baksub(n + 1, aij, piv, col);
                });
            }
        }
        for i in 1..=n {
            for j in n + 1..=nt {
                let v = self.bijc[j * ldb + i];
                self.dij.set(i, j, v);
            }
        }

        // influence of sources on wake velocities
        {
            let me = &*self;
            let f = |o: &mut PsiOut, i: usize| {
                    me.psilin(i, me.x[i], me.y[i], me.nx[i], me.ny[i], false, true, o);
                    let cij = o.dqdg[1..=n].to_vec();
                    let dm1 = o.dqdm[1..=n].to_vec();
                    me.pswlin(i, me.x[i], me.y[i], me.nx[i], me.ny[i], o);
                    let dm2 = o.dqdm[n + 1..=nt].to_vec();
                    (cij, dm1, dm2)
            };
            let rows: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = if me.par_inner {
                (n + 1..=nt).into_par_iter().map_init(PsiOut::new, f).collect()
            } else {
                let mut o = PsiOut::new();
                (n + 1..=nt).map(|i| f(&mut o, i)).collect()
            };
            for (k, (cij, dm1, dm2)) in rows.into_iter().enumerate() {
                let i = n + 1 + k;
                let iw = i - n;
                self.cij.row_mut(iw)[1..=n].copy_from_slice(&cij);
                let r = self.dij.row_mut(i);
                r[1..=n].copy_from_slice(&dm1);
                r[n + 1..=nt].copy_from_slice(&dm2);
            }
        }

        // add effect of all sources on airfoil vorticity which affects wake Qtan
        {
            let cij = &self.cij;
            let bijc = &self.bijc;
            let (top, bot) = self.dij.d.split_at_mut((n + 1) * self.dij.ld);
            let ld = self.dij.ld;
            let work = |(k, row): (usize, &mut [f64])| {
                let iw = k + 1;
                let c = &cij.row(iw)[1..=n];
                // airfoil surface sources: sum_k CIJ(iw,k) * DIJ(k,j)
                let mut acc = vec![0.0; n + 1];
                for (kk, &ck) in c.iter().enumerate() {
                    let drow = &top[(kk + 1) * ld + 1..(kk + 1) * ld + n + 1];
                    for (a, &d) in acc[1..].iter_mut().zip(drow.iter()) {
                        *a += ck * d;
                    }
                }
                for j in 1..=n {
                    row[j] += acc[j];
                }
                // wake sources: sum_k CIJ(iw,k) * BIJ(k,j)
                for j in n + 1..=nt {
                    let col = &bijc[j * ldb + 1..j * ldb + n + 1];
                    let mut sum = 0.0;
                    for (&ck, &b) in c.iter().zip(col.iter()) {
                        sum += ck * b;
                    }
                    row[j] += sum;
                }
            };
            if self.par_inner {
                bot[..nw * ld].par_chunks_mut(ld).enumerate().for_each(work);
            } else {
                bot[..nw * ld].chunks_mut(ld).enumerate().for_each(work);
            }
        }
        // first wake point has same velocity as TE
        for j in 1..=nt {
            let v = self.dij.at(n, j);
            self.dij.set(n + 1, j, v);
        }
        self.lwdij = true;
    }

    /// Fused PSILIN (GEOLIN=SIGLIN=false) for several normals at one point:
    /// returns dPsi/dn (with current GAM) and the alpha=0,90 tangential
    /// velocities QTAN1/QTAN2 for each normal.  Per-output operation order is
    /// identical to `psilin`, so results are bit-identical.
    pub fn psilin_vel<const K: usize>(&self, io: usize, xi: f64, yi: f64, nrm: [(f64, f64); K]) -> ([f64; K], [f64; K], [f64; K]) {
        let n = self.n;
        let x = &self.x;
        let y = &self.y;
        let seps = (self.s[n] - self.s[1]) * 1.0e-5;
        let cosa = self.alfa.cos();
        let sina = self.alfa.sin();
        let mut psi_ni = [0.0f64; K];
        let mut qtan1 = [0.0f64; K];
        let mut qtan2 = [0.0f64; K];
        let (scs, sds) = if self.sharp { (1.0, 0.0) } else { (self.ante / self.dste, self.aste / self.dste) };
        let gamu1 = &self.gamu[1];
        let gamu2 = &self.gamu[2];
        let gam = &self.gam;
        let on_surface = io >= 1 && io <= n;
        for jo in 1..=n {
            let mut jp = jo + 1;
            if jo == n {
                jp = 1;
                if (x[jo] - x[jp]).powi(2) + (y[jo] - y[jp]).powi(2) < seps * seps {
                    break;
                }
            }
            let dso = ((x[jo] - x[jp]).powi(2) + (y[jo] - y[jp]).powi(2)).sqrt();
            if dso == 0.0 {
                continue;
            }
            let dsio = 1.0 / dso;
            let apan = self.apanel[jo];
            let rx1 = xi - x[jo];
            let ry1 = yi - y[jo];
            let rx2 = xi - x[jp];
            let ry2 = yi - y[jp];
            let sx = (x[jp] - x[jo]) * dsio;
            let sy = (y[jp] - y[jo]) * dsio;
            let x1 = sx * rx1 + sy * ry1;
            let x2 = sx * rx2 + sy * ry2;
            let yy = sx * ry1 - sy * rx1;
            let rs1 = rx1 * rx1 + ry1 * ry1;
            let rs2 = rx2 * rx2 + ry2 * ry2;
            let sgn = if on_surface { 1.0 } else { sign1(yy) };
            let (g1, t1) = if io != jo && rs1 > 0.0 {
                (rs1.ln(), (sgn * x1).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI)
            } else {
                (0.0, 0.0)
            };
            let (g2, t2) = if io != jp && rs2 > 0.0 {
                (rs2.ln(), (sgn * x2).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI)
            } else {
                (0.0, 0.0)
            };
            if jo == n {
                let psigx1 = -(t1 - apan);
                let psigx2 = t2 - apan;
                let psigyy = 0.5 * (g1 - g2);
                let pgamx1 = 0.5 * g1;
                let pgamx2 = -0.5 * g2;
                let pgamyy = t1 - t2;
                let sigte1 = 0.5 * scs * (gamu1[jp] - gamu1[jo]);
                let sigte2 = 0.5 * scs * (gamu2[jp] - gamu2[jo]);
                let gamte1 = -0.5 * sds * (gamu1[jp] - gamu1[jo]);
                let gamte2 = -0.5 * sds * (gamu2[jp] - gamu2[jo]);
                let sigte = 0.5 * scs * (gam[jp] - gam[jo]);
                let gamte = -0.5 * sds * (gam[jp] - gam[jo]);
                for k in 0..K {
                    let (nxi, nyi) = nrm[k];
                    let x1i = sx * nxi + sy * nyi;
                    let x2i = sx * nxi + sy * nyi;
                    let yyi = sx * nyi - sy * nxi;
                    let psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
                    let pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;
                    psi_ni[k] += HOPI * (psigni * sigte + pgamni * gamte);
                    qtan1[k] += HOPI * (psigni * sigte1 + pgamni * gamte1);
                    qtan2[k] += HOPI * (psigni * sigte2 + pgamni * gamte2);
                }
                break;
            }
            let dxinv = 1.0 / (x1 - x2);
            let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
            let psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;
            let psx1 = 0.5 * g1;
            let psx2 = -0.5 * g2;
            let psyy = t1 - t2;
            let pdx1 = ((x1 + x2) * psx1 + psis - x1 * g1 - psid) * dxinv;
            let pdx2 = ((x1 + x2) * psx2 + psis + x2 * g2 + psid) * dxinv;
            let pdyy = ((x1 + x2) * psyy - yy * (g1 - g2)) * dxinv;
            let gsum1 = gamu1[jp] + gamu1[jo];
            let gsum2 = gamu2[jp] + gamu2[jo];
            let gdif1 = gamu1[jp] - gamu1[jo];
            let gdif2 = gamu2[jp] - gamu2[jo];
            let gsum = gam[jp] + gam[jo];
            let gdif = gam[jp] - gam[jo];
            for k in 0..K {
                let (nxi, nyi) = nrm[k];
                let x1i = sx * nxi + sy * nyi;
                let x2i = sx * nxi + sy * nyi;
                let yyi = sx * nyi - sy * nxi;
                let psni = psx1 * x1i + psx2 * x2i + psyy * yyi;
                let pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi;
                psi_ni[k] += QOPI * (gsum * psni + gdif * pdni);
                qtan1[k] += QOPI * (gsum1 * psni + gdif1 * pdni);
                qtan2[k] += QOPI * (gsum2 * psni + gdif2 * pdni);
            }
        }
        let qinf = self.qinf;
        for k in 0..K {
            let (nxi, nyi) = nrm[k];
            psi_ni[k] += qinf * (cosa * nyi - sina * nxi);
            qtan1[k] += qinf * nyi;
            qtan2[k] -= qinf * nxi;
        }
        (psi_ni, qtan1, qtan2)
    }

    /// Wake trajectory (XYWAKE).  Also fills the wake QINVU values (what
    /// QWCALC would compute) in the same influence pass.
    pub fn xywake(&mut self) {
        let n = self.n;
        if self.verbose {
            crate::outln!(" Calculating wake trajectory ...");
        }
        self.nw = n / 8 + 2;
        if self.nw > crate::xfoil::IWX {
            self.nw = crate::xfoil::IWX;
        }
        let nw = self.nw;
        let ds1 = 0.5 * (self.s[2] - self.s[1] + self.s[n] - self.s[n - 1]);
        let mut snew = vec![0.0; nw + 1];
        setexp(&mut snew, ds1, self.waklen * self.chord, nw);
        self.xte = 0.5 * (self.x[1] + self.x[n]);
        self.yte = 0.5 * (self.y[1] + self.y[n]);

        let i = n + 1;
        let sx = 0.5 * (self.yp[n] - self.yp[1]);
        let sy = 0.5 * (self.xp[1] - self.xp[n]);
        let smod = (sx * sx + sy * sy).sqrt();
        self.nx[i] = sx / smod;
        self.ny[i] = sy / smod;
        self.x[i] = self.xte - 0.0001 * self.ny[i];
        self.y[i] = self.yte + 0.0001 * self.nx[i];
        self.s[i] = self.s[n];

        let set_normal = |me: &mut XFoil, i: usize, psi_x: f64, psi_y: f64| {
            let m = (psi_x * psi_x + psi_y * psi_y).sqrt();
            me.nx[i + 1] = -psi_x / m;
            me.ny[i + 1] = -psi_y / m;
            me.apanel[i] = psi_y.atan2(psi_x);
        };
        let (p, _, _) = self.psilin_vel(i, self.x[i], self.y[i], [(1.0, 0.0), (0.0, 1.0)]);
        set_normal(self, i, p[0], p[1]);
        for i in n + 2..=n + nw {
            let ds = snew[i - n] - snew[i - n - 1];
            self.x[i] = self.x[i - 1] - ds * self.ny[i];
            self.y[i] = self.y[i - 1] + ds * self.nx[i];
            self.s[i] = self.s[i - 1] + ds;
            if i == n + nw {
                let (_, q1, q2) = self.psilin_vel(i, self.x[i], self.y[i], [(self.nx[i], self.ny[i])]);
                self.qinvu[1][i] = q1[0];
                self.qinvu[2][i] = q2[0];
                continue;
            }
            let (p, q1, q2) = self.psilin_vel(i, self.x[i], self.y[i], [(1.0, 0.0), (0.0, 1.0), (self.nx[i], self.ny[i])]);
            self.qinvu[1][i] = q1[2];
            self.qinvu[2][i] = q2[2];
            set_normal(self, i, p[0], p[1]);
        }
        self.qw_valid = true;
        self.lwake = true;
        self.awake = self.alfa;
        self.lwdij = false;
    }

    /// Stagnation point (STFIND).
    pub fn stfind(&mut self) {
        let n = self.n;
        let mut found = None;
        for i in 1..n {
            if self.gam[i] >= 0.0 && self.gam[i + 1] < 0.0 {
                found = Some(i);
                break;
            }
        }
        let i = match found {
            Some(i) => i,
            None => {
                if self.verbose {
                    crate::outln!(" STFIND: Stagnation point not found. Continuing ...");
                }
                n / 2
            }
        };
        self.ist = i;
        let dgam = self.gam[i + 1] - self.gam[i];
        let ds = self.s[i + 1] - self.s[i];
        if self.gam[i] < -self.gam[i + 1] {
            self.sst = self.s[i] - ds * (self.gam[i] / dgam);
        } else {
            self.sst = self.s[i + 1] - ds * (self.gam[i + 1] / dgam);
        }
        if self.sst <= self.s[i] {
            self.sst = self.s[i] + 1.0e-7;
        }
        if self.sst >= self.s[i + 1] {
            self.sst = self.s[i + 1] - 1.0e-7;
        }
        self.sst_go = (self.sst - self.s[i + 1]) / dgam;
        self.sst_gp = (self.s[i] - self.sst) / dgam;
    }

    /// BL station -> panel pointers (IBLPAN).
    pub fn iblpan(&mut self) {
        let n = self.n;
        let nw = self.nw;
        let mut ibl = 1;
        for i in (1..=self.ist).rev() {
            ibl += 1;
            self.ipan[1][ibl] = i;
            self.vti[1][ibl] = 1.0;
        }
        self.iblte[1] = ibl;
        self.nbl[1] = ibl;
        let mut ibl = 1;
        for i in self.ist + 1..=n {
            ibl += 1;
            self.ipan[2][ibl] = i;
            self.vti[2][ibl] = -1.0;
        }
        self.iblte[2] = ibl;
        for iw in 1..=nw {
            let i = n + iw;
            let ibl = self.iblte[2] + iw;
            self.ipan[2][ibl] = i;
            self.vti[2][ibl] = -1.0;
        }
        self.nbl[2] = self.iblte[2] + nw;
        for iw in 1..=nw {
            self.ipan[1][self.iblte[1] + iw] = self.ipan[2][self.iblte[2] + iw];
            self.vti[1][self.iblte[1] + iw] = 1.0;
        }
        let iblmax = self.iblte[1].max(self.iblte[2]) + nw;
        if iblmax > crate::xfoil::IVX {
            panic!("BL array overflow. Increase IVX to at least {}", iblmax);
        }
        self.lipan = true;
    }

    /// BL arc length arrays (XICALC).
    pub fn xicalc(&mut self) {
        const XFEPS: f64 = 1.0e-7;
        let n = self.n;
        let xeps = XFEPS * (self.s[n] - self.s[1]);
        self.xssi[1][1] = 0.0;
        for ibl in 2..=self.iblte[1] {
            let i = self.ipan[1][ibl];
            self.xssi[1][ibl] = (self.sst - self.s[i]).max(xeps);
        }
        self.xssi[2][1] = 0.0;
        for ibl in 2..=self.iblte[2] {
            let i = self.ipan[2][ibl];
            self.xssi[2][ibl] = (self.s[i] - self.sst).max(xeps);
        }
        let ibl1 = self.iblte[1] + 1;
        self.xssi[1][ibl1] = self.xssi[1][ibl1 - 1];
        let ibl2 = self.iblte[2] + 1;
        self.xssi[2][ibl2] = self.xssi[2][ibl2 - 1];
        for ibl in self.iblte[2] + 2..=self.nbl[2] {
            let i = self.ipan[2][ibl];
            let dxssi = ((self.x[i] - self.x[i - 1]).powi(2) + (self.y[i] - self.y[i - 1]).powi(2)).sqrt();
            let ibl1 = self.iblte[1] + ibl - self.iblte[2];
            let ibl2 = self.iblte[2] + ibl - self.iblte[2];
            self.xssi[1][ibl1] = self.xssi[1][ibl1 - 1] + dxssi;
            self.xssi[2][ibl2] = self.xssi[2][ibl2 - 1] + dxssi;
        }
        let telrat = 2.50;
        let crosp = (self.xp[1] * self.yp[n] - self.yp[1] * self.xp[n])
            / ((self.xp[1].powi(2) + self.yp[1].powi(2)) * (self.xp[n].powi(2) + self.yp[n].powi(2))).sqrt();
        let mut dwdxte = crosp / (1.0 - crosp * crosp).sqrt();
        dwdxte = dwdxte.max(-3.0 / telrat).min(3.0 / telrat);
        let aa = 3.0 + telrat * dwdxte;
        let bb = -2.0 - telrat * dwdxte;
        if self.sharp {
            for iw in 1..=self.nw {
                self.wgap[iw] = 0.0;
            }
        } else {
            for iw in 1..=self.nw {
                let ibl = self.iblte[2] + iw;
                let zn = 1.0 - (self.xssi[2][ibl] - self.xssi[2][self.iblte[2]]) / (telrat * self.ante);
                self.wgap[iw] = 0.0;
                if zn >= 0.0 {
                    self.wgap[iw] = self.ante * (aa + bb * zn) * zn * zn;
                }
            }
        }
    }

    /// Inviscid Ue from panel QINV (UICALC).
    pub fn uicalc(&mut self) {
        for is in 1..=2 {
            self.uinv[is][1] = 0.0;
            self.uinv_a[is][1] = 0.0;
            for ibl in 2..=self.nbl[is] {
                let i = self.ipan[is][ibl];
                self.uinv[is][ibl] = self.vti[is][ibl] * self.qinv[i];
                self.uinv_a[is][ibl] = self.vti[is][ibl] * self.qinv_a[i];
            }
        }
    }

    /// Viscous Ue from QVIS (UECALC).
    pub fn uecalc(&mut self) {
        for is in 1..=2 {
            self.uedg[is][1] = 0.0;
            for ibl in 2..=self.nbl[is] {
                let i = self.ipan[is][ibl];
                self.uedg[is][ibl] = self.vti[is][ibl] * self.qvis[i];
            }
        }
    }

    /// QVIS from UEDG (QVFUE).
    pub fn qvfue(&mut self) {
        for is in 1..=2 {
            for ibl in 2..=self.nbl[is] {
                let i = self.ipan[is][ibl];
                self.qvis[i] = self.vti[is][ibl] * self.uedg[is][ibl];
            }
        }
    }

    /// Inviscid panel tangential velocity for current alpha (QISET).
    pub fn qiset(&mut self) {
        let cosa = self.alfa.cos();
        let sina = self.alfa.sin();
        for i in 1..=self.n + self.nw {
            self.qinv[i] = cosa * self.qinvu[1][i] + sina * self.qinvu[2][i];
            self.qinv_a[i] = -sina * self.qinvu[1][i] + cosa * self.qinvu[2][i];
        }
    }

    /// GAM from QVIS (GAMQV).
    pub fn gamqv(&mut self) {
        for i in 1..=self.n {
            self.gam[i] = self.qvis[i];
            self.gam_a[i] = self.qinv_a[i];
        }
    }

    /// Move stagnation point (STMOVE).
    pub fn stmove(&mut self) {
        let istold = self.ist;
        self.stfind();
        if istold == self.ist {
            self.xicalc();
        } else {
            self.iblpan();
            self.uicalc();
            self.xicalc();
            self.iblsys();
            let ist = self.ist;
            if ist > istold {
                let idif = ist - istold;
                self.itran[1] += idif;
                self.itran[2] -= idif;
                for ibl in (idif + 2..=self.nbl[1]).rev() {
                    self.ctau[1][ibl] = self.ctau[1][ibl - idif];
                    self.thet[1][ibl] = self.thet[1][ibl - idif];
                    self.dstr[1][ibl] = self.dstr[1][ibl - idif];
                    self.uedg[1][ibl] = self.uedg[1][ibl - idif];
                }
                let dudx = self.uedg[1][idif + 2] / self.xssi[1][idif + 2];
                for ibl in (2..=idif + 1).rev() {
                    self.ctau[1][ibl] = self.ctau[1][idif + 2];
                    self.thet[1][ibl] = self.thet[1][idif + 2];
                    self.dstr[1][ibl] = self.dstr[1][idif + 2];
                    self.uedg[1][ibl] = dudx * self.xssi[1][ibl];
                }
                for ibl in 2..=self.nbl[2] {
                    self.ctau[2][ibl] = self.ctau[2][ibl + idif];
                    self.thet[2][ibl] = self.thet[2][ibl + idif];
                    self.dstr[2][ibl] = self.dstr[2][ibl + idif];
                    self.uedg[2][ibl] = self.uedg[2][ibl + idif];
                }
            } else {
                let idif = istold - ist;
                self.itran[1] -= idif;
                self.itran[2] += idif;
                for ibl in (idif + 2..=self.nbl[2]).rev() {
                    self.ctau[2][ibl] = self.ctau[2][ibl - idif];
                    self.thet[2][ibl] = self.thet[2][ibl - idif];
                    self.dstr[2][ibl] = self.dstr[2][ibl - idif];
                    self.uedg[2][ibl] = self.uedg[2][ibl - idif];
                }
                let dudx = self.uedg[2][idif + 2] / self.xssi[2][idif + 2];
                for ibl in (2..=idif + 1).rev() {
                    self.ctau[2][ibl] = self.ctau[2][idif + 2];
                    self.thet[2][ibl] = self.thet[2][idif + 2];
                    self.dstr[2][ibl] = self.dstr[2][idif + 2];
                    self.uedg[2][ibl] = dudx * self.xssi[2][ibl];
                }
                for ibl in 2..=self.nbl[1] {
                    self.ctau[1][ibl] = self.ctau[1][ibl + idif];
                    self.thet[1][ibl] = self.thet[1][ibl + idif];
                    self.dstr[1][ibl] = self.dstr[1][ibl + idif];
                    self.uedg[1][ibl] = self.uedg[1][ibl + idif];
                }
            }
            let ueps = 1.0e-7;
            for is in 1..=2 {
                for ibl in 2..=self.nbl[is] {
                    let i = self.ipan[is][ibl];
                    if self.uedg[is][ibl] <= ueps {
                        self.uedg[is][ibl] = ueps;
                        self.qvis[i] = self.vti[is][ibl] * ueps;
                        self.gam[i] = self.vti[is][ibl] * ueps;
                    }
                }
            }
        }
        for is in 1..=2 {
            for ibl in 2..=self.nbl[is] {
                self.mass[is][ibl] = self.dstr[is][ibl] * self.uedg[is][ibl];
            }
        }
    }

    /// Ue from inviscid Ue plus source influence (UESET).
    pub fn ueset(&mut self) {
        for is in 1..=2 {
            for ibl in 2..=self.nbl[is] {
                let i = self.ipan[is][ibl];
                let drow = self.dij.row(i);
                let mut dui = 0.0;
                for js in 1..=2 {
                    for jbl in 2..=self.nbl[js] {
                        let j = self.ipan[js][jbl];
                        let ue_m = -self.vti[is][ibl] * self.vti[js][jbl] * drow[j];
                        dui += ue_m * self.mass[js][jbl];
                    }
                }
                self.uedg[is][ibl] = self.uinv[is][ibl] + dui;
            }
        }
    }

    /// DSTR from MASS/UEDG (DSSET).
    pub fn dsset(&mut self) {
        for is in 1..=2 {
            for ibl in 2..=self.nbl[is] {
                self.dstr[is][ibl] = self.mass[is][ibl] / self.uedg[is][ibl];
            }
        }
    }
}

/// Normal unit vectors at panel nodes (NCALC).
pub fn ncalc(x: &[f64], y: &[f64], s: &[f64], n: usize, xn: &mut [f64], yn: &mut [f64]) {
    if n <= 1 {
        return;
    }
    segspl(x, xn, s, n);
    segspl(y, yn, s, n);
    for i in 1..=n {
        let sx = yn[i];
        let sy = -xn[i];
        let smod = (sx * sx + sy * sy).sqrt();
        xn[i] = sx / smod;
        yn[i] = sy / smod;
    }
    for i in 1..n {
        if s[i] == s[i + 1] {
            let sx = 0.5 * (xn[i] + xn[i + 1]);
            let sy = 0.5 * (yn[i] + yn[i + 1]);
            let smod = (sx * sx + sy * sy).sqrt();
            xn[i] = sx / smod;
            yn[i] = sy / smod;
            xn[i + 1] = sx / smod;
            yn[i + 1] = sy / smod;
        }
    }
}
