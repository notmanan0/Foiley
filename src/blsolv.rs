//! Cache-blocked, bit-exact reformulation of XFoil's BLSOLV.
//!
//! The original solver eliminates one 3x3 block row per step `iv` and then
//! sweeps *every* remaining row once per step (an O(nsys^3) pass that is
//! L2-bandwidth bound).  Here the steps are grouped into panels of `PB`
//! pivots.  Rows inside the panel (plus the row after it and the wake row
//! coupled via VZ) are processed with the original algorithm; every other
//! ("trailing") row only ever receives "lower VM column" updates from the
//! panel pivots, so all of a panel's updates to a trailing row are applied in
//! one pass while the row is hot in cache.  The per-element sequence of
//! floating point operations is identical to the Fortran, so results are
//! bit-for-bit the same.  Trailing rows are independent and are processed in
//! parallel.
//!
//! The right-hand sides VDEL(k,1:2,iv) are carried as two extra matrix
//! columns (nsys+1, nsys+2): they undergo exactly the same row operations.

use rayon::prelude::*;

const PB: usize = 16;

#[inline(always)]
unsafe fn row(vm: *mut f64, ldm: usize, k: usize, iv: usize) -> *mut f64 {
    vm.add((iv * 3 + k) * ldm)
}

/// One original BLSOLV step `iv` restricted to the given set of "lower" rows.
/// Column `iv` is processed first (it yields the step's scalars); then all of
/// the step's row operations are applied in a single fused pass over columns
/// iv+1..=end.  Per element the operation sequence equals the original.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
unsafe fn step_inside(
    vm: *mut f64,
    ldm: usize,
    va: &mut [[[f64; 2]; 3]],
    vb: &[[[f64; 2]; 3]],
    vz: &[[f64; 2]; 3],
    iv: usize,
    nsys: usize,
    end: usize,
    ivte1: usize,
    ivz: usize,
    vacc: [f64; 3],
    lower: &[usize],
) {
    let ivp = iv + 1;
    let r0 = row(vm, ldm, 0, iv);
    let r1 = row(vm, ldm, 1, iv);
    let r2 = row(vm, ldm, 2, iv);
    // ---- column iv and the step scalars
    let p0 = 1.0 / va[iv][0][0];
    va[iv][0][1] *= p0;
    *r0.add(iv) *= p0;
    let va10 = va[iv][1][0];
    va[iv][1][1] -= va10 * va[iv][0][1];
    *r1.add(iv) -= va10 * *r0.add(iv);
    let va20 = va[iv][2][0];
    va[iv][2][1] -= va20 * va[iv][0][1];
    *r2.add(iv) -= va20 * *r0.add(iv);
    let p1 = 1.0 / va[iv][1][1];
    *r1.add(iv) *= p1;
    let va21 = va[iv][2][1];
    *r2.add(iv) -= va21 * *r1.add(iv);
    let p2 = 1.0 / *r2.add(iv);
    let vt1 = *r0.add(iv);
    let vt2 = *r1.add(iv);
    let vt = va[iv][0][1];
    let sc = Scal { p0, va10, va20, p1, va21, p2, vt1, vt2, vt };

    // ---- targets of this step
    let mut vbr: [(*mut f64, f64, f64, f64); 3] = [(std::ptr::null_mut(), 0.0, 0.0, 0.0); 3];
    let mut nvb = 0;
    let mut vzr: [(*mut f64, f64, f64); 3] = [(std::ptr::null_mut(), 0.0, 0.0); 3];
    let mut nvz = 0;
    let mut low: [(*mut f64, f64); 3 * (PB + 2)] = [(std::ptr::null_mut(), 0.0); 3 * (PB + 2)];
    let mut nlow = 0;
    if iv != nsys {
        for k in 0..3 {
            let d = row(vm, ldm, k, ivp);
            vbr[k] = (d, vb[ivp][k][0], vb[ivp][k][1], *d.add(iv));
        }
        nvb = 3;
        if iv == ivte1 {
            for k in 0..3 {
                vzr[k] = (row(vm, ldm, k, ivz), vz[k][0], vz[k][1]);
            }
            nvz = 3;
        }
        if ivp != nsys {
            for &kv in lower {
                if kv < iv + 2 || kv > nsys {
                    continue;
                }
                for k in 0..3 {
                    let d = row(vm, ldm, k, kv);
                    let c = *d.add(iv);
                    if c.abs() > vacc[k] {
                        low[nlow] = (d, c);
                        nlow += 1;
                    }
                }
            }
        }
    }
    fused_pass(r0, r1, r2, &sc, &vbr[..nvb], &vzr[..nvz], &low[..nlow], ivp, end);
}

#[derive(Clone, Copy)]
struct Scal {
    p0: f64,
    va10: f64,
    va20: f64,
    p1: f64,
    va21: f64,
    p2: f64,
    vt1: f64,
    vt2: f64,
    vt: f64,
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
unsafe fn fused_pass(
    r0: *mut f64,
    r1: *mut f64,
    r2: *mut f64,
    s: &Scal,
    vbr: &[(*mut f64, f64, f64, f64)],
    vzr: &[(*mut f64, f64, f64)],
    low: &[(*mut f64, f64)],
    lo: usize,
    hi: usize,
) {
    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    {
        fused_pass_avx2(r0, r1, r2, s, vbr, vzr, low, lo, hi);
    }
    #[cfg(not(all(target_arch = "x86_64", target_feature = "avx2")))]
    {
        fused_scalar(r0, r1, r2, s, vbr, vzr, low, lo, hi);
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments, dead_code)]
unsafe fn fused_scalar(
    r0: *mut f64,
    r1: *mut f64,
    r2: *mut f64,
    s: &Scal,
    vbr: &[(*mut f64, f64, f64, f64)],
    vzr: &[(*mut f64, f64, f64)],
    low: &[(*mut f64, f64)],
    lo: usize,
    hi: usize,
) {
    for l in lo..=hi {
        let mut a0 = *r0.add(l) * s.p0;
        let mut a1 = *r1.add(l) - s.va10 * a0;
        let mut a2 = *r2.add(l) - s.va20 * a0;
        a1 *= s.p1;
        a2 -= s.va21 * a1;
        a2 *= s.p2;
        a0 -= s.vt1 * a2;
        a1 -= s.vt2 * a2;
        a0 -= s.vt * a1;
        *r0.add(l) = a0;
        *r1.add(l) = a1;
        *r2.add(l) = a2;
        for &(d, b0, b1, b2) in vbr {
            *d.add(l) -= b0 * a0 + b1 * a1 + b2 * a2;
        }
        for &(d, z0, z1) in vzr {
            *d.add(l) -= z0 * a0 + z1 * a1;
        }
        for &(d, c) in low {
            *d.add(l) -= c * a2;
        }
    }
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline(always)]
#[allow(clippy::too_many_arguments)]
unsafe fn fused_pass_avx2(
    r0: *mut f64,
    r1: *mut f64,
    r2: *mut f64,
    s: &Scal,
    vbr: &[(*mut f64, f64, f64, f64)],
    vzr: &[(*mut f64, f64, f64)],
    low: &[(*mut f64, f64)],
    lo: usize,
    hi: usize,
) {
    use std::arch::x86_64::*;
    let p0 = _mm256_set1_pd(s.p0);
    let va10 = _mm256_set1_pd(s.va10);
    let va20 = _mm256_set1_pd(s.va20);
    let p1 = _mm256_set1_pd(s.p1);
    let va21 = _mm256_set1_pd(s.va21);
    let p2 = _mm256_set1_pd(s.p2);
    let vt1 = _mm256_set1_pd(s.vt1);
    let vt2 = _mm256_set1_pd(s.vt2);
    let vt = _mm256_set1_pd(s.vt);
    let mut l = lo;
    while l + 4 <= hi + 1 {
        let a0 = _mm256_mul_pd(_mm256_loadu_pd(r0.add(l)), p0);
        let a1 = _mm256_sub_pd(_mm256_loadu_pd(r1.add(l)), _mm256_mul_pd(va10, a0));
        let a2 = _mm256_sub_pd(_mm256_loadu_pd(r2.add(l)), _mm256_mul_pd(va20, a0));
        let a1 = _mm256_mul_pd(a1, p1);
        let a2 = _mm256_sub_pd(a2, _mm256_mul_pd(va21, a1));
        let a2 = _mm256_mul_pd(a2, p2);
        let a0 = _mm256_sub_pd(a0, _mm256_mul_pd(vt1, a2));
        let a1 = _mm256_sub_pd(a1, _mm256_mul_pd(vt2, a2));
        let a0 = _mm256_sub_pd(a0, _mm256_mul_pd(vt, a1));
        _mm256_storeu_pd(r0.add(l), a0);
        _mm256_storeu_pd(r1.add(l), a1);
        _mm256_storeu_pd(r2.add(l), a2);
        for &(d, b0, b1, b2) in vbr {
            let t = _mm256_add_pd(
                _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(b0), a0), _mm256_mul_pd(_mm256_set1_pd(b1), a1)),
                _mm256_mul_pd(_mm256_set1_pd(b2), a2),
            );
            _mm256_storeu_pd(d.add(l), _mm256_sub_pd(_mm256_loadu_pd(d.add(l)), t));
        }
        for &(d, z0, z1) in vzr {
            let t = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(z0), a0), _mm256_mul_pd(_mm256_set1_pd(z1), a1));
            _mm256_storeu_pd(d.add(l), _mm256_sub_pd(_mm256_loadu_pd(d.add(l)), t));
        }
        for &(d, c) in low {
            let t = _mm256_mul_pd(_mm256_set1_pd(c), a2);
            _mm256_storeu_pd(d.add(l), _mm256_sub_pd(_mm256_loadu_pd(d.add(l)), t));
        }
        l += 4;
    }
    if l <= hi {
        fused_scalar(r0, r1, r2, s, vbr, vzr, low, l, hi);
    }
}

// Apply all lower-VM updates from panel pivots b0..b1 to one trailing row.
/// `d` is the row (columns 0..=end), `piv[j]` the base pointer of pivot row 2 of block b0+j.
#[inline(always)]
unsafe fn trailing_row(d: *mut f64, piv: &[usize], b0: usize, b1: usize, end: usize, ldm: usize, vacc: f64) {
    let np = b1 - b0;
    let mut cs = [0.0f64; PB];
    let mut pp = [std::ptr::null::<f64>(); PB];
    let mut nop = 0;
    if np == PB && b0 + PB <= ldm {
        // 1) panel columns in a fixed-width local buffer.  Columns <= iv are
        //    dead once their decision is made, so full-width updates are safe;
        //    live columns (> iv) see exactly the original operation sequence.
        let mut x = [0.0f64; PB];
        std::ptr::copy_nonoverlapping(d.add(b0), x.as_mut_ptr(), PB);
        for j in 0..PB {
            let c = x[j];
            if c.abs() > vacc {
                let p = piv[j] as *const f64;
                crate::simd::axmy16(&mut x, c, p.add(b0));
                cs[nop] = c;
                pp[nop] = p.add(b1);
                nop += 1;
            }
        }
    } else {
        for iv in b0..b1 {
            let c = *d.add(iv);
            if c.abs() > vacc {
                let p = piv[iv - b0] as *const f64;
                for l in iv + 1..b1 {
                    *d.add(l) -= c * *p.add(l);
                }
                cs[nop] = c;
                pp[nop] = p.add(b1);
                nop += 1;
            }
        }
    }
    if nop == 0 {
        return;
    }
    // 2) trailing columns: same op order per element, row chunk kept in registers
    crate::simd::multi_axmy(d.add(b1), &cs[..nop], &pp[..nop], end + 1 - b1);
}

struct SendPtr(*mut f64);
unsafe impl Send for SendPtr {}
unsafe impl Sync for SendPtr {}
impl SendPtr {
    #[inline(always)]
    fn get(&self) -> *mut f64 {
        self.0
    }
}

/// Blocked BLSOLV.  `vm` rows: VM(k,l,iv) at ((iv*3+k)*ldm + l); ldm >= nsys+3.
#[allow(clippy::too_many_arguments)]
pub fn blsolv_blocked(
    vm: &mut [f64],
    ldm: usize,
    va: &mut [[[f64; 2]; 3]],
    vb: &[[[f64; 2]; 3]],
    vz: &[[f64; 2]; 3],
    vdel: &mut [[[f64; 2]; 3]],
    nsys: usize,
    ivte1: usize,
    ivz: usize,
    vacc: [f64; 3],
) {
    assert!(ldm >= nsys + 3);
    assert!(vm.len() >= (nsys * 3 + 3) * ldm);
    let end = nsys + 2;
    let base = vm.as_mut_ptr();
    unsafe {
        for iv in 1..=nsys {
            for k in 0..3 {
                let r = row(base, ldm, k, iv);
                *r.add(nsys + 1) = vdel[iv][k][0];
                *r.add(nsys + 2) = vdel[iv][k][1];
            }
        }
    }
    let par = crate::par_fine() && nsys > 96;
    let mut b0 = 1;
    let mut lower: Vec<usize> = Vec::with_capacity(PB + 4);
    while b0 <= nsys {
        let b1 = (b0 + PB).min(nsys + 1);
        let ivz_inside = ivte1 >= b0 && ivte1 < b1 && ivz > b1;
        lower.clear();
        for kv in b0..=b1.min(nsys) {
            lower.push(kv);
        }
        if ivz_inside {
            lower.push(ivz);
        }
        unsafe {
            for iv in b0..b1 {
                step_inside(base, ldm, va, vb, vz, iv, nsys, end, ivte1, ivz, vacc, &lower);
            }
        }
        // trailing rows
        if b1 < nsys {
            let mut piv = [0usize; PB];
            for iv in b0..b1 {
                piv[iv - b0] = unsafe { row(base, ldm, 2, iv) as usize };
            }
            let piv = &piv[..b1 - b0];
            let first = b1 + 1;
            let nrows = nsys - b1;
            let skip = if ivz_inside { ivz } else { usize::MAX };
            let bp = SendPtr(base);
            let work = |kv: usize| {
                if kv == skip {
                    return;
                }
                for k in 0..3 {
                    unsafe {
                        let d = row(bp.get(), ldm, k, kv);
                        trailing_row(d, piv, b0, b1, end, ldm, vacc[k]);
                    }
                }
            };
            if par && nrows * (nsys - b1) > 12000 {
                (first..nsys + 1).into_par_iter().with_min_len(8).for_each(|kv| work(kv));
            } else {
                for kv in first..=nsys {
                    work(kv);
                }
            }
        }
        b0 = b1;
    }
    // back substitution (row-oriented, same per-element subtraction order as
    // the original; the 6 independent chains of a block row are interleaved)
    unsafe {
        let mut x3 = vec![[0.0f64; 2]; nsys + 2];
        for kv in (1..=nsys).rev() {
            let r0 = row(base, ldm, 0, kv);
            let r1 = row(base, ldm, 1, kv);
            let r2 = row(base, ldm, 2, kv);
            let mut a = [
                *r0.add(nsys + 1),
                *r0.add(nsys + 2),
                *r1.add(nsys + 1),
                *r1.add(nsys + 2),
                *r2.add(nsys + 1),
                *r2.add(nsys + 2),
            ];
            for iv in (kv + 1..=nsys).rev() {
                let x = *x3.get_unchecked(iv);
                let m0 = *r0.add(iv);
                let m1 = *r1.add(iv);
                let m2 = *r2.add(iv);
                a[0] -= m0 * x[0];
                a[1] -= m0 * x[1];
                a[2] -= m1 * x[0];
                a[3] -= m1 * x[1];
                a[4] -= m2 * x[0];
                a[5] -= m2 * x[1];
            }
            vdel[kv][0] = [a[0], a[1]];
            vdel[kv][1] = [a[2], a[3]];
            vdel[kv][2] = [a[4], a[5]];
            x3[kv] = [a[4], a[5]];
        }
    }
}
