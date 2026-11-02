//! Explicit SIMD kernels (AVX2 with runtime detection, scalar fallback).
//! Multiplication and subtraction are kept as separate instructions (no FMA
//! contraction) so results are bit-identical to the scalar Fortran.

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[inline(always)]
fn has_avx2() -> bool {
    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    {
        true
    }
    #[cfg(all(target_arch = "x86_64", not(target_feature = "avx2")))]
    {
        static F: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
        *F.get_or_init(|| is_x86_feature_detected!("avx2"))
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// d[l] = (((d[l] - c0*p0[l]) - c1*p1[l]) - ...)  for l in 0..n
///
/// # Safety
/// `d` must be valid for n writes, each `ps[o]` valid for n reads, and `d`
/// must not overlap any `ps[o]`.
#[inline]
pub unsafe fn multi_axmy(d: *mut f64, cs: &[f64], ps: &[*const f64], n: usize) {
    #[cfg(target_arch = "x86_64")]
    if has_avx2() {
        return multi_axmy_avx2(d, cs, ps, n);
    }
    multi_axmy_scalar(d, cs, ps, n, 0);
}

#[inline(always)]
unsafe fn multi_axmy_scalar(d: *mut f64, cs: &[f64], ps: &[*const f64], n: usize, start: usize) {
    for l in start..n {
        let mut x = *d.add(l);
        for (c, p) in cs.iter().zip(ps.iter()) {
            x -= *c * *p.add(l);
        }
        *d.add(l) = x;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn multi_axmy_avx2(d: *mut f64, cs: &[f64], ps: &[*const f64], n: usize) {
    let nop = cs.len();
    let mut l = 0;
    while l + 16 <= n {
        let mut x0 = _mm256_loadu_pd(d.add(l));
        let mut x1 = _mm256_loadu_pd(d.add(l + 4));
        let mut x2 = _mm256_loadu_pd(d.add(l + 8));
        let mut x3 = _mm256_loadu_pd(d.add(l + 12));
        for o in 0..nop {
            let c = _mm256_set1_pd(*cs.get_unchecked(o));
            let p = ps.get_unchecked(o).add(l);
            x0 = _mm256_sub_pd(x0, _mm256_mul_pd(c, _mm256_loadu_pd(p)));
            x1 = _mm256_sub_pd(x1, _mm256_mul_pd(c, _mm256_loadu_pd(p.add(4))));
            x2 = _mm256_sub_pd(x2, _mm256_mul_pd(c, _mm256_loadu_pd(p.add(8))));
            x3 = _mm256_sub_pd(x3, _mm256_mul_pd(c, _mm256_loadu_pd(p.add(12))));
        }
        _mm256_storeu_pd(d.add(l), x0);
        _mm256_storeu_pd(d.add(l + 4), x1);
        _mm256_storeu_pd(d.add(l + 8), x2);
        _mm256_storeu_pd(d.add(l + 12), x3);
        l += 16;
    }
    while l + 4 <= n {
        let mut x0 = _mm256_loadu_pd(d.add(l));
        for o in 0..nop {
            let c = _mm256_set1_pd(*cs.get_unchecked(o));
            x0 = _mm256_sub_pd(x0, _mm256_mul_pd(c, _mm256_loadu_pd(ps.get_unchecked(o).add(l))));
        }
        _mm256_storeu_pd(d.add(l), x0);
        l += 4;
    }
    multi_axmy_scalar(d, cs, ps, n, l);
}

/// x[j] -= c * p[j] for j in 0..16 (separate mul and sub).
///
/// # Safety
/// `p` must be valid for 16 reads.
#[inline(always)]
pub unsafe fn axmy16(x: &mut [f64; 16], c: f64, p: *const f64) {
    #[cfg(target_arch = "x86_64")]
    if has_avx2() {
        return axmy16_avx2(x, c, p);
    }
    for j in 0..16 {
        x[j] -= c * *p.add(j);
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn axmy16_avx2(x: &mut [f64; 16], c: f64, p: *const f64) {
    let cv = _mm256_set1_pd(c);
    let xp = x.as_mut_ptr();
    for k in 0..4 {
        let v = _mm256_loadu_pd(xp.add(4 * k));
        let w = _mm256_mul_pd(cv, _mm256_loadu_pd(p.add(4 * k)));
        _mm256_storeu_pd(xp.add(4 * k), _mm256_sub_pd(v, w));
    }
}

macro_rules! dispatch {
    ($avx:ident, $sc:ident, $($a:expr),*) => {{
        #[cfg(target_arch = "x86_64")]
        if has_avx2() {
            return $avx($($a),*);
        }
        $sc($($a),*)
    }};
}

/// d[l] *= c, l in 0..n
///
/// # Safety
/// `d` valid for n reads/writes.
#[inline]
pub unsafe fn scale(d: *mut f64, c: f64, n: usize) {
    dispatch!(scale_avx2, scale_sc, d, c, n)
}
#[inline(always)]
unsafe fn scale_sc(d: *mut f64, c: f64, n: usize) {
    for l in 0..n {
        *d.add(l) *= c;
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn scale_avx2(d: *mut f64, c: f64, n: usize) {
    let cv = _mm256_set1_pd(c);
    let mut l = 0;
    while l + 4 <= n {
        _mm256_storeu_pd(d.add(l), _mm256_mul_pd(_mm256_loadu_pd(d.add(l)), cv));
        l += 4;
    }
    while l < n {
        *d.add(l) *= c;
        l += 1;
    }
}

/// d[l] -= c*s[l], l in 0..n
///
/// # Safety
/// valid pointers, `d` not overlapping `s`.
#[inline]
pub unsafe fn axmy(d: *mut f64, s: *const f64, c: f64, n: usize) {
    dispatch!(axmy_avx2, axmy_sc, d, s, c, n)
}
#[inline(always)]
unsafe fn axmy_sc(d: *mut f64, s: *const f64, c: f64, n: usize) {
    for l in 0..n {
        *d.add(l) -= c * *s.add(l);
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn axmy_avx2(d: *mut f64, s: *const f64, c: f64, n: usize) {
    let cv = _mm256_set1_pd(c);
    let mut l = 0;
    while l + 4 <= n {
        let v = _mm256_sub_pd(_mm256_loadu_pd(d.add(l)), _mm256_mul_pd(cv, _mm256_loadu_pd(s.add(l))));
        _mm256_storeu_pd(d.add(l), v);
        l += 4;
    }
    while l < n {
        *d.add(l) -= c * *s.add(l);
        l += 1;
    }
}

/// d[l] -= (a*s0[l] + b*s1[l] + c*s2[l])   (Fortran left-to-right order)
///
/// # Safety
/// valid pointers, `d` not overlapping sources.
#[inline]
pub unsafe fn axmy3(d: *mut f64, s0: *const f64, s1: *const f64, s2: *const f64, a: f64, b: f64, c: f64, n: usize) {
    dispatch!(axmy3_avx2, axmy3_sc, d, s0, s1, s2, a, b, c, n)
}
#[inline(always)]
unsafe fn axmy3_sc(d: *mut f64, s0: *const f64, s1: *const f64, s2: *const f64, a: f64, b: f64, c: f64, n: usize) {
    for l in 0..n {
        *d.add(l) -= a * *s0.add(l) + b * *s1.add(l) + c * *s2.add(l);
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn axmy3_avx2(d: *mut f64, s0: *const f64, s1: *const f64, s2: *const f64, a: f64, b: f64, c: f64, n: usize) {
    let (av, bv, cv) = (_mm256_set1_pd(a), _mm256_set1_pd(b), _mm256_set1_pd(c));
    let mut l = 0;
    while l + 4 <= n {
        let t = _mm256_add_pd(
            _mm256_add_pd(_mm256_mul_pd(av, _mm256_loadu_pd(s0.add(l))), _mm256_mul_pd(bv, _mm256_loadu_pd(s1.add(l)))),
            _mm256_mul_pd(cv, _mm256_loadu_pd(s2.add(l))),
        );
        _mm256_storeu_pd(d.add(l), _mm256_sub_pd(_mm256_loadu_pd(d.add(l)), t));
        l += 4;
    }
    while l < n {
        *d.add(l) -= a * *s0.add(l) + b * *s1.add(l) + c * *s2.add(l);
        l += 1;
    }
}

/// d[l] -= (a*s0[l] + b*s1[l])
///
/// # Safety
/// valid pointers, `d` not overlapping sources.
#[inline]
pub unsafe fn axmy2(d: *mut f64, s0: *const f64, s1: *const f64, a: f64, b: f64, n: usize) {
    dispatch!(axmy2_avx2, axmy2_sc, d, s0, s1, a, b, n)
}
#[inline(always)]
unsafe fn axmy2_sc(d: *mut f64, s0: *const f64, s1: *const f64, a: f64, b: f64, n: usize) {
    for l in 0..n {
        *d.add(l) -= a * *s0.add(l) + b * *s1.add(l);
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn axmy2_avx2(d: *mut f64, s0: *const f64, s1: *const f64, a: f64, b: f64, n: usize) {
    let (av, bv) = (_mm256_set1_pd(a), _mm256_set1_pd(b));
    let mut l = 0;
    while l + 4 <= n {
        let t = _mm256_add_pd(_mm256_mul_pd(av, _mm256_loadu_pd(s0.add(l))), _mm256_mul_pd(bv, _mm256_loadu_pd(s1.add(l))));
        _mm256_storeu_pd(d.add(l), _mm256_sub_pd(_mm256_loadu_pd(d.add(l)), t));
        l += 4;
    }
    while l < n {
        *d.add(l) -= a * *s0.add(l) + b * *s1.add(l);
        l += 1;
    }
}
