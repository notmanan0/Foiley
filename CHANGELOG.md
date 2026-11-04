# Changelog of reported numbers

Every change to a number in a golden test or the README, per
docs/REMEDIATION_PLAN.md ("Record").

| Item | Quantity | Fixture | Before | After | Why |
|---|---|---|---|---|---|
| F014 | cancel latency during LU factorization (README, "Why it does not freeze") | 14,100-unknown dense factorization, cancel 1 s in, 10 runs (`crates/core/tests/f014.rs`) | "within milliseconds" (untrue: the single faer call could not be interrupted, ~5-7 s) | 17 ms p95 (1.6-16.9 ms) | blocked right-looking LU with cancellation checks between column tiles; the matrix buffer is freed off-thread on cancel |
| F007 | T2 lift/weight with a relaxed wake (`benchmark_review` example, `f007_t2_lift_equals_weight`) | review rectangle, VLM2, 4 free-wake iterations, alpha 8 | 1.000857465 (speed balance from the rigid wake, loads from the relaxed one) | 1.000000000 | every driver goes through `Evaluator::loads`; the balance uses the relaxed wake at unit speed, which does not depend on the speed |
| F007 | T3 relaxed minus rigid CL (`benchmark_review`, `f007_t3_relaxed_differs_from_rigid`) | same | 0 (the glide loop re-solved with the rigid wake), 0 wake lines kept | 0.000566630691, 17 wake lines | the glide loop uses the relaxed wake |
| F007 / F019 | T6 CL with a relaxed wake (`benchmark_review`, `f019_t1_t6_equivalent`) | same, control parameter 0 | rejected (free wake unsupported in T6) | 0.661387325, equal to T1 | T6 and T7 go through the evaluator; free wake is allowed in every VLM driver |
| F022c | triangles from `read_step` of a bilinear B-spline patch, density 4 (`cad::tests::step_spline_surface`) | inline 2 x 2 B_SPLINE_SURFACE_WITH_KNOTS | 32 (fixed 4 x 4 grid) | adaptive, asserted only >= 32 plus area 1.2807 within 2e-3 | STEP faces are tessellated adaptively (chordal / length tolerance) instead of on a fixed parameter grid |
| F015b | the README's "bit-exact with the Fortran" claim (lines 6, 33, 246, 247) | 48-case raw-bit corpus against `build-ref-dp` (`tests/f015.rs::f015_matrix`), 4 November 2026 | "numbers identical to the original", "bit-exact" | 0 of 48 cases bit-identical; largest scaled difference 1.5e-5 (XTR1), CL 4.5e-7. `FOILEY_PIPELINE=1`, `FOILEY_FINE_PAR=1`, the x86-64-v2 build with `FOILEY_SIMD=0` and with runtime AVX2 are all bit-identical to the sequential run on 48 of 48; the parallel sweep (`FOILEY_JOBS=4`) on 15 of 48, otherwise within 1.3e-4 | the claim was never measured; the matrix measures it and the README states what it found |
