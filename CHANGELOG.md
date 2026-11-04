# Changelog of reported numbers

Every change to a number in a golden test or the README, per
docs/REMEDIATION_PLAN.md ("Record").

| Item | Quantity | Fixture | Before | After | Why |
|---|---|---|---|---|---|
| F014 | cancel latency during LU factorization (README, "Why it does not freeze") | 14,100-unknown dense factorization, cancel 1 s in, 10 runs (`crates/core/tests/f014.rs`) | "within milliseconds" (untrue: the single faer call could not be interrupted, ~5-7 s) | 17 ms p95 (1.6-16.9 ms) | blocked right-looking LU with cancellation checks between column tiles; the matrix buffer is freed off-thread on cancel |
| F022c | triangles from `read_step` of a bilinear B-spline patch, density 4 (`cad::tests::step_spline_surface`) | inline 2 x 2 B_SPLINE_SURFACE_WITH_KNOTS | 32 (fixed 4 x 4 grid) | adaptive, asserted only >= 32 plus area 1.2807 within 2e-3 | STEP faces are tessellated adaptively (chordal / length tolerance) instead of on a fixed parameter grid |
