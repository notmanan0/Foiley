# Changelog of reported numbers

Every change to a number in a golden test or the README, per
docs/REMEDIATION_PLAN.md ("Record").

| Item | Quantity | Fixture | Before | After | Why |
|---|---|---|---|---|---|
| F014 | cancel latency during LU factorization (README, "Why it does not freeze") | 14,100-unknown dense factorization, cancel 1 s in, 10 runs (`crates/core/tests/f014.rs`) | "within milliseconds" (untrue: the single faer call could not be interrupted, ~5-7 s) | 17 ms p95 (1.6-16.9 ms) | blocked right-looking LU with cancellation checks between column tiles; the matrix buffer is freed off-thread on cancel |
