# Changelog of reported numbers

Every change to a number in a golden test or the README, per
docs/REMEDIATION_PLAN.md ("Record").

| Item | Quantity | Fixture | Before | After | Why |
|---|---|---|---|---|---|
| F014 | cancel latency during LU factorization (README, "Why it does not freeze") | 14,100-unknown dense factorization, cancel 1 s in, 10 runs (`crates/core/tests/f014.rs`) | "within milliseconds" (untrue: the single faer call could not be interrupted, ~5-7 s) | 17 ms p95 (1.6-16.9 ms) | blocked right-looking LU with cancellation checks between column tiles; the matrix buffer is freed off-thread on cancel |
| F007 | T2 lift/weight with a relaxed wake (`benchmark_review` example, `f007_t2_lift_equals_weight`) | review rectangle, VLM2, 4 free-wake iterations, alpha 8 | 1.000857465 (speed balance from the rigid wake, loads from the relaxed one) | 1.000000000 | every driver goes through `Evaluator::loads`; the balance uses the relaxed wake at unit speed, which does not depend on the speed |
| F007 | T3 relaxed minus rigid CL (`benchmark_review`, `f007_t3_relaxed_differs_from_rigid`) | same | 0 (the glide loop re-solved with the rigid wake), 0 wake lines kept | 0.000566630691, 17 wake lines | the glide loop uses the relaxed wake |
| F007 / F019 | T6 CL with a relaxed wake (`benchmark_review`, `f019_t1_t6_equivalent`) | same, control parameter 0 | rejected (free wake unsupported in T6) | 0.661387325, equal to T1 | T6 and T7 go through the evaluator; free wake is allowed in every VLM driver |
