# Foiley remediation and capability plan

Prepared 4 November 2026 against `7e7fb7f`. Replaces the phased plan of the
same date. Covers every finding (F001–F024) and every gap (G01–G47) from the
4 November 2026 code review, except G33 (marine/hydrofoil) and G42
(GPU/distributed), which are dropped by decision.
The review documents have been removed; each item below states its own
problem and fix.

This plan is a hill climb, not a schedule. There is one scoreboard (§1).
Every merged change must leave every scoreboard metric no worse and at least
one better. There are no phases: the backlog (§4) is a queue of items with
explicit dependencies, and whoever works it, person or agent, takes the
highest-ranked unblocked item. This file is the checklist. Status lives in
the queue table in §4 and is updated in the same commit that lands the item.

## 0. How to work this plan

**Completion condition.** Every row of the queue in §4 is `[x]`, every target
in §1 is met, the perf gate (§3) is green on `master`, and `CHANGELOG.md`
records every changed number. Until then there is always a next item; the
plan never ends in "what would you like next".

**The loop**, once per item:

1. **Pick** the lowest rank in §4 whose `blocked by` items are all `[x]`.
   Ties go to the item that moves more scoreboard metrics. Say in one line
   which item is starting and why that one.
2. **Explore before editing.** Open every file in the item's `touches` list,
   grep every caller of each function the spec names, and run the item's
   proof (it should fail or be ignored) before changing anything. Most review
   findings came from drivers that duplicate each other (T1–T8 speed balance,
   glide loop, controls, trim, derivatives), so assume the symptom appears in
   more places than the spec lists, and fix all of them.
3. **Implement** the spec. Write the proof first if it does not exist yet.
4. **Prove.** The proof passes un-ignored and `cargo test --release
   --workspace` is green.
5. **Measure.** `foiley-bench compare` passes the gate in §3 and the item's
   own `budget` holds. Record the numbers.
6. **Record.** Update the §1 numbers, flip the item's checkbox in §4 with the
   commit hash and the measured numbers, and add any changed numbers to
   `CHANGELOG.md` (§5).
7. **Land.** Branch `item/<ID>-<slug>`, PR to `origin/master`, squash merge.
   A PR that lowers any metric does not merge; it gets fixed or split.
   Then go to 1.

**What "up" means.** A change is accepted when S3 stays inside the gate on
every case and at least one of S1, S2, S4–S10 improves. Infrastructure items
(B1–B4) are accepted on their proof alone. An item that cannot be proven by a
test or a measured number is not finished; find the number.

**Wanted stops.** Stop and say what is blocking only when:
- the item needs something only the user can supply. In this plan that is
  exactly one thing: real `.xfl`/`.fl5` files for the last step of G45. It
  blocks nothing else, so it is not a reason to stop before everything else
  is done;
- an action would delete data, rewrite published history, or force-push;
- the perf gate fails and the only remaining fix is a written exception (§3).
  Exceptions are the user's call; everything else is not.

**Unwanted stops.** Ending a turn with a summary that names the next item
instead of starting it; offering to continue unless told otherwise; listing
the decisions in §2 as questions (they are decided); pausing at a milestone
to report. Status notes go in the commit message and PR body, and the work
carries on. If the spec and the codebase disagree and both readings are
defensible, state both in the PR, pick one, and continue.

**Effort and time.** Default effort is `medium`. Four items are marked
`hard` (F007, F005, F022, G41) because the design has a known trap; run those
at `high`, and go to `xhigh` only if a first attempt at `high` measurably
misses the budget. Each item carries an advisory time budget in hours. Pace
to finish inside it. If an item passes 2× its budget, land what already
passes its proof and append the remainder to §4 as a new item with its own
rank; running out of time is a reason to split, never to stop. When a lead
agent delegates items to subagents, append `elapsed Ns / Ms` (M = the item's
budget in seconds) to each message returned to the subagent.

**Reporting.** One line of intent before the first tool call of each item. A
PR body that states what changed, what was measured, and what the scoreboard
reads now. Nothing else needs narrating.

## 1. Scoreboard

Numbers in "Now" are as of `7e7fb7f`. Every item in §4 names which rows it
moves. Update this table in the commit that moves a number.

| # | Metric | Measured by | Now | Target |
|---|---|---|---|---|
| S1 | Finding proofs passing | `review_findings` tests in `crates/plane/tests`, `crates/foil/tests` and `tests/` (un-ignored and green) | 0 / 24 (the files do not exist; `benchmark_review.rs` only prints 9 of them) | 24 / 24 |
| S2 | Gap proofs passing | each G item's named proof | 0 / 45 | 45 / 45 |
| S3 | Perf suite: per-case median of 5 interleaved runs vs `bench/baseline.json` | `foiley-bench compare` | 1.00 (baseline not yet recorded) | geometric mean ≤ 0.70; no case > 1.02 at any merge |
| S4 | Peak RSS at 14,100 unknowns | perf suite case `trilinear_14100` | ≈ 1.6 GB (N²·8 B, estimated) | ≤ 0.4 GB |
| S5 | Cancel latency during a 14,100-unknown factorization | perf suite case `cancel_14100` | one whole factorization | ≤ 100 ms |
| S6 | 2D configurations with a proven bit-identity status vs XFoil 6.97 | F015 raw-bit matrix | 0 proven, 4 claimed | every configuration in the matrix proven identical or stated with its max ulp |
| S7 | Silent failure paths in 3D drivers (`.ok()`, `if let Ok`, `unwrap_or` on solver results) | grep proof `no_silent_ok` | count taken by B2 | 0 |
| S8 | README performance and accuracy claims backed by a runner result | audit list in `bench/claims.md` (B3) | 0 of all | all |
| S9 | Validation error per track vs reference solvers and experiment | G01/G02 scorecards | none | published per track; no single averaged number |
| S10 | `.foiley` fixtures that round-trip load→save→load | `project_roundtrip` test | not tested | all files under `examples/`, `spitfire/`, `tests/data` |

## 2. Decisions

All decided. None of these is a question for later.

| Topic | Decision |
|---|---|
| Scope | All findings and all gaps except two dropped by decision: G33 (marine/hydrofoil, no free-surface or cavitation work) and G42 (GPU/distributed, CPU only). |
| Overriding constraint | **Extreme efficiency.** Every item has a budget and must pass the perf gate (§3). |
| Objective | The scoreboard (§1). An item is done when its proof passes and its numbers are recorded. |
| Ordering | By rank in §4, subject to `blocked by`. Ranks may be edited; dependencies may not be skipped. |
| Delivery | One branch and one PR per item to `origin/master` (`origin` = `github.com/notmanan0/Foiley`, default branch `master`; the earlier plan's "no remote" was wrong). Squash merge. CI (B4) runs tests and the correctness half of the perf gate; timing numbers count only from the baseline machine. |
| Baseline machine | Intel i7-14650HX, 24 logical CPUs, Windows 11. Power plan, thread pins and `target-cpu=native` are stated in `bench/baseline.json`. |
| Multiple bodies (F004) | Mesh every body, with per-body junction trimming and carry-through. |
| Stability/trim (F005/F006) | Two models selected by a setting: `Inviscid` (default, current numbers) and `FullLoad`. Nonzero-beta reference states supported. |
| Free wake (F007/F019) | One shared wake-aware load evaluator used by every driver, converging on a residual. |
| LU cancellation (F014) | Blocked LU with a cancel check between panels. |
| STEP (F022) | Pure-Rust trimmed B-rep tessellation. No OpenCascade. |
| Project format | `.foiley` v2. v1 loads through a migration; saves are always v2. |
| Changed numbers | No legacy switch. Each PR records before/after in `CHANGELOG.md` and updates golden tests (§5). |
| Bit-exact claim (F015) | Prove with raw-bit tests; keep the claim only for configurations that pass. |
| External solvers (G01) | XFOIL 6.97 + 6.99, AVL, flow5, VSPAERO 7, MachLine, MachUpX, FLOWPanel.jl. |
| Datasets (G02) | Download scripts with SHA-256 checksums, never vendored data. This holds regardless of licence (UIUC data is free to use but carries no redistribution licence; NACA reports are public domain) because it keeps the repo small and makes the licence question moot. |
| Interchange fixtures (G45) | The parser is built against hand-written minimal `.xfl`/`.fl5` files and the flow5 API example already at `refsrc/flow5/API_examples/`. Real files from XFLR5 6.x and flow5 7.58 go in `tests/data/interchange/` when the user adds them; that is the only step in this plan that waits on the user, and it is G45's last sub-step, not its first. |
| G15 (shock-capable 2D) | A time-boxed spike with a scored exit (see G15). The decision to continue is made by the number, not by a meeting. |
| Physics priority after remediation | 1. fixed-wing RC/UAV, 2. powered/unsteady, 3. airfoil design, 4. flexible wings. Expressed as ranks in §4, not as phases. |

## 3. Standing rules

**Perf gate.** B1 adds `crates/bench` (binary `foiley-bench`) with a fixed
suite. It records phase timings and peak RSS as JSON and compares against
the committed `bench/baseline.json`:

| Case | What it stresses |
|---|---|
| `xfoil_8` — the 8 existing 2D cases from `bench/run_bench.py` | BL Newton, SIMD, sweeps |
| `vlm2_t1_100` — VLM2 T1, 100 alphas, ≈2k panels | assembly, unit-solution reuse, `reduce` |
| `vlm2_t2_relaxed_20` — VLM2 T2 relaxed wake, 20 alphas | wake evaluator |
| `panel_t1_6k` — thick quad panels T1, ≈6k panels | assembly + dense LU |
| `trilinear_14100` — TriLinear fine wing, 14,100 unknowns | LU at scale, memory (S4) |
| `cancel_14100` — same mesh, cancel issued 1 s into the factorization | cancel latency (S5) |
| `t7_stability_10` — T7 stability, 10 control steps | trim + derivatives + rebuilds |
| `step_fuselage` — STEP import of the reference fuselage | CAD tessellation |

Rules:

- A PR fails if the median of 5 interleaved runs regresses more than **2%**
  on any case. Exceptions need a written justification in the PR and the
  user's agreement. Correctness fixes that must add work (meshing a second
  body, for example) are judged against the budget stated in their item, not
  against zero.
- Hot loops (per point, per strip, per panel) stay allocation-free. No string
  hashing in them: foil and family names resolve to integer indices once per run.
- Diagnostics cost at most O(N²) and are skipped when not requested. Nothing
  in this plan may add an O(N³) step to the default path.
- Serial and full-machine modes are measured separately (`FOILEY_JOBS=1
  FOILEY_PIPELINE=0 FOILEY_FINE_PAR=0`, plus pinned Rayon/faer threads).
- Prefer algorithmic wins over micro-tuning. The big ones are named in their
  items: speed-invariant free wake (F007), low-rank wake updates (F007),
  mixed-precision LU (G40), H-matrices (G41).

**Tests.** `cargo test --release --workspace` is the correctness gate. Proofs
for F items live in `crates/plane/tests/review_findings.rs`,
`crates/foil/tests/review_findings.rs`, or `tests/raw_bits.rs` (root crate),
named `f0xx_<what>`. Proofs for G items are named in each item. A proof that
is `#[ignore = "F0xx"]` counts as not passing.

**Sizes and budgets.** Time budgets are in hours and advisory (§0). Perf
budgets are measured by the suite and binding.

## 4. Backlog

### Queue

The checklist. Lowest unblocked rank first. `hard` marks the four items to
run at higher effort. Flip the box in the commit that lands the item, and
write the commit hash and measured numbers in the status column.

| Rank | Item | Title | Blocked by | Moves | Budget h | Status |
|---|---|---|---|---|---|---|
| 1 | [ ] B1 | Perf suite and baseline | — | S3, S4, S5 | 6 | |
| 2 | [ ] B2 | Proof skeleton: 24 ignored finding tests, `no_silent_ok` | — | S1, S7 | 6 | |
| 3 | [ ] B3 | `CHANGELOG.md`, `bench/claims.md`, scoreboard protocol | — | S8 | 1 | |
| 4 | [ ] B4 | CI on `origin` | B2 | — | 2 | |
| 5 | [ ] P1 | Free performance wins already found | B1 | S3 | 3 | |
| 6 | [ ] F009 | Input validation and finite output (+G06) | B2 | S1, S2 | 6 | |
| 7 | [ ] F011 | No silent stability failures | B2 | S1, S7 | 3 | |
| 8 | [ ] F017 | Real 2D iteration counts | B2 | S1 | 2 | |
| 9 | [ ] F020 | T3 convergence flag | B2 | S1 | 2 | |
| 10 | [ ] F023 | LLT settings honesty | F009 | S1 | 3 | |
| 11 | [ ] F012 | Viscous envelope status | B2 | S1 | 4 | |
| 12 | [ ] F013 | Numerical diagnostics (+G05 part) | B1 | S1, S2 | 8 | |
| 13 | [ ] F014 | Cancellable blocked LU | B1 | S1, S5 | 8 | |
| 14 | [ ] F001 | Interim 2D bench harness fixes | B1 | S1 | 4 | |
| 15 | [ ] F015 | Raw-bit evidence for the 2D port | B2 | S1, S6, S8 | 12 | |
| 16 | [ ] F002 | Explicit polar families (+G07) | P1 | S1, S2, S3 | 12 | |
| 17 | [ ] F003 | Flap-deflected polars | F002 | S1 | 6 | |
| 18 | [ ] F016 | Provenance and format v2 (+G03) | F002 | S1, S2, S10 | 10 | |
| 19 | [ ] F021 | Coverage-aware auto-polars (+G08) | F003 | S1, S2 | 8 | |
| 20 | [ ] F018 | Node-preserving 2D input (+G09) | B2 | S1, S2 | 6 | |
| 21 | [ ] F007 | Shared wake-aware evaluator (+F019, G29) `hard` | F020, F014 | S1, S2, S3 | 30 | |
| 22 | [ ] F008 | LLT glide balance | F020 | S1 | 3 | |
| 23 | [ ] F005 | Stability models and nonzero beta (+F006, G34) `hard` | F007 | S1, S2 | 20 | |
| 24 | [ ] F004 | Multiple bodies (+G23) | B1 | S1, S2 | 16 | |
| 25 | [ ] F010 | AVL and interchange fidelity (+G45) | F004 | S1, S2 | 16 | |
| 26 | [ ] F024 | GUI CAD import as a job | — | S1 | 3 | |
| 27 | [ ] F022 | Pure-Rust trimmed STEP (+G26) `hard` | B1 | S1, S2 | 40 | |
| 28 | [ ] G25 | Imported lifting meshes | F022 | S2 | 12 | |
| 29 | [ ] G04 | Instrumentation | F016, F017 | S2 | 4 | |
| 30 | [ ] G01 | Benchmark runner | F001, G04 | S2, S8, S9 | 40 | |
| 31 | [ ] G02 | Experimental validation | G01 | S2, S9 | 12 | |
| 32 | [ ] G05 | Refinement studies (+G44) | G01 | S2, S9 | 12 | |
| 33 | [ ] G40 | Memory-scalable dense systems | F013, F014 | S2, S3, S4 | 24 | |
| 34 | [ ] G41 | Hierarchical acceleration `hard` | G40 | S2, S3, S4 | 60 | |
| 35 | [ ] G43 | Adjoints and AD | G40 | S2 | 24 | |
| 36 | [ ] G28 | Freestream-aligned rigid wake | F007 | S2 | 4 | |
| 37 | [ ] G35 | General constrained trim | F005 | S2 | 16 | |
| 38 | [ ] G37 | 3D hinge loads | F003 | S2 | 6 | |
| 39 | [ ] G22 | VLM body model | F004 | S2 | 12 | |
| 40 | [ ] G24 | Mixed thin/thick systems | G22 | S2 | 16 | |
| 41 | [ ] G27 | 3D compressibility | F002 | S2 | 6 | |
| 42 | [ ] G21 | Nonlinear VLM/panel viscous coupling | F002, F007 | S2 | 16 | |
| 43 | [ ] G20 | General nonlinear LLT | F008, F023 | S2 | 24 | |
| 44 | [ ] G36 | Broader flight dynamics | G35 | S2 | 24 | |
| 45 | [ ] G46 | Field output | G41 | S2 | 12 | |
| 46 | [ ] G47 | C and Python bindings | F016 | S2 | 16 | |
| 47 | [ ] G32 | Actuator disks and propellers | G22, G41 | S2 | 40 | |
| 48 | [ ] G30 | UVLM / unsteady panels | G41, F007 | S2 | 60 | |
| 49 | [ ] G31 | Vortex particles | G30 | S2 | 40 | |
| 50 | [ ] G19 | 2D hinge moments | F017 | S2 | 4 | |
| 51 | [ ] G13 | Post-stall extension | F016 | S2 | 4 | |
| 52 | [ ] G11 | 2D ground and cascade | F018 | S2 | 12 | |
| 53 | [ ] G10 | Multi-element 2D | G11 | S2 | 60 | |
| 54 | [ ] G16 | 2D sensitivities | F017 | S2 | 24 | |
| 55 | [ ] G17 | Multipoint optimization | G16 | S2 | 16 | |
| 56 | [ ] G18 | Surrogate screening | G17 | S2 | 24 | |
| 57 | [ ] G12 | Rotational corrections | G13 | S2 | 12 | |
| 58 | [ ] G14 | Time-accurate 2D | G10 | S2 | 40 | |
| 59 | [ ] G15 | Shock-capable 2D (spike with scored exit) | G10 | S2 | 24 | |
| 60 | [ ] G38 | Structural coupling | G43, G41 | S2 | 60 | |
| 61 | [ ] G39 | Frequency-domain DLM | G30 | S2 | 40 | |
| — | dropped | G33 marine/hydrofoil | | | | |
| — | dropped | G42 GPU/distributed | | | | |

Ranks 1–35 are the remediation and the infrastructure that proves it. Ranks
36 onward are new capability, in physics-priority order. Rank is advice;
`blocked by` is law. Anything unblocked may be taken out of rank order when
there is a reason, stated in the PR.

### Item specs

Each spec has: **problem**, **fix**, **proof** (the test that flips the
box), **budget** (perf, binding), **touches** (where to look first). Where
the spec names a function, grep its callers before editing.

---

#### B1 · Perf suite and baseline

**Fix.** Add `crates/bench` with binary `foiley-bench run|compare` and the
suite in §3. Each case reports phase timings (import, mesh, assembly, factor,
RHS, wake, reduce, derivatives, output), wall time and peak RSS as one JSON
line. `compare` takes the median of 5 interleaved runs per case and exits
nonzero on any case above the 2% gate. Commit `bench/baseline.json` recorded
on the baseline machine (§2) with CPU, power plan, thread counts and
`target-cpu=native` stated in the file.
**Proof.** `foiley-bench compare` on an unchanged tree reports every case
inside the gate; `cancel_14100` and `trilinear_14100` report S4 and S5
numbers into §1.
**Touches.** `bench/run_bench.py` (the case list), `crates/plane/examples/*bench*.rs`,
`crates/core/src/progress.rs`.

#### B2 · Proof skeleton

**Fix.** Create `crates/plane/tests/review_findings.rs`,
`crates/foil/tests/review_findings.rs` and `tests/raw_bits.rs` with one test
per finding, each asserting the *correct* behaviour and `#[ignore = "F0xx"]`
until its item lands. The nine observations printed by
[benchmark_review.rs](../crates/plane/examples/benchmark_review.rs) become
the first nine; the other fifteen are written from the specs below. Keep the
example for printing. Add `no_silent_ok`: a test that greps the 3D driver
functions (`solve_point`, `run_polar`, `run_control`, `run_stability`, the
glide loop) for `.ok()`, `if let Ok` and `unwrap_or` applied to solver
results, and asserts the count against a number that only goes down.
**Proof.** `cargo test --release --workspace -- --list` shows 24 finding
tests; `no_silent_ok` passes with its starting count written into §1 S7.

#### B3 · CHANGELOG, claims list, scoreboard protocol

**Fix.** `CHANGELOG.md` with a "Changed numbers" section (§5).
`bench/claims.md` listing every performance or accuracy claim in `README.md`
(the "identical to the original" and "bit-exact" sentences at lines 6 and 33,
the timing table at 225–231, the 2.0×/2.5×/4.2×/5.1× ratios, the "within
milliseconds" cancel sentence at line 50) with a column for the runner result
that backs it. S8 is the count of rows with a result.
**Proof.** Both files exist; `bench/claims.md` lists every claim with
"unbacked" in its result column.

#### B4 · CI

**Fix.** `.github/workflows/ci.yml`: `cargo test --release --workspace`,
`cargo clippy`, and `foiley-bench run --correctness-only` (no timing
comparisons; CI machines are not the baseline). Required for merge.
**Proof.** A PR from a branch shows the checks.

#### P1 · Free performance wins

**Problem.** `Ctx::reduce` computes bending moments in an O(strips²) loop per
wing per point ([analysis.rs:676](../crates/plane/src/analysis.rs)).
`strip_viscous` ([analysis.rs:596](../crates/plane/src/analysis.rs)) does two
`HashMap<String, _>` lookups per strip per point.
**Fix.** Single outboard prefix sum for bending, O(strips). Pre-resolve foil
names to indices in `Ctx::new` (F002 finishes this with family indices).
**Proof.** Bending values bit-identical before and after on the T1 suite
case; `vlm2_t1_100` median improves.
**Budget.** Strip-heavy cases must get faster, not just not slower.

#### F009 · Input validation and finite output (+G06)

**Problem.** A T8 point at zero speed, T6 ranges that cross zero and
nonfinite references all reach the solver and produce NaN rows in `points`.
**Fix.** `PlanePolarDef::validate` expands every derived operating point:
each T8 `(alpha, beta, speed)`, T6 velocity ranges (both ends > 0), T4 range,
and `is_finite` on every scalar, fluid property and reference. One strict
compatibility table, `Method::supports(kind, setting)`, used by validate, the
GUI (greys out options) and the CLI. `solve_point` and every `run_*` send a
point with any nonfinite coefficient to `failed` with a reason.
**Proof.** `f009_t8_zero_speed_rejected`, `g06_no_nonfinite_points`.
**Touches.** `crates/plane/src/analysis.rs` (`Method`, `PlanePolarDef`),
`crates/gui/src/plane_ui.rs`, `crates/cli`.

#### F011 · No silent stability failures

**Problem.** `solve_point`/`run_control` use `.ok()` and `run_stability`
chains `if let Ok`, so a failed trim yields a point with missing stability
and no message.
**Fix.** Warnings on the point plus `stability_error: Option<String>`. The GUI
marks points with missing stability.
**Proof.** `f011_stability_failure_reported`; `no_silent_ok` count drops.

#### F017 · Real 2D iteration counts

**Problem.** `crates/foil/src/analysis.rs:342` writes `iterations: 0`.
**Fix.** Carry the engine's Newton count and recovery attempts (strategy and
substeps) into `FoilPoint`.
**Proof.** `f017_iterations_nonzero` on a converged viscous point.

#### F020 · T3 convergence flag

**Problem.** The glide loop stops after a fixed count and reports whatever it
has.
**Fix.** Return `converged` and the residual `|L_res − W|/W`; an unconverged
point goes to `failed`. F007 rewrites the loop and must keep this proof green.
**Proof.** `f020_t3_unconverged_fails` (a fixture with a stall-limited glide).

#### F023 · LLT settings honesty

**Problem.** LLT accepts beta, ground height, flaps, `fuse_drag`,
`derivatives` and T5–T8, ignores them, and hardcodes CY/roll/yaw to zero.
**Fix.** `Method::supports` rejects those for LLT. CY/roll/yaw become
`None`/absent in the result type instead of 0.
**Proof.** `f023_llt_rejects_unsupported`, `f023_llt_lateral_absent`.

#### F012 · Viscous envelope status

**Problem.** `FoilDb` interpolation returns a bool for "outside", so clamping
in Re is indistinguishable from clamping in CL, and nothing can fail on it.
**Fix.** Return `Envelope { Inside, ClampedRe, ClampedCl, ClampedBoth, Missing }`.
Warnings say "clamped" with counts per kind. `PlanePolarDef.strict_envelope`
fails the point on any clamped strip. Per-strip status is a `u8` in `StripResult`.
**Proof.** `f012_envelope_kinds_reported`, `f012_strict_envelope_fails_point`.
**Budget.** No change on `vlm2_t1_100` (the status write is one byte per strip).

#### F013 · Numerical diagnostics (+G05 part)

**Problem.** `Lu::pivot_ratio` ([linalg.rs:111](../crates/core/src/linalg.rs))
is `dmin/dmax` of the diagonal, not a pivot ratio, and nothing reports
conditioning or residuals.
**Fix.** Rename to `diag_ratio` (serde alias for old files) and document it.
Hager–Higham 1-norm condition estimate using 4–5 solves with the existing
factors: O(N²), off by default, on in `diagnostics=cheap`. Backward error
without a copy of A: re-assemble k sampled rows (default 64 plus all TE/Kutta
rows) and compute ‖r‖/(‖A‖‖x‖+‖b‖) from them, O(kN) influence evaluations.
`diagnostics=full` keeps an f32 copy of A for an exact residual and one step
of iterative refinement. Store near-field vs Trefftz drag and pressure vs wake
lift balance per point (already computed; just keep the difference).
**Proof.** `f013_condition_estimate_within_10x` against a direct SVD on a
200-panel case; `g05_backward_error_reported`.
**Budget.** Default path unchanged on every case; `diagnostics=cheap` ≤ 5%
on `panel_t1_6k`.

#### F014 · Cancellable blocked LU

**Problem.** The single faer `lu_in_place` call cannot be interrupted, so the
README's "within milliseconds" (line 50) is false at scale.
**Fix.** Right-looking blocked LU in `crates/core/src/linalg.rs`: for each
column panel (nb tuned 256–512 on the suite) faer `lu_in_place` on the tall
panel, row swaps, triangular solve for U₁₂, parallel `matmul` for the Schur
update, `prog.check()` between panels. N < ≈2,000 keeps the single call.
Update the README sentence to the measured latency.
**Proof.** `f014_cancel_latency` measured by `cancel_14100` ≤ 100 ms; results
bit-identical to monolithic faer for N < 2,000 and within 1e-13 relative above.
**Budget.** ≤ 3% slower than monolithic faer at N = 14,100 with Rayon on.

#### F001 · Interim 2D bench harness fixes

**Problem.** [run_bench.py](../bench/run_bench.py) keeps the best of N
(`best = min(...)`), never checks the exit code, and compares whatever polar
rows happen to exist.
**Fix.** `subprocess.run(..., timeout=, check=False)`; nonzero exit or timeout
marks the run failed. Delete the polar before each run and require a fresh
mtime. Derive the expected alpha/CL set from `ASEQ`/`CSEQ`; score each point
matched / missing-in-ref / missing-in-foiley / nonfinite. `compare` fails on
an empty common set or coverage below 95% of the reference's converged
points. Keep timing and output of the same repetition; store every sample;
report median and IQR. G01 replaces this script.
**Proof.** `f001_harness_rejects_partial_polar` (a Python test under `bench/`).

#### F015 · Raw-bit evidence for the 2D port

**Problem.** The README claims bit-exact results (lines 6, 33, 246, 247) and
nothing in the repo checks bits.
**Fix.** `RAWDUMP` hook in the instrumented reference build (`refsrc`) writing
CL, CD, CDp, CM, XTR and the BL arrays (UEDG, THET, DSTR, CTAU) as hex `Z16`
per converged point; Foiley writes `f64::to_bits` through the same hook in
`src/` behind a cargo feature. Corpus: the 8 bench cases plus a generated 40
(NACA 4/5 digit, e387, dae11, trips, Mach 0–0.6, Re 5e4–6e6, CSEQ/ASEQ,
inverse-design round trips). Matrix: sequential, `FOILEY_PIPELINE=1`,
`FOILEY_FINE_PAR=1`, SIMD on/off. README restates the claim per configuration
that is bit-identical and gives the measured max ulp for the rest.
**Proof.** `tests/raw_bits.rs::f015_matrix` passes with every configuration
labelled; S6 updated.

#### F002 · Explicit polar families (+G07)

**Problem.** Polars are matched by foil name and Re only; Mach, Ncrit and
trips are implicit, and `reduce` hashes strings per strip per point.
**Fix.** `FamilyKey { foil_geom_hash, polar_type, mach, ncrit, xtr_top,
xtr_bot, flap_deg, mom_ref_x }` with a blake3 hash of normalized coordinates.
`FoilDb` stores `families: Vec<Family>` (polars sorted by Re) plus an index
`HashMap<(foil, FamilyKey), FamilyId>`. `PlanePolarDef.family: FamilySelector
{ mach, ncrit, xtr_top, xtr_bot }`. `Ctx::new` resolves every strip end foil
to a `FamilyId` once and errors on no match or an ambiguous match, listing
candidates. `reduce` indexes `families[id]`: no strings, no hashing.
**Proof.** `f002_polar_order_irrelevant` (the "polar insertion order" print
from the example), `g07_ambiguous_family_errors`.
**Budget.** `vlm2_t1_100` ≤ 0.95× baseline (the lookups were on the hot path).

#### F003 · Flap-deflected polars

**Problem.** Deflected strips use the undeflected section polar.
**Fix.** `Strip.flap_deg` from the mesher. Family lookup brackets on flap
angle (linear between the two nearest deflected families at the same
selector) then on Re. No deflected family: warn by default, fail in strict mode.
**Proof.** `f003_deflected_strip_uses_deflected_polar`.

#### F016 · Provenance and format v2 (+G03)

**Fix.** `RunManifest { foiley_version, git_commit, build_profile,
target_features, thread_config, settings, geometry_hash, mesh_hash,
family_ids+hashes, started_at }` on every 2D and 3D polar; `git_commit` and
features embedded by `build.rs`; hashing once per run. `FoilPolar.settings:
AnalysisSettings` snapshot. `FORMAT_VERSION = 2` in
[project.rs](../crates/app/src/project.rs); v1 migrates (families get
`{mach, ncrit, xtr}` from stored polar fields, a `legacy` marker, hashes
computed on load). Geometry edits invalidate dependent families by hash
mismatch, generalizing the existing `upsert_foil` clearing.
**Proof.** `project_roundtrip` loads, saves, reloads and compares every
`.foiley` under `examples/`, `spitfire/`, `tests/data` (S10);
`f016_manifest_present`.

#### F021 · Coverage-aware auto-polars (+G08)

**Problem.** `crates/app/src/auto.rs` generates polars at hardcoded Mach 0 and
Ncrit 9, checks only for existence, and ignores flap angles.
**Fix.** Required envelope from the actual states: T8 list, T6/T7 ranges,
mass/CG ranges; Re per strip = V·c/ν over min/max chords; CL/alpha span from a
cheap inviscid pre-pass that reuses unit solutions (RHS solves only). Mach,
Ncrit and trips from the selected family. Generate only missing Re stations
and missing flap angles.
**Proof.** `g08_auto_covers_t8_envelope`, `f021_generates_deflected_families`.

#### F018 · Node-preserving 2D input (+G09)

**Fix.** `AnalysisSettings.paneling: Paneling::{Pangen(params), Preserve}`;
`Preserve` takes the `LOAD`/`PCOP` path. Results report the solved node count
and a hash of the node coordinates.
**Proof.** `f018_preserve_keeps_nodes`.

#### F007 · Shared wake-aware evaluator (+F019, G29) `hard`

**Problem.** T2/T3 ignore the relaxed wake (the example prints
`lift/weight` ≠ 1 and `delta_CL` = 0), T1 and T6 give different answers for
equivalent settings, and the wake segment count is hardcoded to 24 at
[analysis.rs:904](../crates/plane/src/analysis.rs).
**Fix.** New `crates/plane/src/evaluator.rs`:

```text
Evaluator { sys: System, wake: WakeModel, cache }
  loads(dir, speed, omega) -> (RawLoads, WakeReport)
WakeModel = Rigid | Relaxed { tol, max_iter, nseg, relax }
```

Every driver (T1–T8 speed balance, glide loop, T6 controls, T7 trim,
derivatives) calls `Evaluator::loads`; direct `sys.loads` calls are removed
and a grep test enforces it. Convergence: max wake-node displacement / c_ref
< `tol` **and** ‖ΔΓ‖∞/‖Γ‖∞ < `tol`, capped at `max_iter`; `WakeReport` holds
iterations and both residuals; not converged → warning, failure in strict
mode. `nseg` and `relax` are exposed in the definition, GUI and CLI.

The design exists for performance:
1. *Speed invariance.* With fixed geometry and inviscid potential flow the
   relaxed wake shape depends on direction and omega/V, not |V|. Cache the
   converged wake per direction; T2/T3 speed iterations cost an O(N) rescale.
2. *Low-rank update.* Moving the wake changes only trailing-edge columns
   (k ≪ N). Each relaxation step solves with the existing LU through
   Sherman–Morrison–Woodbury, O(N²k + k³), refactoring only when the
   capacitance matrix's condition estimate exceeds 1e10. The surface–surface
   block is assembled once.
3. *Warm start.* Each alpha/beta starts from the neighbouring converged wake
   in the existing outward-marching order.
4. *Derivatives* use a frozen converged wake at the reference state by
   default (reusing unit solutions); `relaxed_derivatives` re-relaxes at a
   documented cost.
The panel-method relaxed wake (G29) uses the same evaluator with the panel
system's TE columns.
**Trap.** The Woodbury capacitance matrix goes singular when wake columns
nearly coincide (tip vortices rolling into each other). The condition check
and refactor fallback are not optional.
**Proof.** `f007_t2_lift_equals_weight`, `f019_t1_t6_equivalent`,
`f007_t3_relaxed_differs_from_rigid`, `g29_panel_relaxed_converges`,
`no_direct_sys_loads` (grep).
**Budget.** `vlm2_t2_relaxed_20` ≤ 1.1× T1 relaxed at the same alphas;
relaxed-wake per-iteration cost ≤ 0.2× a full refactorization at 2k panels.

#### F008 · LLT glide balance

**Problem.** The example's constant-CD fixture gives `resultant/weight` ≠ 1.
**Fix.** In `llt.rs`, for T3 iterate V from √(CL²+CD²) with CD including
viscous drag at the current Re (each LLT iteration already re-interpolates).
Converge on |resultant − W|/W < 1e-4 and report failure otherwise. T2 keeps
the CL-only balance.
**Proof.** `f008_llt_glide_resultant_equals_weight`.

#### F005 · Stability models and nonzero beta (+F006, G34) `hard`

**Problem.** Added parasite drag leaves Xu unchanged (example print), and a
nonzero-beta reference state is silently treated as beta = 0.
**Fix.** `StabilityModel::{Inviscid, FullLoad}` and `trim_model` following
it; default `Inviscid` so existing numbers hold. `FullLoad` finite-differences
the `reduce`d loads (viscous strip drag, extra drag, fuselage friction); raw
loads for the 12 perturbations still come from unit solutions, so the cost is
12 `reduce` calls, O(strips) each. For rotation rates `reduce` takes an
optional `omega` and uses per-strip local velocity (V + ω×r) for Re, drag
direction and dynamic pressure. Nonzero beta: stability axes from
`wind_dir(alpha, beta)`; longitudinal and lateral do not decouple, so produce
a coupled 8×8 `[u w q θ v p r φ]` matrix with modes classified by eigenvector
energy fractions. At beta = 0 the coupled matrix is block diagonal and must
reproduce the current 4×4 results exactly. `FullLoad` trim uses the reduced
Cm including viscous section Cm and drag moments.
**Trap.** The 8×8 at beta = 0 reproduces the 4×4 only if the two blocks use
identical perturbation sizes and the same frozen wake; build the 4×4 from the
8×8, not alongside it.
**Proof.** `f005_fullload_xu_moves_with_drag`, `f005_inviscid_xu_unchanged`,
`f006_beta_zero_block_diagonal_matches_4x4`, `g34_phugoid_lanchester` (ζ ≈
1/(√2·L/D) on a clean glider fixture).
**Budget.** `t7_stability_10` unchanged under `Inviscid`; ≤ 1.05× under `FullLoad`.

#### F004 · Multiple bodies (+G23)

**Problem.** `mesh.rs:233` meshes `plane.fuses.first()` only; the example
shows one and two fuselages giving the same body panel count.
**Fix.** Iterate `plane.fuses`. Each wing root is trimmed against the body it
intersects (bounding-box prefilter, then the existing intersection routine);
carry-through segments per wing–body pair; body panels get `body_id`.
Body–body proximity and intersection check (BVH over triangles) with a clear
error on interpenetration.
**Proof.** `f004_two_pods_double_panels_mirror_loads`,
`f004_single_body_bit_identical`, `g23_twin_boom_fixture`.
**Budget.** Mesh time linear in body count; no change for single-body planes.

#### F010 · AVL and interchange fidelity (+G45)

**Problem.** The AVL export writes `YDUPLICATE` for asymmetric wings (the
example shows the right 0012 and left 0009 collapsing to one side).
**Fix.** `YDUPLICATE 0.0` only when the wing is two-sided, centred and
mirror-identical in foils, twist and controls; otherwise separate `SURFACE`s
per side. Offset two-sided wings always get an explicit left half. One
`CONTROL` per actual flap with unique names; ailerons `SgnDup −1`, elevators
and flaps `+1`; gains follow the definition's mixing. Add an AVL geometry
importer (sections, `NACA`/`AFILE`, `YDUPLICATE`, `CONTROL`) and check round
trips numerically (LE, chord, incidence, foil per side within 1e-9).
`.xfl`/`.fl5`: parser against hand-written minimal files and the flow5 API
example; import 3D polar definitions where they map; every import returns a
`ConversionReport` shown in GUI and CLI. Last sub-step: real XFLR5 6.x and
flow5 7.58 files in `tests/data/interchange/` when the user supplies them
(the plan's one wanted stop; land everything else first).
**Proof.** `f010_asymmetric_export_two_surfaces`, `f010_avl_roundtrip_1e9`,
`g45_conversion_report_lists_skips`.

#### F024 · GUI CAD import as a job

**Problem.** [plane_ui.rs:1206](../crates/gui/src/plane_ui.rs) reads and
tessellates STL/STEP on the UI thread and swallows errors.
**Fix.** Move reading and tessellation into the job manager (owned input
bytes, progress, cancel). Notes and errors in the result panel.
**Proof.** `f024_cad_import_is_job` (job manager test with a bad file
producing a visible error).

#### F022 · Pure-Rust trimmed STEP (+G26) `hard`

**Problem.** `cad.rs` ignores trimming, so faces with holes and B-spline
surfaces come out wrong or not at all.
**Fix.** Split into `cad/step_parse.rs`, `cad/surfaces.rs`, `cad/trim.rs`,
`cad/tessellate.rs`, `cad/check.rs`. Surfaces: `PLANE`, `CYLINDRICAL_`,
`CONICAL_`, `SPHERICAL_`, `TOROIDAL_SURFACE`, `SURFACE_OF_REVOLUTION`,
`SURFACE_OF_LINEAR_EXTRUSION`, rational and non-rational `B_SPLINE_SURFACE`
including `_WITH_KNOTS` complex entities. Topology: `ADVANCED_FACE` →
`FACE_OUTER_BOUND`/`FACE_BOUND` → `EDGE_LOOP` → `ORIENTED_EDGE` →
`EDGE_CURVE` (`LINE`, `CIRCLE`, `ELLIPSE`, `B_SPLINE_CURVE`); `PCURVE`s when
present, else invert 3D edge points onto the surface (Newton with a
coarse-grid seed). Tessellation: sample each edge once by chordal deviation
and angle tolerance and share vertices between adjacent faces (watertight by
construction); constrained Delaunay of the trimmed parameter domain with
holes via `spade`, refined by curvature; faces in parallel with Rayon.
Checks: manifold/open-edge count, orientation, signed volume, area vs a
high-density reference, self-intersection (BVH). Unsupported entities listed
in import notes; `strict` fails on any omission. Corpus: existing test plus
FreeCAD exports (box with hole, cylinder fuselage, revolved nacelle, lofted
B-spline fuselage, multi-shell) with FreeCAD-reported area and volume.
**Trap.** Edge inversion without `PCURVE`s fails at surface seams
(cylinders, revolutions); seam handling needs its own test before the corpus.
**Proof.** `f022_corpus_area_volume_within_0_1pct`,
`g26_watertight_and_oriented`.
**Budget.** `step_fuselage` ≤ current time at equal triangle count.

#### G25 · Imported lifting meshes

**Fix.** TE detection on STL/STEP wings (sharp-edge dihedral threshold plus
user override), lifting/nonlifting classification, wake attachment, and a
validation report using the G26 checks.
**Proof.** `g25_imported_wing_matches_native_cl` within 2% on a STEP export
of a native wing.

#### G04 · Instrumentation

**Fix.** Phase timers at phase boundaries only (`Instant`), peak RSS, Newton
and recovery counts from F017, all written into the RunManifest.
**Proof.** `g04_manifest_has_phase_times`.
**Budget.** Not measurable in the perf gate.

#### G01 · Benchmark runner

**Fix.** `foiley-bench run <suite>` as a Rust orchestrator. Cases are
versioned TOML files in `bench/cases/` with geometry, conventions and
per-solver adapter inputs. Every run writes raw output, a status per
requested point (ok / nonconverged / rejected / crash / timeout), timing
phases, machine and solver versions, as JSON lines. Adapters, one module
each, subprocess or file driven:

| Solver | Adapter notes |
|---|---|
| XFOIL 6.97 (sp/dp) | existing builds; stdin scripts |
| XFOIL 6.99 | separate build, labelled separately |
| AVL | official Windows binary; generated `.avl`/`.run`; parses `FT`/`ST` |
| flow5 7.58 | pinned release; C++ API example or batch XML |
| VSPAERO 7 | pinned OpenVSP; `.vspaero` files; thin and thick modes |
| MachLine | pinned commit, gfortran; JSON input |
| MachUpX | pinned commit in a locked Python venv |
| FLOWPanel.jl | pinned Manifest; warm-up excluded from kernel time, kept in end-to-end |

Toolchains pinned in `bench/toolchains.md` (versions, hashes, install
steps); native solvers in `bench/tools/` (gitignored). Scorecards per track:
verification, physical validation, performance (cold single point, warm
20- or 100-point sweep, runtime to a set error, unknowns, memory), workflow.
Each track has its own speed–accuracy curve; never one averaged number.
First milestone: a matched 2D polar, a thin rectangular-wing mesh study and a
thick-wing Cp study across ≥ 5 solvers, published from the runner with raw
outputs saved, and the README claims in `bench/claims.md` replaced by runner
results (S8).
**Proof.** `g01_runner_scores_first_milestone` (runs the three studies in
CI-skippable mode); `bench/claims.md` has no "unbacked" rows.

#### G02 · Experimental validation

**Fix.** UIUC low-speed airfoil data (E387 and others) and the NACA
straight/35° swept wing pressure report, via download scripts with SHA-256
checksums (§2). Geometry, conditions and corrections reconstructed as case
files with stated uncertainty.
**Proof.** `g02_e387_cl_within_experimental_band`.

#### G05 · Refinement studies (+G44)

**Fix.** `foiley plane refine`: 3-level mesh and wake refinement with
Richardson extrapolation and observed order for CL, CDi, Cm, Cp; reports the
grid convergence index. Adaptive refinement from per-panel error indicators
(doublet-gradient jumps).
**Proof.** `g05_observed_order_near_2_on_rectangular_wing`.

#### G40 · Memory-scalable dense systems

**Fix.** Memory budget check before assembly (N²·8 B plus workspace vs
available RAM); refuse or switch strategy with a message, never OOM.
Mixed-precision LU: factor in f32, then f64 iterative refinement against the
f64 operator applied matrix-free (from G41's H-matrix or re-assembled rows),
stopping at backward error < 1e-12, falling back to f64 LU when refinement
stalls (condition > ≈1e6). Right-preconditioned restarted GMRES with a
block-Jacobi near-field LU (per wing or body); the six unit solutions use
block GMRES or deflated restarts.
**Proof.** `g40_mixed_precision_matches_f64_1e10` on the suite;
`g40_refuses_oversize_system`.
**Budget.** `trilinear_14100` ≤ 0.6× time and ≤ 0.55× RSS (S4).

#### G41 · Hierarchical acceleration `hard`

**Fix.** H-matrix with ACA+ on a cluster tree of panel centroids,
admissibility η ≈ 1–2, per-block rank by tolerance: O(N log N) storage and
matvec, used as the GMRES operator and for LU-free large meshes. Replaces the
point far-field approximation with this error-controlled version
(`integration_tol`). Barnes–Hut or FMM for vortex-segment and particle sums
(G30–G32) and off-body velocity (G46). Crossover chosen by measured phase
timings, not fixed N.
**Trap.** ACA on doublet kernels with mixed signs can terminate early on a
zero row; use ACA+ with the partial-pivoting restart and verify every block
against a dense sample during development.
**Proof.** `g41_matches_dense_1e6_cl_at_14100`.
**Budget.** `trilinear_14100` RSS ≤ 25% of dense (S4 target), time ≤ dense.

#### G43 · Adjoints and AD

**Fix.** Discrete 3D adjoint: transpose solve with the existing LU or
H-matrix. Geometry sensitivities as hand-differentiated influence kernels,
checked against finite differences and complex step on small meshes.
**Proof.** `g43_adjoint_gradient_matches_fd_1e6`.

#### G28 · Freestream-aligned rigid wake

**Fix.** Option next to body-x, default in benchmark cases; rebuild per
direction avoided by F007's Woodbury path.
**Proof.** `g28_aligned_wake_cl_vs_avl_within_1pct`.

#### G35 · General constrained trim

**Fix.** Newton on (α, β, controls, V or thrust, bank) with force and moment
targets: level, climb, coordinated turn. AVL-style run cases. Jacobian from
unit solutions and control derivatives.
**Proof.** `g35_coordinated_turn_residual_1e8`.

#### G37 · 3D hinge loads

**Fix.** Integrate panel pressure about the actual hinge line per control
for hinge moment coefficients. O(control panels).
**Proof.** `g37_hinge_moment_sign_and_magnitude_vs_2d`.

#### G22 · VLM body model

**Fix.** AVL-style slender-body source and doublet lines for fuselages and
nacelles; adds O(N_body·N) columns.
**Proof.** `g22_body_lift_curve_vs_avl`.

#### G24 · Mixed thin/thick systems

**Fix.** Thin VLM surfaces and thick panels in one system (fins on a thick
body); one assembly; block structure reused by the G40 preconditioner.
**Proof.** `g24_fin_on_body_matches_all_thick_within_3pct`.

#### G27 · 3D compressibility

**Fix.** Prandtl–Glauert/Göthert transform at low subsonic Mach, limit
M ≤ 0.6 with a warning above 0.5; transformed systems cached per Mach.
**Proof.** `g27_pg_slope_scaling`.

#### G21 · Nonlinear VLM/panel viscous coupling

**Fix.** Decambering / α-coupling iteration: section polar CL(α_eff) drives
an RHS correction; each iteration is RHS solves with the factored LU, no
refactorization.
**Proof.** `g21_stall_rollover_vs_machupx`.

#### G20 · General nonlinear LLT

**Fix.** Phillips-style numerical LLT: multiple surfaces, sweep, dihedral,
winglets, sideslip; Newton on circulation with an analytic Jacobian; dense LU
of the station count. MachUpX is the benchmark.
**Proof.** `g20_swept_wing_vs_machupx_within_2pct`.

#### G36 · Broader flight dynamics

**Fix.** Non-level and nonzero-rate reference states; nonlinear 6-DOF RK4 on
an aero model built from unit solutions plus polar lookup; no per-step solves.
**Proof.** `g36_trimmed_state_is_rk4_fixed_point`.

#### G46 · Field output

**Fix.** Binary appended VTU with Cp, μ, forces; off-body probes and
streamlines through G41's evaluator; per-component force/moment breakdown;
streamlines in parallel.
**Proof.** `g46_vtu_roundtrip_parse`.

#### G47 · C and Python bindings

**Fix.** C ABI (`foiley-capi`, cbindgen) and Python (PyO3 + maturin, numpy
zero-copy). Bindings hold the solver handle so warm unit-solution reuse
survives across calls.
**Proof.** `g47_python_sweep_matches_cli`.

#### G32 · Actuator disks and propellers

**Fix.** Actuator disk as a source/vortex sheet giving a slipstream tube with
a momentum-theory velocity increment, then blade-element props with
prescribed helical wakes; wing–tail slipstream interaction. VSPAERO/DUST are
the benchmarks.
**Proof.** `g32_disk_thrust_momentum_balance`.

#### G30 · UVLM / unsteady panels

**Fix.** Time stepping with wake shedding and convection, added mass, moving
controls, gusts. Rigid motion in body axes keeps the surface matrix
constant: factor once, O(N²) per step; wake–wake and wake–body induction
through G41's tree code. PteraSoftware is the benchmark.
**Proof.** `g30_theodorsen_lift_deficiency`.

#### G31 · Vortex particles

**Fix.** Classic VPM with regularized core and particle strength exchange,
tree-code accelerated, as a wake option for rotors. Separate physics track
from potential flow.
**Proof.** `g31_ring_vortex_self_induction`.

#### G19 · 2D hinge moments

**Fix.** Expose the engine's hinge quantities in the service with sign
conventions. Validated against XFOIL `HINC`.
**Proof.** `g19_hinc_matches_xfoil`.

#### G13 · Post-stall extension

**Fix.** Viterna–Corrigan extrapolation stored as a separate polar segment
with provenance, never mixed silently with solved points.
**Proof.** `g13_extrapolated_segment_flagged`.

#### G11 · 2D ground and cascade

**Fix.** Image panels for ground; periodic kernel for cascades. JavaFoil is
the check.
**Proof.** `g11_ground_effect_cl_increase_vs_javafoil`.

#### G10 · Multi-element 2D

**Fix.** Multiple contours in the panel system, a BL per element, wakes
passing over downstream elements (no confluence model at first, documented).
VIIFLOW/MSES are the benchmarks.
**Proof.** `g10_two_element_cl_vs_mses_within_5pct`.

#### G16 · 2D sensitivities

**Fix.** Discrete adjoint of the existing Newton system (the Jacobian is
already built); transpose solve of the block system.
**Proof.** `g16_adjoint_matches_fd`.

#### G17 · Multipoint optimization

**Fix.** Weighted or constrained multipoint objectives using G16 gradients,
with geometry constraints (thickness, area, TE angle); candidates evaluated
in parallel as today.
**Proof.** `g17_improves_weighted_ld_under_constraints`.

#### G18 · Surrogate screening

**Fix.** Small MLP trained offline on Foiley batch meshes; pure-Rust
inference; domain and uncertainty checks; candidates always re-verified by
the solver.
**Proof.** `g18_surrogate_never_replaces_solver_result`.

#### G12 · Rotational corrections

**Fix.** RFOIL-style rotational BL corrections behind a flag, with documented
closures.
**Proof.** `g12_rotational_flag_changes_post_stall_only`.

#### G14 · Time-accurate 2D

**Fix.** Unsteady panel method with shed wake, plus a Beddoes–Leishman
dynamic-stall model as an engineering option.
**Proof.** `g14_wagner_response`.

#### G15 · Shock-capable 2D (spike with scored exit)

**Fix.** A 24-hour spike: full-potential outer flow coupled to the existing
BL on a single transonic NACA 0012 case. Exit criterion, decided by the
number: continue only if the spike reproduces the MSES shock location within
5% chord at M 0.75 inside the budget. Otherwise record the result in
`CHANGELOG.md` and mark the item `[x] spiked, not pursued`.
**Proof.** The spike report with the measured shock location.

#### G38 · Structural coupling

**Fix.** Timoshenko beam FE per wing (wingbox properties), consistent load
transfer from strips or panels, static aeroelastic equilibrium by Newton;
small deformations reuse the aero factorization through a linearized
geometry update, falling back to rebuild. Then p-k flutter and gust response.
OpenAeroStruct is the benchmark.
**Proof.** `g38_static_deflection_vs_openaerostruct`.

#### G39 · Frequency-domain DLM

**Fix.** Doublet-lattice AICs at reduced frequency and Mach, sharing the VLM
mesh. PanelAero is the benchmark.
**Proof.** `g39_aic_vs_panelaero`.

## 5. Changed-numbers protocol

No legacy switch. When an item changes a result that a golden test or the
README states:

1. `CHANGELOG.md` gets a row: item, quantity, fixture, before, after, reason.
2. The golden test is updated in the same PR.
3. `bench/claims.md` is updated if a README claim moves, and the README
   sentence is rewritten to the measured value.

## 6. Traceability

Every review ID is an item or folded into one. Merged IDs: F019 and G29 into
F007; F006 and G34 into F005; G06 into F009; G07 into F002; G03 into F016;
G08 into F021; G09 into F018; G23 into F004; G45 into F010; G26 into F022;
G05's cheap diagnostics into F013, the rest into G05; G44 into G05. G33
(marine/hydrofoil) and G42 (GPU/distributed) are dropped by the user's
decision; the phased plan omitted G33 without saying so. Everything else
keeps its own ID in the queue.
