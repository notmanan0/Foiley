# Foiley fix list

Revision 3, 4 November 2026, against `5395ce7`. The findings (F001–F024) from
the 4 November 2026 code review, as a checklist. Each item states its own
problem, fix and the tests that prove it.

## How to work this list

Asked to "do N items": take the first N unticked boxes in order, skipping any
whose **after** items are not ticked yet. For each one:

1. **Look first.** Open the files the spec names and grep every caller of
   the functions it names. Many findings come from drivers that duplicate each
   other (T1–T8 speed balance, glide loop, controls, trim, derivatives), so
   expect the same bug in more places than the spec lists, and fix them all.
2. **Write the tests first.** The tests go in one file per item:
   `crates/<crate>/tests/<id>.rs` (for example `crates/plane/tests/f009.rs`),
   `tests/<id>.rs` for the root crate, or `bench/tests/test_<id>.py`.
   Use the test names the spec gives.
3. **Fix it.**
4. **Check.** The item's tests pass and are not `#[ignore]`d.
   `cargo test --release --workspace` and `python -m pytest bench/tests` are
   green.
5. **Record.** If a number in a golden test or the README changed, add a row
   to `CHANGELOG.md`: item, quantity, fixture, before, after, why. Update the
   golden test, and the README sentence if one moved.
6. **Tick the box** in this file, in the same commit.

One branch and one commit (or PR) per item: `item/<ID>-<slug>`. Don't stop
between items to summarize; keep going until N are done or an item truly
needs something only the user can supply (only F010c does).

If the spec and the code disagree and both readings are defensible, pick one,
say which in the commit message, and carry on.

**Code rules.** Hot loops (per point, strip or panel) stay allocation-free,
with no string hashing in them. Diagnostics are off by default and cost at
most O(N²). Nothing adds an O(N³) step to the default path.

**Reference files.** The GPL reference sources and builds (`refsrc/`,
`xfoil-ref/`, `build-ref-*`) stay gitignored and are found through
`FOILEY_REF_ROOT`. A test that needs them and can't find them prints
`SKIPPED(ref)`. With `FOILEY_REQUIRE_REF=1` it fails instead.

## Checklist

- [x] **S1** · Track test fixtures and gitignore the big files
- [x] **F017** · Real 2D iteration counts
- [x] **F020** · T3 convergence flag
- [x] **F001** · 2D bench harness fixes
- [x] **F009** · Input validation and finite output
- [x] **F011** · No silent stability failures
- [x] **F013** · Numerical diagnostics
- [x] **F023** · LLT settings honesty — after F009
- [x] **F012** · Viscous envelope status
- [x] **F014** · Cancellable blocked LU
- [ ] **F008** · LLT glide balance — after F020
- [ ] **F018** · Node-preserving 2D input
- [ ] **F024** · GUI CAD import as a job
- [ ] **F002** · Explicit polar families
- [x] **F007** · One wake-aware evaluator for every driver (+F019) — after F020
- [ ] **F004a** · Multi-body meshing
- [ ] **F010a** · AVL export fidelity and importer
- [ ] **F003** · Flap-deflected polars — after F002, F004a
- [x] **F005a** · `FullLoad` stability model — after F007
- [ ] **F016** · Provenance and project format v2 — after F002
- [ ] **F015a** · Raw-bit hooks and harness
- [ ] **F021** · Coverage-aware auto-polars — after F003
- [ ] **F004b** · GUI editing of every body — after F004a, F024
- [ ] **F010b** · `.xfl`/`.fl5` import and conversion report
- [ ] **F015b** · Raw-bit corpus and configuration matrix — after F015a
- [x] **F005b** · Coupled 8×8 stability and nonzero beta (+F006) — after F005a
- [ ] **F022a** · STEP parse, topology, analytic surfaces, seams
- [ ] **F022b** · B-spline surfaces and curves — after F022a
- [ ] **F022c** · Watertight tessellation — after F022a
- [ ] **F022d** · Mesh checks and CAD corpus — after F022b, F022c
- [ ] **F010c** · Real XFLR5/flow5 files — after F010b; **waits on the user's files**

F004, F005, F010, F015 and F022 are split into lettered parts; the finding is
fixed when all its parts are ticked. F006 is fixed by F005b and F019 by F007.

**Not in this list (for later):** the perf suite and speed gate, the free
performance item (O(strips) bending moments), wake-relaxation speedups
(Woodbury updates, speed-invariant wake cache), the panel-method relaxed
wake, and the capability gaps (G01–G47) such as reference-solver adapters,
H-matrices, UVLM, propellers and aeroelasticity.

---

## Specs

#### S1 · Track test fixtures and gitignore the big files

**Fix.**
- Commit `crates/plane/examples/benchmark_review.rs`,
  `spitfire/Spitfire_MkI_geometry.foiley` and `spitfire/make_spitfire.py`.
- **Before staging anything**, gitignore `spitfire/Spitfire_MkI_6DoF.foiley`
  (844 MB), `spitfire/results/`, `spitfire/run_log.txt`,
  `spitfire/run_err.txt` and `spitfire/spitfire_6dof_job.json`. Check the
  staged sizes before committing.
- Copy `e387.dat` and `dae11.dat` from `xfoil-ref/runs/` to
  `tests/data/foils/` (GPL, same licence as Foiley).
- The flow5 STL/STEP used by
  [validation.rs:263](../crates/plane/tests/validation.rs): copy into
  `tests/data/cad/` if `refsrc/flow5`'s licence allows redistribution,
  otherwise read them through `FOILEY_REF_ROOT`. Either way the test stops
  silently skipping and follows the reference-files rule above.
- Turn the observations printed by `benchmark_review.rs` into `#[ignore]`d
  tests in their items' files, using the test names the specs below give.
  Each item removes the `#[ignore]` when it fixes the bug. Keep the example
  for printing.

**Done when.** A clean clone has every file a test reads, apart from the
`FOILEY_REF_ROOT` set, and `git status` shows none of the big files.

#### F017 · Real 2D iteration counts

**Problem.** [crates/foil/src/analysis.rs:342](../crates/foil/src/analysis.rs)
writes `iterations: 0`.
**Fix.** Carry the engine's Newton count (the `iter` of the `viscal` loop,
[src/xfoil.rs:842](../src/xfoil.rs)) and recovery attempts (strategy and
substeps) into `FoilPoint`, on the sweep path too. Storing a counter must not
change any number.
**Tests.** `f017_iterations_nonzero` on a converged viscous point.

#### F020 · T3 convergence flag

**Problem.** The glide loop
([analysis.rs:911](../crates/plane/src/analysis.rs)) stops after 30 passes and
reports whatever it has.
**Fix.** Return `converged` and the residual `|L_res − W|/W`. An unconverged
point goes to `failed`. F007 rewrites the loop and must keep this test green.
**Tests.** `f020_t3_unconverged_fails` (a fixture with a stall-limited glide).

#### F001 · 2D bench harness fixes

**Problem.** [run_bench.py](../bench/run_bench.py) keeps the best of N runs
(`best = min(...)`, line 60), never checks the exit code, and compares
whatever polar rows happen to exist. It already deletes the polar before each
run (line 53).
**Fix.**
- Run with `subprocess.run(..., timeout=, check=False)`. A nonzero exit or a
  timeout marks the run failed.
- Require that the polar exists after the run and is newer than its start.
- Derive the expected alpha/CL set from `ASEQ`/`CSEQ`, and score each point
  matched, missing-in-ref, missing-in-foiley or nonfinite.
- `compare` fails on an empty common set, or on coverage below 95% of the
  reference's converged points.
- Keep the timing and output of the same repetition. Store every sample and
  report median and IQR.

**Tests.** `bench/tests/test_f001.py::f001_harness_rejects_partial_polar`.

#### F009 · Input validation and finite output

**Problem.** A T8 point at zero speed, T6 ranges that cross zero and
nonfinite references all reach the solver and produce NaN rows in `points`.
**Fix.**
- `PlanePolarDef::validate` expands every derived operating point:
  - each T8 `(alpha, beta, speed)`;
  - T6 velocity ranges (both ends > 0);
  - the T4 range;
  - `is_finite` on every scalar, fluid property and reference.
- One strict compatibility table, `Method::supports(kind, setting)`, used by
  validate, the GUI (which greys out options) and the CLI.
- `solve_point` and every `run_*` send a point with any nonfinite
  coefficient to `failed`, with a reason.

**Tests.** `f009_t8_zero_speed_rejected`, `f009_t6_range_crossing_zero_rejected`,
`g06_no_nonfinite_points`, `f009_cli_and_gui_share_supports_table`.
**Look at.** `crates/plane/src/analysis.rs` (`Method`, `PlanePolarDef`),
`crates/gui/src/plane_ui.rs`, `crates/cli`.

#### F011 · No silent stability failures

**Problem.** `solve_point`/`run_control` use `.ok()` and `run_stability`
chains `if let Ok`, so a failed trim yields a point with missing stability
and no message.
**Fix.**
- Warnings on the point plus `stability_error: Option<String>`. The GUI
  marks points with missing stability.
- Add `crates/plane/tests/no_silent_ok.rs`. It greps the 3D driver functions
  (`solve_point`, `run_polar`, `run_control`, `run_stability`, the glide loop)
  for `.ok()`, `if let Ok` and `unwrap_or` applied to solver results, and
  asserts the count is at most a constant. Set the constant to the count after
  this fix; later items may only lower it.

**Tests.** `f011_stability_failure_reported`, `no_silent_ok`.

#### F013 · Numerical diagnostics

**Problem.** `Lu::pivot_ratio` ([linalg.rs:111](../crates/core/src/linalg.rs))
is `dmin/dmax` of the diagonal, not a pivot ratio, and nothing reports
conditioning or residuals.
**Fix.**
- **Rename.** `pivot_ratio` becomes `diag_ratio`, with a serde alias for old
  files, and is documented.
- **`diagnostics=cheap`.** A Hager–Higham 1-norm condition estimate, using
  4–5 solves with the existing factors: O(N²) and off by default. Also a
  sampled backward error: re-assemble k rows in f64 (default 64, plus all
  TE/Kutta rows) and compute ‖r‖/(‖A‖‖x‖+‖b‖) from them.
- **`diagnostics=full`.** The exact residual r = b − A·x, with A
  **re-assembled in f64 row blocks** (no stored copy). Then one step of
  iterative refinement, kept only if ‖r‖ decreases. Never use a rounded copy
  of A: its residual can read zero for a wrong x.
- **Balances.** Store the near-field vs Trefftz drag and the pressure vs
  wake lift balance per point. They are already computed; keep the
  differences.

**Tests.**
- `f013_condition_estimate_within_10x`: against the exact
  κ₁(A) = ‖A‖₁‖A⁻¹‖₁ from an explicit inverse on a 200-panel case (not an
  SVD, which gives the 2-norm).
- `f013_full_residual_uses_f64_operator`: a case where an f32-rounded A gives
  a zero residual and the f64 residual is nonzero.
- `g05_backward_error_reported`.

#### F023 · LLT settings honesty

**Problem.** LLT accepts beta, ground height, flaps, `fuse_drag`,
`derivatives` and T5–T8, ignores them, and hardcodes CY/roll/yaw to zero.
**Fix.** `Method::supports` (from F009) rejects those for LLT. CY/roll/yaw
become `None`/absent in the result type instead of 0.
**Tests.** `f023_llt_rejects_unsupported`, `f023_llt_lateral_absent`.

#### F012 · Viscous envelope status

**Problem.** `FoilDb` interpolation returns a bool for "outside", so clamping
in Re is indistinguishable from clamping in CL, and nothing can fail on it.
**Fix.** Return `Envelope { Inside, ClampedRe, ClampedCl, ClampedBoth, Missing }`.
Warnings say "clamped" with counts per kind. `PlanePolarDef.strict_envelope`
fails the point on any clamped strip. Per-strip status is a `u8` in
`StripResult`.
**Tests.** `f012_envelope_kinds_reported`, `f012_strict_envelope_fails_point`.

#### F014 · Cancellable blocked LU

**Problem.** The single faer `lu_in_place` call cannot be interrupted, so the
README's "within milliseconds" (line 50) is false at scale.
**Fix.** A right-looking blocked LU in `crates/core/src/linalg.rs`. For each
column panel (nb in 128–256):
1. faer `lu_in_place` on the tall panel;
2. row swaps;
3. a triangular solve for U₁₂;
4. the Schur update **split into column tiles of width t**, each a parallel
   faer `matmul`, with `prog.check()` between tiles.

Pick nb and t so the largest uninterruptible unit takes ≤ 50 ms at
N = 14,100. Checks only between panels can't meet 100 ms: one 256-column
panel's trailing update is ≈ 10¹¹ flops (≈ 0.3 s). Below N ≈ 2,000 keep the
single call. Update the README sentence to the measured latency. Note the
blocked vs monolithic factor time in the commit message.
**Tests.**
- `f014_cancel_latency`: on the TriLinear fine wing (14,100 unknowns),
  cancel 1 s into the factorization; p95 of 10 cancels ≤ 100 ms.
- `f014_backward_error_matches_faer`: above 2,000 unknowns, the backward error
  ‖b − Ax‖/(‖A‖‖x‖) is ≤ 2× that of monolithic faer, and CL/CM agree within
  1e-10 relative.

#### F008 · LLT glide balance

**Problem.** The example's constant-CD fixture gives `resultant/weight` ≠ 1.
**Fix.** In `llt.rs`, for T3 iterate V from √(CL²+CD²), with CD including
viscous drag at the current Re (each LLT iteration already re-interpolates).
Converge on |resultant − W|/W < 1e-4 and report failure otherwise. T2 keeps
the CL-only balance.
**Tests.** `f008_llt_glide_resultant_equals_weight`.

#### F018 · Node-preserving 2D input

**Fix.** `AnalysisSettings.paneling: Paneling::{Pangen(params), Preserve}`;
`Preserve` takes the `LOAD`/`PCOP` path. Results report the solved node count
and a hash of the node coordinates.
**Tests.** `f018_preserve_keeps_nodes`, `g09_node_hash_reported_and_stable`.

#### F024 · GUI CAD import as a job

**Problem.** [plane_ui.rs:1206](../crates/gui/src/plane_ui.rs) reads and
tessellates STL/STEP on the UI thread and swallows errors.
**Fix.** Move reading and tessellation into the job manager (owned input
bytes, progress, cancel). Notes and errors go in the result panel.
**Tests.** `f024_cad_import_is_job` (a job manager test in
`crates/app/tests/f024.rs` where a bad file produces a visible error).

#### F002 · Explicit polar families

**Problem.** Polars are matched by foil name and Re only. Mach, Ncrit and
trips are implicit, and `reduce` hashes strings per strip per point
(`strip_viscous`, [analysis.rs:596](../crates/plane/src/analysis.rs)).
**Fix.**
- **Key.** `FamilyKey { foil_geom_hash, polar_type, mach, ncrit, xtr_top,
  xtr_bot, flap_deg, mom_ref_x }`, with a blake3 hash of the normalized
  coordinates.
- **Storage.** `FoilDb` stores `families: Vec<Family>` (polars sorted by Re)
  plus an index `HashMap<(foil, FamilyKey), FamilyId>`.
- **Selection.** `PlanePolarDef.family: FamilySelector { mach, ncrit,
  xtr_top, xtr_bot }`. `Ctx::new` resolves every strip-end foil to a
  `FamilyId` once, and errors on no match or an ambiguous match, listing
  the candidates.
- **Lookup.** `reduce` indexes `families[id]`: no strings, no hashing.

**Tests.** `f002_polar_order_irrelevant` (the "polar insertion order" print
from the example), `g07_ambiguous_family_errors`.

#### F007 · One wake-aware evaluator for every driver (+F019)

**Problem.** T2/T3 ignore the relaxed wake (the example prints
`lift/weight` ≠ 1 and `delta_CL` = 0). T1 and T6 give different answers for
equivalent settings (F019). The wake segment count is hardcoded to 24
([analysis.rs:904](../crates/plane/src/analysis.rs)).
**Fix.** A new `crates/plane/src/evaluator.rs`:

```text
Evaluator { sys: System, wake: WakeModel, cache }
  loads(dir, speed, omega) -> (RawLoads, WakeReport)
WakeModel = Rigid | Relaxed { tol, max_iter, nseg, relax }
```

- **Routing.** Every driver (T1–T8 speed balance, glide loop, T6 controls,
  T7 trim, derivatives) calls `Evaluator::loads`. Direct `sys.loads` calls
  are removed, and a grep test enforces it.
- **Convergence.** Stop when max wake-node displacement / c_ref < `tol`
  **and** ‖ΔΓ‖∞/‖Γ‖∞ < `tol`, capped at `max_iter`. `WakeReport` holds the
  iteration count and both residuals. Not converged gives a warning, and a
  failure in strict mode.
- **Settings.** `nseg` and `relax` are exposed in the definition, GUI and CLI.
- **Relaxation** refactors each step for now.
- **Derivatives** use a frozen converged wake at the reference state by
  default, reusing unit solutions. `relaxed_derivatives` re-relaxes, at a
  documented cost.

**Trap.** Drivers duplicate each other. Route them all before deleting any
`sys.loads` call, or one driver silently keeps the rigid path.
**Tests.** `f007_t2_lift_equals_weight`, `f019_t1_t6_equivalent`,
`f007_t3_relaxed_differs_from_rigid`, `f007_wake_report_residuals`,
`no_direct_sys_loads` (grep); `f020_t3_unconverged_fails` stays green.

#### F004a · Multi-body meshing

**Problem.** `mesh.rs:233` meshes `plane.fuses.first()` only. The example
shows one and two fuselages giving the same body panel count.
**Fix.**
- Iterate `plane.fuses`.
- Trim each wing root against the body it intersects: a bounding-box
  prefilter, then the existing intersection routine.
- Add carry-through segments per wing–body pair.
- Body panels get a `body_id`.
- Check body–body proximity and intersection with a BVH over triangles, and
  give a clear error on interpenetration.

**Tests.** `f004_two_pods_double_panels_mirror_loads`,
`f004_single_body_bit_identical`, `g23_twin_boom_fixture`.

#### F010a · AVL export fidelity and importer

**Problem.** The AVL export ([export.rs](../crates/plane/src/export.rs))
writes `YDUPLICATE` for asymmetric wings. The example shows the right 0012
and left 0009 collapsing to one side.
**Fix.**
- Write `YDUPLICATE 0.0` only when the wing is two-sided, centred, and
  mirror-identical in foils, twist and controls. Otherwise write a separate
  `SURFACE` per side.
- Offset two-sided wings always get an explicit left half.
- Write one `CONTROL` per actual flap, with unique names. Ailerons get
  `SgnDup −1`, elevators and flaps `+1`, and gains follow the definition's
  mixing.
- Add an AVL geometry importer: sections, `NACA`/`AFILE`, `YDUPLICATE`,
  `CONTROL`.

**Tests.** `f010_asymmetric_export_two_surfaces`; `f010_avl_roundtrip_1e9`
(LE, chord, incidence and foil per side within 1e-9);
`f010_control_signs_and_names`.

#### F003 · Flap-deflected polars

**Problem.** Deflected strips use the undeflected section polar.
**Fix.** The mesher sets `Strip.flap_deg`. Family lookup brackets on flap
angle first (linear between the two nearest deflected families at the same
selector), then on Re. With no deflected family: warn by default, fail in
strict mode.
**Tests.** `f003_deflected_strip_uses_deflected_polar`, `f003_strict_fails_without_family`.

#### F005a · `FullLoad` stability model

**Problem.** Added parasite drag leaves Xu unchanged (example print).
**Fix.**
- **Model.** Add `StabilityModel::{Inviscid, FullLoad}`, with `trim_model`
  following it. The default is `Inviscid`, so existing numbers hold.
- **Derivatives.** `FullLoad` finite-differences the `reduce`d loads: viscous
  strip drag, extra drag, fuselage friction. Raw loads for the 12
  perturbations still come from unit solutions, so the cost is 12 `reduce`
  calls.
- **Rotation rates.** `reduce` takes an optional `omega` and uses per-strip
  local velocity (V + ω×r) for Re, drag direction and dynamic pressure.
- **Trim.** `FullLoad` trim uses the reduced Cm, including viscous section
  Cm and drag moments.

**Tests.** `f005_fullload_xu_moves_with_drag`, `f005_inviscid_xu_unchanged`,
`f005_fullload_trim_includes_viscous_cm`.

#### F016 · Provenance and project format v2

**Fix.**
- **Manifest.** Every 2D and 3D polar gets a `RunManifest { foiley_version,
  git_commit, build_profile, target_features, thread_config, settings,
  geometry_hash, mesh_hash, family_ids+hashes, started_at }`. `build.rs`
  embeds `git_commit` and the features; hashing happens once per run.
- **Settings snapshot.** `FoilPolar.settings: AnalysisSettings`.
- **Format.** `FORMAT_VERSION = 2` in
  [project.rs](../crates/app/src/project.rs). v1 migrates: families get
  `{mach, ncrit, xtr}` from the stored polar fields and a `legacy` marker,
  and hashes are computed on load. Saves are always v2.
- **Invalidation.** Geometry edits invalidate dependent families by hash
  mismatch, generalizing the existing `upsert_foil` clearing.
- **Fixtures.** Create `tests/data/projects/`. Add a generator that saves
  each native example plane as a v2 `.foiley`, plus at least two v1 files
  written by the pre-change code and the tracked Spitfire geometry. List
  them all in `INVENTORY.txt`.

**Tests.**
- `project_roundtrip`: loads, saves, reloads and compares every file in
  `INVENTORY.txt`. It fails if the inventory is empty, under 6 files, or
  names a missing file.
- `f016_manifest_present`, `g03_hashes_stable_across_save_load`,
  `f016_geometry_edit_invalidates_family`.

#### F015a · Raw-bit hooks and harness

**Problem.** The README claims bit-exact results (lines 6, 33, 246, 247), and
nothing in the repo checks bits.
**Fix.**
- **Reference hook.** A `RAWDUMP` hook in the instrumented reference build,
  carried as `bench/xfoil-rawdump.patch` against `xfoil-ref/src` and applied
  by `build_ref.sh` into `build-ref-sp`/`build-ref-dp` (found through
  `FOILEY_REF_ROOT`). Per converged point it writes CL, CD, CDp, CM, XTR and
  the BL arrays (UEDG, THET, DSTR, CTAU) as hex `Z16`.
- **Foiley hook.** Foiley writes `f64::to_bits` in the same format from
  `src/`, behind a cargo feature `rawdump`.
- **SIMD switch.** Add `FOILEY_SIMD=0` to `src/simd.rs`'s runtime detection.
  It forces the scalar path, and is honoured only in builds without
  compile-time AVX2, since `.cargo/config.toml`'s `target-cpu=native`
  hard-enables it.
- **Harness.** A comparison that reports identical, or the max ulp per
  quantity.

**Tests.** `tests/f015.rs::f015_hooks_agree_on_naca0012`: one case,
bit-identical in sequential mode.

#### F021 · Coverage-aware auto-polars

**Problem.** `crates/app/src/auto.rs` generates polars at hardcoded Mach 0 and
Ncrit 9, checks only that a polar exists, and ignores flap angles.
**Fix.**
- **Required envelope** from the actual states: the T8 list, T6/T7 ranges,
  and mass/CG ranges.
- **Re per strip** is V·c/ν over the min/max chords.
- **CL/alpha span** from a cheap inviscid pre-pass that reuses unit
  solutions (RHS solves only).
- **Mach, Ncrit and trips** come from the selected family.
- **Generate** only the missing Re stations and missing flap angles.

**Tests.** `g08_auto_covers_t8_envelope`, `f021_generates_deflected_families`,
`f021_uses_selected_family_mach_ncrit`.

#### F004b · GUI editing of every body

**Problem.** The fuselage panel edits only `p.fuses[0]`
([plane_ui.rs:1201](../crates/gui/src/plane_ui.rs)), so F004a's extra bodies
cannot be edited.
**Fix.** A body list with add, remove, select, and per-body CAD import
(F024's job).
**Tests.** `f004b_gui_state_edits_second_body` (UI-state test, no rendering).

#### F010b · `.xfl`/`.fl5` import and conversion report

**Fix.** A parser built against hand-written minimal files and the flow5 API
example. Import 3D polar definitions where they map. Every import (AVL, xfl,
fl5) returns a `ConversionReport` (mapped, approximated, skipped, with
reasons), shown in the GUI and CLI.
**Tests.** `g45_conversion_report_lists_skips`, `f010b_minimal_xfl_and_fl5_load`.

#### F015b · Raw-bit corpus and configuration matrix

**Fix.**
- **Corpus.** The 8 cases in `bench/run_bench.py` plus 40 generated: NACA 4-
  and 5-digit, e387, dae11, trips, Mach 0–0.6, Re 5e4–6e6, CSEQ/ASEQ, and
  inverse-design round trips.
- **Matrix axes.** Sequential, `FOILEY_PIPELINE=1`, `FOILEY_FINE_PAR=1`,
  `FOILEY_CHUNK_PIPE=1`, and SIMD. The SIMD axis compares the `native`
  build against an `x86-64-v2` build run with `FOILEY_SIMD=0`.
- **README.** Keep the bit-exact claim only for configurations that are
  bit-identical, and give the measured max ulp for the rest.

**Tests.** `tests/f015.rs::f015_matrix` passes with every configuration
labelled.

#### F005b · Coupled 8×8 stability and nonzero beta (+F006)

**Problem.** A nonzero-beta reference state is silently treated as beta = 0.
**Fix.**
- **Axes.** Stability axes come from `wind_dir(alpha, beta)`.
- **Matrix.** Always build the coupled 8×8 `[u w q θ v p r φ]` matrix, and
  classify modes by eigenvector energy fractions.
- **4×4 view.** Derive the 4×4 longitudinal and lateral results as blocks of
  the 8×8, never alongside it.
- **Decoupling.** Only a mirror-symmetric aircraft (geometry, controls,
  mass) in symmetric equilibrium (β = 0, no roll or yaw rate, symmetric
  deflections) decouples. Asymmetric wings and mass distributions give real
  coupling even at β = 0, and that coupling is kept.

**Trap.** The blocks reproduce the old 4×4 only if they use identical
perturbation sizes and the same frozen wake.
**Tests.**
- `f006_symmetric_beta0_block_diagonal_matches_4x4`: off-diagonal blocks
  below 1e-12 relative on a symmetric fixture, and the blocks equal the
  previous 4×4 results.
- `f006_asymmetric_beta0_coupling_matches_fd`: on an asymmetric wing or CG
  fixture, coupling terms are nonzero and match central finite differences
  within 1e-6.
- `f006_nonzero_beta_derivs_match_fd`.
- `g34_phugoid_lanchester`: ζ ≈ 1/(√2·L/D) on a clean glider fixture.

#### F022a · STEP parse, topology, analytic surfaces, seams

**Problem.** `cad.rs` ignores trimming, so faces with holes and B-spline
surfaces come out wrong or not at all.
**Fix.** Pure Rust, no OpenCascade.
- **Module split.** `cad/step_parse.rs`, `cad/surfaces.rs`, `cad/trim.rs`,
  `cad/tessellate.rs`, `cad/check.rs`.
- **Analytic surfaces.** `PLANE`, `CYLINDRICAL_`, `CONICAL_`, `SPHERICAL_`,
  `TOROIDAL_SURFACE`, `SURFACE_OF_REVOLUTION`, `SURFACE_OF_LINEAR_EXTRUSION`.
- **Topology.** `ADVANCED_FACE` → `FACE_OUTER_BOUND`/`FACE_BOUND` →
  `EDGE_LOOP` → `ORIENTED_EDGE` → `EDGE_CURVE`, with `LINE`, `CIRCLE` and
  `ELLIPSE` curves.
- **Parameter-space edges.** Use `PCURVE`s when present. Otherwise invert 3D
  edge points onto the surface by Newton from a coarse-grid seed, with
  explicit **seam handling** (periodic parameter unwrapping) for cylinders,
  cones and revolutions.
- **Unsupported entities** are listed in the import notes; `strict` fails on
  any omission.

**Trap.** Edge inversion without `PCURVE`s fails at surface seams. Seam
handling gets its own test before any corpus work.
**Tests.** `f022_seam_inversion_cylinder_and_revolution`,
`f022_analytic_faces_with_holes_parse` (hand-written STEP fixtures in
`tests/data/cad/`).

#### F022b · B-spline surfaces and curves

**Fix.** Rational and non-rational `B_SPLINE_SURFACE`, including
`_WITH_KNOTS` complex entities, and `B_SPLINE_CURVE` edges. Edge inversion
onto B-spline surfaces uses the same Newton-with-seed path, with seam and
pole handling.
**Tests.** `f022_bspline_surface_eval_matches_reference` (closed-form
checkpoints on a hand-built bicubic patch), `f022_bspline_edge_inversion_residual_1e9`.

#### F022c · Watertight tessellation

**Fix.**
- Sample each edge once, by chordal deviation and angle tolerance.
- Share vertices between adjacent faces, so the mesh is watertight by
  construction.
- Triangulate the trimmed parameter domain, holes included, by constrained
  Delaunay via `spade`, then refine by curvature.
- Process faces in parallel with Rayon.

**Tests.** `f022_shared_edge_vertices_identical` (a box with a hole has zero
open edges), `f022_refinement_meets_chordal_tolerance`.

#### F022d · Mesh checks and CAD corpus

**Fix.**
- **Checks** in `cad/check.rs`: manifold and open-edge count, orientation,
  signed volume, area vs a high-density reference, and self-intersection
  (BVH).
- **Corpus** under `tests/data/cad/`: the flow5 example fuselage, plus
  hand-written STEP files with closed-form area and volume for a box with a
  hole, a cylinder fuselage, a revolved nacelle, a lofted B-spline fuselage
  and a multi-shell part.

**Tests.** `f022_corpus_area_volume_within_0_1pct`, `g26_watertight_and_oriented`.

#### F010c · Real XFLR5/flow5 files

**Waits on the user** to supply real XFLR5 6.x and flow5 7.58 files.
**Fix.** Add them to `tests/data/interchange/` and check that their conversion
reports are complete.
**Tests.** `f010c_real_files_import_with_reports`.
