# XFOIL 6.99 upgrade

Foiley now defaults to the 6.99 wake discretization and surface-integrated
pressure-drag output. The integral boundary-layer correlations and their
baseline calibration are unchanged. The original 6.97 numerical conventions
remain selectable for reproducing earlier runs.

The port was checked against the original MIT archives, including the
intervening 6.98 changes:

| Release | Source | SHA-256 |
|---|---|---|
| 6.97 | [MIT archive](https://web.mit.edu/drela/Public/web/xfoil/xfoil6.97.tar.gz) | `293e28e269bd88a62f8b6ba6bf684f16843f6776736e6b77aaa8b6bebc467a81` |
| 6.99 | [MIT archive](https://web.mit.edu/drela/Public/web/xfoil/xfoil6.99.tgz) | `5c0250643f52ce0e75d7338ae2504ce7907f2d49a30f921826717b8ac12ebe40` |

## Behavior and controls

| Area | Change |
|---|---|
| Wake | `N/12 + 10*INT(WAKLEN)` points, capped by the wake array. A 160-node, one-chord wake uses 23 points instead of 22. Wake length is configurable. |
| Transition | Independent top/bottom Ncrit, defaulting to 9 on both surfaces; fractional transition panel indices are recorded in polars. |
| Lag | The third LAG constant controls wake inverse dissipation length. Two-argument scripts retain the existing third constant without prompting. |
| CDp | 6.99 preserves the integrated surface pressure drag. The iteration display still prints `CD-CDf` as upstream does. These quantities need not be equal. |
| Geometry | Zero-normal fallback, repeated endpoint cleanup, trailing-edge pairing in camber/thickness extraction, and stored-airfoil point-count/state fixes. |
| Design | `GDES NMOV` normal offsets and `CAMB WRTT` thickness/slope export. Studio supports dragging geometry nodes into a preview and dragging inverse target speeds in QDES/MDES. |
| Plots | Native Cp overlay control, x/arc-length BL coordinate selection, and `VPLO XPLO` for script-mode variable dumps. |
| Files | Dual Ncrit headers, transition indices, extended DUMP, detailed 24-column DMP, and backward reading of single-Ncrit polar/defaults files. |

New analysis controls are available in Studio and the modern CLI:

```text
foiley foil polar naca2412 --ncrit 5 --ncrit-bot 11 --wake-length 2.5 --wake-dissipation 1.1
foiley foil polar naca2412 --xfoil-version 6.97
```

For XFOIL scripts, the controls live in `OPER` → `VPAR`:

```text
NT 5
NB 11
WAKE 2.5
LAG 5.6 1.0 1.1
```

`N 9` sets both surfaces together. `VERS 6.97` can be used at top level
before `NACA` to reproduce 6.97's seed geometry as well as its wake and CDp
conventions; `VPAR VERS` changes conventions on the current geometry.
The shared `Foil::naca` library generator keeps its existing 120 points per
side so existing 3D meshes remain stable. XFOIL script-mode NACA follows
each version's array-derived seed count.

Wake length accepts 1–20 chords; very long wakes hit the array cap, so
lengthening the wake is not equivalent to a mesh-refinement study. Changing
the wake invalidates its trajectory, BL mapping, initialization and convergence.

## Saved data and exports

`AnalysisSettings` records solver version, optional bottom Ncrit, wake length,
and wake dissipation. Settings serialized before the version field existed
load with 6.97 conventions; newly created settings default to 6.99. A missing
bottom Ncrit follows the existing scalar value. Polar identities and 3D
family keys distinguish bottom Ncrit and solver/wake calibration. Auto-polar
generation reuses the selected family's saved settings.

6.99 `CDp` is not a substitute for `CD-CDf`. Native results retain actual
friction drag separately. Native text exports include a CDf column. On import,
if a 6.99 file omits CDf, that field is recorded as 0 (unavailable), rather than
inferring friction from the difference of total and integrated pressure drag.
Consumers should inspect column labels instead of assuming seven columns.

Foiley writes consistent DUMP columns on both surface and wake rows and
supports delimiters in its detailed DMP export. It does not reproduce the
MIT archive's inconsistent DUMP header or unconditional diagnostic `.bl`
file writes.

## Validation and limits

Six committed Fortran polar fixtures cover NACA 4/5-digit generation,
160/240 nodes, split Ncrit, a 2.5-chord wake, changed dissipation, integrated
CDp, and transition indices. They run in CI without a Fortran installation.
Tests also cover old/new settings and polar/defaults files, input guards,
normal offsets, export layouts, and incompatible polar-family separation.

The optional 6.99 raw-array audit uses the original 48-case corpus plus
three new-control cases. In the local double-precision run, the new-control
cases had maximum scaled differences of 2.3e-11 (split Ncrit), 2.3e-9
(long wake), and 1.1e-8 (lag dissipation). The full corpus is not bit-identical:
the largest aligned difference was 4.9e-5 in lower transition position, and
the low-Re DAE11 case converged at different sets of points (one point only
in Fortran, four only in Foiley). The audit reports those differences;
it does not turn them into a claim of full numerical equivalence. The
new-control cases have a strict 1e-7 scaled bound and matching BL structure.

To reproduce:

```text
python bench/build_reference.py --source <extracted-Xfoil-6.99> --out <build>
# Set FOILEY_REF_699 to <build>/xfoil[.exe], then:
cargo test --release -p foiley --test f015 f699_reference_audit -- --nocapture
```

The builder stages a copy, adds the existing raw-dump hook, and removes only
the archive's post-solve diagnostic writer, which generates invalid Windows
filenames. It leaves numerical calculations and convergence logic intact.

The old 6.97 polar fixtures and raw-array audit explicitly select 6.97
conventions. Native rendering replaces the Fortran Xplot/PPLOT UI; their
platform build changes and Fortran-only reconstructed velocity-profile
plotting are outside this port. The sharp-trailing-edge viscous influence
correction is a later 6.996 change and is not included here. Studio drag
interactions are compiled with the workspace; manual pointer-interaction
QA is still required.
