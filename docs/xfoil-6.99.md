# XFOIL 6.99 upgrade

Foiley's 2D engine follows XFOIL 6.99: its wake discretization, its
surface-integrated pressure drag, separate top and bottom Ncrit, a
configurable wake length and the wake dissipation constant. The integral
boundary-layer correlations and their baseline calibration are unchanged
from 6.97.

The port was checked against the original MIT archives, including the
intervening 6.98 changes:

| Release | Source | SHA-256 |
|---|---|---|
| 6.97 | [MIT archive](https://web.mit.edu/drela/Public/web/xfoil/xfoil6.97.tar.gz) | `293e28e269bd88a62f8b6ba6bf684f16843f6776736e6b77aaa8b6bebc467a81` |
| 6.99 | [MIT archive](https://web.mit.edu/drela/Public/web/xfoil/xfoil6.99.tgz) | `5c0250643f52ce0e75d7338ae2504ce7907f2d49a30f921826717b8ac12ebe40` |

## Behavior and controls

| Area | Change |
|---|---|
| Wake | `N/12 + 10*INT(WAKLEN)` points, capped by the wake array (`IWX` = 48). A 160-node, one-chord wake uses 23 points instead of 22. |
| Transition | Separate top and bottom Ncrit, 9 on both by default. Polars record the fractional transition panel index of each side (`Top_Itr`, `Bot_Itr`). |
| Lag | The third LAG constant is the wake inverse dissipation-length factor (default 0.9). |
| CDp | 6.99 keeps the integrated surface pressure drag instead of overwriting it with `CD-CDf`. The iteration display still prints `CD-CDf`, as upstream does. |
| Geometry | Zero-normal fallback in NCALC, repeated FIXDUP passes, trailing-edge pairing in camber/thickness extraction, stored-airfoil point counts, and LOAD keeping an existing MSES domain string. |
| Design | `GDES NMOV` (normal offset) and `CAMB WRTT` (thickness and slope export). |
| Plots | `VPLO XPLO` switches the BL abscissa between x and arc length s; `OPER CPIO` toggles the inviscid Cp overlay flag (script mode draws no Cp plot, as with `CPGR`). Studio can plot the BL against arc length. |
| Files | Two Ncrit values in polar, polar-dump, defaults and VPLO dump headers; Ncrit and Itr columns; the 6.99 DUMP layout (H*, P, m, K on surface rows) and the 24-column DMP; the 6.99 polar listing, mismatch tables and VPAR display. Single-Ncrit files still read. |

The controls are available in Studio and the modern CLI:

```text
foiley foil polar naca2412 --ncrit 5 --ncrit-bot 11 --wake-length 2.5 --wake-dissipation 1.1
```

For XFOIL scripts they live in `OPER` → `VPAR`, as in 6.99:

```text
NT 5
NB 11
WAKE 2.5
LAG 5.6 1.0 1.1
```

`N` sets both surfaces, `NT` and `NB` one each. Like XFOIL, script mode does
not range-check these values. One deliberate difference: a two-argument
`LAG` (a 6.97 script) sets the first two constants and keeps the third,
where 6.99 would prompt for all three and consume the next script lines.

`Foil::naca` (Studio, the CLI, 3D meshes) keeps its 120 points per side so
saved projects and meshes do not move; script-mode `NACA` uses 6.99's
array-derived 123. After repaneling the two agree to the printed digits of
the coefficients.

Studio and the modern CLI accept a wake length of 1 to 20 chords. Beyond
about 4 chords (at 160 panels) the wake array is full, so a longer wake
spreads the same points further apart: lengthening it is not a
mesh-refinement study.

## Saved data and compatibility

All new analyses run 6.99. There is no user-facing 6.97 mode: the engine
keeps a `SolverVersion` switch only so the regression tests can still
compare against the 6.97 Fortran.

`AnalysisSettings` gains `ncrit_bot`, `wake_length` and `wake_dissipation`.
At their defaults (none, 1, 0.9) they are left out of the JSON, so existing
settings, polar identities and family hashes are unchanged. A bottom Ncrit
equal to the top one is treated as no split.

Polars computed before this change (6.97 CDp = `CD-CDf`) and new ones share a
family when their Mach, Ncrit and trips agree; only a split bottom Ncrit forms
a separate family, selected by the 3D analysis's bottom Ncrit. Wake length
and dissipation are analysis settings like paneling, not family keys;
auto-polar generation carries a family's wake settings to the polars it adds.

The project format is now 3. Format 2 projects load unchanged; builds from
before this change refuse format 3 files instead of silently dropping the
new settings.

6.99 `CDp` is not a substitute for `CD-CDf`. Native results keep friction drag
separately, and native text exports include a `CDf` column. Importing a 6.99
file (recognized by its version line) without `CDf` records friction as 0,
unavailable, rather than inferring it from CD and the integrated CDp.

## Validation and limits

Six committed MIT 6.99 polars cover NACA 4/5-digit generation, 160 and 240
nodes, split Ncrit, a 2.5-chord wake, changed dissipation, integrated CDp and
transition indices. Script mode reproduces them at printed precision; the
library path (`analyze`) is checked against the split-Ncrit polar. They run
in CI without a Fortran installation.

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

The 6.97 polar fixtures and raw-array audit select the engine's 6.97
conventions. Native rendering replaces the Fortran Xplot/PPLOT UI; their
platform build changes and Fortran-only reconstructed velocity-profile
plotting are outside this port. The sharp-trailing-edge viscous influence
correction is a later 6.996 change and is not included here.
