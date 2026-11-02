# Foiley porting conventions

Foiley is a Rust port of XFoil 6.97 (`xfoil-ref/src/*.f`, the RobotLocomotion
fork).  The core solver (panel method, BL, Newton coupling, OPER ALFA/CL/ASEQ/
CSEQ, polar accumulation, NACA, LOAD, PANE/PPAR, SAVE) is done and produces
results bit-identical to the double-precision Fortran build.

## Layout
| file | contents |
|---|---|
| `src/xfoil.rs` | `XFoil` state struct (= XFOIL.INC), MRCL, COMSET, CPCALC, CLCALC, CDCALC, TECALC, SPECAL, SPECCL, VISCAL, FCPMIN |
| `src/panel.rs` | PSILIN (with GEOLIN/SIGLIN), PSWLIN, GGCALC, QDCALC, XYWAKE, STFIND, IBLPAN, XICALC, UICALC, QVFUE, QISET, GAMQV, STMOVE, UESET, NCALC, APCALC, LUDCMP/BAKSUB |
| `src/bl.rs`, `src/blsys.rs`, `src/blsolv.rs` | BL equations, SETBL, MRCHUE/MRCHDU, UPDATE, BLSOLV |
| `src/geom.rs` | NACA4/5, LEFIND, SOPPS, NORM, GEOPAR, AECALC, TCCALC, CANG, AREAD, LOAD, ABCOPY, PANGEN, ATANC, GETFLT |
| `src/spline.rs` | SPLINE, SPLIND, SPLINA, TRISOL, SEVAL, DEVAL, D2VAL, CURV, CURVS, SINVRT, SCALC, SEGSPL, SEGSPLD, SETEXP |
| `src/cli.rs` | `Session` (menus), `Input` (ASKC/ASKS/ASKR/ASKI/ASKL), top-level menu, OPER, VPAR, PPAR |
| `src/polar.rs` | polar storage + POLWRIT-format writer |
| `src/fmt.rs` | Fortran `F`, `E`, `G` edit-descriptor formatting (`f_fmt`, `e_fmt`, `g_fmt`, `fixed`) |
| `src/extra.rs` | dispatch to `gdes`, `qdes`, `mdes`, `opextra` |

## Conventions
* **Arrays are 1-based** exactly like the Fortran (`Vec` of size n+1, index 0
  unused).  Field names are the Fortran names in lower case
  (`xb, yb, xbp, ybp, sb, nb` = buffer airfoil; `x, y, xp, yp, s, n` = current
  paneled airfoil; `gam, gamu[1..2], qinv, qvis, cpi, cpv, ...`).
* 2-D Fortran arrays `A(I,IS)` with IS=1,2 are `[Vec<f64>; 3]` used as `a[is][i]`.
* Translate routines **faithfully**, keeping the same operation order so that
  results match the Fortran.  Fortran `X**2` → `x*x` (or `powi(2)`), real
  exponents → `powf`.
* Console output: use `crate::outln!(...)` / `crate::out!(...)` (buffered),
  never `println!`.  Reproduce XFoil's messages/formats where reasonable.
* Input: `s.input.askc(prompt)` returns `(COMAND, COMARG)` (command upper-cased,
  ≤4 chars); `s.input.ask(prompt)` reads a raw line (ASKS); `askr`, `aski`,
  `askl` for numbers / yes-no.  `crate::cli::args(&comarg)` returns
  `(reals, ints)` like GETFLT/GETINT.  Respect `s.input.eof` (treat as exit).
* **No graphics**: plotting commands are accepted and ignored.  Commands that
  need cursor input should print a short "cursor input not available" note
  and return.
* Persistent per-menu state goes in your module's `*State` struct, which is a
  field of `Session` (`s.gdes_state`, `s.qdes_state`, `s.mdes_state`,
  `s.extra_state`).  Avoid adding fields to `XFoil` unless the Fortran
  variable is genuinely shared with the core solver.
* After any geometry change of the *current* airfoil call the same
  invalidation the Fortran does (see `XFoil::invalidate_geometry`,
  `finish_current_airfoil`, `abcopy`, `pangen`).

## Validation
Reference executables (graphics stubbed): `build-ref-dp/xfoil.exe` (double
precision) and `build-ref-sp/xfoil.exe` (single precision, the usual build).
Feed both programs the same command script on stdin, e.g.

```
PLOP
G F

NACA 2412
...
QUIT
```

(the `PLOP / G F / blank` prefix disables graphics in the reference).  The
Fortran truncates file names to 64 characters and some to 32 — run from the
output directory and use short relative file names.  Compare output files
(coordinates, polars, dumps) numerically; the double-precision build should
match to the printed digits.

Build: `cargo build --release` (binary `target/release/foiley.exe`).
