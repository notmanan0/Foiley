# Foiley

Airfoil, wing and plane design and analysis in Rust. It covers the
xflr5 / flow5 workflow:

* **2D:** XFoil 6.97, ported line by line, with numbers identical to the original.
* **3D:** LLT, VLM1/VLM2 and a 3D panel method.
* **Flight dynamics:** stability analysis.

There are two front ends:

* **Foiley Studio** (`foiley-studio`): a desktop GUI that never blocks.
  Every analysis, mesh, save and load runs in the background, with live
  progress, stage text, point counts, ETA and a Cancel button.
* **`foiley`**: a headless CLI for scripts and servers. It shows progress
  bars or JSON progress events, writes tables, CSV or JSON, and runs JSON job
  scripts. It also stays a drop-in XFoil: pipe an XFoil script on stdin.

```
cargo build --release
target/release/foiley-studio            # GUI
target/release/foiley --help            # CLI
```

`.cargo/config.toml` builds with `-C target-cpu=native`. Remove that line to
build portable binaries; the SIMD kernels detect CPU features at run time.
Use `cargo build --profile dist` for a fully optimised (fat LTO) release.

## Layout

| crate | contents |
|---|---|
| `foiley` (repo root, `src/`) | the XFoil 6.97 engine: panel method, boundary layer, Newton solver, GDES/QDES/MDES. Results are bit-exact with the Fortran. |
| `crates/core` | `Vec3`, progress and cancellation handles, dense LU and eigenvalues (faer, multithreaded), interpolation, ISA atmosphere |
| `crates/foil` | `Foil` geometry and file I/O (Selig and Lednicer .dat), NACA 4/5 digit, every XFoil design operation as a function, inverse design (QDES/MDES), the optimiser (MOPSO/GA), the 2D analysis service with convergence recovery, polar types 1/2/3/4/6, batch polar meshes |
| `crates/plane` | plane/wing/fuselage model (elliptic, NURBS, STL/STEP), reference dimensions, inertia, meshing (thin and thick, quads or triangles, flaps, tips, fuselage), LLT, VLM1, VLM2 with optional free wake, source/doublet panel methods, polar types T1 to T8, trim, stability and control derivatives, modes, time responses, XML/AVL/STL exchange |
| `crates/app` | project file (`.foiley`, JSON), imports (including xflr5 `.xfl`, flow5 `.fl5`, `.plr`), automatic 2D polar meshes for 3D runs, the background **job manager** |
| `crates/cli` | the `foiley` binary |
| `crates/gui` | the `foiley-studio` binary (egui) |

## Why it does not freeze (and other design decisions)

The xflr5 / flow5 sources show why those programs freeze or crash. The
problems found and how Foiley avoids them:

| flow5 / xflr5 | Foiley |
|---|---|
| Streamlines, meshing, CAD import and project saving run on the GUI thread, and the GUI thread waits on worker futures, so Windows marks the window "Not Responding". | The UI thread never computes. Jobs run on a dedicated pool that leaves one core free for the interface. Meshing, saving and loading are background jobs too. |
| Some task paths exit without posting their final message, and the dialog then waits forever. | Every job ends in exactly one terminal state (Done, Failed or Cancelled), even if the solver panics. Panics are caught per job (the build uses `panic = "unwind"`). |
| Cancel flags are global, non-atomic statics, so starting one task un-cancels the others. | Each job has its own atomic cancel handle. Cancellation is checked between points and matrix blocks, and above 2,000 unknowns the LU factorization is blocked and checks between column tiles, so it takes effect within milliseconds: 17 ms p95 over 10 cancels of a 14,100-unknown factorization (i7-14650HX). |
| Workers write into polars the GUI is drawing, and deleting a foil during a run can crash. | A job owns immutable copies of its inputs. Results stream through a mutex-protected buffer and are merged on the UI thread, and finished points appear live in the plots. |
| Unconverged XFoil points are skipped, with no retry. | Recovery strategy: re-initialise the boundary layer with more iterations, then march in sub-steps from the last converged point or from easier neighbouring points. Points that still fail are listed and marked in the plots, never silently lost. Sweeps march outwards from zero lift so each point warm-starts on the attached-flow branch. |
| Analysis settings are process-global. | Every setting is stored with the polar, so results are reproducible from the project file alone. |
| Lednicer .dat files are misread, and the clockwise check is a no-op. | Selig and Lednicer are both parsed, the orientation is detected from the data, and errors are reported clearly. |
| Project save is not atomic. | Saves write to a temp file and rename it. There is an autosave every 3 minutes when modified, and undo/redo covers 40 steps. |

All settings have hover help that explains the parameter and its typical
values (for example Ncrit: 9 is an average wind tunnel, 11–12 quiet free air).

## Command line

```
foiley foil info naca2412 my_foil.dat
foiley foil polar "NACA 4412" --re 3e5 --alpha=-4:14:0.5 -f csv -o polar.csv
foiley foil polar naca0012 --type 2 --re 2e5 --cl 0.1:1.2:0.1      # fixed-lift polar
foiley foil polar naca0012 --re-sweep 1e5:1e6:1e5 --at-alpha 4      # type 4
foiley foil polar naca2412 --flap-sweep=-10:20:2 --hinge 0.75,0.5   # type 6
foiley foil point naca2412 --alpha 5 -o cp_and_bl.csv               # Cp and boundary layer
foiley foil batch naca2412 e387.dat --re 5e4:1e6:8 --ncrit 9,11 --out-dir polars --save mesh.foiley
foiley foil design naca0012 --flap 0.7,0.5,10 --thickness 0.10 --te-gap 0.002 -o new.dat
foiley foil inverse naca2412 --alpha 2 --segment upper,0.3,0.7 --scale 0.97 -o inv.dat   # QDES
foiley foil inverse naca2412 --mode full --segment lower,0.2,0.7 --bump 0.03 -o inv.dat  # MDES
foiley foil optimize naca0012 --re 5e5 --alpha 3 --objectives "cl/cd:max,cm:0" -o best.dat

foiley plane info glider.xml                     # xflr5 / flow5 XML or a .foiley project
foiley plane analyze glider.xml --type t2 --method vlm2 --range=-2:10:1 --strips span.csv
foiley plane analyze glider.xml --type t7 --gain 1,0,5 --gain 1,1,5 --response doublet:0.5 --response-out r.csv
foiley plane analyze glider.xml --method tripanel --inviscid --speed 20 -f json -o out.json
foiley plane analyze glider.xml --method trilinear --inviscid --range=0:6:2  # Galerkin linear doublets
foiley plane analyze glider.xml --free-wake 4 --range 2:8:2            # relaxed VLM wake
foiley plane import-fuselage glider.xml --file body.step --units auto --mass 0.3 --save g.foiley
foiley plane export g.foiley --avl g.avl --mass g.mass --stl g.stl --xml g.xml
foiley project old_design.xfl --save converted.foiley   # xflr5 .xfl, flow5 .fl5, .plr polars
foiley run docs/examples/glider_study.json       # scripted study (see the example)

foiley < script.xfoil                            # XFoil-compatible menu mode
```

* Missing 2D polars for viscous or LLT 3D runs are computed automatically
  over the Reynolds range the analysis needs (`--auto-polars N`).
* NACA foils referenced by an imported plane are generated automatically.
* Ctrl-C cancels the running job; a second press quits.
* `--json-progress` emits progress and log events as JSON lines on stderr.

## GUI

There are five workspaces: **Foil design**, **2D analysis**, **Plane design**,
**3D analysis** and **Stability**. A project explorer sits on the left, and
the Jobs and Log panel is always visible at the bottom.

* **Foil design.**
  * Every geometric operation is previewed as an overlay before you apply it
    or save it as a new foil: flap, thickness and camber, max-point positions,
    LE radius, TE gap, repanelling, normalisation, interpolation.
  * **Inverse design** (mixed QDES / full MDES): plot the surface speed,
    scale, bump or smooth a segment, see the target, and solve.
  * **Optimiser** (MOPSO particle swarm or genetic algorithm) with
    Hicks–Henne bumps and an optional flap, on Cl, Cd, Cl/Cd, Cpmin and Cm
    (minimise, maximise or target). The best design and the score history
    update live.
* **2D analysis.** Polar types 1/2/3/4/6 and single points.
  * Four configurable polar plots.
  * Click a point to see its Cp and boundary-layer distributions.
  * Batch polar meshes.
  * Unconverged points are marked.
* **Plane design.**
  * Wing, section and point-mass tables, with a live 3D view meshed in the
    background.
  * Fuselages: elliptic frames, NURBS frames of control points, or an
    imported STL/STEP surface (with unit choice).
  * A summary shows span, area, MAC, aspect ratio, wing loading, tail
    volumes, CoG and inertia.
  * Export to AVL (.avl + foil files, .mass), STL and xflr5 XML.
* **3D analysis.**
  * Polar types T1 to T8 with the methods LLT, VLM1, VLM2, quad panels and
    triangle panels.
  * Optional free (relaxed) wake for VLM, drawn in the 3D view.
  * Viscous drag from the 2D polars; ground effect; altitude or water.
  * Mass overrides, flap presets, extra drag, T6 ranges and T7 gains.
  * Plots update live while a job runs.
  * Span distributions, and a 3D view coloured by Cp or loading with lift vectors.
* **Stability.** Trim, static margin and neutral point; a modes table;
  dimensional and non-dimensional derivatives; state matrices; root locus.
  Modal time responses, and forced responses to step, pulse, doublet or
  sine control inputs (T7 gains give the control derivatives).

`foiley-studio project.foiley --module 3d` opens a project in a given
workspace. File › Import also reads xflr5 `.xfl`, flow5 `.fl5` and `.plr`
files.

## Feature coverage compared with xflr5 / flow5

| area | status |
|---|---|
| XFoil direct analysis, polar types 1/2/3/4/6, Cp and boundary layer, batch polar meshes | ✔ |
| Foil design: NACA, flap, TE gap, LE radius, thickness/camber, max-point positions, repanel, derotate/normalise, interpolation | ✔ |
| Mixed and full inverse design (QDES / MDES) | ✔ GUI and CLI (all other XFoil menus via `foiley xfoil`) |
| Foil optimiser: MOPSO and genetic algorithm, Hicks–Henne bumps, flap, multi-objective with Pareto archive | ✔ |
| Planes: multi-section wings, dihedral, twist, sweep, fins, flaps (control surfaces), point masses, inertia | ✔ |
| Fuselages: elliptic frames, NURBS (xflr5/flow5 frames), STL and STEP surfaces | ✔ (STEP: B-spline/NURBS surfaces and planar faces; other analytic surface types are skipped with a note) |
| LLT, VLM1, VLM2, 3D panels on quads and triangles, Galerkin linear-doublet triangles, ground effect | ✔ |
| Free (relaxed) wake for VLM | ✔ |
| Polar types T1 fixed speed, T2 fixed lift, T3 glide, T4 fixed α, T5 sideslip, T6 control, T7 stability, T8 list | ✔ |
| Stability derivatives, control derivatives, state matrices, modes, root locus, modal and forced time responses | ✔ |
| Import: xflr5 `.xfl` and flow5 `.fl5` projects (foils, polars, planes), `.plr` polars, plane/wing XML, `.dat`, XFoil polar text | ✔ (projects' 3D polars and operating points are not imported; the readers were tested against byte streams built from the C++ write code, not real project files) |
| Export: AVL `.avl` and `.mass`, STL mesh, xflr5 XML, `.dat`, XFoil polar text, CSV, JSON | ✔ |

LLT treats the wing as straight and planar (sweep and dihedral are ignored), as
in xflr5 and flow5; use VLM or panels for swept or dihedral wings.

## Validation

`cargo test --release --workspace` runs all of the checks below.

**XFoil engine**

* Polars are compared line for line with the double-precision Fortran
  XFoil reference polars, for the sequential, pipelined and parallel sweeps.

**2D service**

* The NACA 0012 drag and lift slope match XFoil.
* Fixed-lift polars satisfy Re = Re√CL/√CL.
* Recovery: a NACA 2412 point at Re 4.5e5 fails from a cold start in XFoil
  and is recovered.
* Cancellation and invalid inputs return errors and never hang.
* Inverse design: a mixed inverse changes only the targeted segment.
* Optimiser: NACA 0012 Cl/Cd improves (53 → 73 in a short run).

**3D panel kernels**

* The doublet kernel reproduces the solid angle: −1 inside a closed cube,
  0 outside, +½ just above a panel.
* The source potential matches brute-force numerical integration, and its
  analytic velocity matches the finite-difference gradient.

**3D aerodynamics** (`crates/plane/tests/validation.rs`)

* An aspect-ratio-8 rectangular wing gives CLα ≈ 4.6/rad and span
  efficiency 0.98, consistent with lifting-surface theory. VLM1 and VLM2 agree.
* In potential flow, a 12% thick wing lifts more than a flat plate; a
  symmetric section at 0° gives exactly zero lift.
* Triangle panels agree with quad panels within 3%, with no spurious yaw.
* With a fuselage, the induced drag stays positive and scales with CL², because
  circulation is carried across the fuselage.
* Ground effect reduces induced drag and increases lift.
* Galerkin linear-doublet triangles (`crates/plane/tests/trilinear.rs`) pass
  thin-wing lift-slope, trailing-edge pressure, mesh and quadrature refinement,
  symmetry, camber, ground, fuselage and flap checks. See
  [the formulation and validation limits](docs/TRILINEAR.md).
* The free wake stays within 5% of the rigid-wake lift and trails downwards.
* The sideslip response is antisymmetric.
* LLT and VLM agree to within 12% in the linear range.
* T2 lift equals the weight, and T3 aerodynamic force equals the weight.
* The elevator flap moves Cm the right way.
* A NURBS fuselage matches the elliptic one's wetted area. flow5's STEP and
  STL example fuselages import, the STEP file with its foot units resolved.

**Stability**

* The default glider trims with a positive static margin and
  Cmα < 0, Cmq < 0, Clp < 0, Cnβ > 0.
* Zu = −2mg/u₀.
* The phugoid frequency matches the short-period-residualised approximation
  within 15%.
* All five modes are found: short period, phugoid, Dutch roll, roll damping
  and spiral.
* With elevator gains, the control derivative is nose-down for trailing edge
  down, and a step input pitches the aircraft nose-down.

## Performance

Measured on an i7-14650HX, release build.

| case | time |
|---|---|
| NACA 4412, 20-point viscous polar (Re 3e5, transition and stall) | < 0.1 s |
| 2D polar mesh: 2 foils × 6 Re × 57 points (696 points) | 2.7 s |
| T2 VLM2 polar of the default glider (634 panels) | 0.03 s |
| T7 trim, derivatives and modes | 0.2 s |
| 3D panel method, 6,120 unknowns, 15 operating points (single factorisation, unit-solution superposition) | 1.1 s |
| complete scripted study: XFoil polar, batch, XML import, T2, T7, T6 | 2.0 s |

**Engine speed against Fortran XFoil.** The XFoil engine keeps its own
optimisations: a cache-blocked BLSOLV, AVX2 kernels, a pipelined wake
computation and parallel sweeps. Whole-run benchmarks against the
gfortran-built Fortran XFoil are in `bench/run_bench.py`:
* exact mode is about 2.0× faster than double precision and 2.5× faster
  than single precision;
* default mode is about 4.2× and 5.1× faster.

**Parallelism controls** (environment variables, XFoil mode):

| variable | default | effect |
|---|---|---|
| `FOILEY_JOBS` | min(8, logical CPUs / 2) | parallel sweep chunks; `1` gives exact XFoil semantics (bit-identical output) |
| `FOILEY_PIPELINE` | `1` | `0` disables the bit-exact helper-thread pipelining of ASEQ |
| `FOILEY_FINE_PAR` | `0` | `1` enables intra-iteration parallelism |

**3D method details.** The panel method takes forces from the far field
(Kutta–Joukowski on the shed circulation) and moments from the integrated
pressures, as flow5 does. Far-field forces converge much faster with
chordwise panel count. In the default 3D methods all operating points of a
polar share one LU factorisation. `docs/PORTING.md` documents the XFoil port
conventions.

`trilinear` uses three independent doublet values per triangle, linear
Galerkin test functions, and a wake tied to the actual trailing-edge values.
Select it in the GUI's method list or use `"method": "TriLinear"` in JSON jobs.
It needs more memory and computation than `tripanel`; start with fewer panels.

## Reference sources

The specs used for this work were taken from:

* **flow5:** github.com/techwinder/flow5, cloned into `refsrc/` (git-ignored).
* **xflr5 v6.61:** github.com/polmes/xflr5-ubuntu, a source mirror; xflr5
  itself is published only in SVN on SourceForge.
