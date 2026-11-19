# Capturing the 3D viewport

Build a development executable with `cargo build -p foiley-cli`. Use the
default `snap` feature or run `cargo snap`. The examples below run from
the repo root at 1400 × 900 and scale 1 unless the preset sets a different
size. On Windows the executable ends in `.exe`.

```sh
cargo snap --preset views
cargo snap --preset orbit --crop viewport
cargo snap --preset nav --crop viewport
cargo snap --preset analysis
cargo snap --preset narrow --script docs/snap-review/narrow.snap
```

The manifest lists frames, camera state, picks, overlay labels and solver
jobs. The canvas coordinates are in egui points. Script positions without
the `vp:` prefix are image pixels. A `vp:` position uses fractions of
the 3D canvas. For example, `vp:0.5,0.5` is its center. The narrow preset
checks the pick against the CPU image at that position.

## Review cases

Generate three small projects under `target/snap/fixtures` using Python
and the built executable.

```sh
python docs/snap-review/fixtures.py target/debug/foiley
```

Use `--script docs/snap-review/<script>` with the arguments below. You can
set an explicit output directory with `--out`.

| Script | Other arguments | What to look for |
|---|---|---|
| `ctrl-mmb.snap` | `--module plane --nav cad --crop viewport` | Horizontal movement changes neither target nor distance. Vertical movement changes distance while the target stays fixed. |
| `fit-interrupt.snap` | `--module plane --crop viewport` | The orientation at the end of Fit matches the one at its start. |
| `rear-view.snap` | `--module plane` | Shift+1 moves yaw to 0° and pitch to 0°. |
| `plot-click.snap` | `--module 3d target/snap/fixtures/control.foiley --crop viewport` | When control changes from 0 to 1, the old mesh and pressure colors disappear. The new mesh and pressure colors appear together. |
| `rerun.snap` | `--module 3d target/snap/fixtures/rerun.foiley` | The same-name solve moves the seeded CoG marker back to the origin. The fixture differs only in the stored CoG so it can exercise the result cache. |
| `autorun.snap` | `--module 3d --run` | The finished run has a Cp legend and a CoG marker. |
| `llt.snap` | `--module 3d target/snap/fixtures/llt.foiley` | The main wing shows lift vectors at the LLT stations and a CoG marker. |
| `rename.snap` | `--module plane` | The plane name and mesh key change. The zoom and target stay the same. |
| `narrow.snap` | `--preset narrow` | The hover readout identifies the visible panel under the canvas center at 900 × 600. |

The PNGs show the CPU painter. Face sorting, line antialiasing, and wake
and lift-line occlusion can differ from the OpenGL renderer.
