#!/bin/bash
# Build the reference Fortran XFoil (single + double precision) with
# gfortran -O2 into build-ref-sp / build-ref-dp.
#
# The sources are a staged copy of xfoil-ref/src (build-ref-src, line ends
# normalised) with bench/xfoil-rawdump.patch applied: a RAWDUMP hook in
# PLRADD that writes every stored point's coefficients and BL arrays as hex
# bits when the RAWDUMP environment variable names a file, and does nothing
# otherwise.  xfoil-ref itself is not modified.
#
# The GPL sources and builds are gitignored; they are found through
# FOILEY_REF_ROOT (default: this checkout).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${FOILEY_REF_ROOT:-$HERE}"
for d in /ucrt64/bin /c/msys64/ucrt64/bin /mingw64/bin /c/msys64/mingw64/bin; do
  if [ -x "$d/gfortran" ] || [ -x "$d/gfortran.exe" ]; then export PATH="$d:/usr/bin:$PATH"; break; fi
done
command -v gfortran >/dev/null || { echo "gfortran not found (looked in PATH, /ucrt64, /c/msys64/ucrt64, mingw64)"; exit 1; }
SRC0="$ROOT/xfoil-ref/src"; PLT="$ROOT/xfoil-ref/plotlib"
[ -d "$SRC0" ] || { echo "no reference sources at $SRC0 (set FOILEY_REF_ROOT)"; exit 1; }
STUB="$HERE/refstub/nullgfx.c"; [ -f "$STUB" ] || STUB="$ROOT/refstub/nullgfx.c"
SRC="$ROOT/build-ref-src"
rm -rf "$SRC"; mkdir -p "$SRC"
for f in "$SRC0"/*; do [ -f "$f" ] && sed 's/\r$//' "$f" > "$SRC/$(basename "$f")"; done
(cd "$SRC" && sed 's/\r$//' "$HERE/bench/xfoil-rawdump.patch" | patch -p1 --quiet) || { echo "FAIL applying bench/xfoil-rawdump.patch"; exit 1; }
XS="xfoil.f xpanel.f xoper.f xtcam.f xgdes.f xqdes.f xmdes.f xsolve.f xbl.f xblsys.f xpol.f xplots.f pntops.f xgeom.f xutils.f modify.f blplot.f polplt.f aread.f naca.f spline.f plutil.f iopol.f gui.f sort.f dplot.f profil.f userio.f frplot0.f"
PS="plt_base.f plt_font.f plt_util.f plt_color.f set_subs.f gw_subs.f ps_subs.f"
for v in sp dp; do
  OUT="$ROOT/build-ref-$v"; mkdir -p "$OUT"; cd "$OUT" || exit 1
  FL="-O2 -std=legacy -w"; [ $v = dp ] && FL="$FL -fdefault-real-8"
  rm -f ./*.o xfoil.exe
  for f in $PS; do gfortran $FL -I"$PLT" -c "$PLT/$f" || echo "FAIL $f"; done
  gcc -O2 -c "$STUB" || echo "FAIL nullgfx"
  for f in $XS; do gfortran $FL -I"$SRC" -c "$SRC/$f" || echo "FAIL $f"; done
  gfortran -o xfoil.exe ./*.o -static || echo LINKFAIL
  ls -la xfoil.exe
done
