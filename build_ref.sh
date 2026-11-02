#!/bin/bash
# Build reference Fortran XFoil (single + double precision) with gfortran -O2
export PATH=/ucrt64/bin:/usr/bin:$PATH
ROOT="/c/Users/Manan Malavia/Foiley"
SRC="$ROOT/xfoil-ref/src"; PLT="$ROOT/xfoil-ref/plotlib"
XS="xfoil.f xpanel.f xoper.f xtcam.f xgdes.f xqdes.f xmdes.f xsolve.f xbl.f xblsys.f xpol.f xplots.f pntops.f xgeom.f xutils.f modify.f blplot.f polplt.f aread.f naca.f spline.f plutil.f iopol.f gui.f sort.f dplot.f profil.f userio.f frplot0.f"
PS="plt_base.f plt_font.f plt_util.f plt_color.f set_subs.f gw_subs.f ps_subs.f"
for v in sp dp; do
  OUT="$ROOT/build-ref-$v"; mkdir -p "$OUT"; cd "$OUT"
  FL="-O2 -std=legacy -w"; [ $v = dp ] && FL="$FL -fdefault-real-8"
  rm -f *.o; for f in $PS; do gfortran $FL -I"$PLT" -c "$PLT/$f" || echo FAIL $f; done
  gcc -O2 -c "$ROOT/refstub/nullgfx.c" || echo FAIL nullgfx
  for f in $XS; do gfortran $FL -I"$SRC" -c "$SRC/$f" || echo FAIL $f; done
  gfortran -o xfoil.exe *.o -static || echo LINKFAIL
  ls -la xfoil.exe
done
