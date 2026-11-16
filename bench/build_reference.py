"""Build a headless, raw-dump-instrumented MIT XFOIL source tree.

Example: python bench/build_reference.py --source /path/to/Xfoil --out work/ref699
No source files are changed. Requires gcc and gfortran from the same toolchain.
"""
import argparse
import os
from pathlib import Path
import shutil
import subprocess


def build(source, output, compiler):
    source, output = source.resolve(), output.resolve()
    src, plot = source / "src", source / "plotlib"
    if not (src / "xfoil.f").is_file() or not (plot / "plt_base.f").is_file():
        raise ValueError("--source must be an extracted Xfoil directory")
    if output == source or output in source.parents:
        raise ValueError("--out must be separate from the source tree")
    output.mkdir(parents=True, exist_ok=True)
    staged = output / "src"
    staged.mkdir(exist_ok=True)
    for f in src.iterdir():
        if f.is_file() and f.suffix.lower() in (".f", ".inc"):
            (staged / f.name).write_text(f.read_text())
    xpol = staged / "xpol.f"
    text = xpol.read_text()
    start = text.index("      SUBROUTINE PLRADD")
    hook = text.index("      LOGICAL ERROR", start) + len("      LOGICAL ERROR")
    patch = (Path(__file__).parent / "xfoil-rawdump.patch").read_text()
    raw = patch[patch.index("+      SUBROUTINE RAWDMP"):]
    raw = "\n".join(l[1:] for l in raw.splitlines() if l.startswith("+"))
    xpol.write_text(text[:hook] + "\n      CALL RAWDMP" + text[hook:] + "\n" + raw + "\n")
    # The 6.99 tarball contains a post-solve diagnostic that writes names
    # such as ':00.bl', invalid on Windows. Remove only that diagnostic;
    # all aerodynamic calculations and convergence logic remain untouched.
    xoper = staged / "xoper.f"
    text = xoper.read_text()
    diagnostic = text.find("\n        is = 1", text.index("      SUBROUTINE VISCAL"))
    if diagnostic >= 0:
        end = text.index("\n      RETURN", diagnostic)
        xoper.write_text(text[:diagnostic] + text[end:])
    fc = Path(compiler).resolve() if Path(compiler).is_file() else Path(shutil.which(compiler) or compiler)
    gcc = fc.with_name("gcc.exe" if os.name == "nt" else "gcc")
    env = dict(os.environ)
    env["PATH"] = str(fc.parent) + os.pathsep + env.get("PATH", "")
    flags = ["-O2", "-std=legacy", "-w", "-fdefault-real-8", "-fallow-argument-mismatch"]
    names = "xfoil xpanel xoper xtcam xgdes xqdes xmdes xsolve xbl xblsys xpol xplots pntops xgeom xutils modify blplot polplt aread naca spline plutil iopol gui sort dplot profil userio frplot0".split()
    plots = "plt_base plt_font plt_util plt_color set_subs gw_subs ps_subs".split()
    objects = []
    for directory, units in [(staged, names), (plot, plots)]:
        for name in units:
            obj = output / (name + ".o")
            subprocess.run([str(fc), *flags, "-I" + str(staged), "-I" + str(plot), "-c", str(directory / (name + ".f")), "-o", str(obj)], check=True, env=env)
            objects.append(str(obj))
    stub = Path(__file__).resolve().parent.parent / "refstub/nullgfx.c"
    obj = output / "nullgfx.o"
    subprocess.run([str(gcc), "-O2", "-c", str(stub), "-o", str(obj)], check=True, env=env)
    objects.append(str(obj))
    exe = output / ("xfoil.exe" if os.name == "nt" else "xfoil")
    subprocess.run([str(fc), "-o", str(exe), *objects, "-static"], check=True, env=env)
    print(exe, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--compiler", default="gfortran")
    args = parser.parse_args()
    build(args.source, args.out, args.compiler)
