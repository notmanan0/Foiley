"""Benchmark + validation harness: runs reference Fortran XFoil and Foiley on
identical command scripts, compares the polar files and reports wall times.

usage: python bench/run_bench.py [case-filter] [--reps N]
"""
import os, subprocess, sys, time, re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "bench", "out")
RUNS = os.path.join(ROOT, "xfoil-ref", "runs")
EXES = {
    "xfoil-sp": os.path.join(ROOT, "build-ref-sp", "xfoil.exe"),
    "xfoil-dp": os.path.join(ROOT, "build-ref-dp", "xfoil.exe"),
    "foiley":   os.path.join(ROOT, "target", "release", "foiley.exe"),
}

def geom(g):
    if g.startswith("NACA"):
        return f"{g}\n"
    import shutil
    shutil.copy(os.path.join(RUNS, g), os.path.join(OUT, g))
    return f"LOAD {g}\nPANE\n"

def script(g, oper, polfile):
    return ("PLOP\nG F\n\n" + geom(g) + "OPER\n" + oper.replace("@POL", polfile) + "\n\nQUIT\n")

CASES = {
 "n0012_re1e6":   ("NACA 0012", "VISC 1e6\nITER 100\nPACC\n@POL\n\nASEQ 0 16 0.5"),
 "n2412_re3e5":   ("NACA 2412", "VISC 3e5\nITER 100\nPACC\n@POL\n\nASEQ -5 14 0.25"),
 "n4415_re6e6_m3":("NACA 4415", "VISC 6e6\nMACH 0.3\nITER 100\nPACC\n@POL\n\nASEQ -4 12 0.5"),
 "n23012_cseq":   ("NACA 23012", "VISC 2e6\nITER 100\nPACC\n@POL\n\nCSEQ 0 1.2 0.1"),
 "n0012_trip":    ("NACA 0012", "VPAR\nXTR 0.1 0.2\nN 7\n\nVISC 5e5\nITER 100\nPACC\n@POL\n\nASEQ 0 10 1"),
 "e387_re2e5":    ("e387.dat", "VISC 2e5\nITER 100\nPACC\n@POL\n\nASEQ -2 12 0.5"),
 "dae11_re5e5":   ("dae11.dat", "VISC 5e5\nITER 100\nPACC\n@POL\n\nASEQ 0 10 0.5"),
 "n0012_inv":     ("NACA 0012", "PACC\n@POL\n\nASEQ -5 15 1"),
}

def read_polar(fn):
    rows = []
    if not os.path.exists(fn): return rows
    ok = False
    for l in open(fn):
        if l.strip().startswith("---"): ok = True; continue
        if ok and l.strip():
            rows.append([float(t) for t in l.split()])
    return rows

def run(exe, case, reps):
    g, oper = CASES[case]
    polrel = f"{case}.{exe}.pol"; pol = os.path.join(OUT, polrel)
    best = 1e9
    for r in range(reps):
        if os.path.exists(pol): os.remove(pol)
        inp = script(g, oper, polrel)
        env = dict(os.environ)
        if exe == "foiley" and JOBS is not None:
            env["FOILEY_JOBS"] = str(JOBS)
        t0 = time.perf_counter()
        subprocess.run([EXES[exe]], input=inp.encode(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=OUT, env=env)
        best = min(best, time.perf_counter() - t0)
    return best, read_polar(pol)

def compare(a, b):
    da = {round(r[0], 3): r for r in a}; db = {round(r[0], 3): r for r in b}
    common = sorted(set(da) & set(db))
    mx = [0.0]*7
    for k in common:
        for j in range(1, 7):
            mx[j] = max(mx[j], abs(da[k][j]-db[k][j]))
    return len(a), len(b), len(common), mx

JOBS = None

if __name__ == "__main__":
    reps = 3
    if "--jobs" in sys.argv: JOBS = int(sys.argv[sys.argv.index("--jobs")+1])
    if "--reps" in sys.argv: reps = int(sys.argv[sys.argv.index("--reps")+1])
    filt = [a for i,a in enumerate(sys.argv[1:]) if not a.startswith("--") and not a.isdigit()]
    exes = ["xfoil-dp", "xfoil-sp", "foiley"]
    tot = {e: 0.0 for e in exes}
    for case in CASES:
        if filt and not any(f in case for f in filt): continue
        res = {e: run(e, case, reps) for e in exes}
        for e in exes: tot[e] += res[e][0]
        na, nb, nc, mx = compare(res["xfoil-dp"][1], res["foiley"][1])
        print(f"{case:16s} dp {res['xfoil-dp'][0]*1e3:7.1f}ms  sp {res['xfoil-sp'][0]*1e3:7.1f}ms  rs {res['foiley'][0]*1e3:7.1f}ms "
              f"| pts dp/rs/common {na}/{nb}/{nc}  max|dCL| {mx[1]:.4f} |dCD| {mx[2]:.5f} |dCM| {mx[4]:.4f}")
    print(f"{'TOTAL':16s} dp {tot['xfoil-dp']*1e3:7.1f}ms  sp {tot['xfoil-sp']*1e3:7.1f}ms  rs {tot['foiley']*1e3:7.1f}ms  speedup vs dp {tot['xfoil-dp']/tot['foiley']:.2f}x  vs sp {tot['xfoil-sp']/tot['foiley']:.2f}x")
