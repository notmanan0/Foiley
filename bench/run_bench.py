"""Benchmark + validation harness: runs reference Fortran XFoil and Foiley on
identical command scripts, compares the polar files and reports wall times.

usage: python bench/run_bench.py [case-filter] [--reps N] [--jobs N] [--timeout S]

Every repetition is checked: a nonzero exit, a timeout, or a polar file that
is missing or older than the run marks that repetition failed.  Timings are
reported as median and IQR over the successful repetitions, and the polar
compared is the one written by the median-time repetition.  Each expected
operating point (from the ASEQ/CSEQ command) is scored as matched,
missing-in-ref, missing-in-foiley or nonfinite; a case fails when no point is
in common or coverage is below 95% of the reference's converged points.
The process exits nonzero if any case fails.
"""
import math, os, re, shutil, statistics, subprocess, sys, time
from dataclasses import dataclass, field

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "bench", "out")
# foil files: the reference tree if present, else the tracked copies
RUNS = [os.path.join(ROOT, "xfoil-ref", "runs"), os.path.join(ROOT, "tests", "data", "foils")]
EXES = {
    "xfoil-sp": [os.path.join(ROOT, "build-ref-sp", "xfoil.exe")],
    "xfoil-dp": [os.path.join(ROOT, "build-ref-dp", "xfoil.exe")],
    "foiley":   [os.path.join(ROOT, "target", "release", "foiley.exe")],
}
TIMEOUT = 300.0
MIN_COVERAGE = 0.95
JOBS = None
LABEL = {"xfoil-dp": "dp", "xfoil-sp": "sp", "foiley": "rs"}

def geom(g, out):
    if g.startswith("NACA"):
        return f"{g}\n"
    for d in RUNS:
        src = os.path.join(d, g)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(out, g))
            break
    return f"LOAD {g}\nPANE\n"

def script(g, oper, polfile, out=OUT):
    return ("PLOP\nG F\n\n" + geom(g, out) + "OPER\n" + oper.replace("@POL", polfile) + "\n\nQUIT\n")

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

# ---------------------------------------------------------------- polars

def _num(t):
    try:
        return float(t)
    except ValueError:  # Fortran overflow '********' and the like
        return math.nan

def read_polar(fn):
    rows = []
    if not os.path.exists(fn):
        return rows
    ok = False
    for l in open(fn):
        if l.strip().startswith("---"):
            ok = True
            continue
        if ok and l.strip():
            rows.append([_num(t) for t in l.split()])
    return rows

def expected_points(oper):
    """(column, values) of the operating points a command script asks for:
    column 0 (alpha) for ASEQ, 1 (CL) for CSEQ.  Follows XFoil's own count
    NP = int((last - first)/step + 0.5) + 1."""
    col, vals = None, []
    for l in oper.splitlines():
        t = l.split()
        if not t or t[0].upper() not in ("ASEQ", "CSEQ"):
            continue
        col = 0 if t[0].upper() == "ASEQ" else 1
        a, b, d = (float(x) for x in t[1:4])
        d = -abs(d) if b < a else abs(d)
        n = int((b - a) / d + 0.5) + 1 if d != 0 else 1
        vals += [a + d * i for i in range(max(n, 0))]
    return col, vals

def _key(v):
    return round(v + 0.0, 3)

@dataclass
class Score:
    expected: int = 0
    matched: list = field(default_factory=list)
    missing_in_ref: list = field(default_factory=list)
    missing_in_foiley: list = field(default_factory=list)
    nonfinite: list = field(default_factory=list)
    ref_converged: int = 0
    max_diff: list = field(default_factory=lambda: [0.0] * 7)
    failures: list = field(default_factory=list)

    @property
    def coverage(self):
        return len(self.matched) / self.ref_converged if self.ref_converged else 0.0

    @property
    def ok(self):
        return not self.failures

def _index(rows, col):
    """Rows keyed by the swept column, finite rows only; and the keys of rows
    with any nonfinite value."""
    good, bad = {}, set()
    for r in rows:
        if len(r) <= col or not math.isfinite(r[col]):
            continue
        k = _key(r[col])
        if all(math.isfinite(x) for x in r[:7]):
            good[k] = r
        else:
            bad.add(k)
    return good, bad

def compare(ref, foiley, oper):
    """Score the Foiley polar against the reference on the expected points."""
    col, want = expected_points(oper)
    s = Score(expected=len(want))
    if col is None:
        s.failures.append("no ASEQ/CSEQ in the command script")
        return s
    # CSEQ rows land on the requested CL only to the printed precision
    tol = 0.0 if col == 0 else 0.0015
    def find(idx, v):
        k = _key(v)
        if k in idx:
            return idx[k]
        if tol:
            for kk, r in idx.items():
                if abs(kk - v) <= tol:
                    return r
        return None
    ri, rbad = _index(ref, col)
    fi, fbad = _index(foiley, col)
    s.ref_converged = sum(1 for v in want if find(ri, v) is not None)
    for v in want:
        a, b = find(ri, v), find(fi, v)
        if _key(v) in rbad or _key(v) in fbad:
            s.nonfinite.append(v)
        elif a is None and b is None:
            s.missing_in_ref.append(v)
            s.missing_in_foiley.append(v)
        elif a is None:
            s.missing_in_ref.append(v)
        elif b is None:
            s.missing_in_foiley.append(v)
        else:
            s.matched.append(v)
            for j in range(1, 7):
                s.max_diff[j] = max(s.max_diff[j], abs(a[j] - b[j]))
    if not s.matched:
        s.failures.append("no operating point in common")
    elif s.coverage < MIN_COVERAGE:
        s.failures.append(f"coverage {s.coverage:.1%} of the reference's {s.ref_converged} converged points is below {MIN_COVERAGE:.0%}")
    if s.nonfinite:
        s.failures.append(f"{len(s.nonfinite)} nonfinite point(s)")
    return s

# ---------------------------------------------------------------- runs

@dataclass
class Sample:
    seconds: float
    ok: bool
    reason: str = ""
    rows: list = field(default_factory=list)

@dataclass
class RunResult:
    samples: list

    @property
    def good(self):
        return [s for s in self.samples if s.ok]

    @property
    def ok(self):
        return bool(self.good)

    def median_sample(self):
        """The successful repetition with the median time (lower median)."""
        g = sorted(self.good, key=lambda s: s.seconds)
        return g[(len(g) - 1) // 2] if g else None

    def median(self):
        m = self.median_sample()
        return m.seconds if m else math.nan

    def iqr(self):
        t = sorted(s.seconds for s in self.good)
        if len(t) < 2:
            return 0.0
        q = statistics.quantiles(t, n=4, method="inclusive")
        return q[2] - q[0]

    def rows(self):
        m = self.median_sample()
        return m.rows if m else []

    def reasons(self):
        return sorted({s.reason for s in self.samples if not s.ok})

def run_once(cmd, inp, pol, cwd, env=None, timeout=TIMEOUT):
    if os.path.exists(pol):
        os.remove(pol)
    wall0 = time.time()
    t0 = time.perf_counter()
    try:
        p = subprocess.run(cmd, input=inp.encode(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           cwd=cwd, env=env, timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        return Sample(time.perf_counter() - t0, False, f"timeout after {timeout:g} s")
    dt = time.perf_counter() - t0
    if p.returncode != 0:
        return Sample(dt, False, f"exit code {p.returncode}")
    if not os.path.exists(pol):
        return Sample(dt, False, "no polar written")
    # allow for coarse filesystem timestamps
    if os.path.getmtime(pol) < wall0 - 2.0:
        return Sample(dt, False, "polar older than the run")
    return Sample(dt, True, "", read_polar(pol))

def run(exe, case, reps, out=OUT, timeout=TIMEOUT):
    g, oper = CASES[case]
    polrel = f"{case}.{exe}.pol"
    pol = os.path.join(out, polrel)
    env = dict(os.environ)
    if exe == "foiley" and JOBS is not None:
        env["FOILEY_JOBS"] = str(JOBS)
    inp = script(g, oper, polrel, out)
    return RunResult([run_once(EXES[exe], inp, pol, out, env, timeout) for _ in range(reps)])

def main(argv):
    global JOBS
    reps, timeout = 3, TIMEOUT
    if "--jobs" in argv: JOBS = int(argv[argv.index("--jobs") + 1])
    if "--reps" in argv: reps = int(argv[argv.index("--reps") + 1])
    if "--timeout" in argv: timeout = float(argv[argv.index("--timeout") + 1])
    vals = {argv[i + 1] for i, a in enumerate(argv[:-1]) if a in ("--jobs", "--reps", "--timeout")}
    filt = [a for a in argv if not a.startswith("--") and a not in vals]
    os.makedirs(OUT, exist_ok=True)
    exes = ["xfoil-dp", "xfoil-sp", "foiley"]
    tot = {e: 0.0 for e in exes}
    failed = []
    for case in CASES:
        if filt and not any(f in case for f in filt): continue
        res = {e: run(e, case, reps, OUT, timeout) for e in exes}
        for e in exes:
            tot[e] += res[e].median() if res[e].ok else 0.0
        t = "  ".join(f"{LABEL[e]} {res[e].median()*1e3:7.1f}ms (IQR {res[e].iqr()*1e3:5.1f})" for e in exes)
        bad = [f"{e}: {', '.join(res[e].reasons())}" for e in exes if res[e].reasons()]
        if not (res["xfoil-dp"].ok and res["foiley"].ok):
            print(f"{case:16s} {t} | FAILED run: {'; '.join(bad)}")
            failed.append(case)
            continue
        s = compare(res["xfoil-dp"].rows(), res["foiley"].rows(), CASES[case][1])
        print(f"{case:16s} {t} | expected {s.expected} matched {len(s.matched)} "
              f"miss-ref {len(s.missing_in_ref)} miss-rs {len(s.missing_in_foiley)} nonfinite {len(s.nonfinite)} "
              f"cov {s.coverage:.0%} max|dCL| {s.max_diff[1]:.4f} |dCD| {s.max_diff[2]:.5f} |dCM| {s.max_diff[4]:.4f}"
              + (f" | FAILED: {'; '.join(s.failures)}" if s.failures else "")
              + (f" | partial runs: {'; '.join(bad)}" if bad else ""))
        if s.failures:
            failed.append(case)
    if tot["foiley"] > 0:
        print(f"{'TOTAL (median)':16s} dp {tot['xfoil-dp']*1e3:7.1f}ms  sp {tot['xfoil-sp']*1e3:7.1f}ms  rs {tot['foiley']*1e3:7.1f}ms  "
              f"speedup vs dp {tot['xfoil-dp']/tot['foiley']:.2f}x  vs sp {tot['xfoil-sp']/tot['foiley']:.2f}x")
    if failed:
        print(f"FAILED cases: {', '.join(failed)}")
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
