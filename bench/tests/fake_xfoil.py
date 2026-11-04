"""Stand-in for an XFoil executable in harness tests.

usage: fake_xfoil.py MODE   (reads an XFoil command script on stdin)
MODE: full | partial | nonfinite | crash | nopolar | hang | stale
"""
import sys, time

mode = sys.argv[1]
lines = sys.stdin.read().splitlines()
pol = lines[lines.index("PACC") + 1]
seq = next(l.split() for l in lines if l.split()[:1] in (["ASEQ"], ["CSEQ"]))
a, b, d = (float(x) for x in seq[1:4])
n = int((b - a) / d + 0.5) + 1
alphas = [a + d * i for i in range(n)]
if mode == "hang":
    time.sleep(60)
if mode == "partial":
    alphas = alphas[: len(alphas) // 2]
if mode != "nopolar":
    with open(pol, "w") as f:
        f.write(" header\n   alpha    CL        CD       CDp       CM     Top_Xtr  Bot_Xtr\n  ------ --------\n")
        for i, al in enumerate(alphas):
            cd = "NaN" if (mode == "nonfinite" and i == 1) else "0.00600"
            f.write(f"  {al:6.3f}  {0.11*al:7.4f}   {cd}   0.00100  -0.0100   0.5000   0.5000\n")
if mode == "stale":
    import os
    os.utime(pol, (0, 0))
sys.exit(3 if mode == "crash" else 0)
