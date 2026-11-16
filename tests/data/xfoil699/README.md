These six polars were generated from MIT's XFOIL 6.99 archive:
https://web.mit.edu/drela/Public/web/xfoil/xfoil6.99.tgz
SHA-256: 5c0250643f52ce0e75d7338ae2504ce7907f2d49a30f921826717b8ac12ebe40

Build with `python bench/build_reference.py --source <extracted-Xfoil> --out <build>`.
The build uses gfortran -O2 -fdefault-real-8 and the headless graphics stub.
It adds the existing RAWDUMP hook and removes the archive's unconditional
post-solve diagnostic file writer (invalid Windows filenames); numerical
calculations and convergence logic are preserved.

Each matching .in file contains the complete input; replace @POL with the
output filename and run it in a fresh directory. Cases exercise default
160-node and 240-node wakes, NACA 4/5-digit generation, separate top/bottom
Ncrit, a 2.5-chord wake, the third LAG parameter, CDp and transition indices.

Trailing whitespace in the Fortran headers is stripped; data values and columns are unchanged.
