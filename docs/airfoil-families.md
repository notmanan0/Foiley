# Airfoil families

Foiley 1.2 generates foils from these families. Every generated foil records
how it was made (its *recipe*), so Studio can reopen the parameters later
with **Create → Edit parameters**. Any operation that changes the shape
clears the recipe: a flap, a TE gap, scaling, the optimiser and the inverse
methods. Repaneling keeps it.

| Family | Studio | Command line | Library |
|---|---|---|---|
| NACA (all below) | Create → NACA | `foiley foil gen naca "<designation>"`, or any `"NACA …"` foil argument | `Foil::naca`, `gen::naca::parse` |
| CST (Kulfan) | Create → CST, *Fit to selected foil* | `foiley foil gen cst`, `foiley foil fit-cst` | `gen::cst::Cst` |
| PARSEC | Create → PARSEC | `foiley foil gen parsec` | `gen::parsec::Parsec` |
| Kármán–Trefftz | Create → Karman-Trefftz | `foiley foil gen kt [--exact a1,a2]` | `gen::kt::KarmanTrefftz` |
| Cambered plate | Create → cambered plate | `foiley foil gen plate` | `gen::plate::Plate` |

## NACA designations

The "NACA" prefix, spaces and letter case are optional.

| Family | Examples | Meaning |
|---|---|---|
| 4-digit | `2412`, `0012` | camber %, position of max camber (tenths), t/c % |
| 5-digit | `23012`, `43012` | design CL = 0.15 × 1st digit, max camber at 2nd digit / 20, 3rd digit 0, t/c % |
| 5-digit reflexed | `23112`, `25112` | 3rd digit 1: reflexed mean line with zero moment about c/4 (positions 221 to 251) |
| 4-digit modified | `0012-64`, `2412-63` | 4-digit camber and t/c, then LE radius index I (0 to 9) and position of max thickness J/10 (J = 2 to 6) |
| 16-series | `16-009`, `16-212` | design CL in tenths, t/c %: thickness `00tt-45`, a = 1.0 mean line |
| 6-series | `63-206`, `64(2)-415`, `65_2-415`, `652-415`, `63-412 a=0.5` | thickness family 63 to 67, low-drag-range subscript (name only), design CL in tenths, t/c %, mean line a = 1.0 unless `a=` is given |
| 6A-series | `64A010`, `63A210` | 63A, 64A or 65A thickness; modified a = 0.8 mean line; straight surfaces aft of x/c = 0.82 |

Add `closed TE` to a 4- or 5-digit code (`2412 closed TE`) to use −0.1036
for the x⁴ thickness coefficient instead of −0.1015. This closes the
trailing edge, which is useful for 3D meshes and STL export.

NACA defined the 6A-series only for 63A, 64A and 65A. The 6-series thickness
data run from 63 to 67.

### How the shapes are built

Each family is built as its NACA report defines it. The thickness is laid
off **perpendicular to the mean line** (Abbott & von Doenhoff, eq. 6.1):

    xu = x − yt sin θ,  yu = yc + yt cos θ,  xl = x + yt sin θ,  yl = yc − yt cos θ

* **4-digit modified and 16-series thickness.** Abbott & von Doenhoff
  p. 117:
  * LE radius r = 1.1019 (t·I/6)².
  * The trailing-edge slope d1 comes from the table for the position of
    maximum thickness.
  * The a and d polynomial coefficients are solved from the continuity
    conditions at that position.
* **5-digit mean lines.** NASA TM-4741:
  * The standard line uses the tabulated r and k1 for CL = 0.3, scaled by
    the design CL.
  * The reflexed line takes k2/k1 from the zero-moment condition of NACA
    Report 537. The value printed in TM-4741 is wrong, as NACA456 notes.
* **a-type mean lines.** Abbott & von Doenhoff eq. 4-26/4-27.
  * The 6A-series uses NACA Report 903's modified a = 0.8 line: the
    ordinates are scaled by 0.97948 and run straight at slope −0.24521 CL
    aft of x/c = 0.86.
* **6- and 6A-series thickness.** Conformal mapping of a circle through the
  NASA Langley ε(φ) and ψ(φ) functions, 201 points per family, taken from
  Ralph Carmichael's public-domain NACA456 (PDAS).
  * NACA456 starts from its fitted scale factor. Foiley then iterates the
    scale until the maximum thickness equals t/c exactly, as NASA's
    original NACA6 program did.
  * For the 6A-series, both surfaces aft of x/c = 0.82 are replaced by
    straight lines to the trailing edge, as NACA6 does.

The codes XFoil's own `NACA` command accepts (4-digit, and 5-digit 210xx to
250xx) still come from XFoil's routine, so they are unchanged from Foiley
1.1 and match XFoil and AVL. That routine adds the thickness vertically,
so a cambered section differs slightly near the nose from the perpendicular
construction. Saved projects store coordinates, so they never change.

The new families use 239 points with 120 stations per side, bunched at the
nose only (half-cosine spacing, about 0.013 c at the trailing edge, close to
NACA's own tables). See "Near-cusped trailing edges" below for why.

### Validation against NASA's programs

`crates/foil/tests/naca_reference.rs` compares the generators with the
ordinates printed by NASA Langley's original NACA4 and NACA6 programs.
Those outputs are distributed with NACA456 and are public domain. They are
5-decimal tables (`crates/foil/tests/data/naca`). Largest differences:

| Section | aft of x/c 0.0125 | nose (x/c < 0.0125) |
|---|---|---|
| 23012 | 5.5e-6 | 5.0e-6 |
| 23112 (reflexed) | 6.4e-6 | 5.0e-6 |
| 0010-64 | 1.6e-5 | 1.4e-5 |
| 16-009 | 1.2e-5 | 1.0e-5 |
| 63-006 | 8.0e-6 | 4.3e-5 |
| 63A010 | 1.1e-5 | 3.1e-5 |
| 63-206 | 8.1e-6 | 6.3e-5 |
| 63A409 | 1.1e-5 | 1.9e-4 |

* **4-digit modified sections.** The NASA runs were given the LE radius
  rounded to three figures (0.0110 where the formula gives 0.011019), which
  accounts for the 4-digit-modified differences.
* **The nose.** NASA's program splines the mapped thickness and evaluates
  the logarithmic mean-line slope differently there, which accounts for the
  larger differences at x/c < 0.0125.
* **16-212.** The 16-212 run in that package used an a = 0.8 (6A) mean line
  rather than the a = 1.0 line that defines the 16-series, so it is not
  used.

### Near-cusped trailing edges

The 6-series thickness forms close in a near-cusp: the included angle is
under 1°. The a = 1.0 mean line's slope also runs to −∞ (logarithmically)
at the trailing edge.

**Where XFoil fails.** For thin sections with high design CL and a = 1.0,
XFoil's sharp-trailing-edge Kutta condition can put its control point
outside the body. The inviscid solution is then wrong: CL of order 10 or
more. NACA 63-406 shows this, even with NASA's own ordinates. NACA 63-206,
63-412, 64(2)-415, 65-610 and 66-209 solve normally.

**Workarounds.**
* Add a small trailing-edge gap: Modify → Leading / trailing edge → TE gap
  0.0005, or `foiley foil design --te-gap 0.0005`.
* Use an a < 1 mean line, e.g. `63-406 a=0.8`.

## CST (Kulfan)

Each surface is defined by

    y(x) = x^N1 (1−x)^N2 · Σ Aᵢ Kᵢ xⁱ (1−x)ⁿ⁻ⁱ + x·Δz

with N1 = 0.5 and N2 = 1 for a round nose and a sharp tail. The weights
control the shape:

* **A₀** sets the LE radius (r = A₀²/2).
* **Aₙ** sets the trailing-edge wedge angle.
* **Δz** is the trailing-edge height of that surface.

**Fitting.** *Fit* (`Cst::fit`, `foiley foil fit-cst`) normalises the foil
and solves for the weights of each surface by least squares (QR). It
reports the RMS and largest vertical error.
* An order-8 fit of NACA 2412 has an RMS error of 4e-5 and a maximum of
  2e-4 of chord.
* `foiley foil fit-cst FOIL --json` prints the recipe, which
  `foiley foil gen cst --json FILE` reads back.

**In the optimiser**, select *CST weights* in Studio, or
`foiley foil optimize --shape cst --cst-order 5 --cst-delta 0.05`. The genes
are changes of the weights, added to the base foil as
`dy = x^0.5 (1−x) Σ ΔAᵢ Bᵢ(x)`. Zero genes therefore reproduce the base
foil exactly, and the fit error never enters the first design. Each gene
moves thickness, camber or the LE radius over the whole chord, without the
local waviness of Hicks–Henne bumps.

## PARSEC

Sobieczky's 11 parameters:
* the LE radius;
* the upper and lower crest positions, heights and curvatures;
* the TE height and thickness;
* the TE direction (degrees below the chord line);
* the wedge angle.

Each surface is z = Σ aₙ x^(n−½) (n = 1…6). Its six coefficients come from
one 6×6 linear solve. The defaults are close to NACA 0012. Foiley rejects
inconsistent inputs and surfaces that cross. PARSEC cannot represent
reflexed or strongly aft-loaded sections well.

## Kármán–Trefftz

The section is the map

    z = n ((ζ+1)ⁿ + (ζ−1)ⁿ) / ((ζ+1)ⁿ − (ζ−1)ⁿ),  n = 2 − τ/π

of a circle through ζ = 1, centred at (−x, y). Here τ is the trailing-edge
angle, and τ = 0 gives the Joukowsky section.

**Exact solution.** The potential flow about the circle is known, so
`KarmanTrefftz::exact(alpha)` returns the exact inviscid CL, Cm about c/4
and Cp. CL is lift relative to the chord of the normalised section; Cm and
Cp are integrated from the exact surface speed. Studio shows CL and Cm next
to the parameters, and `foiley foil gen kt --exact 0,4` prints them.

**Check against XFoil.** `crates/foil/tests/karman_trefftz.rs` checks
XFoil's inviscid solution against the exact one on three sections:
* CL agrees within 0.3%;
* Cm agrees within 4e-4.

## Cambered plates

A plate of constant thickness bent to a circular-arc or parabolic camber
line (the NACA 4-digit mean line, with its position of maximum camber).
The leading edge is a semicircle of diameter equal to the thickness. This
family is meant for very low Reynolds numbers: indoor models, micro air
vehicles, and sail or wing-section baselines.

The trailing edge can be **square** (the default) or **round**:

* **Square:** cut square, a blunt edge as thick as the plate. XFoil
  analyses this one; the default 2%, 4%-camber plate converges at Re 1e5
  from 0° to 6°.
* **Round:** a semicircle like the leading edge. Use it for geometry, such
  as 3D meshes and STL. XFoil's viscous solution does not converge on it at
  any Reynolds number tried (1e5 to 5e5), because the boundary layer cannot
  follow a cap that closes to a point. A tapered trailing edge fails too.
  The inviscid solution is fine.

## Sources

* I. H. Abbott and A. E. von Doenhoff, *Theory of Wing Sections*, Dover 1959.
* NACA Reports 460, 537, 824 and 903.
* C. L. Ladson, C. W. Brooks, A. S. Hill and D. W. Sproles, *Computer Program
  to Obtain Ordinates for NACA Airfoils*, NASA TM-4741, 1996.
* R. L. Carmichael, *Algorithm for Calculating Coordinates of Cambered NACA
  Airfoils at Specified Chord Locations*, AIAA 2001-5235. NACA456 source,
  ε/ψ tables and NASA reference outputs:
  [pdas.com/naca456.html](https://www.pdas.com/naca456.html), public domain.
* B. M. Kulfan, *Universal Parametric Geometry Representation Method*,
  J. Aircraft 45(1), 2008.
* H. Sobieczky, *Parametric Airfoils and Wings*, Notes on Numerical Fluid
  Mechanics 68, 1998.
