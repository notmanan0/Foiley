"""Build a Foiley project of the Spitfire Mk I and a job script that fills it with a 6-DoF aerodynamic database.

    python spitfire/make_spitfire.py            # writes Spitfire_MkI_geometry.foiley and spitfire_6dof_job.json
    target/release/foiley run spitfire/spitfire_6dof_job.json    # runs every case -> Spitfire_MkI_6DoF.foiley
    python spitfire/make_spitfire.py --tables   # after the run: derivative tables from the saved project

How it works
------------
1. Geometry. The wing, tailplane, fin, fuselage and spinner come from the Individual Project's
   `04 Aero Model/cad_geometry.json` (the Spitfire CAD repository's published dimensions, the same geometry flow5_run.py
   gives flow5). The control surfaces come from the CAD repository's `spitfire_cad/aero_config.py` (copied below
   with its values, so this script does not import the CAD package). Inches in, metres out.
   Foiley axes: x aft, y starboard, z up, origin at the CAD datum on the thrust line.

2. Control surfaces. In Foiley a wing panel is a flap when the foils at both of its sections carry a hinge, and every
   such panel is a separate flap, numbered per wing from the left tip to the right tip. So each control surface is a
   run of hinged sections, and the job deflects all of its panels together:
     * ailerons: s = 127 to 209.5 in, hinge straight in plan view from 78.5 % to 66 % chord; 5 panels per side.
       Antisymmetric: +delta_a = port (left) trailing edge down, starboard up, which rolls the aircraft to starboard
       (the CAD convention).
     * elevator: z = 4.5 in to the tip, hinge at x = 237 in; 7 panels per side. +delta_e = trailing edge down.
     * rudder: fin root (y = 12 in, the top of the tail cone) to the hinge's top end (y = 42 in), along the CAD's
       inclined hinge line; 7 panels. +delta_r = trailing edge to port (the CAD convention); the sign of Foiley's fin
       flap is checked against this (FIN_FLAP_SIGN).
   Horn balances are not deflected (the parts of the elevator and rudder ahead of their hinges stay fixed), the rudder
   below the tail cone is left out (the fin starts at the top of the cone, as in the CAD's fin table), and the split
   flaps are not modelled (a lattice cannot represent a split flap, and they are not a 6-DoF control).

3. Cases (all VLM2 ring-vortex lattice, viscous drag from XFoil polars of each section, V = 94 m/s, sea level),
   0.5 deg steps throughout unless stated:
     * "grid clean": alpha -10 to 20 deg x beta -15 to 15 deg, every control at zero (3,721 points);
     * "grid de-30.0" ... "grid dr+30.0": the same grid at every deflection of each control over its full travel, one
       control at a time (110 elevator, 100 aileron, 120 rudder polars);
     * "combo de.. da.. dr..": all three controls deflected together, every combination at 5 deg steps over the
       full travel (1,716 polars), on a 5 deg alpha x beta grid. These measure how far the control increments
       add: a 0.5 deg three-way set would be 1.36 million separate lattices;
     * "speed grid clean": alpha x speed 40 to 160 m/s (speed changes only the viscous drag, through Re);
     * "rates alpha sweep": the stability derivatives at every alpha (rate derivatives p, q, r, and u, alpha, beta);
     * T7 trimmed polars: an elevator trim curve (0.05 deg steps, the range that trims), and the aileron and
       rudder control derivatives at trim.
   A 6-DoF table is usually built up this way: the clean grid gives the base forces and moments, each control's
   grid gives its increment against the clean grid at the same (alpha, beta), the combination set shows what is
   lost by adding increments, and the rate derivatives add the damping terms. Span distributions are not kept
   (keep_strips false): at 55 kB a point they would make the project about 70 GB.

Why VLM2: a lattice handles sideslip, dihedral and deflected panels, and is fast enough for every combination; the
ring lattice is the flow5 default the 3-DOF model used. Foiley's lattices deflect a control by turning the flap
panels' normals about the hinge (the small-disturbance model of AVL), not by moving the panels: with the panels
moved, VLM2's flap increments fell with alpha (symmetric ailerons +10 deg: 0.151 at alpha -6 deg, 0.026 at +10 deg),
because the ring sides at a flap's side edge separated into a close vortex pair (fixed 2 Nov,
MeshOptions::flap_normals). It is inviscid in lift: nothing here stalls, so results beyond about 12 deg alpha (the
plan's strip-stall edge) are the linear extrapolation, not the aircraft.

Sign conventions in the output. The polar tables' Cl_roll and Cn are moments about Foiley's x-aft and z-up wind
axes, so they have the opposite sign to the usual flight-mechanics Cl (+ right wing down) and Cn (+ nose right);
CY (+ starboard) and Cm (+ nose up) are the usual ones. The T7 derivatives and control derivatives are in the usual
stability axes (x forward, z down): Clp < 0, Cnb > 0, aileron L > 0, rudder N < 0.
"""
import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
GEOMETRY = pathlib.Path(r"C:\Users\Manan Malavia\Individual Project\04 Aero Model\cad_geometry.json")
IN = 0.0254

PLANE_NAME = "Spitfire Mk I"
GEOMETRY_FILE = "Spitfire_MkI_geometry.foiley"
RESULT_FILE = "Spitfire_MkI_6DoF.foiley"
JOB_FILE = "spitfire_6dof_job.json"

# --- control surfaces, from Spitfire CAD spitfire_cad/aero_config.py (inches, degrees)
AILERON = dict(s_in=127.0, s_out=209.5, frac_in=0.785, frac_out=0.66, travel=(-25.0, 25.0))
ELEVATOR = dict(hinge_x=237.0, z_root=4.5, travel=(-30.0, 25.0))
RUDDER = dict(hinge_top=(235.25, 41.98), hinge_bottom=(246.03, -15.30), travel=(-30.0, 30.0))
FIN_TC = (0.09, 0.07)

# --- stations (in). Fewer than the CAD's 68 wing stations, because each hinged panel is its own flap.
WING_S = [0.0, 15.5, 31.0, 50.0, 69.0, 88.0, 107.0, 127.0, 143.5, 160.0, 176.5, 193.0, 209.5, 214.0, 218.0, 220.5, 222.4]
TAIL_Z = [0.0, 3.23, 4.4, 4.5, 12.9, 25.8, 38.8, 48.4, 54.9, 59.8, 63.0]
FIN_Y = [12.03, 13.85, 15.7, 20.0, 26.0, 32.0, 37.0, 41.98, 44.5, 46.7, 47.9]
WING_NX, TAIL_NX, FIN_NX = 12, 10, 10
WING_DY, TAIL_DY, FIN_DY = 6.0, 3.0, 3.0      # spanwise panel size (in)

# --- flight condition and references (the 3-DOF model's, 04 Aero Model/aero_model.json)
V = 94.0
RHO, NU, A_SOUND = 1.225, 1.46e-5, 340.294
REF = dict(area=22.36418496323411, span=11.244655731293165, chord=1.9888723583591803)

# --- the analysis grid (deg), 0.5 deg steps throughout
def steps(a, b, d):
    return [round(a + i * d, 3) for i in range(int(round((b - a) / d)) + 1)]


STEP = 0.5
ALPHAS = steps(-10.0, 20.0, STEP)                 # 61
BETAS = steps(-15.0, 15.0, STEP)                  # 61
DE_SET = [d for d in steps(*ELEVATOR["travel"], STEP) if d != 0.0]    # 110, full travel
DA_SET = [d for d in steps(*AILERON["travel"], STEP) if d != 0.0]     # 100
DR_SET = [d for d in steps(*RUDDER["travel"], STEP) if d != 0.0]      # 120
# Every control at once: each combination is its own lattice (new normals, new factorisation), so the three-way
# combinations use 5 deg steps over the full travel (12 x 11 x 13 = 1716) on a 5 deg alpha x beta grid; at 0.5 deg
# they would be 111 x 101 x 121 = 1.36 million lattices.
COMBO_STEP = 5.0
COMBO_ALPHAS = steps(-10.0, 20.0, COMBO_STEP)
COMBO_BETAS = steps(-15.0, 15.0, COMBO_STEP)
SPEEDS = steps(40.0, 160.0, 10.0)                 # clean alpha x speed grid: speed enters only the viscous drag
# Foiley's fin flap is + trailing edge "down" in the fin's own axes. The fin is the right half of a one-sided wing
# rotated +90 deg about x, so its local "down" (-z) becomes +y: + fin flap = trailing edge to starboard, the
# opposite of the CAD's rudder sign. Checked by the run: + delta_r must give a nose-left (negative) Cn.
FIN_FLAP_SIGN = -1.0
# Elevator preset for the aileron and rudder T7 polars, chosen from the elevator trim curve so that they trim near
# 94 m/s (see the T7 elevator polar).
DE_TRIM = 0.985           # from the trim curve: delta_e +0.9 trims at 87.4 m/s, +1.0 at 95.2 m/s (VLM2)


def naca4(m, p, t, n=81):
    """NACA four-digit foil, Selig order (TE upper -> LE -> TE lower), cosine spacing, closed trailing edge."""
    b = np.linspace(0.0, math.pi, n)
    x = 0.5 * (1.0 - np.cos(b))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)
    if m > 0:
        yc = np.where(x < p, m / p ** 2 * (2 * p * x - x ** 2), m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2))
        dy = np.where(x < p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    else:
        yc, dy = np.zeros_like(x), np.zeros_like(x)
    th = np.arctan(dy)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    xu[-1] = xl[-1] = 1.0
    yu[-1] = yl[-1] = 0.0
    xs = np.r_[xu[::-1], xl[1:]]
    ys = np.r_[yu[::-1], yl[1:]]
    return [round(float(v), 7) for v in xs], [round(float(v), 7) for v in ys]


def foil(name, m, p, t, hinge=None, note=""):
    x, y = naca4(m, p, t)
    return {"name": name, "x": x, "y": y, "hinge": None if hinge is None else [round(hinge, 5), 0.5], "note": note}


def section(y, chord, offset, foil_name, nx, ny, dihedral=0.0, twist=0.0):
    return {"y": y, "chord": chord, "offset": offset, "dihedral": dihedral, "twist": twist,
            "foil_left": foil_name, "foil_right": foil_name, "nx": nx, "ny": ny, "xdist": "Cosine", "ydist": "Uniform"}


def vec(x, y, z):
    return {"x": x, "y": y, "z": z}


def ny_for(stations, i, dy):
    return 1 if i + 1 >= len(stations) else max(1, round((stations[i + 1] - stations[i]) / dy))


def build_plane(geo):
    w, t, f, fu = geo["wing"], geo["tailplane"], geo["fin"], geo["fuselage"]
    foils, flaps = [], {}

    # ---- wing: NACA 22xx (2 % camber at 20 % chord), t/c from the CAD table; ailerons hinged straight in plan view
    s_all = np.array(w["s"])
    at = lambda key, s: float(np.interp(s, s_all, w[key]))
    le0 = at("le_x", 0.0)
    hx_in = at("le_x", AILERON["s_in"]) + AILERON["frac_in"] * at("chord", AILERON["s_in"])
    hx_out = at("le_x", AILERON["s_out"]) + AILERON["frac_out"] * at("chord", AILERON["s_out"])
    secs, n_ail = [], 0
    for i, s in enumerate(WING_S):
        tc, c, le = at("tc", s), at("chord", s), at("le_x", s)
        hinged = AILERON["s_in"] - 1e-6 <= s <= AILERON["s_out"] + 1e-6
        if hinged:
            hx = hx_in + (hx_out - hx_in) * (s - AILERON["s_in"]) / (AILERON["s_out"] - AILERON["s_in"])
            frac = (hx - le) / c
            name = f"Wing s{s:05.1f} aileron"
            foils.append(foil(name, 0.02, 0.2, tc, frac, f"NACA 22xx t/c {tc:.4f}, aileron hinge {frac:.3f} c"))
            n_ail += 1
        else:
            name = f"Wing s{s:05.1f}"
            foils.append(foil(name, 0.02, 0.2, tc, None, f"NACA 22xx t/c {tc:.4f}"))
        secs.append(section(s * IN, c * IN, (le - le0) * IN, name, WING_NX, ny_for(WING_S, i, WING_DY),
                            dihedral=0.0 if s < w["dihedral_s"] - 1e-6 else w["dihedral_deg"],
                            twist=at("incidence_deg", s)))
    wing = {"name": "Main wing", "kind": "Main", "sections": secs,
            "position": vec(le0 * IN, 0.0, w["chord_plane_y_at_front_spar"] * IN),
            "rx": 0.0, "ry": 0.0, "two_sided": True, "symmetric": True, "close_inner": False, "mass": 0.0,
            "point_masses": []}
    k = n_ail - 1                         # aileron panels per side
    flaps["aileron_left"] = list(range(0, k))           # left side is numbered tip -> root
    flaps["aileron_right"] = list(range(k, 2 * k))      # right side root -> tip

    # ---- tailplane: symmetric 9 % -> 7.5 %, elevator aft of x = 237 in outboard of z = 4.5 in
    z_all = np.array(t["z"])
    tle0 = float(np.interp(0.0, z_all, t["le_x"]))
    tsecs, n_el = [], 0
    for i, z in enumerate(TAIL_Z):
        le, te = float(np.interp(z, z_all, t["le_x"])), float(np.interp(z, z_all, t["te_x"]))
        tc = t["tc_root"] + (t["tc_tip"] - t["tc_root"]) * min(z / z_all[-1], 1.0)
        if z >= ELEVATOR["z_root"] - 1e-6:
            frac = (ELEVATOR["hinge_x"] - le) / (te - le)
            name = f"Tail z{z:04.1f} elevator"
            foils.append(foil(name, 0, 0, tc, frac, f"NACA 00xx t/c {tc:.4f}, elevator hinge {frac:.3f} c"))
            n_el += 1
        else:
            name = f"Tail z{z:04.1f}"
            foils.append(foil(name, 0, 0, tc, None, f"NACA 00xx t/c {tc:.4f}"))
        tsecs.append(section(z * IN, (te - le) * IN, (le - tle0) * IN, name, TAIL_NX, ny_for(TAIL_Z, i, TAIL_DY),
                             twist=t["incidence_deg"]))
    tail = {"name": "Tailplane", "kind": "Elevator", "sections": tsecs,
            "position": vec(tle0 * IN, 0.0, t["chord_plane_y"] * IN),
            "rx": 0.0, "ry": 0.0, "two_sided": True, "symmetric": True, "close_inner": False, "mass": 0.0,
            "point_masses": []}
    k = n_el - 1
    flaps["elevator_left"] = list(range(0, k))
    flaps["elevator_right"] = list(range(k, 2 * k))

    # ---- fin: one-sided wing rotated +90 deg about x; rudder on the CAD's inclined hinge line up to its top end
    y_all = np.array(f["y"])
    fle0 = float(np.interp(FIN_Y[0], y_all, f["le_x"]))
    (xt, yt), (xb, yb) = RUDDER["hinge_top"], RUDDER["hinge_bottom"]
    fsecs, n_ru = [], 0
    for i, y in enumerate(FIN_Y):
        le, te = float(np.interp(y, y_all, f["le_x"])), float(np.interp(y, y_all, f["te_x"]))
        tc = FIN_TC[0] + (FIN_TC[1] - FIN_TC[0]) * (y - FIN_Y[0]) / (FIN_Y[-1] - FIN_Y[0])
        if y <= yt + 1e-6:
            hx = xb + (xt - xb) * (y - yb) / (yt - yb)
            frac = (hx - le) / (te - le)
            name = f"Fin y{y:04.1f} rudder"
            foils.append(foil(name, 0, 0, tc, frac, f"NACA 00xx t/c {tc:.4f}, rudder hinge {frac:.3f} c"))
            n_ru += 1
        else:
            name = f"Fin y{y:04.1f}"
            foils.append(foil(name, 0, 0, tc, None, f"NACA 00xx t/c {tc:.4f} (above the rudder hinge: horn, fixed)"))
        fsecs.append(section((y - FIN_Y[0]) * IN, (te - le) * IN, (le - fle0) * IN, name, FIN_NX,
                             ny_for(FIN_Y, i, FIN_DY)))
    fin = {"name": "Fin", "kind": "Fin", "sections": fsecs, "position": vec(fle0 * IN, 0.0, FIN_Y[0] * IN),
           "rx": 90.0, "ry": 0.0, "two_sided": False, "symmetric": True, "close_inner": True, "mass": 0.0,
           "point_masses": []}
    flaps["rudder"] = list(range(0, n_ru - 1))

    # ---- fuselage and spinner, elliptic frames (used by the panel methods and the 3D view; VLM ignores it)
    frames = [dict(x=x * IN, z=0.0, half_width=r * IN, half_height=r * IN)
              for x, r in zip(fu["spinner_x"], fu["spinner_r"]) if x < fu["x"][0] - 0.5]
    for x, top, bot, hw in list(zip(fu["x"], fu["y_top"], fu["y_bot"], fu["half_width"]))[::3] + \
            [(fu["x"][-1], fu["y_top"][-1], fu["y_bot"][-1], fu["half_width"][-1])]:
        frames.append(dict(x=x * IN, z=0.5 * (top + bot) * IN, half_width=max(hw, 0.05) * IN,
                           half_height=max(0.5 * (top - bot), 0.05) * IN))
    fuse = {"name": "Fuselage", "position": vec(0.0, 0.0, 0.0), "frames": frames, "exponent": 2.0,
            "nx": 40, "nh": 10, "mass": 0.0, "nurbs": None, "triangles": []}

    plane = {"name": PLANE_NAME,
             "description": "Spitfire Mk I from the Spitfire CAD repository (cad_geometry.json, aero_config.py): "
                            "wing with ailerons, tailplane with elevator, fin with rudder, fuselage. "
                            "Built by spitfire/make_spitfire.py.",
             "wings": [wing, tail, fin], "fuses": [fuse], "point_masses": []}
    return plane, foils, flaps


def inertia(geo):
    """CAD mass, CG and inertia in Foiley axes. The CAD body tensor (x fwd, y stbd, z down) has off-diagonals equal
    to minus the products (CAD mass_properties.md); Foiley stores the products +integral(x z dm) in x aft, z up, so
    Ixz keeps the CAD product's sign and Ixy, Iyz change sign. The CG's 0.9 in starboard offset is dropped so the
    clean aircraft is symmetric (it would put a rolling moment of the weight into every point)."""
    m = geo["mass"]
    tb = m["tensor_body_kgm2"]
    cg = m["cg_model_inches"]
    return {"mass": m["mass_kg"], "cog": vec(cg["X_aft_of_datum"] * IN, 0.0, cg["Y_above_thrust_line"] * IN),
            "ixx": tb[0][0], "iyy": tb[1][1], "izz": tb[2][2],
            "ixz": -tb[0][2], "ixy": tb[0][1], "iyz": -tb[1][2]}


def deflections(flaps, control, d):
    """Per-panel Foiley flap angles for a control deflection d (deg) in the CAD convention."""
    out = []
    if control == "elevator":
        for k in ("elevator_left", "elevator_right"):
            out += [[[1, i], d] for i in flaps[k]]
    elif control == "aileron":
        out += [[[0, i], d] for i in flaps["aileron_left"]] + [[[0, i], -d] for i in flaps["aileron_right"]]
    else:
        out += [[[2, i], FIN_FLAP_SIGN * d] for i in flaps["rudder"]]
    return out


def gains(flaps, control):
    """T7 gains: degrees of each panel per unit control parameter, so the control derivatives are per degree."""
    return [[{"Flap": {"wing": w, "index": i}}, angle] for (w, i), angle in deflections(flaps, control, 1.0)]


def job(flaps, mass):
    mu = RHO * NU
    base = {"method": "Vlm2", "velocity": V, "viscous": True, "viscous_from_cl": True,
            "fluid": {"rho": RHO, "mu": mu, "a": A_SOUND}, "reference": {"Custom": REF}, "inertia": mass,
            "include_fuse": False, "keep_cp": False, "keep_strips": False}
    grid = [[a, b, V] for a in ALPHAS for b in BETAS]
    out = [{"plane_polar": {"def": dict(base, name="grid clean", kind="List", list=grid)}},
           {"plane_polar": {"def": dict(base, name="speed grid clean", kind="List",
                                        list=[[a, 0.0, v] for v in SPEEDS for a in ALPHAS])}}]
    for control, tag, values in (("elevator", "de", DE_SET), ("aileron", "da", DA_SET), ("rudder", "dr", DR_SET)):
        for d in values:
            out.append({"plane_polar": {"def": dict(base, name=f"grid {tag}{d:+.1f}", kind="List", list=grid,
                                                    flaps=deflections(flaps, control, d))}})
    combo_grid = [[a, b, V] for a in COMBO_ALPHAS for b in COMBO_BETAS]
    for de in steps(*ELEVATOR["travel"], COMBO_STEP):
        for da in steps(*AILERON["travel"], COMBO_STEP):
            for dr in steps(*RUDDER["travel"], COMBO_STEP):
                f = deflections(flaps, "elevator", de) + deflections(flaps, "aileron", da) +                     deflections(flaps, "rudder", dr)
                out.append({"plane_polar": {"def": dict(base, name=f"combo de{de:+.0f} da{da:+.0f} dr{dr:+.0f}",
                                                        kind="List", list=combo_grid, flaps=f)}})
    out.append({"plane_polar": {"def": dict(base, name="rates alpha sweep", kind="FixedSpeed", derivatives=True,
                                            range={"from": ALPHAS[0], "to": ALPHAS[-1], "step": STEP})}})
    out.append({"plane_polar": {"def": dict(base, name="T7 elevator trim", kind="Stability",
                                            gains=gains(flaps, "elevator"),
                                            range={"from": -0.2, "to": 1.6, "step": 0.05})}})
    for control in ("aileron", "rudder"):
        out.append({"plane_polar": {"def": dict(base, name=f"T7 {control}", kind="Stability",
                                                  gains=gains(flaps, control),
                                                  flaps=deflections(flaps, "elevator", DE_TRIM),
                                                  range={"from": 0.0, "to": 0.0, "step": 1.0})}})
    return {"project": GEOMETRY_FILE, "save": RESULT_FILE, "output_dir": "results", "steps": out}


def tables():
    """Derivative tables from the saved project: rates_alpha_sweep.csv and t7_control_derivatives.csv."""
    proj = json.loads((HERE / RESULT_FILE).read_text())
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    for pol in proj["plane_polars"]:
        pts = [p for p in pol["points"] if p.get("stability")]
        if not pts:
            continue
        keys = list(pts[0]["stability"]["nondim"].keys())
        rows = ["alpha,V,CL,Cm,static_margin," + ",".join(keys) + ",ctrl_X,ctrl_Y,ctrl_Z,ctrl_L,ctrl_M,ctrl_N"]
        for p in pts:
            s = p["stability"]
            c = s.get("control") or [float("nan")] * 6
            rows.append(",".join(f"{v:.6g}" for v in [p["alpha"], p["velocity"], p["cl"], p["cm"], s["static_margin"]]
                                 + [s["nondim"][k] for k in keys] + list(c)))
        (out / f"{pol['name'].replace(' ', '_')}_derivatives.csv").write_text("\n".join(rows) + "\n")
        print("wrote", out / f"{pol['name'].replace(' ', '_')}_derivatives.csv")


def main():
    if "--tables" in sys.argv:
        tables()
        return
    geo = json.loads(GEOMETRY.read_text())
    plane, foils, flaps = build_plane(geo)
    project = {"format": 1, "name": "Spitfire Mk I 6-DoF", "foils": foils, "foil_polars": [], "planes": [plane],
               "plane_polars": []}
    (HERE / GEOMETRY_FILE).write_text(json.dumps(project, indent=1))
    (HERE / JOB_FILE).write_text(json.dumps(job(flaps, inertia(geo)), separators=(",", ":")))
    print(f"{GEOMETRY_FILE}: {len(foils)} foils; flaps {json.dumps(flaps)}")
    print(f"{JOB_FILE}: {len(job(flaps, inertia(geo))['steps'])} polars")


if __name__ == "__main__":
    main()
