"""Write the minimal xflr5 (.xfl) and flow5 (.fl5) project files of this folder.

The files are written field by field in the order of the C++ load code
(xflr5 v6.61 `MainFrame::serializeProjectXFL`, `WPolar::serializeWPlrXFL`;
flow5 `FileIO::serializeProjectFl5`, `serial::serializePlanePolarFl5v750`),
as Qt QDataStreams: big-endian, doubles for C++ floats, QString = u32 byte
count + UTF-16BE.  Each file holds one small plane, one foil, plane polar
definitions (one for a plane that is not in the file) and the objects the
conversion report must list as approximated or skipped.

    python make_minimal.py      (rewrites minimal_v6.xfl and minimal_v7.fl5)
"""
import math
import os
import struct


class W:
    def __init__(self):
        self.b = bytearray()

    def i(self, *v):
        for x in v:
            self.b += struct.pack(">i", x)
        return self

    def bo(self, *v):
        for x in v:
            self.b.append(1 if x else 0)
        return self

    def f(self, *v):
        for x in v:
            self.b += struct.pack(">d", x)
        return self

    def fs(self, n, v=0.0):
        return self.f(*([v] * n))

    def is_(self, n, v=0):
        return self.i(*([v] * n))

    def bos(self, n, v=False):
        return self.bo(*([v] * n))

    def s(self, text):
        u = text.encode("utf-16-be")
        self.b += struct.pack(">I", len(u)) + u
        return self

    def color(self):
        self.b.append(1)
        for c in (255, 40, 80, 160, 0):
            self.b += struct.pack(">H", c)
        return self


def foil_points():
    """21-point NACA 0012-like foil, Selig order."""
    xs, ys = [], []
    for k in range(21):
        t = math.pi * k / 20.0
        x = 0.5 * (1.0 + math.cos(t))
        yt = 0.6 * (0.2969 * math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1036 * x ** 4) * 0.2
        xs.append(x)
        ys.append(yt if k <= 10 else -yt)
    return xs, ys


# ---------------------------------------------------------------- xflr5 v6


def x_style(w):  # LineStyle::serializeXfl: stipple, width, symbol, colour, visible
    w.i(0, 1, 0).color().bo(True)


def x_style_fl5(w):  # LineStyle::serializeFl5
    w.i(500001, 0, 1, 0).color().bo(True).s("")


def x_wpolar(w, plane, name, method, ptype, flags, height, rho, nu, refdim, ref, auto_inertia, mass, cog, inertia, gains, wake, qinf, alpha, beta, npts, extra):
    w.i(200014).s(plane).s(name).f(*ref)  # area, chord, span
    x_style(w)
    w.i(method, ptype)
    w.bo(*flags)  # VLM1, thin surfaces, tilted geometry, Dirichlet, viscous, ignore body, ground
    w.f(height, rho, nu).i(refdim).bo(auto_inertia).f(mass, *cog).f(*inertia)
    w.i(len(gains)).f(*gains)
    w.i(wake[0]).f(wake[1], wake[2])
    w.f(qinf, alpha, beta)
    w.i(npts).fs(36 * npts, 0.1)
    w.is_(19).i(0).fs(35)
    w.f(extra[0], 0.0, 0.0, 0.0).f(extra[1], 0.0, 0.0, 0.0)  # extra drag areas, coefficients
    w.fs(7)  # inertia gains


def x_wing(w, name, secs, mass, kind):
    w.i(100001).s(name).s("").color().bo(True).i(len(secs))
    for (y, chord, offset, dihedral, twist, foil) in secs:
        w.s(foil).s(foil).f(chord, y, offset, dihedral, twist).i(7, 9, 1, 0)
    w.f(mass).i(0)  # no point masses
    w.is_(19).i(kind).fs(50)


def write_xfl(path):
    w = W()
    w.i(200002).is_(6)  # project format, unit indices
    # the polar dialog's defaults (not a polar of the project)
    x_wpolar(w, "", "", 2, 1, [False] * 7, 0.0, 1.225, 1.5e-5, 1, (0, 0, 0), True, 0, (0, 0, 0), (0, 0, 0, 0), [], (1, 100, 1.1), 10, 0, 0, 0, (0, 0))
    # one plane: main wing and elevator, no body
    w.i(1)
    w.i(100002).s("Mini").s("minimal xflr5 plane")
    x_style_fl5(w)
    x_wing(w, "Main", [(0.0, 0.25, 0.0, 3.0, 0.0, "M1"), (1.0, 0.15, 0.05, 3.0, -2.0, "M1")], 0.4, 0)
    x_wing(w, "", [(0.0, 0.2, 0.0, 0.0, 0.0, "M1"), (0.5, 0.1, 0.0, 0.0, 0.0, "M1")], 0.0, 1)
    x_wing(w, "Stab", [(0.0, 0.12, 0.0, 0.0, 0.0, "M1"), (0.3, 0.08, 0.02, 0.0, 0.0, "M1")], 0.05, 2)
    x_wing(w, "Fin", [(0.0, 0.15, 0.0, 0.0, 0.0, "M1"), (0.2, 0.1, 0.05, 0.0, 0.0, "M1")], 0.03, 3)
    w.bo(False, True, False, False, False, False)  # biplane, stab, fin, double fin, sym fin, -
    for (p, tilt) in [((0.0, 0.0, 0.0), 1.0), ((0.0, 0.0, 0.0), 0.0), ((0.8, 0.0, 0.05), -1.5), ((0.9, 0.0, 0.0), 0.0)]:
        w.f(*p).f(tilt)
    w.bo(False).f(0.0, 0.0)  # no body (position still stored)
    w.i(1).f(0.3, -0.1, 0.0, 0.0).s("battery")
    w.is_(20).fs(50)
    # plane polars
    w.i(2)
    x_wpolar(w, "Mini", "T2 trimmed", 2, 2, [False, True, True, True, True, False, True], 0.5, 1.2, 1.6e-5, 3, (0.42, 0.2, 2.1), False, 1.1, (0.08, 0.0, 0.01), (0.05, 0.06, 0.1, 0.001), [0.5, 0.0], (5, 50.0, 1.1), 12.0, 3.0, 0.0, 1, (0.01, 0.5))
    x_wpolar(w, "Ghost", "orphan", 1, 1, [False] * 7, 0.0, 1.225, 1.5e-5, 2, (0, 0, 0), True, 0, (0, 0, 0), (0, 0, 0, 0), [], (1, 100, 1.1), 10, 0, 0, 0, (0, 0))
    w.i(0)  # plane operating points
    # one foil with a deflected trailing-edge flap
    xs, ys = foil_points()
    w.i(1)
    w.i(100007).s("M1").s("minimal foil")
    x_style(w)
    w.bo(False, False, True).f(0.0, 20.0, 50.0).f(5.0, 75.0, 50.0).i(len(xs))
    for x, y in zip(xs, ys):
        w.f(x, y)
    # one 2D polar
    w.i(1)
    w.i(100005).s("M1").s("T1_Re0.200_M0.00_N9.0")
    x_style(w)
    w.i(1, 1, 1).f(200000.0, 0.0, 0.0, 1.0, 1.0, 9.0)
    pts = [(-2.0, 0.009, 0.004, -0.2), (0.0, 0.008, 0.003, 0.0), (2.0, 0.0085, 0.0035, 0.2)]
    w.i(len(pts))
    for (a, cd, cdp, cl) in pts:
        w.f(a, cd, cdp, cl, 0.0, 0.6, 0.9, 0.0, -0.8, 200000.0, 0.25)
    w.is_(19).fs(50)
    w.i(0)  # foil operating points (the reader stops before them)
    open(path, "wb").write(bytes(w.b))


# ---------------------------------------------------------------- flow5 v7.54


def f_style(w):  # serial::serializeLineStyleFl5
    w.i(500756, 0, 1, 0).color().bo(True).s("")


def f_spline(w):
    w.i(500001)
    f_style(w)
    w.bo(False).i(79).i(2).f(0.0, 0.0, 1.0, 0.0).i(2).f(1.0, 1.0).bo(False, False).f(0.0).i(0, 0).i(0, 0)


def f_bspline(w):
    f_spline(w)
    w.i(3, 0, 0)


def f_polar2d(w, foil, name):
    w.i(500750).s(foil).s(name)
    f_style(w)
    w.i(1, 0, 1, 1).f(300000.0, 0.0, 0.0, 1.0, 1.0, 9.0)
    w.i(0, 12, 0)  # no points
    w.i(0, 1).f(0.0)


def f_planepolar(w, plane, name, method, ptype, ignore_body, extra, ground, rho, nu, beta, phi, auto_inertia, mass, cog, inertia, viscous, on_the_fly, wake, visc_from_cl, qinf, alpha, thin, refdim, ref, flaps, ranges, npts):
    # Polar3d v7.50
    w.i(500750).s(name)
    f_style(w)
    w.i(method).bo(True, True, ignore_body).i(ptype)  # Trefftz, Dirichlet
    w.i(len(extra))
    for (tag, a, c) in extra:
        w.s(tag).f(a, c)
    w.bo(ground, False).f(ground and 0.4 or 0.0, 0.0, rho, nu, beta, phi)
    w.bo(auto_inertia).f(mass, *cog).f(inertia[0], inertia[1], inertia[2], 0.0, inertia[3], 0.0)
    w.bo(viscous, on_the_fly).f(9.0, 1.0, 1.0)
    w.bo(False).i(wake[0]).f(wake[1], wake[2]).bo(visc_from_cl).fs(4).i(10)
    w.bos(10).is_(20).fs(20)
    # PlanePolar v7.50
    w.s(plane).f(qinf, alpha).bo(thin, False).i(refdim).f(*ref)
    w.i(len(flaps))
    for vals in flaps:
        w.i(500750).s("flaps").i(len(vals)).f(*vals)
    w.i(0)  # AVL-type controls
    w.bo(False)  # adjusted velocity
    w.i(1).i(len(ranges))
    for (n, lo, hi) in ranges:
        w.s(n).f(lo, hi)
    w.i(0).i(0)  # operating range, inertia range
    w.bo(False, False).i(1)  # fuselage inertia, fuselage drag, drag method
    w.i(0).i(npts)
    for _ in range(npts):
        w.i(500750).fs(27, 0.1).fs(7 + 16, 0.2)
    w.bo(False)
    f_bspline(w)
    w.bo(True, True).f(0.0).bos(10).is_(20).fs(20)


def f_inertia(w, mass, pms):
    w.i(500001).f(mass).fs(9)
    w.i(len(pms))
    for (tag, m, p) in pms:
        w.f(m, *p).s(tag)
    w.i(0, 0)


def f_part(w, name, mass, pms):
    w.i(500757).s(name).s("")
    f_style(w)
    f_inertia(w, mass, pms)
    w.bo(True, False).f(0.01, 0.1).i(23).f(0.01, 0.1).i(23).i(0).i(0)


def f_wing(w, name, secs, two_sided, kind, pos, rot):
    f_part(w, name, 0.4, [])
    w.i(500001).bo(True).i(len(secs))
    for (y, chord, offset, foil) in secs:
        w.i(500001).s(foil).s(foil).f(chord, y, offset, 2.0, 0.0).i(13, 11, 1, 0).i(0)
    w.bo(two_sided, not two_sided, False).i(kind).f(*pos).f(*rot).i(0, 0)


def write_fl5(path):
    w = W()
    w.i(500754)
    # project meta data: the dialogs' default 2D and plane polars
    f_polar2d(w, "", "")
    f_planepolar(w, "", "", 3, 1, False, [], False, 1.225, 1.5e-5, 0, 0, True, 0, (0, 0, 0), (0, 0, 0, 0), True, False, (1, 100.0, 1.1), True, 10.0, 0.0, False, 2, (0, 0, 0), [], [], 0)
    w.s("")
    w.i(500001)  # spline foil
    f_style(w)
    f_bspline(w)
    f_bspline(w)
    w.bos(5).i(0, 0)
    f_bspline(w)
    f_bspline(w)
    f_bspline(w)  # B-spline and cubic spline
    f_spline(w)
    w.i(0, 0)
    w.i(0, 0)  # meta data spares
    # 2D objects
    w.i(500750).i(1)
    xs, ys = foil_points()
    w.i(500753).s("N1").s("minimal fl5 foil").i(0, 10, 0).color().bo(True)
    w.bo(False, False, False).f(0.0, 0.2, 0.5, 0.0, 0.7, 0.5).f(0.0).i(3).i(len(xs))
    for x, y in zip(xs, ys):
        w.f(x, y)
    w.i(0)  # 2D polars
    w.i(0)  # foil operating points
    # one xfl-type plane: a two-sided wing and a fin, no fuselage
    w.i(1).i(0)
    w.i(500003).i(2).s("Flow mini").s("minimal flow5 plane").i(0, 1, 0).color()
    f_wing(w, "Wing", [(0.0, 0.3, 0.0, "N1"), (1.2, 0.2, 0.05, "N1")], True, 0, (0.0, 0.0, 0.0), (0.0, 1.5, 0.0))
    f_wing(w, "Fin", [(0.0, 0.2, 0.0, "N1"), (0.3, 0.1, 0.05, "N1")], False, 3, (1.0, 0.0, 0.0), (-90.0, 0.0, 0.0))
    w.i(0)  # fuselages
    w.bo(True)
    f_inertia(w, 0.0, [("nose", 0.2, (-0.2, 0.0, 0.0))])
    w.i(0, 0)
    # plane polars
    w.i(2)
    f_planepolar(w, "Flow mini", "T6 flap sweep", 6, 6, True, [("wheels", 0.01, 0.5)], True, 1.2, 1.6e-5, 0.0, 5.0, False, 1.3, (0.1, 0.0, 0.0), (0.04, 0.05, 0.09, 0.002), True, True, (6, 40.0, 1.2), False, 15.0, 2.0, True, 1, (0.5, 0.25, 2.4), [[0.0, 5.0], [0.0]], [("flap 1", 0.0, 10.0)], 2)
    f_planepolar(w, "Nope", "orphan", 3, 1, False, [], False, 1.225, 1.5e-5, 0, 0, True, 0, (0, 0, 0), (0, 0, 0, 0), True, False, (1, 100.0, 1.1), True, 10.0, 0.0, False, 2, (0, 0, 0), [], [], 0)
    w.i(0)  # plane operating points (not read)
    open(path, "wb").write(bytes(w.b))


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    write_xfl(os.path.join(here, "minimal_v6.xfl"))
    write_fl5(os.path.join(here, "minimal_v7.fl5"))
