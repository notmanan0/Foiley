"""Writes the hand-made STEP (AP214 B-rep) test fixtures in this folder.

Each part has a closed-form area and volume (see crates/plane/tests/f022.rs).
Run from anywhere:  python tests/data/cad/make_fixtures.py
The output is deterministic; the generated .stp files are committed.
"""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
TAU = 2.0 * math.pi
X, Y, Z = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)


def fmt(x):
    x = float(x)
    if x == 0.0:
        return "0."
    s = "%.17g" % x
    if "e" in s:
        m, e = s.split("e")
        if "." not in m:
            m += "."
        return m + "E" + e
    if "." not in s:
        s += "."
    return s


def norm(v):
    n = math.sqrt(sum(c * c for c in v))
    return tuple(c / n for c in v)


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


class Step:
    def __init__(self, scale=1.0):
        self.s = scale  # file units per metre
        self.ents = []
        self.vtx = {}
        self.edges = {}
        self.items = []
        self.ctx2 = None

    def add(self, txt):
        self.ents.append(txt)
        return len(self.ents)

    # ---- geometry
    def pt(self, p):
        return self.add("CARTESIAN_POINT('',(%s))" % ",".join(fmt(c * self.s) for c in p))

    def pt2(self, p):  # parameter-space point (not scaled)
        return self.add("CARTESIAN_POINT('',(%s))" % ",".join(fmt(c) for c in p))

    def dirn(self, d):
        return self.add("DIRECTION('',(%s))" % ",".join(fmt(c) for c in norm(d)))

    def ax2(self, o, z, x):
        return self.add("AXIS2_PLACEMENT_3D('',#%d,#%d,#%d)" % (self.pt(o), self.dirn(z), self.dirn(x)))

    def ax1(self, o, z):
        return self.add("AXIS1_PLACEMENT('',#%d,#%d)" % (self.pt(o), self.dirn(z)))

    def vec(self, d):
        mag = math.sqrt(sum(c * c for c in d))
        return self.add("VECTOR('',#%d,%s)" % (self.dirn(d), fmt(mag * self.s)))

    def line(self, p, d):
        return self.add("LINE('',#%d,#%d)" % (self.pt(p), self.vec(d)))

    def line2(self, p, d):
        dd = self.add("DIRECTION('',(%s))" % ",".join(fmt(c) for c in norm(d)))
        mag = math.sqrt(sum(c * c for c in d))
        v = self.add("VECTOR('',#%d,%s)" % (dd, fmt(mag)))
        return self.add("LINE('',#%d,#%d)" % (self.pt2(p), v))

    def circle(self, o, z, x, r):
        return self.add("CIRCLE('',#%d,%s)" % (self.ax2(o, z, x), fmt(r * self.s)))

    def ellipse(self, o, z, x, a, b):
        return self.add("ELLIPSE('',#%d,%s,%s)" % (self.ax2(o, z, x), fmt(a * self.s), fmt(b * self.s)))

    def bspline_curve(self, deg, ctrl, mults, knots, weights=None):
        cp = ",".join("#%d" % self.pt(p) for p in ctrl)
        m = ",".join(str(k) for k in mults)
        k = ",".join(fmt(x) for x in knots)
        if weights is None:
            return self.add("B_SPLINE_CURVE_WITH_KNOTS('',%d,(%s),.UNSPECIFIED.,.F.,.F.,(%s),(%s),.UNSPECIFIED.)" % (deg, cp, m, k))
        w = ",".join(fmt(x) for x in weights)
        return self.add(
            "(BOUNDED_CURVE() B_SPLINE_CURVE(%d,(%s),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS((%s),(%s),.UNSPECIFIED.) "
            "CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((%s)) REPRESENTATION_ITEM(''))" % (deg, cp, m, k, w)
        )

    def bspline_surface(self, du, dv, net, mu, mv, ku, kv, weights=None):
        cp = ",".join("(" + ",".join("#%d" % self.pt(p) for p in row) + ")" for row in net)
        fm = lambda l: ",".join(str(x) for x in l)
        fk = lambda l: ",".join(fmt(x) for x in l)
        if weights is None:
            return self.add(
                "B_SPLINE_SURFACE_WITH_KNOTS('',%d,%d,(%s),.UNSPECIFIED.,.F.,.F.,.F.,(%s),(%s),(%s),(%s),.UNSPECIFIED.)" % (du, dv, cp, fm(mu), fm(mv), fk(ku), fk(kv))
            )
        w = ",".join("(" + ",".join(fmt(x) for x in row) + ")" for row in weights)
        return self.add(
            "(BOUNDED_SURFACE() B_SPLINE_SURFACE(%d,%d,(%s),.UNSPECIFIED.,.F.,.F.,.F.) B_SPLINE_SURFACE_WITH_KNOTS((%s),(%s),(%s),(%s),.UNSPECIFIED.) "
            "GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE((%s)) REPRESENTATION_ITEM('') SURFACE())" % (du, dv, cp, fm(mu), fm(mv), fk(ku), fk(kv), w)
        )

    # ---- surfaces
    def plane(self, o, z, x):
        return self.add("PLANE('',#%d)" % self.ax2(o, z, x))

    def cyl(self, o, z, x, r):
        return self.add("CYLINDRICAL_SURFACE('',#%d,%s)" % (self.ax2(o, z, x), fmt(r * self.s)))

    def cone(self, o, z, x, r, ang):
        return self.add("CONICAL_SURFACE('',#%d,%s,%s)" % (self.ax2(o, z, x), fmt(r * self.s), fmt(ang)))

    def sphere(self, o, z, x, r):
        return self.add("SPHERICAL_SURFACE('',#%d,%s)" % (self.ax2(o, z, x), fmt(r * self.s)))

    def torus(self, o, z, x, R, r):
        return self.add("TOROIDAL_SURFACE('',#%d,%s,%s)" % (self.ax2(o, z, x), fmt(R * self.s), fmt(r * self.s)))

    def revolution(self, curve, o, z):
        return self.add("SURFACE_OF_REVOLUTION('',#%d,#%d)" % (curve, self.ax1(o, z)))

    def extrusion(self, curve, d):
        return self.add("SURFACE_OF_LINEAR_EXTRUSION('',#%d,#%d)" % (curve, self.vec(d)))

    def pcurve(self, surf, curve2d):
        if self.ctx2 is None:
            self.ctx2 = self.add("( GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('2D SPACE','') )")
        d = self.add("DEFINITIONAL_REPRESENTATION('',(#%d),#%d)" % (curve2d, self.ctx2))
        return self.add("PCURVE('',#%d,#%d)" % (surf, d))

    def surface_curve(self, c3, pcurves, seam=False):
        kind = "SEAM_CURVE" if seam else "SURFACE_CURVE"
        return self.add("%s('',#%d,(%s),.CURVE_3D.)" % (kind, c3, ",".join("#%d" % p for p in pcurves)))

    # ---- topology
    def vertex(self, p):
        key = tuple(round(c, 12) for c in p)
        if key not in self.vtx:
            self.vtx[key] = self.add("VERTEX_POINT('',#%d)" % self.pt(p))
        return self.vtx[key]

    def edge(self, p0, p1, curve, same=True):
        return self.add("EDGE_CURVE('',#%d,#%d,#%d,%s)" % (self.vertex(p0), self.vertex(p1), curve, ".T." if same else ".F."))

    def oe(self, e, fwd=True):
        return self.add("ORIENTED_EDGE('',*,*,#%d,%s)" % (e, ".T." if fwd else ".F."))

    def loop(self, oes):
        return self.add("EDGE_LOOP('',(%s))" % ",".join("#%d" % o for o in oes))

    def face(self, loops, surf, same=True):
        """loops: lists of (edge, forward); the first one is the outer bound"""
        bs = []
        for i, l in enumerate(loops):
            lp = self.loop([self.oe(e, f) for e, f in l])
            bs.append(self.add("%s('',#%d,.T.)" % ("FACE_OUTER_BOUND" if i == 0 else "FACE_BOUND", lp)))
        return self.add("ADVANCED_FACE('',(%s),#%d,%s)" % (",".join("#%d" % b for b in bs), surf, ".T." if same else ".F."))

    def seg(self, p0, p1):
        """straight edge between two points, shared between faces"""
        k0, k1 = self.vertex(p0), self.vertex(p1)
        if (k1, k0) in self.edges:
            return self.edges[(k1, k0)], False
        if (k0, k1) not in self.edges:
            self.edges[(k0, k1)] = self.edge(p0, p1, self.line(p0, sub(p1, p0)))
        return self.edges[(k0, k1)], True

    def poly_face(self, pts, holes=()):
        """planar face from a polygon counter-clockwise seen from outside"""
        n = (0.0, 0.0, 0.0)
        for i in range(len(pts)):
            a, b = pts[i], pts[(i + 1) % len(pts)]
            n = (n[0] + (a[1] - b[1]) * (a[2] + b[2]), n[1] + (a[2] - b[2]) * (a[0] + b[0]), n[2] + (a[0] - b[0]) * (a[1] + b[1]))
        outer = [self.seg(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
        surf = self.plane(pts[0], n, sub(pts[1], pts[0]))
        return self.face([outer] + list(holes), surf)

    def solid(self, faces, closed=True):
        sh = self.add("%s('',(%s))" % ("CLOSED_SHELL" if closed else "OPEN_SHELL", ",".join("#%d" % f for f in faces)))
        if closed:
            self.items.append(self.add("MANIFOLD_SOLID_BREP('',#%d)" % sh))
        else:
            self.items.append(self.add("SHELL_BASED_SURFACE_MODEL('',(#%d))" % sh))

    def write(self, name, prefix="$"):
        u = ".%s." % prefix if prefix != "$" else "$"
        lines = ["#%d=%s;" % (i + 1, e) for i, e in enumerate(self.ents)]
        k = len(self.ents)
        lu, au, su, un, ctx, rep = k + 1, k + 2, k + 3, k + 4, k + 5, k + 6
        lines.append("#%d=( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(%s,.METRE.) );" % (lu, u))
        lines.append("#%d=( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );" % au)
        lines.append("#%d=( NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT() );" % su)
        lines.append("#%d=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#%d,'distance_accuracy_value','');" % (un, lu))
        lines.append(
            "#%d=( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#%d)) GLOBAL_UNIT_ASSIGNED_CONTEXT((#%d,#%d,#%d)) REPRESENTATION_CONTEXT('','') );"
            % (ctx, un, lu, au, su)
        )
        lines.append("#%d=ADVANCED_BREP_SHAPE_REPRESENTATION('%s',(%s),#%d);" % (rep, name, ",".join("#%d" % i for i in self.items), ctx))
        txt = (
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('Foiley test fixture %s'),'2;1');\n"
            "FILE_NAME('%s.stp','2026-11-04T00:00:00',('Foiley'),(''),'make_fixtures.py','','');\n"
            "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n%s\nENDSEC;\nEND-ISO-10303-21;\n"
        ) % (name, name, "\n".join(lines))
        with open(os.path.join(HERE, name + ".stp"), "w", newline="\n") as f:
            f.write(txt)


def box_hole(st, o=(0.0, 0.0, 0.0), L=2.0, W=1.0, H=0.5, r=0.2):
    """box [0,L]x[0,W]x[0,H] (offset by o) with a round hole through it along z"""
    p = lambda x, y, z: (o[0] + x, o[1] + y, o[2] + z)
    A, B, C, D = p(0, 0, 0), p(L, 0, 0), p(L, W, 0), p(0, W, 0)
    A1, B1, C1, D1 = p(0, 0, H), p(L, 0, H), p(L, W, H), p(0, W, H)
    c0, c1 = p(L / 2, W / 2, 0), p(L / 2, W / 2, H)
    h0, h1 = p(L / 2 + r, W / 2, 0), p(L / 2 + r, W / 2, H)
    eb = st.edge(h0, h0, st.circle(c0, Z, X, r))
    et = st.edge(h1, h1, st.circle(c1, Z, X, r))
    es = st.edge(h0, h1, st.line(h0, (0, 0, H)))
    faces = [
        st.poly_face([A, D, C, B], [[(eb, True)]]),
        st.poly_face([A1, B1, C1, D1], [[(et, False)]]),
        st.poly_face([A, B, B1, A1]),
        st.poly_face([B, C, C1, B1]),
        st.poly_face([C, D, D1, C1]),
        st.poly_face([D, A, A1, D1]),
        st.face([[(es, True), (et, True), (es, False), (eb, False)]], st.cyl(c0, Z, X, r), same=False),
    ]
    st.solid(faces)


def cylinder(st, R=0.5, L=2.0, pcurves=False):
    """closed cylinder along +x with a tilted seam, planar caps"""
    z, x = X, (0.0, 0.6, 0.8)
    v0 = (0.0, R * 0.6, R * 0.8)
    v1 = (L, R * 0.6, R * 0.8)
    surf = st.cyl((0, 0, 0), z, x, R)
    cap0 = st.plane((0, 0, 0), z, x)
    cap1 = st.plane((L, 0, 0), z, x)
    c0, c1 = st.circle((0, 0, 0), z, x, R), st.circle((L, 0, 0), z, x, R)
    sl = st.line(v0, (L, 0, 0))
    if pcurves:
        c0 = st.surface_curve(c0, [st.pcurve(surf, st.line2((0, 0), (1, 0)))])
        c1 = st.surface_curve(c1, [st.pcurve(surf, st.line2((0, L), (1, 0)))])
        sl = st.surface_curve(sl, [st.pcurve(surf, st.line2((0, 0), (0, L))), st.pcurve(surf, st.line2((TAU, 0), (0, L)))], seam=True)
    e0, e1, es = st.edge(v0, v0, c0), st.edge(v1, v1, c1), st.edge(v0, v1, sl)
    st.solid([
        st.face([[(e0, True), (es, True), (e1, False), (es, False)]], surf),
        st.face([[(e0, False)]], cap0, same=False),
        st.face([[(e1, True)]], cap1),
    ])


def frustum(st, r0=0.3, r1=0.5, H=1.0):
    """surface of revolution of a straight generatrix about +z, planar caps"""
    v0, v1 = (r0, 0.0, 0.0), (r1, 0.0, H)
    gen = st.line(v0, (r1 - r0, 0.0, H))
    eb = st.edge(v0, v0, st.circle((0, 0, 0), Z, X, r0))
    et = st.edge(v1, v1, st.circle((0, 0, H), Z, X, r1))
    es = st.edge(v0, v1, gen)
    st.solid([
        st.face([[(eb, True), (es, True), (et, False), (es, False)]], st.revolution(gen, (0, 0, 0), Z)),
        st.face([[(eb, False)]], st.plane((0, 0, 0), Z, X), same=False),
        st.face([[(et, True)]], st.plane((0, 0, H), Z, X)),
    ])


def cyl_partial(st):
    """open cylindrical patch spanning -60..+60 degrees (crosses the u = 0 seam)"""
    a = math.pi / 3
    P = lambda t, z: (math.cos(t), math.sin(t), z)
    va, vb, vc, vd = P(-a, 0), P(a, 0), P(a, 1), P(-a, 1)
    ab = st.edge(va, vb, st.circle((0, 0, 0), Z, X, 1.0))
    at = st.edge(vd, vc, st.circle((0, 0, 1), Z, X, 1.0))
    lr = st.edge(vb, vc, st.line(vb, Z))
    ll = st.edge(va, vd, st.line(va, Z))
    st.solid([st.face([[(ab, True), (lr, True), (at, False), (ll, False)]], st.cyl((0, 0, 0), Z, X, 1.0))], closed=False)


def cyl_annulus(st):
    """cylindrical band bounded by two full circles and no seam edge"""
    v0, v1 = (1.0, 0.0, 0.0), (1.0, 0.0, 1.0)
    eb = st.edge(v0, v0, st.circle((0, 0, 0), Z, X, 1.0))
    et = st.edge(v1, v1, st.circle((0, 0, 1), Z, X, 1.0))
    st.solid([st.face([[(eb, True)], [(et, False)]], st.cyl((0, 0, 0), Z, X, 1.0))], closed=False)


def cone_apex(st):
    """cone from its apex (a singular point of the surface) to a planar base"""
    ang = math.atan(0.3)
    apex, vb = (0.0, 0.0, -2.0 / 3.0), (0.5, 0.0, 1.0)
    es = st.edge(apex, vb, st.line(apex, sub(vb, apex)))
    eb = st.edge(vb, vb, st.circle((0, 0, 1), Z, X, 0.5))
    st.solid([
        st.face([[(es, True), (eb, False), (es, False)]], st.cone((0, 0, 0), Z, X, 0.2, ang)),
        st.face([[(eb, True)]], st.plane((0, 0, 1), Z, X)),
    ])


def sphere(st, R=0.5, c=(0.0, 0.0, 0.0)):
    """full sphere: one face bounded by a meridian seam between the poles"""
    s, n = (c[0], c[1], c[2] - R), (c[0], c[1], c[2] + R)
    es = st.edge(s, n, st.circle(c, (0, -1, 0), X, R))
    st.solid([st.face([[(es, True), (es, False)]], st.sphere(c, Z, X, R))])


def torus(st, R=0.5, r=0.15):
    """full torus: one face, a u seam and a v seam through one vertex"""
    v = (R + r, 0.0, 0.0)
    ee = st.edge(v, v, st.circle((0, 0, 0), Z, X, R + r))
    em = st.edge(v, v, st.circle((R, 0, 0), (0, -1, 0), X, r))
    st.solid([st.face([[(ee, True), (em, True), (ee, False), (em, False)]], st.torus((0, 0, 0), Z, X, R, r))])


def ellipse_prism(st, a=0.4, b=0.2, H=1.0):
    """elliptic prism: linear extrusion of an ellipse, planar caps"""
    v0, v1 = (a, 0.0, 0.0), (a, 0.0, H)
    el0 = st.ellipse((0, 0, 0), Z, X, a, b)
    e0 = st.edge(v0, v0, el0)
    e1 = st.edge(v1, v1, st.ellipse((0, 0, H), Z, X, a, b))
    es = st.edge(v0, v1, st.line(v0, (0, 0, H)))
    st.solid([
        st.face([[(e0, True), (es, True), (e1, False), (es, False)]], st.extrusion(el0, (0, 0, H))),
        st.face([[(e0, False)]], st.plane((0, 0, 0), Z, X), same=False),
        st.face([[(e1, True)]], st.plane((0, 0, H), Z, X)),
    ])


# ---------------------------------------------------------------- B-splines (F022b)

W45 = math.sqrt(0.5)
CIRCLE_W = [1.0, W45, 1.0, W45, 1.0, W45, 1.0, W45, 1.0]
CIRCLE_MULTS = [3, 2, 2, 2, 3]
CIRCLE_KNOTS = [0.0, 0.25, 0.5, 0.75, 1.0]


def circle_net(c, r, axis="x"):
    """control points of a rational quadratic full circle about `axis` through c"""
    sq = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (1, 0)]
    if axis == "x":
        return [(c[0], c[1] + r * a, c[2] + r * b) for a, b in sq]
    return [(c[0] + r * a, c[1] + r * b, c[2]) for a, b in sq]


def bezier_u2v_net():
    """bicubic Bezier patch with S(u, v) = (u, v, u^2 v) exactly"""
    cu = [0.0, 0.0, 1.0 / 3.0, 1.0]  # Bernstein coefficients of u^2
    return [[(i / 3.0, j / 3.0, cu[i] * j / 3.0) for j in range(4)] for i in range(4)]


def bspline_patch(st):
    """two faces of the patch S = (u, v, u^2 v) split along the diagonal
    curve S(t, t) = (t, t, t^3); no PCURVEs, so the diagonal must be inverted"""
    surf = st.bspline_surface(3, 3, bezier_u2v_net(), [4, 4], [4, 4], [0.0, 1.0], [0.0, 1.0])
    p00, p10, p11, p01 = (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0.0, 1.0, 0.0)
    diag = st.edge(p00, p11, st.bspline_curve(3, [(0, 0, 0), (1 / 3.0, 1 / 3.0, 0), (2 / 3.0, 2 / 3.0, 0), (1, 1, 1)], [4, 4], [0.0, 1.0]))
    bottom = st.edge(p00, p10, st.line(p00, (1, 0, 0)))  # S(u, 0)
    right = st.edge(p10, p11, st.line(p10, (0, 1, 1)))  # S(1, v)
    top = st.edge(p01, p11, st.bspline_curve(2, [(0, 1, 0), (0.5, 1, 0), (1, 1, 1)], [3, 3], [0.0, 1.0]))  # S(u, 1)
    left = st.edge(p00, p01, st.line(p00, (0, 1, 0)))  # S(0, v)
    st.solid([
        st.face([[(bottom, True), (right, True), (diag, False)]], surf),
        st.face([[(diag, True), (top, False), (left, False)]], surf),
    ], closed=False)


def bspline_plate_hole(st, r=0.2):
    """bilinear B-spline plate [0,1]^2 (S = (u, v, 0)) with a hole bounded by a
    rational B-spline circle"""
    surf = st.bspline_surface(1, 1, [[(0, 0, 0), (0, 1, 0)], [(1, 0, 0), (1, 1, 0)]], [2, 2], [2, 2], [0.0, 1.0], [0.0, 1.0])
    A, B, C, D = (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)
    outer = [st.seg(A, B), st.seg(B, C), st.seg(C, D), st.seg(D, A)]
    h = (0.5 + r, 0.5, 0.0)
    hole = st.edge(h, h, st.bspline_curve(2, circle_net((0.5, 0.5, 0.0), r, axis="z"), CIRCLE_MULTS, CIRCLE_KNOTS, CIRCLE_W))
    st.solid([st.face([outer, [(hole, False)]], surf)], closed=False)


LOFT = [(0.0, 0.0), (0.3, 0.1), (1.0, 0.1), (1.4, 0.0)]  # (x, radius) stations


def loft_fuselage(st):
    """closed lofted B-spline fuselage: degree 1 along x through the stations
    (nose and tail collapse to points), rational circle around; exact cones
    and a cylinder.  One face bounded by its seam, like a sphere."""
    net = [circle_net((x, 0.0, 0.0), r) for x, r in LOFT]
    w = [CIRCLE_W[:] for _ in LOFT]
    surf = st.bspline_surface(1, 2, net, [2, 1, 1, 2], CIRCLE_MULTS, [0.0, 1 / 3.0, 2 / 3.0, 1.0], CIRCLE_KNOTS, w)
    nose, tail = (0.0, 0.0, 0.0), (1.4, 0.0, 0.0)
    seam = st.edge(nose, tail, st.bspline_curve(1, [(x, r, 0.0) for x, r in LOFT], [2, 1, 1, 2], [0.0, 1 / 3.0, 2 / 3.0, 1.0]))
    # dS/du x dS/dv points inward (u along +x, v counter-clockwise about +x)
    st.solid([st.face([[(seam, True), (seam, False)]], surf, same=False)])


def one(name, build, scale=1.0, prefix="$", **kw):
    st = Step(scale)
    build(st, **kw)
    st.write(name, prefix=prefix)


if __name__ == "__main__":
    one("box_hole", box_hole, scale=1000.0, prefix="MILLI")
    one("cyl_seam", cylinder)
    one("cyl_pcurve", cylinder, pcurves=True)
    one("frustum_rev", frustum)
    one("cyl_partial", cyl_partial)
    one("cyl_annulus", cyl_annulus)
    one("cone_apex", cone_apex)
    one("sphere", sphere)
    one("torus", torus)
    one("ellipse_prism", ellipse_prism)
    one("bspline_patch", bspline_patch)
    one("bspline_plate_hole", bspline_plate_hole)
    one("loft_fuselage", loft_fuselage)
