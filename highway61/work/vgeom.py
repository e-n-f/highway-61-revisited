"""Vector back end for the 1961 constructions.

The construction modules (dim1961, dimE1..dimE6) build letters from half-planes, disks, sectors and polygons
combined with & | ~ on the raster Grid.  VGrid offers the same interface but every mask is an exact shapely
region, so the same construction code yields outlines instead of bitmaps.  Circles are tessellated very finely;
outline.py later recovers the circular arcs from the vertices and writes them as Bezier arcs.
"""
import numpy as np
from shapely.geometry import Polygon, Point, box, MultiPolygon, GeometryCollection
from shapely import affinity
from shapely.ops import unary_union

QS = 2048                                   # quarter-circle segments: 0.044 deg steps, sagitta < 1e-4 unit at r = 700
BIG = 1e5
BOX = box(-3000, -3000, 4000, 4000)
GRIDBOX = box(-60, -80, 780, 760)


def _clean(g):
    if isinstance(g, (Polygon, MultiPolygon)):
        return g
    if isinstance(g, GeometryCollection):
        polys = [p for p in g.geoms if isinstance(p, (Polygon, MultiPolygon))]
        return unary_union(polys) if polys else Polygon()
    return Polygon()


class Region:
    __array_ufunc__ = None

    def __init__(s, geom):
        s.geom = _clean(geom)

    def __and__(s, o): return Region(s.geom.intersection(o.geom))
    def __or__(s, o): return Region(s.geom.union(o.geom))
    def __invert__(s): return Region(BOX.difference(s.geom))
    __rand__ = __and__
    __ror__ = __or__


class Lin:
    """a*X + b*Y + c"""
    __array_ufunc__ = None

    def __init__(s, a, b, c=0.0):
        s.a, s.b, s.c = float(a), float(b), float(c)

    def _l(s, o): return o if isinstance(o, Lin) else Lin(0, 0, float(o))
    def __add__(s, o): o = s._l(o); return Lin(s.a + o.a, s.b + o.b, s.c + o.c)
    __radd__ = __add__
    def __neg__(s): return Lin(-s.a, -s.b, -s.c)
    def __sub__(s, o): return s + (-s._l(o))
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, k): k = float(k); return Lin(s.a * k, s.b * k, s.c * k)
    __rmul__ = __mul__
    def __truediv__(s, k): return s * (1.0 / float(k))

    def _half(s, sign):
        """region where sign*(aX+bY+c) >= 0"""
        a, b, c = s.a * sign, s.b * sign, s.c * sign
        n = np.hypot(a, b)
        if n == 0:
            return Region(BOX if c >= 0 else Polygon())
        nv = np.array([a, b]) / n; p0 = -c / n * nv; t = np.array([-nv[1], nv[0]])
        pts = [p0 - t * BIG, p0 + t * BIG, p0 + t * BIG + nv * BIG, p0 - t * BIG + nv * BIG]
        return Region(Polygon(pts).intersection(BOX))

    def __ge__(s, k): return (s - k)._half(1)
    def __gt__(s, k): return (s - k)._half(1)
    def __le__(s, k): return (s - k)._half(-1)
    def __lt__(s, k): return (s - k)._half(-1)
    def mean(s): return 1.0


def _sector_poly(c, a0, a1):
    a0 %= 360; a1 %= 360
    span = (a1 - a0) % 360
    if span == 0:
        return Polygon()
    ts = np.radians(a0 + np.linspace(0, span, max(3, int(span) + 2)))
    R = 1e4
    pts = [tuple(c)] + [(c[0] + R * np.cos(t), c[1] + R * np.sin(t)) for t in ts]
    return Polygon(pts)


class VGrid:
    def __init__(s, *a, **k):
        s.X = Lin(1, 0); s.Y = Lin(0, 1)

    def disk(s, c, r): return Region(Point(float(c[0]), float(c[1])).buffer(float(r), quad_segs=QS))
    def sector(s, c, a0, a1): return Region(_sector_poly(np.asarray(c, float), float(a0), float(a1)))
    def bitmap(s, m): return (m.geom.intersection(GRIDBOX), 0, 0)   # the raster Grid's extent


def vraster(polys, *a, **k):
    g = Polygon([tuple(map(float, p)) for p in polys[0]])
    for h in polys[1:]:
        g = g.difference(Polygon([tuple(map(float, p)) for p in h]))
    return g, 0, 0


def vpoly_mask(g, poly):
    return Region(Polygon([tuple(map(float, p)) for p in poly]).buffer(0))


def fuse(geom):
    """fuse pieces that meet along shared radii / tangent points but are split by float noise"""
    return geom.buffer(1e-6, join_style='mitre').buffer(-1e-6, join_style='mitre')


def series(S):
    """{char: fused shapely geometry} for a constructed series"""
    if S == 'E':
        return install()
    import importlib
    reg = importlib.import_module('series%s1961' % S).REG
    return {ch: fuse(fn()) for ch, fn in reg.items()}


def install():
    """patch the construction modules to use the vector back end; returns REG of char -> shapely geometry"""
    import dim1961 as p5, dimE1 as e1, dimE2 as e2, dimE3 as e3, dimE4 as e4, dimE6 as e6
    for m in (p5, e1, e2, e3, e4, e6):
        m.Grid = VGrid
    p5.raster = vraster
    e6.raster = vraster
    e6.poly_mask = vpoly_mask

    def letter_K():                     # dimE2.letter_K with the sign-test band written as a strip
        g = VGrid(); D = e2.DIM; A, CAP = e2.A, e2.CAP
        O, M, N, J = D['O'], D['M'], D['N'], D['J']
        P1 = np.array([A, J]); P2 = np.array([O - N, CAP])
        u = (P2 - P1) / np.linalg.norm(P2 - P1); nu = np.array([-u[1], u[0]])
        yc = CAP - M; crotch = P1 + u * ((yc - P1[1]) / u[1])
        v = crotch - np.array([O, 0.0]); v /= np.linalg.norm(v); nv = np.array([-v[1], v[0]])
        if nv[0] > 0: nv = -nv
        stem = (g.X >= 0) & (g.X <= A) & (g.Y >= 0) & (g.Y <= CAP)
        arm = e3.strip(g, P1, u, A, P1 + nu) & (g.X >= A - 1) & (g.Y <= CAP)
        leg = e3.strip(g, np.array([O, 0.0]), v, A, np.array([O, 0.0]) + nv) & (g.Y >= 0) & (g.Y <= yc)
        return g.bitmap(stem | arm | leg)

    def fig9():                         # the 6 turned through 180 degrees about the centre of its box
        geom, _, _ = e6.fig6()[0]
        return (affinity.rotate(geom, 180, origin=(e6.DIM['B'] / 2, e6.CAP / 2)), 0, 0)

    import seriesE1961 as S
    REG = dict(S.REG)
    REG['K'] = letter_K
    REG['9'] = fig9
    out = {}
    for ch, fn in REG.items():
        r = fn()
        geom = r[0] if isinstance(r[0], (Polygon, MultiPolygon)) else r[0][0]
        out[ch] = fuse(geom)
    return out
