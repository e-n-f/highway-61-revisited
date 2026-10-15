"""Exact outlines from the vector constructions: recover the straight edges and circular arcs from a shapely
polygon's vertices and draw them as lines and cubic Bezier arcs (at most 45 degrees per segment)."""
import numpy as np
from shapely.geometry.polygon import orient

TOL_IN = 2e-4          # interior arc vertices are on their circle to float precision (tessellation is exact)
MAX_STEP = np.radians(4)  # largest chord angle of a tessellated arc
TOL_END = 0.02         # run end points are intersections, which lie on a chord, so allow a little more


def _circ(a, b, c):
    d = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]))
    if abs(d) < 1e-12:
        return None, None
    sa, sb, sc = a @ a, b @ b, c @ c
    ux = (sa * (b[1] - c[1]) + sb * (c[1] - a[1]) + sc * (a[1] - b[1])) / d
    uy = (sa * (c[0] - b[0]) + sb * (a[0] - c[0]) + sc * (b[0] - a[0])) / d
    cen = np.array([ux, uy])
    return cen, float(np.linalg.norm(a - cen))


def _clean(pts):
    pts = [p for i, p in enumerate(pts) if np.linalg.norm(p - pts[i - 1]) > 1e-4] if len(pts) > 1 else pts
    changed = True
    while changed and len(pts) > 3:
        changed = False
        out = []
        n = len(pts)
        for i in range(n):
            a, b, c = pts[i - 1], pts[i], pts[(i + 1) % n]
            e1, e2 = b - a, c - b
            l1, l2 = np.linalg.norm(e1), np.linalg.norm(e2)
            cr = (e1[0] * e2[1] - e1[1] * e2[0]) / (l1 * l2)
            dt = (e1 @ e2) / (l1 * l2)
            if (abs(cr) < 1e-7 and dt > 0) or dt < -0.999999:     # collinear vertex, or a hairline spike
                changed = True
                continue
            out.append(b)
        pts = out
    return pts


def segments(ring):
    """ring: closed coordinate sequence -> list of ('L', p0, p1) and ('A', p0, p1, centre, r, sweep_radians)"""
    pts = _clean([np.array(p, float) for p in list(ring)[:-1]])
    n = len(pts)
    el = [np.linalg.norm(pts[(i + 1) % n] - pts[i]) for i in range(n)]
    s0 = int(np.argmax(el))                     # start at the longest edge's start: a corner, not mid-arc
    P = [pts[(s0 + i) % n] for i in range(n)] + [pts[s0]]
    segs = []
    i = 0
    while i < n:
        arc = None
        if i + 4 <= n:
            c, r = _circ(P[i + 1], P[i + 2], P[i + 3])
            # tessellated arcs have short chords; this rejects e.g. rectangle corners, which are concyclic too
            fine = lambda k: np.linalg.norm(P[k + 1] - P[k]) < r * MAX_STEP
            if c is not None and r < 2e4 and fine(i + 1) and fine(i + 2):
                j = i + 3
                while j + 1 <= n and abs(np.linalg.norm(P[j + 1] - c) - r) < TOL_IN + r * 1e-7 and fine(j):
                    j += 1
                if j - (i + 1) >= 3:
                    a0 = i if abs(np.linalg.norm(P[i] - c) - r) < TOL_END and fine(i) else i + 1
                    a1 = j + 1 if j + 1 <= n and abs(np.linalg.norm(P[j + 1] - c) - r) < TOL_END and fine(j) else j
                    m = (a0 + a1) // 2
                    c, r = _circ(P[a0 + 1] if a0 + 1 < m else P[a0], P[m], P[a1 - 1] if a1 - 1 > m else P[a1])
                    arc = (a0, a1, c, r)
        if arc is None:
            segs.append(('L', P[i], P[i + 1])); i += 1
            continue
        a0, a1, c, r = arc
        if a0 > i:
            segs.append(('L', P[i], P[a0]))
        sw = 0.0
        for k in range(a0, a1):                  # signed sweep, summed step by step
            u, v = P[k] - c, P[k + 1] - c
            sw += np.arctan2(u[0] * v[1] - u[1] * v[0], u @ v)
        segs.append(('A', P[a0], P[a1], c, r, sw))
        i = a1
    # merge neighbouring arcs on the same circle
    out = []
    for s in segs:
        if out and s[0] == 'A' and out[-1][0] == 'A' and np.linalg.norm(s[3] - out[-1][3]) < 0.01 and abs(s[4] - out[-1][4]) < 0.01:
            p = out[-1]; out[-1] = ('A', p[1], s[2], p[3], p[4], p[5] + s[5])
        else:
            out.append(s)
    return out


def contours(geom):
    geom = orient(geom, sign=-1.0) if geom.geom_type == 'Polygon' else type(geom)([orient(g, sign=-1.0) for g in geom.geoms])
    polys = [geom] if geom.geom_type == 'Polygon' else list(geom.geoms)
    rings = []
    for p in polys:
        rings.append(p.exterior.coords)
        rings.extend(h.coords for h in p.interiors)
    return [segments(r) for r in rings]


def beziers(seg):
    """an arc segment as a list of cubic (p0, c1, c2, p1), at most 45 degrees each (radial error < 4e-6 r)"""
    _, p0, p1, c, r, sw = seg
    k = max(1, int(np.ceil(abs(sw) / (np.pi / 4) - 1e-9)))
    t0 = np.arctan2(p0[1] - c[1], p0[0] - c[0])
    d = sw / k
    h = 4 / 3 * np.tan(d / 4) * r
    out = []
    prev = p0
    for i in range(k):
        ta, tb = t0 + d * i, t0 + d * (i + 1)
        q = p1 if i == k - 1 else c + r * np.array([np.cos(tb), np.sin(tb)])
        ua = np.array([-np.sin(ta), np.cos(ta)]); ub = np.array([-np.sin(tb), np.cos(tb)])
        out.append((prev, prev + h * ua, q - h * ub, q))
        prev = q
    return out


def draw(geom, pen, dx=0.0, dy=0.0):
    T = lambda p: (float(p[0] + dx), float(p[1] + dy))
    for segs in contours(geom):
        pen.moveTo(T(segs[0][1]))
        for k, s in enumerate(segs):
            if s[0] == 'L':
                if k < len(segs) - 1:            # the closing edge is implied
                    pen.lineTo(T(s[2]))
            else:
                for b in beziers(s):
                    pen.curveTo(T(b[1]), T(b[2]), T(b[3]))
        pen.closePath()


def describe(geom):
    """count of lines and arcs per contour, and the distinct arc radii"""
    rep = []
    for segs in contours(geom):
        rep.append((sum(s[0] == 'L' for s in segs), sum(s[0] == 'A' for s in segs),
                    sorted({round(s[4], 1) for s in segs if s[0] == 'A'})))
    return rep
