"""The dollar sign: the series' lower-case s, as it is, centred on the capitals' height, with a short stub of a vertical
stroke above it and below it and none through the middle, as the FHWA 2000 fonts and the grid drawing of the $ keep
it open.

The lower case's s is symmetric, as the capitals' S is not.  The series' own s (lccondense's, squared arcs and all;
E(M)'s as constructed) is raised by half the difference between the cap height and the x-height, (690 - 517.5)/2 =
86.25 (an inch at 8 inches), so that it is centred on the I.  Each stub is the s's stroke wide, centred on the s's ink,
and runs from the s's outer edge to the cap line (above) or the baseline (below), so the $ is as tall as the I and the
stubs stand about 0.16 of the s's height beyond it, as in the drawing.  The stubs are spliced into the s's outline
where their sides cross its outer edge, so the outline stays the s's exact segments plus four straight lines.

make(S) -> strokeop segment contours."""
import numpy as np
import strokeop, lcE, lccondense
CAP, X = 690.0, 517.5
RAISE = (CAP - X) / 2                          # the s centred on the capitals' height
def lc_s(S):
    """the series' lower-case s (segment contours) and its stroke"""
    if S == 'EM':
        g, pr = lcE.series()['s']; return strokeop.from_lc(g, pr), lccondense.V_EM
    return lccondense.series(S, only=['s'])['s'], lccondense.stroke(S)
def _move(C, dy):
    out = []
    for segs in C:
        cs = []
        for s in segs:
            M = lambda p: np.array([p[0], p[1] + dy])
            if s[0] == 'L': cs.append(('L', M(s[1]), M(s[2])))
            elif s[0] == 'B': cs.append(('B', M(s[1]), M(s[2]), M(s[3]), M(s[4])))
            else:
                cx, cy, a, b = s[3]; cs.append(('E', M(s[1]), M(s[2]), (cx, cy + dy, a, b), s[4], s[5]))
        out.append(cs)
    return out
def _bez(s, t):
    m = 1 - t; return m ** 3 * s[1] + 3 * m * m * t * s[3] + 3 * m * t * t * s[4] + t ** 3 * s[2]
def _crossings(s, x):
    """(parameter, point) where segment s crosses the vertical line at x, inside it"""
    out = []
    if s[0] == 'L':
        a, b = s[1], s[2]
        if (a[0] - x) * (b[0] - x) < 0:
            t = (x - a[0]) / (b[0] - a[0]); out.append((t, a + (b - a) * t))
    elif s[0] == 'B':
        p0, p3, p1, p2 = s[1][0], s[2][0], s[3][0], s[4][0]
        co = [-p0 + 3 * p1 - 3 * p2 + p3, 3 * p0 - 6 * p1 + 3 * p2, -3 * p0 + 3 * p1, p0 - x]
        for r in np.roots(co):
            if abs(r.imag) < 1e-9 and 1e-9 < r.real < 1 - 1e-9: out.append((r.real, _bez(s, r.real)))
    else:
        cx, cy, a, b = s[3]; t0, sw = s[4], s[5]
        if abs(x - cx) < a:
            base = np.arccos((x - cx) / a)
            for r in (base, -base):
                for k in range(-3, 4):
                    u = (r + 2 * np.pi * k - t0) / sw
                    if 1e-9 < u < 1 - 1e-9: out.append((u, np.array([x, cy + b * np.sin(t0 + sw * u)])))
    return out
def _split(s, u, p):
    if s[0] == 'L': return ('L', s[1], p), ('L', p, s[2])
    if s[0] == 'E':
        t0, sw = s[4], s[5]
        return ('E', s[1], p, s[3], t0, sw * u), ('E', p, s[2], s[3], t0 + sw * u, sw * (1 - u))
    P0, P1, P2, P3 = s[1], s[3], s[4], s[2]                 # de Casteljau
    a, b, c = P0 + (P1 - P0) * u, P1 + (P2 - P1) * u, P2 + (P3 - P2) * u
    d, e = a + (b - a) * u, b + (c - b) * u
    return ('B', P0, p, a, d), ('B', p, P3, e, c)
def _cut(c, x, top):
    """c with a joint added where the vertical line at x crosses it highest (top) or lowest; returns (c, index of the
    segment that starts at the joint)"""
    hits = [(i, u, p) for i, s in enumerate(c) for u, p in _crossings(s, x)]
    i, u, p = max(hits, key=lambda h: h[2][1] if top else -h[2][1])
    a, b = _split(c[i], u, p)
    return c[:i] + [a, b] + c[i + 1:], i + 1
def _stub(c, xl, xr, y, top):
    """splice in a stub between xl and xr from the contour's outer edge to the line y"""
    c, _ = _cut(c, xl, top); c, _ = _cut(c, xr, top)
    n = len(c); pt = lambda i: c[i % n][1]
    jl = next(i for i in range(n) if abs(pt(i)[0] - xl) < 1e-9 and (pt(i)[1] > RAISE + X / 2) == top)
    jr = next(i for i in range(n) if abs(pt(i)[0] - xr) < 1e-9 and (pt(i)[1] > RAISE + X / 2) == top)
    # the stretch of outline between the two joints that holds the extreme (topmost or lowest) point goes
    ext = lambda i0, i1: [k % n for k in range(i0, i1 if i1 > i0 else i1 + n)]
    ys = lambda ks: max(c[k][2][1] for k in ks) if top else -min(c[k][2][1] for k in ks)
    run = max((ext(jl, jr), (jl, jr)), (ext(jr, jl), (jr, jl)), key=lambda r: ys(r[0]))
    ks, (j0, j1) = run
    A, B = pt(j0), pt(j1)
    new = [('L', A, np.array([A[0], y])), ('L', np.array([A[0], y]), np.array([B[0], y])), ('L', np.array([B[0], y]), B)]
    keep = [c[k % n] for k in range(j1, j1 + n) if k % n not in ks]
    return keep + new
def make(S):
    C, V = lc_s(S); C = _move(C, RAISE)
    assert len(C) == 1, 'the s is one contour'
    P = np.vstack(strokeop.to_polygons(C)); cx = (P[:, 0].min() + P[:, 0].max()) / 2
    c = _stub(C[0], cx - V / 2, cx + V / 2, CAP, True)
    c = _stub(c, cx - V / 2, cx + V / 2, 0.0, False)
    return [c]
