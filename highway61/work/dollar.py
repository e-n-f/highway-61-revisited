"""The dollar sign: the series' lower-case s, widened, at its own height and centred on the capitals' height, with a
vertical bar through its middle.

The lower case's s is symmetric, as the capitals' S is not.  The series' own s (as lccondense makes it, before its
arcs are squared: lines and ellipses; E(M)'s as constructed) keeps its height and is raised by half the difference
between the cap height and the x-height, (690 - 517.5)/2 = 86.25 (an inch at 8 inches), so that it is centred on the
I.  It is widened (its width less the stroke scaled) by WIDER, and every edge then moved along its normal so that
every stroke is the capitals' stroke T (strokeop's direction-dependent offset, as condense2 makes it; in E(M) this
also brings the s's 134.8 to the capitals' 138).  The bar is T wide (the alphabet has no thin strokes), centred on
the body, and stands EXT = 86.25 beyond it, so it runs from the baseline to the cap line, as tall as the I.  WIDER is
solved (dollarroom.py) so that the bar leaves ROOM strokes clear beside it in the body's two counters.  The outline
stays exact: the bar is united on a fine flattening whose edges lcoutline puts back on the body's lines and ellipses.

make(S) -> strokeop segment contours."""
import numpy as np
from shapely.geometry import Polygon, box
import strokeop, lcE, lccondense, lcoutline, vgeom
CAP, X = 690.0, 517.5
RAISE = (CAP - X) / 2                          # the s centred on the capitals' height
EXT = RAISE                                    # the bar beyond the s: from the baseline to the cap line
ROOM = 0.40                                    # clear space beside the bar in the counters, in strokes
WIDER = {'A': 1.246, 'B': 1.375, 'C': 1.058, 'D': 1.014, 'E': 1.002, 'EM': 1.113}   # the s's ink width x this, to
                                               # leave ROOM beside the bar (dollarroom.py); F has room as it is
def lc_s(S):
    """the series' lower-case s (segment contours, unsquared) and its stroke"""
    if S == 'EM':
        g, pr = lcE.series()['s']; return strokeop.from_lc(g, pr), lccondense.V_EM
    old = lccondense.SQUARE; lccondense.SQUARE = False
    try: return lccondense.series(S, only=['s'])['s'], lccondense.stroke(S)
    finally: lccondense.SQUARE = old
def _affine(C, fx, c, b):
    out = []
    for segs in C:
        cs = []
        for s in segs:
            M = lambda p: np.array([p[0] * fx, p[1] * c + b])
            if s[0] == 'L': cs.append(('L', M(s[1]), M(s[2])))
            else:
                cx, cy, a, bb = s[3]; cs.append(('E', M(s[1]), M(s[2]), (cx * fx, cy * c + b, a * fx, bb * c), s[4], s[5]))
        out.append(cs)
    return out
_S = {}
def body(S, w, T):
    """the s raised by RAISE, its ink width multiplied by w, every stroke T"""
    if S not in _S: _S[S] = lc_s(S)
    C0, V = _S[S]; P = np.vstack(strokeop.to_polygons(C0)); W0 = P[:, 0].max() - P[:, 0].min()
    fx = (w * W0 - T) / (W0 - V)               # the width less the stroke scales
    amt = lambda n: -(V * np.sqrt((fx * n[0]) ** 2 + n[1] ** 2) - T) / 2
    out = []
    for cc in _affine(C0, fx, 1.0, RAISE):
        cc = strokeop.clean(cc); tg = [amt] * len(cc)
        cc, tg = strokeop.offset_contour(cc, tg, tags=tg); out.append(strokeop.clean(cc, tags=tg)[0])
    return out
def stroke(S):
    """the capitals' stroke"""
    return 690 / 5 if S == 'EM' else lccondense.A[S]
def poly(C, step=0.01):
    P = Polygon()
    for r in strokeop.to_polygons(C, step=step): P = P.symmetric_difference(Polygon(r).buffer(0))
    return P
def make(S):
    T = stroke(S); C = body(S, WIDER.get(S, 1.0), T)
    P = poly(C, step=np.radians(0.25)); x0, y0, x1, y1 = P.bounds; cx = (x0 + x1) / 2
    P = vgeom.fuse(P.union(box(cx - T / 2, RAISE - EXT, cx + T / 2, RAISE + X + EXT)))
    prims = [np.array(sg[3], float) for c in C for sg in c if sg[0] == 'E']
    joints = [sg[1] for c in C for sg in c]
    out = []
    for c in lcoutline.contours(P, prims):
        # where an arc meets another arc or the spine, the classifier can leave a speck of line between them: the
        # two meet again at the body's own joint
        changed = True
        while changed:
            changed = False; n = len(c)
            for i in range(n):
                a, l, b = c[i - 1], c[i], c[(i + 1) % n]
                if l[0] != 'L' or 'E' not in (a[0], b[0]) or np.hypot(*(l[2] - l[1])) >= 3: continue
                q = min(joints, key=lambda p: min(np.hypot(*(p - l[1])), np.hypot(*(p - l[2]))))
                if min(np.hypot(*(q - l[1])), np.hypot(*(q - l[2]))) >= 3: continue
                c[i - 1] = strokeop._rebuild(a, a[1], q); c[(i + 1) % n] = strokeop._rebuild(b, q, b[2])
                c = [sg for k, sg in enumerate(c) if k != i]; changed = True; break
        out.append(c)
    return out
