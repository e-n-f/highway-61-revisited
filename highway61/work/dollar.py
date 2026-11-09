"""The dollar sign: the series' lower-case s enlarged to the figures' height, with a vertical bar through its middle.

The lower case's s is symmetric, as the capitals' S is not.  The series' own s (as lccondense makes it, before its
arcs are squared: lines and ellipses) is scaled to the capital S's ink width and to the figures' height with their
overshoots; the width less the stroke and the height less the stroke scale (about 1.2-1.6 and 1.4: nearly uniform), and
every edge is then moved along its normal so that every stroke is the capitals' stroke T again (strokeop's
direction-dependent offset, as condense2 makes it).  The bar is T wide (the alphabet has no thin strokes), centred on
the body, and stands 1/8 cap height (1 inch at 8 inches) beyond the cap line and the baseline, as on signs
("MINIMUM FINE $250").  Where the bar would leave less than ROOM strokes clear beside it in the body's two counters
(Series B), the body is made wider by WIDER (dollarroom.py solves it).  The outline stays exact: the bar is united
on a fine flattening whose edges lcoutline puts back on the body's lines and ellipses.

make(S) -> strokeop segment contours."""
import numpy as np
from shapely.geometry import Polygon, box
import strokeop, lcE, lccondense, lcoutline, vgeom
CAP, OV = 690.0, 690 * (1 / 8) / 8            # the figures' height and overshoot (1/8 in at 8 in, 10.8 units)
EXT = 690 / 8                                  # the bar beyond the cap line and the baseline
ROOM = 0.40                                    # least clear space beside the bar in either counter, in strokes
WIDER = {'B': 1.142}                             # body width over the capital S's (dollarroom.py)
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
def body(S, W, T):
    """the s at the figures' height and ink width W, every stroke T"""
    if S not in _S: _S[S] = lc_s(S)
    C0, V = _S[S]; P = np.vstack(strokeop.to_polygons(C0)); lo, hi = P.min(0), P.max(0)
    fx = (W - T) / (hi[0] - lo[0] - V); c = (CAP + 2 * OV - T) / (hi[1] - lo[1] - V)
    b = -OV - c * lo[1] - (c * V - T) / 2      # the bottom overshoot lands on -OV after the offset
    amt = lambda n: -(V * np.sqrt((fx * n[0]) ** 2 + (c * n[1]) ** 2) - T) / 2
    out = []
    for cc in _affine(C0, fx, c, b):
        cc = strokeop.clean(cc); tg = [amt] * len(cc)
        cc, tg = strokeop.offset_contour(cc, tg, tags=tg); out.append(strokeop.clean(cc, tags=tg)[0])
    return out
def width_stroke(S):
    """the capital S's ink width and the capitals' stroke (E(M): Series E's S widened by its offset, 2 x 9.70)"""
    g = vgeom.series('E' if S == 'EM' else S)['S']; W = g.bounds[2] - g.bounds[0]
    return (W + (690 / 5 - lccondense.A['E']), 690 / 5) if S == 'EM' else (W, lccondense.A[S])
def poly(C, step=0.01):
    P = Polygon()
    for r in strokeop.to_polygons(C, step=step): P = P.symmetric_difference(Polygon(r).buffer(0))
    return P
def make(S):
    W, T = width_stroke(S); C = body(S, W * WIDER.get(S, 1.0), T)
    P = poly(C, step=np.radians(0.25)); x0, _, x1, _ = P.bounds; cx = (x0 + x1) / 2
    P = vgeom.fuse(P.union(box(cx - T / 2, -EXT, cx + T / 2, CAP + EXT)))
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
