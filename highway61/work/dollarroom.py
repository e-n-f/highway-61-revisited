"""The $'s widening in the narrow series (build1961.DOLLAR_F).  The $ is the S with a full-stroke bar through its
middle; the bar stands in the round counters under the S's top arc and over its bottom arc.  For each series this
measures the clear space either side of a centred bar of the capitals' stroke T in those two counters (the widest
horizontal run through the centre line inside each counter, less T, halved, in strokes; the smaller of the two), and
solves the widening f (build1961.widen: x scaled by f, strokes kept at T) that gives Series D's clearance.
python3 dollarroom.py [series ...]       (default B C; copy the factors into build1961.DOLLAR_F)"""
import sys, numpy as np
from shapely.geometry import LineString, Polygon
from scipy.optimize import brentq
import vgeom, strokeop, lccondense
CAP = 690.0
def _poly(C):
    P = Polygon()
    for r in strokeop.to_polygons(C, step=0.01): P = P.symmetric_difference(Polygon(r).buffer(0))
    return P
def widened_S(S, f):
    T = lccondense.A[S]; g = vgeom.series(S)['S']
    if f == 1: return g, T
    C = strokeop.condense2(strokeop.from_caps(g), f, (T - T * f) / 2, 0.0, CAP, V=T, proportional=False)
    return _poly(C), T
def clearance(P, T):
    x0, _, x1, _ = P.bounds; cx = (x0 + x1) / 2; out = []
    for ys in (np.linspace(CAP / 2, CAP, 150), np.linspace(0, CAP / 2, 150)):
        best = 0.0
        for y in ys:
            g = LineString([(x0 - 1, y), (x1 + 1, y)]).difference(P)
            for s in getattr(g, 'geoms', [g]):
                if s.is_empty: continue
                a, b = s.bounds[0], s.bounds[2]
                if a < cx < b and a > x0 + 0.5 and b < x1 - 0.5: best = max(best, b - a)
        out.append((best - T) / 2 / T)
    return min(out)
if __name__ == '__main__':
    target = clearance(*widened_S('D', 1))
    print('Series D clearance %.3f strokes' % target)
    for S in sys.argv[1:] or ['B', 'C']:
        c0 = clearance(*widened_S(S, 1))
        f = 1.0 if c0 >= target else brentq(lambda f: clearance(*widened_S(S, f)) - target, 1.0, 2.5, xtol=1e-4)
        print('%s: clearance %.3f as drawn; widened by %.3f' % (S, c0, f))
