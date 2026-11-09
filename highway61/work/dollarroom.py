"""The dollar sign's widening (dollar.WIDER).  The $ (dollar.py) is the lower-case s, widened, centred on the
capitals' height, with a full-stroke bar through its middle; the bar stands in the s's two counters.  For each series
this measures the clear space either side of the bar in those counters (the widest horizontal run through the centre
line inside each counter, less the bar, halved, in strokes; the smaller of the two) and solves the widening of the
s's ink width that gives dollar.ROOM.
python3 dollarroom.py [series ...]       (default all; copy the factors into dollar.WIDER)"""
import sys, numpy as np
from shapely.geometry import LineString
from scipy.optimize import brentq
import dollar
def clearance(P, T):
    x0, _, x1, _ = P.bounds; cx = (x0 + x1) / 2; out = []
    for ys in (np.linspace(dollar.CAP / 2, dollar.CAP, 150), np.linspace(0, dollar.CAP / 2, 150)):
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
    out = {}
    for S in sys.argv[1:] or ['A', 'B', 'C', 'D', 'E', 'EM', 'F']:
        T = dollar.stroke(S); room = lambda w: clearance(dollar.poly(dollar.body(S, w, T)), T)
        c0 = room(1.0); w = 1.0 if c0 >= dollar.ROOM else brentq(lambda w: room(w) - dollar.ROOM, 1.0, 2.0, xtol=1e-4)
        out[S] = round(w, 3); print('%-2s clearance %.3f as drawn; ink width x %.3f' % (S, c0, w))
    print('WIDER =', out)
