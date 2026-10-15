"""fit a construction (shapely geometry, font units) onto its 1961 drawing: the pages draw letters as outlines, so the
fit minimises the distance from the construction's boundary to the nearest drawn line.  Reports that distance in font
units and writes an overlay (scan grey, construction red)."""
import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt
from scipy.optimize import minimize
SCRATCH = '/tmp/claude-0/-home-user-standard-alphabets/d2c3cfcf-66d0-5f41-8524-bbf06168e396/scratchpad/'
_pages = {}
def page(p):
    if p not in _pages: _pages[p] = np.array(Image.open('p1961/stdalphabets1961-%s.png' % p).convert('L')) < 128
    return _pages[p]
def boundary_pts(g, step=2.0):
    polys = [g] if g.geom_type == 'Polygon' else list(g.geoms)
    pts = []
    for p in polys:
        for r in [p.exterior] + list(p.interiors):
            n = max(8, int(r.length / step)); pts += [r.interpolate(t, normalized=True).coords[0] for t in np.linspace(0, 1, n, endpoint=False)]
    return np.array(pts)
def fit(g, p, box, scale=0.858, out=None):
    x0, y0, x1, y1 = box
    dark = page(p)[y0:y1, x0:x1]
    dt = distance_transform_edt(~dark)
    H, W = dt.shape
    P = boundary_pts(g)
    def cost(v, clip=12.0):
        s, ox, oy = v
        X = np.round(ox + s * P[:, 0]).astype(int); Y = np.round(oy - s * P[:, 1]).astype(int)
        ok = (X >= 0) & (X < W) & (Y >= 0) & (Y < H)
        d = np.full(len(P), clip); d[ok] = np.minimum(dt[Y[ok], X[ok]], clip)
        return d.mean()
    gx0, gy0, gx1, gy1 = g.bounds
    best = None
    for s in scale * np.array([0.98, 1.0, 1.02]):
        for ox in np.arange(-s * gx0, W - s * gx1, 4):
            for oy in np.arange(s * gy1, H + s * gy0, 4):
                c = cost((s, ox, oy))
                if best is None or c < best[0]: best = (c, (s, ox, oy))
    r = minimize(lambda v: cost(v, 6.0), best[1], method='Nelder-Mead', options={'xatol': 0.05, 'fatol': 1e-4})
    s, ox, oy = r.x
    X = ox + s * P[:, 0]; Y = oy - s * P[:, 1]
    ok = (X >= 0) & (X < W - 1) & (Y >= 0) & (Y < H - 1)
    d = dt[np.round(Y[ok]).astype(int), np.round(X[ok]).astype(int)] / s
    if out:
        img = np.stack([np.where(dark, 150, 255)] * 3, -1).astype(np.uint8)
        for x, y in zip(np.round(X[ok]).astype(int), np.round(Y[ok]).astype(int)):
            img[max(0, y - 1):y + 2, max(0, x - 1):x + 2] = (220, 0, 0)
        Image.fromarray(img).save(SCRATCH + out)
    return {'scale': round(s, 4), 'median': round(float(np.median(d)), 2), 'p95': round(float(np.percentile(d, 95)), 2), 'max': round(float(d.max()), 1)}
