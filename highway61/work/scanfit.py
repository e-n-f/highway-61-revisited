"""fit a construction (shapely geometry, font units) onto its 1961 drawing: the pages draw letters as outlines, so the
fit minimises the distance from the construction's boundary to the nearest drawn line.  Reports that distance in font
units and writes an overlay (scan grey, construction red)."""
import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt
from scipy.optimize import minimize
SCRATCH = '/tmp/claude-0/-home-user-standard-alphabets/d2c3cfcf-66d0-5f41-8524-bbf06168e396/scratchpad/'
_pages = {}
_last = {}
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
    _last['v'] = r.x
    X = ox + s * P[:, 0]; Y = oy - s * P[:, 1]
    ok = (X >= 0) & (X < W - 1) & (Y >= 0) & (Y < H - 1)
    d = dt[np.round(Y[ok]).astype(int), np.round(X[ok]).astype(int)] / s
    if out:
        img = np.stack([np.where(dark, 150, 255)] * 3, -1).astype(np.uint8)
        for x, y in zip(np.round(X[ok]).astype(int), np.round(Y[ok]).astype(int)):
            img[max(0, y - 1):y + 2, max(0, x - 1):x + 2] = (220, 0, 0)
        Image.fromarray(img).save(SCRATCH + out)
    return {'scale': round(s, 4), 'median': round(float(np.median(d)), 2), 'p95': round(float(np.percentile(d, 95)), 2), 'max': round(float(d.max()), 1)}


def edge_samples(g, p, box, scale=0.858, reach=12):
    """fit g to the drawing, then snap each boundary point of g to the nearest drawn pixel (within reach px); returns
    (model points, snapped points) in font units, for fitting the drawing's own arcs"""
    from scipy.ndimage import distance_transform_edt
    res = fit(g, p, box, scale)
    x0, y0, x1, y1 = box
    dark = page(p)[y0:y1, x0:x1]
    # recover the transform the fit chose
    P = boundary_pts(g, 1.0)
    s = res['scale']
    best = _last['v']
    s, ox, oy = best
    dt, (iy, ix) = distance_transform_edt(~dark, return_indices=True)
    X = np.round(ox + s * P[:, 0]).astype(int); Y = np.round(oy - s * P[:, 1]).astype(int)
    ok = (X >= 0) & (X < dark.shape[1]) & (Y >= 0) & (Y < dark.shape[0])
    ok[ok] &= dt[Y[ok], X[ok]] <= reach
    sx = ix[Y[ok], X[ok]]; sy = iy[Y[ok], X[ok]]
    S = np.c_[(sx - ox) / s, (oy - sy) / s]
    return P[ok], S


def outline(p, box, width_units):
    """the drawn letter's outer contour (background flood-filled from the box's corner), in font units: x from the
    leftmost point, y from the lowest point, scaled so the contour is width_units wide"""
    import cv2
    x0, y0, x1, y1 = box
    dark = page(p)[y0:y1, x0:x1].astype(np.uint8)
    dark = cv2.dilate(dark, np.ones((3, 3), np.uint8))
    ff = (1 - dark) * 255
    mask = np.zeros((ff.shape[0] + 2, ff.shape[1] + 2), np.uint8)
    cv2.floodFill(ff, mask, (0, 0), 128)
    inside = (ff != 128).astype(np.uint8)
    cs, _ = cv2.findContours(inside, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    c = max(cs, key=cv2.contourArea)[:, 0, :].astype(float)
    s = (c[:, 0].max() - c[:, 0].min()) / width_units
    return np.c_[(c[:, 0] - c[:, 0].min()) / s, (c[:, 1].max() - c[:, 1]) / s], s


def circle_fit(P):
    x, y = P[:, 0], P[:, 1]; M = np.c_[2 * x, 2 * y, np.ones(len(x))]; b = x * x + y * y
    c = np.linalg.lstsq(M, b, rcond=None)[0]; return c[:2], float(np.sqrt(c[2] + c[0] ** 2 + c[1] ** 2))
