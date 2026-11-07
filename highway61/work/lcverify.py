"""The E(M) lower case over the 1966 sheets: python3 lcverify.py -> ../outputs/verify-EM-lc.png

Each letter's outline, as built in otf/Highway61Revisited-EM.otf, drawn in red on the scan of its 1966 grid sheet
(lc1966/p33-p35.jpg, ink lightened to grey so the line shows on it).  The outline is placed by the scan fit (the
construction's offset ox, lc1966/fit.json) and mapped onto the page through its fitted grid lines (lcgrid), so
it sits where the fit puts it on the drawing, including the page's own distortions.  Under each letter: the median
and 95th-percentile distance between the outline and the scan's contour (font units; the fit allows lcmodel.SPREAD
for printing spread, not applied here).  SCALE enlarges the scan for detail."""
import os, sys, string, numpy as np
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import RecordingPen
from scipy.spatial import cKDTree
from PIL import Image, ImageDraw, ImageFont
import lcscan, lcview, lcE
W = os.path.dirname(os.path.abspath(__file__)); DIST = os.path.join(W, '..', '..'); OUT = os.path.join(W, '..', 'outputs')
SQ = 690 / 8; SCALE = 2; SS = 3; COLS = 7; MARGIN = 45
LABEL = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 26)


def rings(gs, name):
    """the glyph's contours as point arrays (curves flattened)"""
    rp = RecordingPen(); gs[name].draw(rp); out = []; cur = []
    for op, a in rp.value:
        if op == 'moveTo': cur = [np.array(a[0], float)]
        elif op == 'lineTo': cur.append(np.array(a[0], float))
        elif op == 'curveTo':
            p0 = cur[-1]; p1, p2, p3 = map(np.array, a)
            for t in np.linspace(0, 1, 25)[1:]: cur.append((1-t)**3*p0 + 3*(1-t)**2*t*p1 + 3*(1-t)*t*t*p2 + t**3*p3)
        elif op == 'closePath': out.append(np.array(cur))
    return out


def to_px(page, gx, gy):
    """grid squares (gx right, gy down, from the page's first lines) -> scan pixels: lcgrid.Page.to_grid inverted"""
    V = np.array(page.V); H = np.array(page.H)
    def lerp(g, A, B, t):                      # position of fractional line g at coordinate t (lines c = a + b t)
        i = np.clip(np.floor(g).astype(int), 0, len(A) - 2); f = g - i
        return (A[i] + B[i] * t) * (1 - f) + (A[i + 1] + B[i + 1] * t) * f
    y = lerp(gy, H[:, 0], H[:, 1], np.zeros_like(gy) + 1000.0); x = lerp(gx, V[:, 0], V[:, 1], y)
    for _ in range(6): y = lerp(gy, H[:, 0], H[:, 1], x); x = lerp(gx, V[:, 0], V[:, 1], y)
    return x, y


def origin(ch):
    """lcview.units' origin for the letter (squares, page grid, y up): the grid line left of it and the baseline"""
    P = np.vstack(lcscan.letters()[ch]['contours'])
    if ch in lcview.DESC: y0 = np.round(P[:, 1].max()) - (8 if ch == 'j' else 6)
    else: y0 = np.round(P[:, 1].min())
    return np.floor(P[:, 0].min() + 0.2), y0


def tile(ch, gs, cm, dims, series):
    L = lcscan.letters()[ch]; page = lcscan.page(L['page']); x0, y0 = origin(ch)
    R = rings(gs, cm[ord(ch)])
    g = series[ch][0]; ox = dims.get(ch, {}).get('ox', 0.0) * SQ
    allR = np.vstack(R); dx = (g.bounds[0] + g.bounds[2]) / 2 + ox - (allR[:, 0].min() + allR[:, 0].max()) / 2
    R = [r + [dx, 0] for r in R]                                       # font units, the scan fit's frame
    scan = lcview.units(ch); S = np.vstack(scan)
    M = np.vstack([np.vstack([np.linspace(a, b, max(2, int(np.hypot(*(b - a)) / 2)), endpoint=False)
                              for a, b in zip(r, np.roll(r, -1, axis=0))]) for r in R])
    d = np.r_[cKDTree(M).query(S)[0], cKDTree(S).query(M)[0]]
    px = lambda r: to_px(page, r[:, 0] / SQ + x0, -(r[:, 1] / SQ + y0))
    allpx = [np.c_[px(r)] for r in R] + [np.c_[px(c)] for c in scan]
    Pp = np.vstack([np.c_[p[0], p[1]] if p.shape[0] == 2 else p for p in allpx])
    bx0, by0 = np.floor(Pp.min(0)) - MARGIN; bx1, by1 = np.ceil(Pp.max(0)) + MARGIN
    crop = Image.fromarray(page.im).crop((int(bx0), int(by0), int(bx1), int(by1)))
    crop = crop.resize((crop.width * SCALE, crop.height * SCALE), Image.LANCZOS)
    base = np.asarray(crop, float); base = 255 - (255 - base) * 0.45            # ink to grey
    big = Image.new('L', (crop.width * SS, crop.height * SS), 0); dr = ImageDraw.Draw(big)
    for r in R:
        x, y = px(r); pts = [((a - bx0) * SCALE * SS, (b - by0) * SCALE * SS) for a, b in zip(x, y)]
        dr.line(pts + pts[:1], fill=255, width=SS + 1)
    m = np.asarray(big.resize(crop.size, Image.BOX), float)[..., None] / 255
    rgb = np.repeat(base[..., None], 3, axis=2) * (1 - m) + np.array([220, 0, 0]) * m
    im = Image.new('RGB', (crop.width, crop.height + 44), 'white'); im.paste(Image.fromarray(rgb.astype(np.uint8)), (0, 0))
    ImageDraw.Draw(im).text((8, crop.height + 8), '%s   median %.1f  p95 %.1f' % (ch, np.median(d), np.percentile(d, 95)),
                            font=LABEL, fill=(0, 0, 0))
    return im, d


def main():
    f = TTFont(os.path.join(DIST, 'otf', 'Highway61Revisited-EM.otf')); gs = f.getGlyphSet(); cm = f.getBestCmap()
    dims = lcE.dims(); series = lcE.series(); tiles = []; stats = []
    for ch in string.ascii_lowercase:
        im, d = tile(ch, gs, cm, dims, series); tiles.append(im); stats.append((ch, np.median(d), np.percentile(d, 95)))
    rows = [tiles[i:i + COLS] for i in range(0, len(tiles), COLS)]
    Wd = max(sum(t.width for t in r) + 20 * (len(r) - 1) for r in rows); Ht = sum(max(t.height for t in r) for r in rows) + 30 * len(rows)
    out = Image.new('RGB', (Wd + 40, Ht + 40), 'white'); y = 20
    for r in rows:
        x = 20
        for t in r: out.paste(t, (x, y)); x += t.width + 20
        y += max(t.height for t in r) + 30
    out.save(os.path.join(OUT, 'verify-EM-lc.png'))
    for ch, md, p95 in stats: print(ch, 'median %.1f p95 %.1f' % (md, p95))
    print('size', out.size)


if __name__ == '__main__': main()
