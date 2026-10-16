"""Highway 61 Revisited <series>: capitals and figures as exact outlines from the 1961 dimensioned drawings (vgeom + outline),
lowercase and punctuation taken over from the matching Highway Even 1966 font.  Writes an OTF (cubic arcs,
fractional coordinates) and a TTF (quadratic, integer coordinates).

usage: python3 build1961.py [series]      (B C D E F, and EM)

E: 1961 Series E capitals and figures with the lower case narrowed from E(M) weight to Series E weight (strokeop).
EM: the 1966 recommendation for capitals to go with the E(M) lower case, "Series E uppercase, but with the stroke
width widened to approximately one-fifth of the letter height": the Series E outlines offset by (690/5 - A)/2 and
squeezed back to the cap height (strokeop), with the E(M) lower case as constructed.
"""
import sys, os
from fontTools.ttLib import TTFont
from fontTools.fontBuilder import FontBuilder
from fontTools.pens.t2CharStringPen import T2CharStringPen
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.pens.cu2quPen import Cu2QuPen
from fontTools.pens.boundsPen import BoundsPen
from fontTools.pens.reverseContourPen import ReverseContourPen
import vgeom, outline, lcoutline, strokeop

S = sys.argv[1] if len(sys.argv) > 1 else 'E'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'outputs')
BASE = TTFont(os.path.join(OUT, 'Highway1966-EM-Constructed.ttf' if S == 'EM' else 'HighwayEven1966-%s.ttf' % S))
FIG = dict(zip('0123456789', ['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine']))
NAME = {c: c for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'}; NAME.update(FIG)
FAMILY = 'Highway 61 Revisited %s' % ('E(M)' if S == 'EM' else S)
PS = 'Highway61Revisited-%s' % S

geoms = vgeom.series('E' if S == 'EM' else S)
order = BASE.getGlyphOrder()
gs = BASE.getGlyphSet()
hmtx = BASE['hmtx']


class SegShifted:
    """draws strokeop segment contours, shifted to keep the base glyph's side bearings"""
    def __init__(s, contours, name):
        import numpy as np
        P = np.vstack(strokeop.to_polygons(contours)); x0, x1 = P[:, 0].min(), P[:, 0].max()
        bp = BoundsPen(gs); gs[name].draw(bp)
        lsb = bp.bounds[0]; rsb = hmtx[name][0] - bp.bounds[2]
        s.c, s.dx = contours, lsb - x0
        s.width = int(round(lsb + (x1 - x0) + rsb))

    def draw(s, pen): strokeop.draw(s.c, pen, dx=s.dx)


class Shifted:
    """draws a construction shifted so it keeps the 1966 glyph's side bearings"""
    def __init__(s, geom, name, prims=None):
        s.prims = prims
        x0, y0, x1, y1 = geom.bounds
        bp = BoundsPen(gs); gs[name].draw(bp)
        lsb = bp.bounds[0]; rsb = hmtx[name][0] - bp.bounds[2]
        s.geom, s.dx = geom, lsb - x0
        s.width = int(round(lsb + (x1 - x0) + rsb))

    def draw(s, pen):
        if s.prims is None: outline.draw(s.geom, pen, dx=s.dx)
        else: lcoutline.draw(s.geom, pen, s.prims, dx=s.dx)


import lcE
A_E = 690 * (3 + 3 / 32) / 18                     # Series E stroke (table A, 3 3/32 in at 18 in)
D_EM = (690 / 5 - A_E) / 2                       # offset that takes Series E strokes to one-fifth of the height
if S == 'EM':
    k = 690 / (690 + 2 * D_EM)                   # squeeze back to baseline and cap line
    new = {}
    for c, g in geoms.items():
        C = strokeop.vmap(strokeop.offset(strokeop.from_caps(g), D_EM), [(-1e9, 1e9, k, 345 * (1 - k))])
        y0, y1 = g.bounds[1], g.bounds[3]
        P = strokeop.to_polygons(C)
        import numpy as np
        lo = min(p[:, 1].min() for p in P); hi = max(p[:, 1].max() for p in P)
        if lo < y0 - 2: C = strokeop.clip_y(C, y0, below=True)     # an acute corner pushed through the baseline (M)
        if hi > y1 + 2: C = strokeop.clip_y(C, y1, below=False)
        new[NAME[c]] = SegShifted(C, NAME[c])
    new.update({c: Shifted(g, c, prims) for c, (g, prims) in lcE.series().items()})
else:
    new = {NAME[c]: Shifted(g, NAME[c]) for c, g in geoms.items()}
if S == 'E':                                     # lower case narrowed by the same amount (the opposite change)
    X = 517.5; d = D_EM
    bands = [(-1e9, d, 1.0, -d), (d, X - d, X / (X - 2 * d), -d * X / (X - 2 * d)), (X - d, 1e9, 1.0, d)]
    for c, (g, prims) in lcE.series().items():
        C = strokeop.vmap(strokeop.offset(strokeop.from_lc(g, prims), -d), bands)
        new[c] = SegShifted(C, c)
width = {n: (new[n].width if n in new else hmtx[n][0]) for n in order}
src = lambda n: new[n] if n in new else gs[n]


DESC = {'E': '1961 Standard Alphabets for Highway Signs, Series E: capitals and figures constructed from the dimensioned '
             'drawings (straight lines and circular arcs); lower case constructed from the 1966 E(M) sheets and narrowed '
             'to Series E stroke weight; punctuation from Highway Even 1966 E.',
        'EM': 'Series E(M): 1961 Series E capitals and figures with the stroke widened to one-fifth of the letter height '
              '(1966 recommendation); lower case constructed from the 1966 E(M) sheets (lines and quarter-ellipses); '
              'punctuation from Highway 1966 E(M).'}.get(S, '1961 Standard Alphabets for Highway Signs, Series %s: '
              'capitals and figures constructed from the dimensioned drawings (straight lines and circular arcs); '
              'lowercase and punctuation from Highway Even 1966 %s.' % (S, S))


def setup(fb, outlines_kind):
    fb.setupGlyphOrder(order)
    fb.setupCharacterMap(BASE.getBestCmap())
    if outlines_kind == 'cff':
        cs = {}
        for n in order:
            pen = T2CharStringPen(width[n], None, roundTolerance=0)
            src(n).draw(ReverseContourPen(pen)); cs[n] = pen.getCharString()   # CFF winds outer contours anticlockwise
        fb.setupCFF(PS + '-Regular', {'FullName': FAMILY + ' Regular', 'FamilyName': FAMILY}, cs, {})
    else:
        glyphs = {}
        for n in order:
            pen = TTGlyphPen(None)
            src(n).draw(Cu2QuPen(pen, max_err=0.05, reverse_direction=False))
            glyphs[n] = pen.glyph()
        fb.setupGlyf(glyphs)
    lsb = {}
    for n in order:
        bp = BoundsPen(None); src(n).draw(bp); lsb[n] = int(round(bp.bounds[0])) if bp.bounds else 0
    fb.setupHorizontalMetrics({n: (width[n], lsb[n]) for n in order})
    hh = BASE['hhea']
    fb.setupHorizontalHeader(ascent=hh.ascent, descent=hh.descent)
    os2 = BASE['OS/2']
    fb.setupOS2(sTypoAscender=os2.sTypoAscender, sTypoDescender=os2.sTypoDescender, sTypoLineGap=os2.sTypoLineGap,
                usWinAscent=os2.usWinAscent, usWinDescent=os2.usWinDescent, sxHeight=os2.sxHeight,
                sCapHeight=os2.sCapHeight, achVendID='NONE', fsType=0)
    fb.setupNameTable({'familyName': FAMILY, 'styleName': 'Regular', 'uniqueFontIdentifier': PS + '-Regular-1.000',
                       'fullName': FAMILY + ' Regular', 'psName': PS + '-Regular', 'version': 'Version 1.000',
                       'description': DESC})
    fb.setupPost(keepGlyphNames=True)


fb = FontBuilder(1000, isTTF=False); setup(fb, 'cff'); fb.save(os.path.join(OUT, PS + '.otf'))
fb = FontBuilder(1000, isTTF=True); setup(fb, 'glyf'); fb.save(os.path.join(OUT, PS + '.ttf'))
print('wrote', PS + '.otf', PS + '.ttf')
