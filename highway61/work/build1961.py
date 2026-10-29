"""Highway 61 Revisited <series>: capitals and figures as exact outlines from the 1961 dimensioned drawings (vgeom + outline),
lowercase and punctuation taken over from the matching Highway Even 1966 font.  Writes an OTF (cubic arcs,
fractional coordinates) and a TTF (quadratic, integer coordinates).

usage: python3 build1961.py [series]      (A B C D E F, and EM)

E: 1961 Series E capitals and figures with the lower case narrowed from E(M) weight to Series E weight (lccondense).
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
DIST = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')   # the fonts go to the repository's top level
BASE = TTFont(os.path.join(OUT, 'Highway1966-EM-Constructed.ttf' if S == 'EM' else 'HighwayEven1966-%s.ttf' % ('B' if S == 'A' else S)))
# Series A has no 1966 font: its spacing is Series B's scaled.  Across B..F the capitals' and figures' mean side bearing
# is linear in the stroke (table A); extrapolated to Series A's stroke it is 0.766 of B's, the space 0.814 of B's.
SB_SCALE = 0.766 if S == 'A' else 1.0
SPACE_SCALE = 0.814 if S == 'A' else 1.0
FIG = dict(zip('0123456789', ['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine']))
NAME = {c: c for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'}; NAME.update(FIG)
FAMILY = 'Highway 61 Revisited %s' % ('E(M)' if S == 'EM' else S)
VERSION = '1.011'          # bump for every release: Font Book compares it to decide whether a font is a new version
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
        lsb = bp.bounds[0] * SB_SCALE; rsb = (hmtx[name][0] - bp.bounds[2]) * SB_SCALE
        s.c, s.dx = contours, lsb - x0
        s.width = int(round(lsb + (x1 - x0) + rsb))

    def draw(s, pen): strokeop.draw(s.c, pen, dx=s.dx)


class Shifted:
    """draws a construction shifted so it keeps the 1966 glyph's side bearings"""
    def __init__(s, geom, name, prims=None):
        s.prims = prims
        x0, y0, x1, y1 = geom.bounds
        bp = BoundsPen(gs); gs[name].draw(bp)
        lsb = bp.bounds[0] * SB_SCALE; rsb = (hmtx[name][0] - bp.bounds[2]) * SB_SCALE
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
if S in 'ABCDEF':                                # condensed / expanded lower case (lccondense), E narrowed to its stroke
    import lccondense
    for c, C in lccondense.series(S).items():
        new[c] = SegShifted(C, c)
W_LIGHT = 0.08                                   # the w's strokes thinner by this fraction of the series stroke
def _lighter_w(C, T):
    """offset inward by W_LIGHT*T/2, then baseline and x-height restored (four strokes in one letter read heavy)"""
    X = 517.5; d = W_LIGHT * T / 2
    bands = [(-1e9, d, 1.0, -d), (d, X - d, X / (X - 2 * d), -d * X / (X - 2 * d)), (X - d, 1e9, 1.0, d)]
    return strokeop.vmap(strokeop.offset(C, -d), bands)
if S in 'ABCDEF':
    new['w'] = SegShifted(_lighter_w(new['w'].c, lccondense.stroke(S)), 'w')
else:
    g, prims = lcE.series()['w']
    new['w'] = SegShifted(_lighter_w(strokeop.from_lc(g, prims), (1 + 9 / 16) * 690 / 8), 'w')
import punct
from shapely.geometry import Polygon as _Poly
def spacing_stroke(S):
    import spacingcheck; return spacingcheck.STROKE[S]
def _glyph_poly(n):
    o = new[n]
    if isinstance(o, SegShifted):
        g = _Poly()
        for r in strokeop.to_polygons(o.c, step=0.004):
            if len(r) > 2: g = g.symmetric_difference(_Poly(r).buffer(0))
        return g
    return o.geom
for n, C in punct.make(_glyph_poly).items():
    new[n] = SegShifted(C, n)
width = {n: (new[n].width if n in new else int(round(hmtx[n][0] * (SPACE_SCALE if n == 'space' else 1)))) for n in order}


class Points:
    """draws one closed polygon (corner points, anticlockwise)"""
    def __init__(s, pts): s.p = [(round(x, 3), round(y, 3)) for x, y in pts]

    def draw(s, pen):
        pen.moveTo(s.p[0])
        for q in s.p[1:]: pen.lineTo(q)
        pen.closePath()


# the slash (not in the base fonts): the 1958 Interstate manual's, as thick as the capitals' stroke, spaced code 1
import lccondense
_T = 690 / 5 if S == 'EM' else lccondense.A[S]
_pts, width['slash'] = punct.slash(_T, lccondense.F['E' if S == 'EM' else S], 1.5 * spacing_stroke(S))
new['slash'] = Points(_pts); order = order + ['slash']

# spacing by the standard tables (spacing.py): side bearings, and class kerning for the exact pair spaces
import spacing, string as _string
_gn = lambda ch: spacing.NAME.get(ch, ch)
_cur, _box = {}, {}
for ch in _string.ascii_letters + _string.digits:
    _bp = BoundsPen(None); new[_gn(ch)].draw(_bp); x0, _, x1, _ = _bp.bounds
    _box[ch] = (x0, x1); _cur[ch] = (x0, width[_gn(ch)] - x1)
_geo = spacing.geometry(lambda ch: _glyph_poly(_gn(ch)), _string.ascii_letters + _string.digits)
_sb, KERN = spacing.solve(S, _cur, geo=_geo)
for ch, (l, r) in _sb.items():
    n = _gn(ch); x0, x1 = _box[ch]
    new[n].dx += l - x0; width[n] = int(round(l + (x1 - x0) + r))
FEA = spacing.fea(KERN, _gn)
src = lambda n: new[n] if n in new else gs[n]


DESC = ('1961 Standard Alphabets for Highway Signs, Series %s: capitals and figures constructed from the dimensioned '
        'drawings (straight lines and circular arcs); lower case from the 1966 E(M) constructions (lines and '
        'quarter-ellipses), its width and stroke set for the series; punctuation made from the letters.' % S)
if S == 'A':
    DESC += ' Series A was dropped in 1966 and never had a lower case: this one follows the rule of the others; spacing ' \
            'is Series B\'s scaled to Series A\'s stroke.'
if S == 'EM':
    DESC = ('Series E(M): 1961 Series E capitals and figures with the stroke widened to one-fifth of the letter height '
            '(1966 recommendation); lower case constructed from the 1966 E(M) sheets (lines and quarter-ellipses); '
            'punctuation made from the letters.')


def setup(fb, outlines_kind):
    fb.setupGlyphOrder(order)
    fb.setupCharacterMap({**BASE.getBestCmap(), ord('/'): 'slash'})
    if outlines_kind == 'cff':
        cs = {}
        for n in order:
            pen = T2CharStringPen(width[n], None, roundTolerance=0)
            src(n).draw(ReverseContourPen(pen)); cs[n] = pen.getCharString()   # CFF winds outer contours anticlockwise
        fb.setupCFF(PS + '-Regular', {'FullName': FAMILY + ' Regular', 'FamilyName': FAMILY, 'version': VERSION}, cs, {})
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
    fb.setupNameTable({'familyName': FAMILY, 'styleName': 'Regular', 'uniqueFontIdentifier': PS + '-Regular-' + VERSION,
                       'fullName': FAMILY + ' Regular', 'psName': PS + '-Regular', 'version': 'Version ' + VERSION,
                       'description': DESC})
    fb.font['head'].fontRevision = float(VERSION)
    fb.setupPost(keepGlyphNames=True)
    from fontTools.feaLib.builder import addOpenTypeFeaturesFromString
    addOpenTypeFeaturesFromString(fb.font, FEA)


if __name__ == '__main__':                        # lcspace.py imports this for the spacing inputs only
    fb = FontBuilder(1000, isTTF=False); setup(fb, 'cff'); fb.save(os.path.join(DIST, 'otf', PS + '.otf'))
    fb = FontBuilder(1000, isTTF=True); setup(fb, 'glyf'); fb.save(os.path.join(DIST, 'ttf', PS + '.ttf'))
    print('wrote', PS + '.otf', PS + '.ttf')
