"""Highway 61 Revisited <series>: capitals and figures as exact outlines from the 1961 dimensioned drawings (vgeom + outline),
lowercase and punctuation taken over from the matching Highway Even 1966 font.  Writes an OTF (cubic arcs,
fractional coordinates) and a TTF (quadratic, integer coordinates).

usage: python3 build1961.py [series]      (E and D are constructed so far)
"""
import sys, os
from fontTools.ttLib import TTFont
from fontTools.fontBuilder import FontBuilder
from fontTools.pens.t2CharStringPen import T2CharStringPen
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.pens.cu2quPen import Cu2QuPen
from fontTools.pens.boundsPen import BoundsPen
from fontTools.pens.reverseContourPen import ReverseContourPen
import vgeom, outline, lcoutline

S = sys.argv[1] if len(sys.argv) > 1 else 'E'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'outputs')
BASE = TTFont(os.path.join(OUT, 'HighwayEven1966-%s.ttf' % S))
FIG = dict(zip('0123456789', ['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine']))
NAME = {c: c for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'}; NAME.update(FIG)
FAMILY = 'Highway 61 Revisited %s' % S
PS = 'Highway61Revisited-%s' % S

geoms = vgeom.series(S)
order = BASE.getGlyphOrder()
gs = BASE.getGlyphSet()
hmtx = BASE['hmtx']


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


new = {NAME[c]: Shifted(g, NAME[c]) for c, g in geoms.items()}
if S == 'E':                                     # lower case constructed so far (lcE: quarter-ellipses and lines)
    import lcE
    new.update({c: Shifted(g, c, prims) for c, (g, prims) in lcE.series().items()})
width = {n: (new[n].width if n in new else hmtx[n][0]) for n in order}
src = lambda n: new[n] if n in new else gs[n]


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
                       'description': '1961 Standard Alphabets for Highway Signs, Series %s: capitals and figures constructed '
                                      'from the dimensioned drawings (straight lines and circular arcs); lowercase and '
                                      'punctuation from Highway Even 1966 %s.' % (S, S)})
    fb.setupPost(keepGlyphNames=True)


fb = FontBuilder(1000, isTTF=False); setup(fb, 'cff'); fb.save(os.path.join(OUT, PS + '.otf'))
fb = FontBuilder(1000, isTTF=True); setup(fb, 'glyf'); fb.save(os.path.join(OUT, PS + '.ttf'))
print('wrote', PS + '.otf', PS + '.ttf')
