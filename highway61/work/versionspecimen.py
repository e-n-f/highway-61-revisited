"""The versions of each series' lower case, one under another: python3 versionspecimen.py
-> ../../specimen/specimen-lc-versions.png

Main (lccondense.F_LC, spacing.LC_TIGHT: the a-z as wide as FHWA 2000's), Matched (one width factor and one spacing
factor per series making the a-z's ink widths and advances, s and w left out, add up to FHWA's; lcmatch.py), the FHWA
2000 font itself where there is one (B-F, and Series E Modified beside E(M); fhwa2otf.py, skipped if its fonts are
not there); see variant.py.  The Matched fonts are not in the repository at present: build them first
(build1961.py <S> matched), or their rows are left out.  All at equal cap height.  Words are set with each font's
advances; the FHWA subsets have no space, so they get the main font's space of the same series, in cap heights (as
in fhwaspecimen.py)."""
import os
from PIL import Image, ImageDraw, ImageFont
from fhwaspecimen import info
W = os.path.dirname(os.path.abspath(__file__)); DIST = os.path.join(W, '..', '..'); FH = os.path.join(W, '..', 'fhwa2000')
CAP = 97; LINE = 150; TEXTS = ['abcdefghijklmnopqrstuvwxyz', 'Bayshore Fwy Business District']


def versions(S):
    """[(label, font file, space in cap heights)] for series S"""
    ours = lambda v: os.path.join(DIST, 'otf', 'Highway61Revisited%s-%s.otf' % (v, S))
    sp = info(ours(''))[1]
    out = [('Main', ours(''), sp)]
    if os.path.exists(ours('Matched')):                 # (built only on request: build1961.py <S> matched)
        out.append(('Matched', ours('Matched'), info(ours('Matched'))[1]))
    fh = os.path.join(FH, 'Series%s2000-Normal.otf' % ('ModE' if S == 'EM' else S))
    if os.path.exists(fh): out.append(('FHWA 2000', fh, sp))
    return out


def setline(d, x, y, fn, text, space):
    cap, _ = info(fn); f = ImageFont.truetype(fn, int(round(CAP / cap)))
    for w in text.split(' '):
        d.text((x, y), w, font=f, fill='black', anchor='ls'); x += d.textlength(w, font=f) + space * CAP


def main():
    SER = ['A', 'B', 'C', 'D', 'E', 'EM', 'F']; V = {S: versions(S) for S in SER}
    heights = {S: 90 + len(V[S]) * 2 * LINE + 40 for S in SER}
    im = Image.new('RGB', (3000, sum(heights.values()) + 20), 'white'); d = ImageDraw.Draw(im)
    lab = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 36)
    small = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 26)
    y0 = 0
    for S in SER:
        y = y0 + 60
        d.text((20, y), 'Series ' + ('E(M)' if S == 'EM' else S), font=lab, fill=(80, 80, 80), anchor='ls')
        for j, (name, fn, sp) in enumerate(V[S]):
            yy = y + 30 + j * 2 * LINE
            d.text((20, yy + 70), name, font=small, fill=(150, 150, 150), anchor='ls')
            for k, t in enumerate(TEXTS): setline(d, 230, yy + 110 + k * LINE, fn, t, sp)
        y0 += heights[S]
    im.save(os.path.join(DIST, 'specimen', 'specimen-lc-versions.png'))


if __name__ == '__main__': main()
